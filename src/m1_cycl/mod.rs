//! Precomputation for Method 1 where each point set is a cyclic subgroup of the evaluation domain, with blst optimizations
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::Zero;
use core::ops::Deref;

use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use ark_std::vec::Vec;
use core::ops::Mul;
use merlin::Transcript;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::method1::{M1NoPrecomp, Proof};
use crate::poly_ops::{ev_points, SplitEvalDomain};
use crate::traits::{Committer, MSMEngine, PolyMultiProof};
use crate::{
    cfg_iter, check_opening_sizes, check_verify_sizes, gen_powers, get_challenge, get_field_size,
    linear_combination, transcribe_points_and_evals, Commitment, Error,
};

/// Method 1 with blst optimization and precomputed lagrange polynomials/vanishing polys
#[derive(Clone)]
pub struct M1CyclPrecomp<E: Pairing, M: MSMEngine<E = E>> {
    /// The inner method 1 object without precomputation
    pub inner: M1NoPrecomp<E, M>,
    split_domain: SplitEvalDomain<E::ScalarField>,
    point_set_groups: Vec<Radix2EvaluationDomain<E::ScalarField>>,
    num_point_sets: usize,
    base_size: usize,
    g2_zeros: Vec<E::G2Affine>,
}

fn is_power_of_two(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
}

impl<E: Pairing, M: MSMEngine<E = E>> M1CyclPrecomp<E, M> {
    /// Make a precompute-optimized version of a method 1 object for the given sets of points
    pub fn from_inner(
        inner: M1NoPrecomp<E, M>,
        base_size: usize,
        num_point_sets: usize,
    ) -> Result<Self, Error> {
        if !is_power_of_two(base_size) || !is_power_of_two(num_point_sets) {
            todo!() // return error
        }
        if inner.powers_of_g1.len() < base_size {
            return Err(Error::DomainConstructionFailed(base_size));
        }
        let split_domain = SplitEvalDomain::<E::ScalarField>::new(base_size, num_point_sets)
            .ok_or(Error::DomainConstructionFailed(inner.powers_of_g1.len()))?;
        let point_set_groups = split_domain.subgroups();
        let vanishing_polys: Vec<_> = cfg_iter!(point_set_groups)
            .map(|(_, sg)| sg.vanishing_polynomial())
            .collect();
        let g2_zeros = cfg_iter!(vanishing_polys)
            .map(|(_, p)| {
                let coeffs = p.deref();
                let mut accum = E::G2::zero();
                for (i0, p0) in coeffs {
                    accum += inner
                        .powers_of_g2
                        .get(*i0)
                        .ok_or(Error::TooManyScalars {
                            n_coeffs: inner.powers_of_g1.len(),
                            expected_max: i0 + 1,
                        })?
                        .mul(p0);
                }
                Ok(accum.into_affine())
            })
            .collect::<Result<Vec<_>, Error>>()?;

        Ok(Self {
            inner,
            base_size,
            split_domain,
            point_set_groups,
            num_point_sets,
            g2_zeros,
        })
    }

    /// Returns the SplitEvalDomain beign used for the multiproof scheme.
    /// In order to figure out which points map to which evaluation index, you should use this
    /// object
    pub fn point_sets(&self) -> &SplitEvalDomain<E::ScalarField> {
        &self.split_domain
    }
}

impl<E: Pairing, M: MSMEngine<E = E>> Committer<E> for M1CyclPrecomp<E, M> {
    fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error> {
        self.inner.commit(poly)
    }
}

impl<E: Pairing, M: MSMEngine<E = E>> PolyMultiProof<E> for M1CyclPrecomp<E, M> {
    type Proof = Proof<E>;

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        point_set_index: usize,
    ) -> Result<Self::Proof, Error> {
        check_opening_sizes(evals, polys, self.base_size / self.num_point_sets)?;

        // Commit the evals and the points to the transcript
        // TODO: better error
        let subgroup = self
            .point_set_groups
            .get(point_set_index)
            .ok_or(Error::NoPointsGiven)?;
        let points = ev_points(subgroup);
        let field_size_bytes = get_field_size::<E::ScalarField>();
        transcribe_points_and_evals(transcript, &points, evals, field_size_bytes)?;

        // Read the challenge
        let gamma = get_challenge::<E::ScalarField>(transcript, b"open gamma", field_size_bytes);
        // Make the gamma powers
        let gammas = gen_powers::<E::ScalarField>(gamma, self.inner.powers_of_g1.len());
        // Take a linear combo of gammas with the polynomials
        let fsum = linear_combination::<E::ScalarField>(polys, &gammas)
            .ok_or(Error::NoPolynomialsGiven)?;

        // Polynomial divide, the remained would contain the gamma * ri_s,
        // The result is the correct quotient
        let (q, _) = DensePolynomial::from_coefficients_vec(fsum)
            .divide_by_vanishing_poly(*subgroup)
            .expect("This always succeeds");
        // Open to the resulting polynomial
        Ok(Proof(
            M::multi_scalar_mul_g1(&self.inner.g1_precomp, &q.deref())?.into_affine(),
        ))
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        point_set_index: usize,
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error> {
        check_verify_sizes(commits, evals, self.base_size / self.num_point_sets)?;

        let field_size_bytes = get_field_size::<E::ScalarField>();
        // TODO: better error
        let subgroup = self
            .point_set_groups
            .get(point_set_index)
            .ok_or(Error::NoPointsGiven)?;
        let points = ev_points(subgroup);
        transcribe_points_and_evals(transcript, &points, evals, field_size_bytes)?;
        let gamma = get_challenge(transcript, b"open gamma", field_size_bytes);
        // Aggregate the r_is and then do a single msm of just the ri's and gammas
        let gammas = gen_powers(gamma, evals.len());

        // We first get the values of sum_i gamma^i-1 r_i,j (z_j)
        let mut gamma_ris = linear_combination(&evals, &gammas).expect("TODO");
        // Then we find the coefficients
        subgroup.ifft_in_place(&mut gamma_ris);
        let gamma_ris_pt = M::multi_scalar_mul_g1(&self.inner.g1_precomp, &gamma_ris)?;

        // Then do a single msm of the gammas and commitments
        let cms_prep = M::prepare_g1(commits.into_iter().map(|i| i.0).collect());
        let gamma_cm_pt = M::multi_scalar_mul_g1(&cms_prep, &gammas)?;

        let g2 = self.inner.powers_of_g2[0];

        let lhsg1 = (gamma_cm_pt - gamma_ris_pt).into_affine();
        let lhsg2 = g2;
        let rhsg1 = proof.0;
        let rhsg2 = self.g2_zeros[point_set_index];
        Ok(M::pairing_eq_check(lhsg1, lhsg2, rhsg1, rhsg2))
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
    use merlin::Transcript;

    use super::M1CyclPrecomp;
    use crate::{
        method1::M1NoPrecomp,
        msm::ArkMSMEngine,
        poly_ops::ev_points,
        test_rng,
        testing::test_basic_precomp,
        traits::{Committer, PolyMultiProof},
    };

    #[test]
    fn test_basic_open_works() {
        let s = <M1NoPrecomp<Bls12_381, ArkMSMEngine<Bls12_381>>>::new(256, 256, &mut test_rng());
        let s = M1CyclPrecomp::from_inner(s, 256, 2).expect("Failed to construct");
        let points = ev_points(&s.point_set_groups[0]);
        test_basic_precomp(&s, &points);
    }

    #[test]
    fn test_complex_open_works() {
        let s = <M1NoPrecomp<Bls12_381, ArkMSMEngine<Bls12_381>>>::new(256, 256, &mut test_rng());
        let s = M1CyclPrecomp::from_inner(s, 256, 2).expect("Failed to construct");
        let polys = (0..2)
            .map(|_| DensePolynomial::<Fr>::rand(255, &mut test_rng()).coeffs)
            .collect::<Vec<_>>();
        let commits = polys
            .iter()
            .map(|p| s.commit(p).expect("Commit failed"))
            .collect::<Vec<_>>();
        let points = ev_points(s.split_domain.base());
        let naive_evals = polys
            .iter()
            .map(|poly| {
                points
                    .iter()
                    .map(|p| DensePolynomial::from_coefficients_vec(poly.clone()).evaluate(p))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let evals = polys
            .iter()
            .map(|p| s.split_domain.base().fft(p))
            .collect::<Vec<_>>();
        assert_eq!(evals, naive_evals);
        for gi in 0..s.point_set_groups.len() {
            let trimmed_evals: Vec<_> = naive_evals
                .iter()
                .map(|ev| {
                    s.split_domain
                        .take_subgroup_indices(gi, ev.clone())
                        .unwrap()
                })
                .collect();
            let proof = s
                .open(&mut Transcript::new(b"test"), &trimmed_evals, &polys, gi)
                .expect("Failed to open");
            assert_eq!(
                Ok(true),
                s.verify(
                    &mut Transcript::new(b"test"),
                    &commits,
                    gi,
                    &trimmed_evals,
                    &proof
                )
            );
        }
    }
}
