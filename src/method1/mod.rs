//! # BDFG Method 1
//! This contains a pure ark implementation of BDFG21 method 1
use crate::{
    check_opening_sizes, check_verify_sizes,
    lagrange::LagrangeInterpContext,
    traits::{Committer, PolyMultiProofNoPrecomp},
};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{vec::Vec, UniformRand};
use merlin::Transcript;

use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::rand::RngCore;

use crate::{get_challenge, get_field_size, transcribe_points_and_evals, Commitment};

use super::{
    gen_curve_powers, gen_powers, linear_combination, poly_div_q_r, vanishing_polynomial, Error,
};

pub mod precompute;

/// A method 1 proof scheme with no precomputation of lagrange polynomials
#[derive(Clone, Debug)]
pub struct M1NoPrecomp<E: Pairing> {
    /// The given powers tau in G1
    pub powers_of_g1: Vec<E::G1Affine>,
    /// The given powers tau in G2
    pub powers_of_g2: Vec<E::G2Affine>,
}

/// A method 1 proof
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: Pairing>(pub E::G1Affine);

impl<E: Pairing> M1NoPrecomp<E> {
    /// Make a new random scheme
    pub fn new(max_coeffs: usize, max_pts: usize, rng: &mut impl RngCore) -> Self {
        let x = E::ScalarField::rand(rng);
        let g1 = E::G1::rand(rng);
        let g2 = E::G2::rand(rng);
        Self::new_from_scalar(x, g1, g2, max_coeffs, max_pts)
    }

    /// Make a new scheme from a given secret scalar
    pub fn new_from_scalar(
        x: E::ScalarField,
        g1: E::G1,
        g2: E::G2,
        max_coeffs: usize,
        max_pts: usize,
    ) -> Self {
        let n_g2_powers = max_pts + 1;
        let x_powers = gen_powers(x, core::cmp::max(max_coeffs, n_g2_powers));

        let powers_of_g1 = gen_curve_powers::<E::G1>(x_powers.as_ref(), g1);
        let powers_of_g2 = gen_curve_powers::<E::G2>(x_powers[..n_g2_powers].as_ref(), g2);

        Self::new_from_affine(powers_of_g1, powers_of_g2)
    }

    /// Make a new scheme from the given projective powers
    pub fn new_from_powers(powers_of_g1: &[E::G1], powers_of_g2: &[E::G2]) -> Self {
        Self {
            powers_of_g1: powers_of_g1.iter().map(|s| s.into_affine()).collect(),
            powers_of_g2: powers_of_g2.iter().map(|s| s.into_affine()).collect(),
        }
    }

    /// Make a new scheme from the given powers in affine form
    pub fn new_from_affine(powers_of_g1: Vec<E::G1Affine>, powers_of_g2: Vec<E::G2Affine>) -> Self {
        Self {
            powers_of_g1,
            powers_of_g2,
        }
    }

    fn open_with_vanishing_poly(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
        vp: &DensePolynomial<E::ScalarField>,
    ) -> Result<Proof<E>, Error> {
        // Check sizes
        check_opening_sizes(evals, polys, points.len())?;
        // Commit the evals and the points to the transcript
        let field_size_bytes = get_field_size::<E::ScalarField>();
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;

        // Read the challenge
        let gamma = get_challenge::<E::ScalarField>(transcript, b"open gamma", field_size_bytes);
        // Make the gamma powers
        let gammas = gen_powers::<E::ScalarField>(gamma, self.powers_of_g1.len());
        // Take a linear combo of gammas with the polynomials
        let fsum = linear_combination::<E::ScalarField>(polys, &gammas)
            .ok_or(Error::NoPolynomialsGiven)?;

        // Polynomial divide, the remained would contain the gamma * ri_s,
        // The result is the correct quotient
        let (q, _) = poly_div_q_r(
            DensePolynomial::from_coefficients_vec(fsum).into(),
            vp.into(),
        )?;
        // Open to the resulting polynomial
        Ok(Proof(
            super::curve_msm::<E::G1>(&self.powers_of_g1, &q)?.into_affine(),
        ))
    }

    fn verify_with_lag_ctx_g2_zeros(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        points: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
        lag_ctx: &LagrangeInterpContext<E::ScalarField>,
        g2_zeros: &E::G2,
    ) -> Result<bool, Error> {
        check_verify_sizes(commits, evals, points.len())?;

        let field_size_bytes = get_field_size::<E::ScalarField>();
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;
        let gamma = get_challenge(transcript, b"open gamma", field_size_bytes);
        // Aggregate the r_is and then do a single msm of just the ri's and gammas
        let gammas = gen_powers(gamma, evals.len());

        // Get the gamma^i r_i polynomials with lagrange interp. This does both the lagrange interp
        // and the gamma mul in one step so we can just lagrange interp once.
        let gamma_ris = lag_ctx.lagrange_interp_linear_combo(evals, &gammas)?.coeffs;
        let gamma_ris_pt = super::curve_msm::<E::G1>(&self.powers_of_g1, gamma_ris.as_ref())?;

        // Then do a single msm of the gammas and commitments
        let cms = commits.iter().map(|i| i.0).collect::<Vec<_>>();
        let gamma_cm_pt = super::curve_msm::<E::G1>(&cms, gammas.as_ref())?;

        let g2 = self.powers_of_g2[0];

        Ok(E::pairing(gamma_cm_pt - gamma_ris_pt, g2) == E::pairing(proof.0, g2_zeros))
    }
}

impl<E: Pairing> Committer<E> for M1NoPrecomp<E> {
    fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error> {
        let res = super::curve_msm::<E::G1>(&self.powers_of_g1, poly.as_ref())?;
        Ok(Commitment(res.into_affine()))
    }
}

impl<E: Pairing> PolyMultiProofNoPrecomp<E> for M1NoPrecomp<E> {
    type Proof = Proof<E>;

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
    ) -> Result<Proof<E>, Error> {
        let vp = vanishing_polynomial(points.as_ref());
        self.open_with_vanishing_poly(transcript, evals, polys, points, &vp)
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        points: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        let vp = vanishing_polynomial(points);
        let g2_zeros = super::curve_msm::<E::G2>(&self.powers_of_g2, &vp)?;
        let lag_ctx = LagrangeInterpContext::new_from_points(points)?;
        self.verify_with_lag_ctx_g2_zeros(
            transcript, commits, points, evals, proof, &lag_ctx, &g2_zeros,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::M1NoPrecomp;
    use crate::{
        test_rng,
        testing::{test_basic_no_precomp, test_size_errors},
    };
    use ark_bls12_381::Bls12_381;

    #[test]
    fn test_basic_open_works() {
        let s = M1NoPrecomp::<Bls12_381>::new(256, 30, &mut test_rng());
        test_basic_no_precomp(&s);
        test_size_errors(&s);
    }
}
