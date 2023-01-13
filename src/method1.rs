use ark_poly::univariate::DensePolynomial;
use ark_std::UniformRand;
use merlin::Transcript;
use std::usize;

use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::rand::RngCore;

use crate::{get_challenge, get_field_size, transcribe_points_and_evals, MultiOpenKzg};

use super::{
    gen_curve_powers, gen_powers, lagrange_interp, linear_combination, poly_div_q_r,
    vanishing_polynomial, Error,
};

pub struct Setup<E: Pairing> {
    powers_of_g1: Vec<E::G1Affine>,
    powers_of_g2: Vec<E::G2Affine>,
}

#[derive(Debug)]
pub struct Commitment<E: Pairing>(E::G1Affine);
#[derive(Debug)]
pub struct Proof<E: Pairing>(E::G1Affine);

impl<E: Pairing> MultiOpenKzg<E> for Setup<E> {
    type Commitment = Commitment<E>;
    type Proof = Proof<E>;

    fn new(max_degree: usize, max_pts: usize, rng: &mut impl RngCore) -> Setup<E> {
        let num_scalars = max_degree + 1;

        let x = E::ScalarField::rand(rng);
        let x_powers = gen_powers(x, num_scalars);

        let powers_of_g1 = gen_curve_powers::<E::G1>(x_powers.as_ref(), rng);
        let powers_of_g2 = gen_curve_powers::<E::G2>(x_powers[..max_pts + 1].as_ref(), rng);

        Setup {
            powers_of_g1,
            powers_of_g2,
        }
    }

    fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error> {
        let res = super::curve_msm::<E::G1>(&self.powers_of_g1, poly.as_ref())?;
        Ok(Commitment(res.into_affine()))
    }

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
    ) -> Result<Proof<E>, Error> {
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

        // Compute Z_s
        let z_s = vanishing_polynomial(points.as_ref());
        // Polynomial divide, the remained would contain the gamma * ri_s,
        // The result is the correct quotient
        let (q, _) = poly_div_q_r(DensePolynomial { coeffs: fsum }.into(), z_s.into())?;
        // Open to the resulting polynomial
        Ok(Proof(
            super::curve_msm::<E::G1>(&self.powers_of_g1, &q)?.into_affine(),
        ))
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        points: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        let zeros = vanishing_polynomial(points);
        let zeros = super::curve_msm::<E::G2>(&self.powers_of_g2, &zeros)?;

        // Get the r_i polynomials with lagrange interp. These could be precomputed.
        let ri_s = lagrange_interp(evals, points);

        let field_size_bytes = get_field_size::<E::ScalarField>();
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;
        let gamma = get_challenge(transcript, b"open gamma", field_size_bytes);
        // Aggregate the r_is and then do a single msm of just the ri's and gammas
        let gammas = gen_powers(gamma, evals.len());
        let gamma_ris =
            linear_combination(&ri_s.iter().map(|i| &i.coeffs).collect::<Vec<_>>(), &gammas)
                .ok_or(Error::NoPolynomialsGiven)?;
        let gamma_ris_pt = super::curve_msm::<E::G1>(&self.powers_of_g1, gamma_ris.as_ref())?;

        // Then do a single msm of the gammas and commitments
        let cms = commits.iter().map(|i| i.0).collect::<Vec<_>>();
        let gamma_cm_pt = super::curve_msm::<E::G1>(&cms, gammas.as_ref())?;

        let g2 = self.powers_of_g2[0];

        Ok(E::pairing(gamma_cm_pt - gamma_ris_pt, g2) == E::pairing(proof.0, zeros))
    }
}

#[cfg(test)]
mod tests {
    use super::Setup;
    use crate::test_rng;
    use crate::MultiOpenKzg;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::UniformRand;
    use merlin::Transcript;

    #[test]
    fn test_basic_open_works() {
        let s = Setup::<Bls12_381>::new(256, 32, &mut test_rng());
        let points = (0..30)
            .map(|_| Fr::rand(&mut test_rng()))
            .collect::<Vec<_>>();
        let polys = (0..20)
            .map(|_| DensePolynomial::<Fr>::rand(50, &mut test_rng()))
            .collect::<Vec<_>>();
        let evals: Vec<Vec<_>> = polys
            .iter()
            .map(|p| points.iter().map(|x| p.evaluate(x)).collect())
            .collect();
        let coeffs = polys.iter().map(|p| p.coeffs.clone()).collect::<Vec<_>>();
        let commits = coeffs
            .iter()
            .map(|p| s.commit(p).expect("Commit failed"))
            .collect::<Vec<_>>();
        let mut transcript = Transcript::new(b"testing");
        let open = s
            .open(&mut transcript, &evals, &coeffs, &points)
            .expect("Open failed");
        let mut transcript = Transcript::new(b"testing");
        assert_eq!(
            Ok(true),
            s.verify(&mut transcript, &commits, &points, &evals, &open)
        );
    }
}
