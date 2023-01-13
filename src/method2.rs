use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::{One, UniformRand};
use merlin::Transcript;
use std::{
    ops::{Div, Mul, Sub},
    usize,
};

use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_std::rand::RngCore;

use crate::{
    get_challenge, get_field_size, transcribe_generic, transcribe_points_and_evals, MultiOpenKzg,
};

use super::{
    gen_curve_powers, gen_powers, lagrange_interp, linear_combination, poly_div_q_r,
    vanishing_polynomial, Error,
};

pub struct Setup<E: Pairing> {
    powers_of_g1: Vec<E::G1Affine>,
    g2: E::G2Affine,
    g2x: E::G2Affine,
}

#[derive(Debug)]
pub struct Commitment<E: Pairing>(E::G1Affine);
#[derive(Debug)]
pub struct Proof<E: Pairing>(E::G1Affine, E::G1Affine);

impl<E: Pairing> MultiOpenKzg<E> for Setup<E> {
    type Commitment = Commitment<E>;
    type Proof = Proof<E>;

    fn new(max_degree: usize, _max_pts: usize, rng: &mut impl RngCore) -> Setup<E> {
        let num_scalars = max_degree + 1;

        let x = E::ScalarField::rand(rng);
        let x_powers = gen_powers(x, num_scalars);

        let powers_of_g1 = gen_curve_powers::<E::G1>(x_powers.as_ref(), rng);
        let g2 = E::G2::rand(rng).into_affine();
        let g2x = (g2 * x).into_affine();

        Setup {
            powers_of_g1,
            g2,
            g2x,
        }
    }

    fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error> {
        let res = super::curve_msm::<E::G1>(&self.powers_of_g1, poly.as_ref())?;
        Ok(Commitment(res.into_affine()))
    }

    fn open(
        &self,
        transcript: &mut Transcript,
        polys: &[impl AsRef<[E::ScalarField]>],
        evals: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
    ) -> Result<Proof<E>, Error> {
        let field_size_bytes = get_field_size::<E::ScalarField>();
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;

        let gamma = get_challenge(transcript, b"open gamma", field_size_bytes);

        let gammas = gen_powers::<E::ScalarField>(gamma, self.powers_of_g1.len());
        let gamma_fis = linear_combination::<E::ScalarField>(polys, &gammas)
            .ok_or(Error::NoPolynomialsGiven)?;
        let gamma_fis_poly = DensePolynomial::from_coefficients_vec(gamma_fis);

        let z_s = vanishing_polynomial(points.as_ref());
        let (h, gamma_ris_over_zs) = poly_div_q_r((&gamma_fis_poly).into(), (&z_s).into())?;

        let w_1 = super::curve_msm::<E::G1>(&self.powers_of_g1, &h)?.into_affine();

        transcribe_generic(transcript, b"open W", &w_1)?;
        let chal_z = get_challenge(transcript, b"open z", field_size_bytes);

        let gamma_ri_z = DensePolynomial::from_coefficients_vec(gamma_ris_over_zs)
            .mul(&z_s)
            .evaluate(&chal_z);

        let f_z = gamma_fis_poly.sub(&DensePolynomial::from_coefficients_vec(vec![gamma_ri_z])); // XXX
        let l = f_z.sub(&DensePolynomial::from_coefficients_vec(h).mul(z_s.evaluate(&chal_z)));

        let x_minus_z =
            DensePolynomial::from_coefficients_vec(vec![-chal_z, E::ScalarField::one()]);
        let l_quotient = l.div(&x_minus_z);

        let w_2 = super::curve_msm::<E::G1>(&self.powers_of_g1, &l_quotient)?.into_affine();
        Ok(Proof(w_1, w_2))
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        points: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        let field_size_bytes = get_field_size::<E::ScalarField>();
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;

        let gamma = get_challenge(transcript, b"open gamma", field_size_bytes);
        transcribe_generic(transcript, b"open W", &proof.0)?;
        let chal_z = get_challenge(transcript, b"open z", field_size_bytes);

        let zeros = vanishing_polynomial(points);
        let zeros_z = zeros.evaluate(&chal_z);

        // Get the r_i polynomials with lagrange interp. These could be precomputed.
        let gammas = gen_powers(gamma, evals.len());
        let ri_s = lagrange_interp(evals, points);

        // Aggregate the r_is and then evaluate at chal_z
        let gamma_ris =
            linear_combination(&ri_s.iter().map(|i| &i.coeffs).collect::<Vec<_>>(), &gammas)
                .ok_or(Error::NoPolynomialsGiven)?;
        let gamma_ris_z = DensePolynomial::from_coefficients_vec(gamma_ris).evaluate(&chal_z);
        let gamma_ris_z_pt = self.powers_of_g1[0].mul(gamma_ris_z);

        // Then do a single msm of the gammas and commitments
        let cms = commits.iter().map(|i| i.0).collect::<Vec<_>>();
        let gamma_cm_pt = super::curve_msm::<E::G1>(&cms, gammas.as_ref())?;

        let f = gamma_cm_pt - gamma_ris_z_pt - proof.0.mul(zeros_z);

        let x_minus_z = self.g2x.into_group() - self.g2.into_group().mul(&chal_z);
        Ok(E::pairing(f, self.g2) == E::pairing(proof.1, x_minus_z))
    }
}

#[cfg(test)]
mod tests {
    use super::Setup;
    use crate::{test_rng, MultiOpenKzg};
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::UniformRand;
    use merlin::Transcript;

    #[test]
    fn test_basic_open_works() {
        let s = Setup::<Bls12_381>::new(256, 0, &mut test_rng());
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
        let mut open_transcript = Transcript::new(b"testing");
        let open = s
            .open(&mut open_transcript, &coeffs, &evals, &points)
            .expect("Open failed");

        let mut verify_transcript = Transcript::new(b"testing");
        assert_eq!(
            Ok(true),
            s.verify(&mut verify_transcript, &commits, &points, &evals, &open)
        );
    }
}
