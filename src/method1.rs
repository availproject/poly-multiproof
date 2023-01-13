use ark_poly::univariate::DensePolynomial;
use ark_std::UniformRand;
use std::usize;

use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::rand::RngCore;

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

impl<E: Pairing> Setup<E> {
    pub fn new(max_degree: usize, max_pts: usize, rng: &mut impl RngCore) -> Setup<E> {
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

    pub fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error> {
        let res = super::curve_msm::<E::G1>(&self.powers_of_g1, poly.as_ref())?;
        Ok(Commitment(res.into_affine()))
    }

    pub fn open(
        &self,
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
        challenge: E::ScalarField,
    ) -> Result<Proof<E>, Error> {
        let gammas = gen_powers::<E::ScalarField>(challenge, self.powers_of_g1.len());
        let fsum = linear_combination::<E::ScalarField>(polys, &gammas)
            .ok_or(Error::NoPolynomialsGiven)?;

        let z_s = vanishing_polynomial(points.as_ref());
        let (q, _) = poly_div_q_r(DensePolynomial { coeffs: fsum }.into(), z_s.into())?;
        Ok(Proof(self.commit(q)?.0))
    }

    pub fn verify(
        &self,
        commits: &[Commitment<E>],
        pts: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
        challenge: E::ScalarField,
    ) -> Result<bool, Error> {
        let zeros = vanishing_polynomial(pts);
        let zeros = super::curve_msm::<E::G2>(&self.powers_of_g2, &zeros)?;

        // Get the r_i polynomials with lagrange interp. These could be precomputed.
        let ri_s = lagrange_interp(evals, pts);

        // Aggregate the r_is and then do a single msm of just the ri's and gammas
        let gammas = gen_powers(challenge, evals.len());
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
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{test_rng, UniformRand};

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
        let challenge = Fr::rand(&mut test_rng());
        let open = s.open(&coeffs, &points, challenge).expect("Open failed");
        assert_eq!(
            Ok(true),
            s.verify(&commits, &points, &evals, &open, challenge)
        );
    }
}
