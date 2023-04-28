use ark_std::{vec, vec::Vec};
use core::ops::Mul;

use ark_bls12_381::{Bls12_381, Fr};
use ark_ec::{pairing::Pairing, AffineRepr};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

use crate::{
    poly_div_q_r,
    traits::{Committer, KZGProof},
};

use super::{fast_msm, M1NoPrecomp};

impl KZGProof<Bls12_381> for M1NoPrecomp {
    type Proof = super::Proof;

    fn open(&self, witness_poly: Vec<Fr>) -> Result<Self::Proof, crate::Error> {
        self.commit(witness_poly)
            .map(|commit| super::Proof { 0: commit.0 })
    }

    fn verify(
        &self,
        commit: &crate::Commitment<Bls12_381>,
        point: Fr,
        value: Fr,
        proof: &Self::Proof,
    ) -> Result<bool, crate::Error> {
        let g1 = self.powers_of_g1[0];
        let g2 = self.powers_of_g2[0];
        let lhsg1 = commit.0.into_group() - g1.mul(value);
        let lhsg2 = g2;

        let rhsg1 = proof.0;
        let rhsg2 = self.powers_of_g2[1] - g2.mul(point);

        Ok(fast_msm::check_pairings_equal(
            lhsg1.into(),
            lhsg2.into(),
            rhsg1,
            rhsg2.into(),
        ))
    }

    fn compute_witness_polynomial(
        &self,
        poly: Vec<Fr>,
        point: <Bls12_381 as Pairing>::ScalarField,
    ) -> Result<Vec<Fr>, crate::Error> {
        use ark_ff::One;
        let poly = DensePolynomial::from_coefficients_vec(poly);
        let divisor = DensePolynomial::from_coefficients_vec(vec![-point, Fr::one()]);
        let (q, _r) = poly_div_q_r(poly.into(), divisor.into())?;
        Ok(q)
    }
}

#[cfg(test)]
mod tests {
    use crate::testing::test_kzg;

    #[test]
    fn test_m1_blst_kzg() {
        let srs = super::M1NoPrecomp::new(256, 32, &mut crate::test_rng());
        test_kzg(&srs);
    }
}
