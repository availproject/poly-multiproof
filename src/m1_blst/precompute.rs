//! Precomputation for Method 1 with blst optimizations
use ark_bls12_381::{Bls12_381, Fr, G2Projective as G2};
use ark_poly::univariate::DensePolynomial;
use ark_std::vec::Vec;
use merlin::Transcript;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::{fast_msm, vanishing_polynomial, Error, Proof};
use crate::lagrange::LagrangeInterpContext;
use crate::traits::{Committer, PolyMultiProof};
use crate::{cfg_iter, Commitment};

/// Method 1 with blst optimization and precomputed lagrange polynomials/vanishing polys
pub struct M1Precomp {
    /// The inner method 1 object without precomputation
    pub inner: super::M1NoPrecomp,
    point_sets: Vec<Vec<Fr>>,
    vanishing_polys: Vec<DensePolynomial<Fr>>,
    g2_zeros: Vec<G2>,
    lagrange_ctxs: Vec<LagrangeInterpContext<Fr>>,
}

impl M1Precomp {
    /// Make a precompute-optimized version of a method 1 object for the given sets of points
    pub fn from_inner(inner: super::M1NoPrecomp, point_sets: Vec<Vec<Fr>>) -> Result<Self, Error> {
        let vanishing_polys: Vec<_> = cfg_iter!(point_sets)
            .map(|(_, ps)| vanishing_polynomial(ps))
            .collect();
        let g2_zeros = cfg_iter!(vanishing_polys)
            .map(|(_, p)| fast_msm::g2_msm(&inner.prepped_g2s, &p, inner.powers_of_g2.len()))
            .collect::<Result<Vec<_>, Error>>()?;
        let lagrange_ctxs = cfg_iter!(point_sets)
            .map(|(_, ps)| LagrangeInterpContext::new_from_points(ps.as_ref()))
            .collect::<Result<Vec<_>, Error>>()?;

        Ok(M1Precomp {
            inner,
            point_sets,
            vanishing_polys,
            g2_zeros,
            lagrange_ctxs,
        })
    }
}

impl Committer<Bls12_381> for M1Precomp {
    fn commit(&self, poly: impl AsRef<[Fr]>) -> Result<Commitment<Bls12_381>, Error> {
        self.inner.commit(poly)
    }
}

impl PolyMultiProof<Bls12_381> for M1Precomp {
    type Proof = Proof;

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[Fr]>],
        polys: &[impl AsRef<[Fr]>],
        point_set_index: usize,
    ) -> Result<Self::Proof, Error> {
        self.inner.open_with_vanishing_poly(
            transcript,
            evals,
            polys,
            &self.point_sets[point_set_index],
            &self.vanishing_polys[point_set_index],
        )
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<Bls12_381>],
        point_set_index: usize,
        evals: &[impl AsRef<[Fr]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error> {
        self.inner.verify_with_lag_ctx_g2_zeros(
            transcript,
            commits,
            &self.point_sets[point_set_index],
            evals,
            proof,
            &self.lagrange_ctxs[point_set_index],
            &self.g2_zeros[point_set_index],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::M1Precomp;
    use crate::{m1_blst::M1NoPrecomp, test_rng, testing::test_basic_precomp};
    use ark_bls12_381::Fr;
    use ark_std::{vec, vec::Vec, UniformRand};

    #[test]
    fn test_basic_open_works() {
        let points = (0..30)
            .map(|_| Fr::rand(&mut test_rng()))
            .collect::<Vec<_>>();
        let s = M1NoPrecomp::new(256, 30, &mut test_rng());
        let s = M1Precomp::from_inner(s, vec![points.clone()]).expect("Failed to construct");
        test_basic_precomp(&s, &points);
    }
}
