use crate::utils::smallest_power_of_2_greater_than;
use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use core::ops::Mul;

fn poly<F: Field>(p: Vec<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(p)
}

/// Compute the inverse of the monic polynomial $u$ mod $x^l$.
fn inv_modl<F: FftField>(u: &DensePolynomial<F>, n: usize) -> DensePolynomial<F> {
    debug_assert!(!u[0].is_zero());

    let u0 = u.coeffs[0].clone();
    let mut v: Vec<F> = vec![F::one() / u0];
    for _i in 0..=smallest_power_of_2_greater_than(n) {
        let vpoly: DensePolynomial<F> = poly(v.clone());
        let doublev: DensePolynomial<F> = vpoly.mul(F::from(2u8));
        let uv: DensePolynomial<F> = u.mul(&vpoly);
        let uv2: DensePolynomial<F> = uv.mul(&vpoly);
        let mut vp = &doublev - &uv2;
        vp.coeffs.truncate(n + 1);
        v = vp.coeffs
    }
    poly(v)
}

pub fn truncate_poly<F: Field>(mut p: DensePolynomial<F>, max_coeffs: usize) -> DensePolynomial<F> {
    p.coeffs.truncate(max_coeffs);
    // We do this to prevent training zeros from messing computations up
    poly(p.coeffs)
}

pub(crate) fn rev_poly<F: Field>(mut p: DensePolynomial<F>) -> DensePolynomial<F> {
    p.coeffs.reverse();
    // We do this to prevent training zeros from messing computations up
    poly(p.coeffs)
}

/// Context for performing polynomial division with a fixed denominator in near-linear time
pub struct FastDivisionContext<F: Field> {
    denom_rev_inv: DensePolynomial<F>,
    max_num_poly_deg: usize,
    denom_degree: usize,
}

impl<F: FftField> FastDivisionContext<F> {
    pub fn new(denom_poly: DensePolynomial<F>, max_num_poly_deg: usize) -> Self {
        let denom_degree = denom_poly.degree();
        let l = max_num_poly_deg - denom_degree + 1;
        let rev_denom = rev_poly(denom_poly);
        Self {
            denom_rev_inv: inv_modl(&rev_denom, l),
            max_num_poly_deg,
            denom_degree,
        }
    }

    pub fn fast_div(&self, num_poly: DensePolynomial<F>) -> Result<DensePolynomial<F>, ()> {
        //TODO: Figure out what degrees are ok to use and error otherwise
        if num_poly.degree() > self.max_num_poly_deg {
            return Err(());
        }
        let l = num_poly.degree() - self.denom_degree + 1;
        let num_rev = rev_poly(num_poly);
        let qrev = num_rev.mul(&self.denom_rev_inv);
        let qrev = truncate_poly(qrev, l);
        Ok(rev_poly(qrev))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly_div_q_r;
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use rand::thread_rng;
    use rayon::prelude::*;

    fn tostr(p: &Vec<Fr>) -> String {
        let a = p
            .iter()
            .map(|a| {
                if a.is_zero() {
                    "0"
                } else if a.is_one() {
                    "1 "
                } else {
                    "N"
                }
            })
            .collect::<Vec<_>>()
            .join(" ");
        format!("len: {}, {}", p.len(), a)
    }

    #[test]
    fn test_modl() {
        let p1 = DensePolynomial::<Fr>::rand(128, &mut thread_rng());
        let num_degree = 256;
        let l = num_degree - p1.degree() + 1;
        let v = inv_modl(&p1, num_degree - p1.degree() + 1);
        let res = &p1 * &v;
        assert_eq!(res[0], One::one());
        for i in 1..=l {
            assert_eq!(res[i], Zero::zero(), "x^{}", i);
        }
        println!("{}", tostr(&res.coeffs));
    }

    #[test]
    fn test_rev_poly() {
        (1..=32).into_par_iter().for_each(|denom_size| {
            let denom = DensePolynomial::<Fr>::rand(denom_size, &mut thread_rng());
            let ctx = FastDivisionContext::new(denom.clone(), 128);
            (denom_size..=128).into_par_iter().for_each(|num_size| {
                let num = DensePolynomial::<Fr>::rand(num_size, &mut thread_rng());
                let q = ctx.fast_div(num.clone()).unwrap();

                let (naive_q, _) = poly_div_q_r(num.clone().into(), denom.clone().into()).unwrap();
                assert_eq!(
                    q.len(),
                    naive_q.len(),
                    "num: {}, denom: {}",
                    num_size,
                    denom_size
                );
                assert_eq!(
                    q.coeffs, naive_q,
                    "num: {}, denom: {}",
                    num_size, denom_size
                );
            })
        });
    }
}
