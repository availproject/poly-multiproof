use ark_bls12_381::{Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::One;
use ark_ff::UniformRand;
use divan::Bencher;
use rayon::prelude::*;
#[divan::bench_group]
mod pairing_benches {
    use super::*;

    #[derive(Clone)]
    struct Inputs {
        g1: G1Projective,
        g2: G2Projective,
        g1aff: G1Affine,
        g2aff: G2Affine,
        h1aff: G1Affine,
        h2aff: G2Affine,
    }

    lazy_static::lazy_static! {

        static ref INPUTS: Inputs = {
            let mut rng = rand::thread_rng();
            let g1 = ark_bls12_381::G1Projective::rand(&mut rng);
            let g2 = ark_bls12_381::G2Projective::rand(&mut rng);
            let z = Fr::rand(&mut rng);
            let h1 = ark_bls12_381::G1Projective::rand(&mut rng) * (Fr::one() / z);
            let h2 = ark_bls12_381::G2Projective::rand(&mut rng) * z;
            let g1aff = g1.into_affine();
            let g2aff = g2.into_affine();
            let h1aff = h1.into_affine();
            let h2aff = h2.into_affine();
            Inputs { g1, g2, g1aff, g2aff, h1aff, h2aff }
        };
    }

    #[divan::bench]
    fn ark_pairing_bench(bencher: Bencher) {
        bencher
            .with_inputs(|| INPUTS.clone())
            .bench_refs(|i| ark_bls12_381::Bls12_381::pairing(i.g1, i.g2));
    }

    #[divan::bench]
    fn blst_pairing_bench(bencher: Bencher) {
        bencher
            .with_inputs(|| INPUTS.clone())
            .bench_refs(|i| poly_multiproof::m1_blst::fast_msm::pairing(i.g1aff, i.g2aff))
    }

    #[divan::bench]
    fn blst_slow_check(bencher: Bencher) {
        bencher.with_inputs(|| INPUTS.clone()).bench_refs(|i| {
            let a1 = poly_multiproof::m1_blst::fast_msm::pairing(i.g1aff, i.g2aff);
            let a2 = poly_multiproof::m1_blst::fast_msm::pairing(i.h1aff, i.h2aff);
            a1 == a2
        });
    }
    #[divan::bench]
    fn blst_fast_check(bencher: Bencher) {
        bencher.with_inputs(|| INPUTS.clone()).bench_refs(|i| {
            poly_multiproof::m1_blst::fast_msm::check_pairings_equal(
                i.g1aff, i.g2aff, i.h1aff, i.h2aff,
            )
        });
    }
}

#[divan::bench_group(sample_count = 5, sample_size = 5)]
mod polydiv_benches {

    use ark_poly::{
        univariate::{DenseOrSparsePolynomial, DensePolynomial},
        DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
    };
    use core::fmt;
    use poly_multiproof::poly_ops::{truncate_poly, FastDivisionContext};
    use rand::thread_rng;
    use std::{ops::Deref, time::Instant};

    use super::*;

    struct Input {
        num: usize,
        denom: usize,
        div_ctx: FastDivisionContext<Fr>,
    }

    impl fmt::Display for Input {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "num: {}, denom: {}", self.num, self.denom,)
        }
    }

    impl Input {
        fn num_poly(&self) -> DensePolynomial<Fr> {
            truncate_poly(NUM.clone(), self.num)
        }
        fn denom_poly(&self) -> DensePolynomial<Fr> {
            truncate_poly(DENOM.clone(), self.denom)
        }
    }
    const MAX_COEFFS: usize = 1024;
    const STEP: usize = 256;

    fn inputs() -> Vec<Input> {
        println!("poly setup");
        let now = Instant::now();
        let denoms: Vec<_> = (STEP..=MAX_COEFFS).step_by(STEP).collect();

        let ret: Vec<_> = denoms
            .into_par_iter()
            .flat_map(|denom| {
                //let nums: Vec<_> = (denom..=MAX_COEFFS).step_by(STEP).collect();
                let nums = vec![MAX_COEFFS];
                nums.into_par_iter()
                    .map(|num| {
                        let denom_poly = truncate_poly(DENOM.clone(), denom);
                        let div_ctx = FastDivisionContext::new(denom_poly, num);
                        println!("num: {}, denom: {}", num, denom);
                        Input {
                            num,
                            denom,
                            div_ctx,
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        println!("poly setup done: {:#?}", now.elapsed());
        ret
    }

    lazy_static::lazy_static! {
        static ref NUM: DensePolynomial<Fr> = {
            DensePolynomial::rand(MAX_COEFFS, &mut thread_rng())
        };
        static ref DENOM: DensePolynomial<Fr> = {
            DensePolynomial::rand(MAX_COEFFS, &mut thread_rng())
        };

        static ref INPUTS: Vec<Input> = {
            inputs()
        };
    }

    #[divan::bench(args = INPUTS.deref())]
    fn naive_div(inp: &Input) {
        let num: DenseOrSparsePolynomial<_> = inp.num_poly().into();
        let denom: DenseOrSparsePolynomial<_> = inp.denom_poly().into();
        num.divide_with_q_and_r(&denom).unwrap();
    }

    #[divan::bench(args = INPUTS.deref())]
    fn fast_rev_div(inp: &Input) {
        inp.div_ctx.fast_div(inp.num_poly()).unwrap();
    }

    #[divan::bench]
    fn div_by_ev_vanishing(bencher: Bencher) {
        bencher
            .with_inputs(|| {
                let ev = Radix2EvaluationDomain::<Fr>::new(256).unwrap();
                let poly = DensePolynomial::<Fr>::rand(1024, &mut thread_rng());
                (ev, poly)
            })
            .bench_values(|(ev, poly)| poly.divide_by_vanishing_poly(ev));
    }
}

fn main() {
    divan::main()
}
