use ark_bls12_381::{Fr, G1Affine, G2Affine};
use ark_ff::One;
use ark_ff::UniformRand;
use divan::Bencher;
use rayon::prelude::*;
#[divan::bench_group]
mod pairing_benches {
    use ark_bls12_381::Bls12_381;
    use poly_multiproof::{
        msm::{blst::BlstMSMEngine, ArkMSMEngine},
        traits::MSMEngine,
    };

    use super::*;

    #[derive(Clone)]
    struct Inputs {
        g1aff: G1Affine,
        g2aff: G2Affine,
        h1aff: G1Affine,
        h2aff: G2Affine,
    }

    lazy_static::lazy_static! {

        static ref INPUTS: Inputs = {
            let mut rng = rand::thread_rng();
            let g1aff = ark_bls12_381::G1Affine::rand(&mut rng);
            let g2aff = ark_bls12_381::G2Affine::rand(&mut rng);
            let z = Fr::rand(&mut rng);
            let h1aff = ark_bls12_381::G1Affine::rand(&mut rng) * (Fr::one() / z);
            let h2aff = ark_bls12_381::G2Affine::rand(&mut rng) * z;
            Inputs { g1aff, g2aff, h1aff: h1aff.into(), h2aff: h2aff.into() }
        };
    }

    #[divan::bench]
    fn ark_pairing_bench(bencher: Bencher) {
        bencher
            .with_inputs(|| INPUTS.clone())
            .bench_refs(|i| <ArkMSMEngine<Bls12_381>>::pairing(i.g1aff, i.g2aff));
    }

    #[divan::bench]
    fn blst_pairing_bench(bencher: Bencher) {
        bencher
            .with_inputs(|| INPUTS.clone())
            .bench_refs(|i| BlstMSMEngine::pairing(i.g1aff, i.g2aff))
    }

    #[divan::bench]
    fn blst_slow_check(bencher: Bencher) {
        bencher.with_inputs(|| INPUTS.clone()).bench_refs(|i| {
            let a1 = BlstMSMEngine::pairing(i.g1aff, i.g2aff);
            let a2 = BlstMSMEngine::pairing(i.h1aff, i.h2aff);
            a1 == a2
        });
    }
    #[divan::bench]
    fn blst_fast_check(bencher: Bencher) {
        bencher
            .with_inputs(|| INPUTS.clone())
            .bench_refs(|i| BlstMSMEngine::pairing_eq_check(i.g1aff, i.g2aff, i.h1aff, i.h2aff));
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
                let poly = DensePolynomial::<Fr>::rand(32768, &mut thread_rng());
                (ev, poly)
            })
            .bench_values(|(ev, poly)| poly.divide_by_vanishing_poly(ev));
    }
}

#[divan::bench_group(max_time = 0.3)]
mod msm {
    use ark_bls12_381::{Fr, G1Affine};
    use ark_ff::UniformRand;
    use divan::Bencher;
    use poly_multiproof::{msm::blst::BlstMSMEngine, traits::MSMEngine};
    use rand::thread_rng;

    fn inputs(to: usize) -> Vec<usize> {
        (1..)
            .map(|x| 2_usize.pow(x))
            .take_while(|&x| x <= to)
            .collect()
    }
    #[divan::bench(args = inputs(2usize.pow(15)))]
    fn bench_msm(bencher: Bencher, size: usize) {
        bencher
            .with_inputs(|| {
                (
                    BlstMSMEngine::prepare_g1(
                        (0..size)
                            .map(|_| G1Affine::rand(&mut thread_rng()))
                            .collect::<Vec<_>>(),
                    ),
                    (0..size)
                        .map(|_| Fr::rand(&mut thread_rng()))
                        .collect::<Vec<_>>(),
                )
            })
            .bench_values(|(g1s, frs)| BlstMSMEngine::multi_scalar_mul_g1(&g1s, &frs));
    }
}

fn main() {
    divan::main()
}
