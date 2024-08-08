use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use ark_std::log2;
use merlin::Transcript;
use poly_multiproof::m1_blst::cyclic_precomp::M1CyclPrecomp;
use poly_multiproof::traits::{Committer, PolyMultiProof};
use poly_multiproof::{method1, method2, Commitment};
use rand::thread_rng;
use rayon::prelude::*;
use std::{
    fmt::{self, Display, Formatter},
    ops::Deref,
    time::Instant,
};

type M1 = method1::M1NoPrecomp<Bls12_381>;
type M1Pc = method1::precompute::M1Precomp<Bls12_381>;
type M2 = method2::M2NoPrecomp<Bls12_381>;
type M2Pc = method2::precompute::M2Precomp<Bls12_381>;
type M1Blst = poly_multiproof::m1_blst::M1NoPrecomp;
type M1BlstPc = poly_multiproof::m1_blst::precompute::M1Precomp;
type M1BlstCyclPc = poly_multiproof::m1_blst::cyclic_precomp::M1CyclPrecomp;

const WIDTH: usize = 32768; // 256 mb block
const HEIGHT: usize = 64;
// This won't let the block width go above 128
const MIN_WIDTH: usize = 256;
const MAX_WIDTH: usize = WIDTH;
#[derive(Clone)]
struct TestGrid<F: Clone> {
    coeffs: Vec<Vec<F>>,
    evals: Vec<Vec<F>>,
    points: Vec<F>,
}

impl<F: PrimeField> TestGrid<F> {
    pub(crate) fn gen_grid(width: usize, height: usize) -> Self {
        let degree = width - 1;
        let ev = Radix2EvaluationDomain::<F>::new(width).unwrap();
        let points = ev.elements().into_iter().collect::<Vec<_>>();
        assert_eq!(points.len(), width);
        let coeffs = (0..height)
            .map(|_| DensePolynomial::<F>::rand(degree, &mut thread_rng()).coeffs)
            .collect::<Vec<_>>();
        let evals: Vec<Vec<_>> = coeffs.iter().map(|p| ev.fft(&p)).collect();
        Self {
            coeffs,
            points,
            evals,
        }
    }
    pub(crate) fn trim_pts(&self, n_pts: usize) -> Self {
        Self {
            coeffs: self.coeffs.clone(),
            evals: self.evals.iter().map(|ev| ev[..n_pts].to_vec()).collect(),
            points: self.points[..n_pts].to_vec(),
        }
    }
}

lazy_static::lazy_static! {
    static ref TEST_GRID: TestGrid<Fr> = {
        TestGrid::<Fr>::gen_grid(WIDTH, HEIGHT)
    };

    static ref M1_PMP: M1 = {
        M1::new(WIDTH, WIDTH, &mut thread_rng())
    };

    static ref M2_PMP: M2 = {
        M2::new_from_affine(M1_PMP.powers_of_g1.clone(), M1BLST_PMP.powers_of_g2[0].into(), M1BLST_PMP.powers_of_g2[1].into())
    };

    static ref M1BLST_PMP: M1Blst = {
        M1Blst::new_from_affine(&M1_PMP.powers_of_g1, &M1_PMP.powers_of_g2)
    };

    static ref COMMITS: Vec<Commitment<Bls12_381>> = {
        TEST_GRID
            .coeffs
            .par_iter()
            .map(|c| M1_PMP.commit(c).unwrap())
            .collect::<Vec<_>>()
    };
}

trait MethodStr: PolyMultiProof<Bls12_381> {
    fn str() -> String;
}

impl MethodStr for M1Pc {
    fn str() -> String {
        "m1".to_string()
    }
}

impl MethodStr for M1BlstPc {
    fn str() -> String {
        "m1blst".to_string()
    }
}

impl MethodStr for M1BlstCyclPc {
    fn str() -> String {
        "m1blst_cycl".to_string()
    }
}

impl MethodStr for M2Pc {
    fn str() -> String {
        "m2".to_string()
    }
}

type EvalSelector<P> = fn(&P, &TestGrid<Fr>, usize, usize) -> Vec<Vec<Fr>>;

struct PmpCase<P: PolyMultiProof<Bls12_381>> {
    width: usize,
    height: usize,
    opening: P::Proof,
    backend: P,
    grid: &'static TestGrid<Fr>,
    eval_selector: EvalSelector<P>,
}

fn open_with_pmp<P: PolyMultiProof<Bls12_381>>(
    pmp: &P,
    grid: &TestGrid<Fr>,
    width: usize,
    height: usize,
    eval_selector: EvalSelector<P>,
) -> P::Proof {
    pmp.open(
        &mut Transcript::new(b"bench"),
        &eval_selector(&pmp, grid, width, height),
        &grid.coeffs[..height],
        0,
    )
    .unwrap()
}

fn verify_with_pmp<P: PolyMultiProof<Bls12_381>>(
    pmp: &P,
    grid: &TestGrid<Fr>,
    commits: &[Commitment<Bls12_381>],
    open: &P::Proof,
    width: usize,
    height: usize,
    eval_selector: EvalSelector<P>,
) {
    let mut transcript = Transcript::new(b"bench");
    assert_eq!(
        pmp.verify(
            &mut transcript,
            &commits[..height],
            0,
            &eval_selector(&pmp, grid, width, height),
            &open
        ),
        Ok(true)
    );
}

impl<P: MethodStr> Display for PmpCase<P> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "method: {}, w: {}, h: {}",
            P::str(),
            self.width,
            self.height
        )
    }
}

impl<P: MethodStr + Send + Sync> Arg for PmpCase<P>
where
    P::Proof: Send + Sync,
{
    fn open(&self) {
        open_with_pmp(
            &self.backend,
            &self.grid,
            self.width,
            self.height,
            self.eval_selector,
        );
    }

    fn verify(&self) {
        verify_with_pmp(
            &self.backend,
            &self.grid,
            &COMMITS,
            &self.opening,
            self.width,
            self.height,
            self.eval_selector,
        );
    }
}

trait Arg: ToString + Send + Sync + Display {
    fn open(&self);
    fn verify(&self);
}

fn powers_of_two_up_to(from: usize, to: usize) -> Vec<usize> {
    let start = log2(from);
    (start..)
        .map(|x| 2_usize.pow(x))
        .take_while(|&x| x <= to)
        .collect()
}

fn input_args() -> Vec<Box<dyn Arg>> {
    println!("Performing Setup");
    let start = Instant::now();
    let widths = powers_of_two_up_to(MIN_WIDTH, MAX_WIDTH);
    let heights = powers_of_two_up_to(1, HEIGHT);

    let args: Vec<Box<dyn Arg>> = widths
        .par_iter()
        .flat_map(|&width| {
            let start = Instant::now();
            //let subgrid = TEST_GRID.trim_pts(width);
            //let point_sets = vec![TEST_GRID.points[..width].to_vec()];
            //let m1blst_pc = M1BlstPc::from_inner(M1BLST_PMP.clone(), point_sets.clone()).unwrap();
            let _basic_eval_selector = |grid: &TestGrid<Fr>, width: usize, height: usize| {
                grid.evals[..height]
                    .iter()
                    .map(|ev| ev[..width].to_vec())
                    .collect::<Vec<_>>()
            };

            let m1blst_cycl_pc =
                M1CyclPrecomp::from_inner(M1BLST_PMP.clone(), WIDTH, WIDTH / width).unwrap();
            let cycl_eval_selector =
                |pmp: &M1CyclPrecomp, grid: &TestGrid<Fr>, _width: usize, height: usize| {
                    grid.evals[..height]
                        .iter()
                        .map(|ev| pmp.point_sets().take_subgroup_indices(0, ev).unwrap())
                        .collect()
                };
            println!("Width {} pc took {:#?}", width, start.elapsed());
            heights
                .par_iter()
                .flat_map(|&height| -> Vec<Box<dyn Arg>> {
                    vec![
                        //Box::new(PmpCase {
                        //    width,
                        //    height,
                        //    opening: open_with_pmp(&m1blst_pc, &subgrid, height),
                        //    backend: m1blst_pc.clone(),
                        //    grid: subgrid.clone(),
                        //}),
                        Box::new(PmpCase {
                            width,
                            height,
                            opening: open_with_pmp(
                                &m1blst_cycl_pc,
                                &TEST_GRID,
                                width,
                                height,
                                cycl_eval_selector,
                            ),
                            backend: m1blst_cycl_pc.clone(),
                            grid: &TEST_GRID,
                            eval_selector: cycl_eval_selector,
                        }),
                    ]
                })
                .collect::<Vec<Box<dyn Arg>>>()
        })
        .collect();
    let _ = COMMITS.deref(); // Load the commits
    println!("Setup took {:#?}", start.elapsed());
    args
}

lazy_static::lazy_static! {
    static ref INPUT_ARGS: Vec<Box<dyn Arg>> = {
        input_args()
    };
}

#[divan::bench_group(sample_size = 3, sample_count = 3)]
mod pmp_benches {

    use super::*;
    #[divan::bench(args = INPUT_ARGS.deref())]
    fn open_bench(arg: &Box<dyn Arg>) {
        arg.open()
    }

    #[divan::bench(args = INPUT_ARGS.deref())]
    fn verify_bench(arg: &Box<dyn Arg>) {
        arg.verify()
    }
}

fn main() {
    divan::main()
}
