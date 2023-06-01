# `poly-multiproof`
`poly-multiproof` is a library for generating BDFG21-esque proofs against standard KZG
commitments. It supports both the method 1 and method 2 proofs from the paper, and contains
an assembly-optimized implementation of method 1 for the BLS12-381 curve. When the points that
will be committed to are known beforehand, separate `precompute` modules can be used which
pre-compute lagrange polynomials and vanishing polynomials for the points, which can speed up
proof generation by a significant amount, especially for larger proof sizes.

### Features
* `blst` enables a specific `bls12-381` implementation which uses `blst` for curve msm.
* `parallel` enables parallel computation for
  * PMP setup generation
  * operations in the `data_availability_grid` example
* `print-trace` enables some tracing that shows the time certain things take to execute

See [the `poly-multiproof` documentation](https://docs.rs/poly-multiproof) for more details.

### Examples

An example of using pmp for a grid data availability scheme with 1d erasure encoding is in `examples/data_availability_grid.rs`. To run it with a nice timer, do
```bash
cargo run --example data_availability_grid --release --features print-trace,blst,parallel
```

which will print out something like this example for a 256x256 grid
```
Start:   create pmp
End:     create pmp ........................................330.880ms
Start:   create grid
··Start:   erasure encoding columns
··End:     erasure encoding columns ........................32.143ms
··Start:   computing polynomials from evals
··End:     computing polynomials from evals ................22.562ms
··Start:   computing commitments
··End:     computing commitments ...........................707.127ms
End:     create grid .......................................771.776ms
Start:   opening to grid
End:     opening to grid ...................................244.446ms
Start:   verifying grid
End:     verifying grid ....................................223.557ms
```

There are nice constants in the top of the file to play with.

### Benchmarks

To run benchmarks with `arkworks-rs` asm optimizations on x86 machines, run
```bash
RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly criterion --features asm
```
or to run with the goal of plotting, run
```bash
RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly criterion --features asm --plotting-backend disabled -- --quick --quiet &> bench_out.txt
```
The logs in `bench_out.txt` can then be parsed and plotted in `Plot Benches.ipynb`.
Using `--quick` is nice since there are many many inputs benchmarked and it will still take an hour or so to run with `--quick`.

