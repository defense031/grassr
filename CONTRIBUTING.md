# Contributing to grassr

Two kinds of contribution help this project: code and compute.

## Code, bugs, and documentation

File bug reports at https://github.com/defense031/grassr/issues with a
minimal reproducible example (the rating matrix or the
`gen_logitnormal()` call that produces it, plus `sessionInfo()`).
Pull requests are welcome for fixes and documentation. Every PR runs
the full test suite on the five-platform GitHub Actions matrix; new
behavior needs a test.

## Compute: extending the calibration

Every number grassr reports is read off a simulated reference, and
the reference is only as fine as the compute that built it. The
bundled calibration spent 133 million simulated panels. Its remaining
bounds are compute limits, not method limits, and the pipeline
that removes them is in this repository under `simulation/`. Every
stage runs from a declared seed and partitions by cell, so
independent machines can extend the calibration. The shipped surfaces
were themselves computed on two machines and merged by cell.

### Open programs

| Program | What it buys | Rough scale |
|---|---|---|
| Subject-prevalence profile compendium | The ICC reference covers 52 profiles (a 48-point logit-normal grid plus four mixtures); a panel whose true profile sits elsewhere carries a misspecification cost on the order of 20 percentile points at small designs. Skew-normal, heavy-tailed, and bimodal families would convert ICC's disclosed marker into a calibrated reading. | Largest open program; scales with profiles x 44,616 cells |
| Prevalence-stratified `delta_hat` null | The shipped null pools five calibration prevalences, and the realized flag size runs to roughly twice nominal at extreme prevalence. Stratifying the null removes the distortion. | ~5x the shipped null program (~110M draws) |
| Denser (k, N) lattice | Designs off the calibrated grid snap to the nearest cell and the card discloses the snap. Intermediate k and N values shrink the snap distance. | Scales with added cells x 2,000 replications |
| Null tail top-ups | 21 of the 385 shipped null cells carry unstable-extreme-tail transparency flags at the current draw counts. Deeper draws resolve them. | Small; 21 cells |
| Power map extension | The shipped power table covers asymmetry A in {0.05, 0.10, 0.15, 0.20, 0.30}. Finer alternatives and off-grid designs sharpen prospective design guidance. | Moderate |

### Protocol: from the installed package

The pipeline is part of the installed package. The null programs
above are divided into 2,166 seeded blocks, and three functions cover
the whole path:

1. **Size it.** `grass_contribute(dir = tempdir(), hours = 5,
   dry_run = TRUE)` benchmarks your machine on two small anchor cells
   and returns the blocks a five-hour budget buys, with estimated wall
   time. Nothing is written.
2. **Claim (optional but kind).** Block selection is randomized so
   uncoordinated runs rarely collide, and duplicate runs cross-verify
   rather than waste. For a large budget, open an issue titled
   `calibration-run: blocks <first>-<last>` and pass the claimed range
   as `blocks = first:last`.
3. **Run it.** The same call without `dry_run` runs each block from
   its manifest seed (`20300000 + block_id`, disjoint from every
   shipped program seed) and writes the bundle to `dir`: one result
   file per block plus `bundle_manifest.csv` (checksums, versions) and
   `SESSIONINFO.txt`. A block's result file holds the per-draw
   statistic vectors (`delta_hat` and the three implied qualities),
   a few hundred KB per block. Raw rating panels are never stored or
   submitted (`**/per_cell/` simulation output stays gitignored).
4. **Check it.** `grass_verify_contribution(dir)` validates checksums,
   completeness, and a seed replay of each block's leading draws.
5. **Submit.** Fork the repository, add your bundle directory under
   `contrib/contributions/<your-name>-<date>/` on the
   `calibration-contrib` branch, and open a pull request against that
   branch (never `main`). Maintainers verify by re-executing your
   seeds, which reproduce bit for bit at the same package version;
   runs that do not reproduce are returned with the diff, not merged.
6. **Shipping.** Verified blocks enter the next release's bundled
   reference through the `data-raw/` rebuild scripts, and the run is
   credited in the package NEWS.

The wider programs (the reference surfaces and the ICC
subject-prevalence compendium) are not yet pre-chunked into package
blocks; they run from the stage scripts under
`simulation/v070_program/` (the scripts take `SPLIT_MOD`/`SPLIT_REM`
environment variables to partition cells across machines, and declare
their seed constants at the top). Open an issue first and we will cut
a block range.

Long runs on macOS: the pipeline's chain scripts pin `caffeinate` so
the machine does not sleep mid-stage. Budget accordingly; the shipped
null regeneration was 44 million draws over roughly 19 hours across
two machines.
