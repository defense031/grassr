# calibration-contrib — community calibration runs

This branch is the landing zone for contributed calibration compute.
It never merges into `main` wholesale; verified blocks enter the
package's bundled reference through the `data-raw/` rebuild scripts at
each release.

## What lands here

grassr's reference surfaces and null distributions are Monte Carlo
objects, and their remaining bounds are compute limits. The installed
package (0.7.2 and later) includes the contributor pipeline:

```r
grass_calibration_manifest()                  # the open blocks
grass_contribute(dir, hours = 5, dry_run = TRUE)   # what 5 hours buys here
grass_contribute(dir, hours = 5)              # run it, write a bundle
grass_verify_contribution(dir)                # check before submitting
```

A bundle is one directory: one result file per block (the per-draw
statistic vectors, a few hundred KB each), `bundle_manifest.csv`
(checksums, seeds, versions), and `SESSIONINFO.txt`. Raw rating panels
are never stored or submitted.

## How to submit

1. Fork the repository and check out this branch.
2. Copy your bundle directory to
   `contrib/contributions/<your-github-name>-<YYYYMMDD>/`.
3. Open a pull request against `calibration-contrib` (never `main`).
   For large claims, open an issue first
   (`calibration-run: blocks <first>-<last>`) so runs do not collide.
   Duplicate runs are not wasted; they cross-verify.

## Intake (maintainers)

1. `grass_verify_contribution(dir)` on the submitted bundle: checksums,
   completeness, and a leading-draw seed replay. The replay re-runs only
   the first few hundred draws from the block's seed at the recorded
   package version. That is enough to prove the bundle came from the real
   pipeline: a contributor cannot reproduce those draws without having run
   it. It costs seconds, never the block's full runtime.
2. Sanity-check the submitted draws against expectation: finite rate,
   `delta_hat` quantiles, implied qualities centered on the cell's `q`, a
   balanced prevalence sweep, and consistency with neighboring cells. A
   block that fails is returned with the diff, not merged. Do not
   re-execute whole blocks at intake. Re-running a block reproduces the
   contributor's compute and defeats the point of distributed calibration.
   A full re-run is only ever a rare, random spot-audit of one small
   block; routine trust comes from the seed replay above plus the
   cross-verification of duplicate runs over time.
3. Merged blocks are marked `claimed` in `contrib/calibration_manifest.csv`
   on this branch (the authoritative live manifest; the package bundles
   a snapshot). Edit the one status cell by hand so the diff is a single
   line; a `write.csv` round-trip reformats the whole file.
4. At the next release, verified draws enter the bundled reference via
   `data-raw/`, and the contribution is credited in NEWS.md.

## Manifest

`contrib/calibration_manifest.csv` on this branch is the live copy of
the manifest. Programs:

- `tail_topup` (21 blocks) — deeper draws for the shipped null cells
  whose extreme tails carry a transparency flag.
- `prev_strata` (1,925 blocks) — a prevalence-stratified delta_hat
  null; the shipped null pools five prevalences and its realized flag
  size runs to roughly twice nominal at extreme prevalence.
- `lattice_k` (220 blocks) — intermediate rater counts between the
  shipped grid points, cutting the design-snap distance.

Every block's seed is `20300000 + block_id`, disjoint from every
shipped program seed, so contributed draws can never replay shipped
draws. Block ids are frozen; the manifest only ever changes `status`
or appends new blocks.
