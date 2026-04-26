# TODO

## Outstanding from 2026-04-25 session

- [ ] **poles() — pole-picking math unclear**. Need to derive on paper.

  Goal: pick 2nd pole = point with max distance AND most perpendicular
  to prev dims. Currently uses fastmap residual:

      g2(r1, r2) = distx(r1, r2)^2 - sum_over_prev_dims( (x1-x2)^2 )

  Claim (unverified): residual already mixes both — big d^2 = far,
  parallel to prev dim shrinks residual via (x1-x2)^2 subtraction,
  so perpendicular+far stays large. Need to convince self this is
  what "perpendicular" really means here.

  Geometry note (capture before forgotten):
  - gap-dist  = distx(d, r1, r2)
  - x, gap    = project(dim, d, r)   # proj returns x; gap = dim.gap
  - need dw   because y = sqrt(max(0, de^2 - x^2))
  - "max dist*angle is y*y/x" — meaning unclear, derive on paper

  Open questions:
  1. Is pairwise residual `g2 = d^2 - sum (x1-x2)^2` right scoring fn,
     or want per-point y-based score?
  2. If y-based: each r has y_k = sqrt(de_k^2 - x_k^2) per prev dim k.
     Combine how? min, prod, avg?
  3. Should `proj` return `(x, gap)` so callers don't reach dim.gap?
  4. Algorithm:
     - O(n^2) all-pairs (current).
     - O(n) linear scan: anchor=t[0], find e=farthest, then w=farthest
       from e. Standard fastmap heuristic.
     - Cached projections: each r projects on each prev dim once
       upfront, pairwise lookup O(prev). Faster.

  Action: sketch math on paper, decide scoring fn, then pick algorithm
  (linear scan + cached proj likely winner once scoring locked).

## Outstanding from 2026-04-23 session

- [ ] **Auto-scale eps by dimensionality**. One-liner in `dims()`:
      `epsy = parentVar * the.eps / max(1, len(d.cols.xs))`
      Removes per-dataset tuning. Tested rule: `eps ≈ 1/d` works.
      Verify across auto93 (4 xs), SS-M (17), simplex (5).

- [ ] **Why only 2 dims on auto93?** Even at eps=0.01 we stop at 2.
      Instrument `newDim` to print whether gap<1e-9 or cuts empty.
      Suspect: after 2 dims residual variance collapses on small budget.

- [ ] **Zero-variance y pathology**. Simplex with Goal-=0 shatters
      (all ys=0 → parentVar=0 → gain=0 ≥ 0 fires every k).
      Guard in `binsOf`: `if parentVar < 1e-12: return []`.

- [ ] **Simplex test data**. Current generator gives uniform Goal=0.
      Add variant where Goal = vertex index (0..5). Then `test__simplex`
      becomes a real separation test — expect ≈6 clusters.

- [ ] **Pyright noise**. `aha` / `_numAdd` Atom-union warnings.
      Narrow via `assert isinstance(v, (int, float))` or inline `type: ignore`.
      Harmless but clutters diagnostics.

- [ ] **Bins-per-dim constant 7 on simplex** is suspicious. May indicate
      Welford drift when `sub` runs through zero-variance stretches.
      Add sanity: `assert R.sd >= 0` after each sub.

## Done in session

- Deleted `nump`, `nest`, `_ingest`, `_symAdd`, `_numAdd`, `symp`.
- Added `isa = isinstance`, type aliases `Atom/Row/Rows`.
- Renamed `summarize` → `add`, made polymorphic Num/Sym/Data.
- Added `sub(i,v)` via weighted Welford.
- Rewrote `binsOf` single-pass with incremental stats (O(n log n)).
- Switched cut test from range to variance-gain vs global parent.
- Added `test__simplex`, `make eps`, `make eps_simplex`.
