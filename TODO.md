# TODO

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
