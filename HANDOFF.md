# HANDOFF: Dim stopping criterion for rash.py

## Goal

Make `dims()` in `rash.py` stop discovering new dimensions when they no longer add useful structure. Specifically: 7D data with only 4 live dims should find ~4 dims, not 6-7.

## What's done

- **Test data**: `etc/spread4of7.py` generates 7D data (4 live dims σ=0.08, 3 dead dims σ=0.001, 5 clusters along a gradient, goal = sum of live center values + noise). Output: `etc/spread4of7.csv` (200 rows).
- **Test**: `test__4d` in `rash.py` loads this data, runs `dims()`, asserts ≤5 dims found.
- **Equal-width bins**: Replaced greedy `binsOf` (variance-reduction cuts) with simple equal-width bins. `--bins=6` param added. `newDim` now divides projection range into `the.bins` equal parts. Removed `binsOf` function entirely.
- **epsx gate**: `newDim` rejects dims where `gap < epsx` (where `epsx = the.eps * sd(sampled pairwise x-distances)`). Works on synthetic data but doesn't fully block dead dims on real data because `distx` mixes all columns and residual gaps stay above threshold.
- **Removed `--min` param** (was used by old `binsOf`). Needs to be re-added for the new stopping criterion.

## What worked

- Equal-width bins are simpler and eliminate the epsy/binsOf threshold problem that was causing 0 dims on real config data.
- epsx gate: good idea in principle but insufficient alone — `distx` uses all dims so residual gaps don't collapse enough on dead dims.
- auto93.csv: `--dim` gives 6 dims/20 clusters with current code (was 2 dims with old greedy binsOf).
- Config files now get dims (were getting 0 with old variance-threshold binsOf).

## What didn't work

- **epsy in binsOf**: `parentVar * the.eps` as variance-reduction threshold was too strict with 30 rows. Most config files got 0 dims. Pre-existing problem, not caused by our changes.
- **Pole y-separation gate**: Poles chosen for max x-distance naturally land in different y-regions regardless of dim quality. Too weak a filter.
- **Residual epsx** (recomputing xSd after each dim from residual distances): Sampling variance too noisy, didn't reliably collapse for dead dims.
- **Budget increase to 100**: Broke original code on auto93 (2 dims → 1 dim). Reverted to 30.

## What's left — THE NEXT STEP

Implement **lives-based stopping**:

1. Add params: `--lives=2`, `--min=2`
2. In `dims()`, after each new dim:
   - Compute `clusters(d, rs, out)` using all dims found so far
   - Count clusters with `len(rows) >= the.min`
   - If count didn't increase vs previous iteration: `lives -= 1`
   - If `lives == 0`: stop, remove last dim(s) that didn't help
3. This naturally handles dead dims: equal-width bins on a dead dim put all rows in 1-2 bins → cluster count doesn't grow → burns a life

The user was about to add `--lives=2` and `--min=2` to the options docstring when they stopped. The `dims()` function needs the lives loop logic. `clusters()` already exists and returns `dict[tuple[int,...], Rows]`.

## Key context

- User uses `i` for `self` in classes (project convention).
- `the` is the global config namespace parsed from docstring defaults.
- Tests run via `python3 rash.py --testname`.
- Real-world test: `for f in ~/gits/moot/optimize/config/*.csv; do ./rash.py --file=$f --dim; done`
- User strongly prefers keeping `distx` as-is (uses all dims equally — supported by literature).
- Caveman mode active (terse responses, no filler).
