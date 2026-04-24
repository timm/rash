#!/usr/bin/env python3 -B
# rash.py: recursive axis-split clustering
# (c) 2026 Tim Menzies, timm@ieee.org, MIT license
"""
Options:

    --seed=1          random seed
    --p=2             Minkowski p
    --eps=0.35        margin multiplier
    --budget=30       label budget
    --bins=6          bins per dim
    --few=32          pole sample size
    --dims=6          max dims
    --file=auto93.csv   input csv
"""
from __future__ import annotations
import re, sys, random
from random import sample, shuffle
from math import exp
from types import SimpleNamespace as S
from typing import Iterable, Iterator

Atom = int | float | str
Row  = list[Atom]
Rows = list[Row]
isa  = isinstance

# ---- 0. Utilities ----
def o(x) -> str:
  """Format recursively with 2dp floats."""
  if isa(x, float): return f"{x:.2f}"
  if isa(x, dict):
    return "{" + ", ".join(f"{k}={o(v)}"
           for k, v in sorted(x.items())) + "}"
  if isa(x, list):  return "{" + ", ".join(map(o, x)) + "}"
  if isa(x, S):     return "S" + o(x.__dict__)
  return str(x)

def thing(txt: str) -> Atom:
  """Coerce str to int/float/bool/str."""
  txt = txt.strip()
  b = lambda s: {"true": 1, "false": 0}.get(s.lower(), s)
  for f in [int, float, b]:
    try: return f(txt)
    except ValueError: pass
  return txt

def csv(f: str) -> Iterator[Row]:
  """Yield typed rows from file."""
  with open(f, encoding="utf-8") as file:
    for s in file:
      r = s.partition("#")[0].split(",")
      if any(x.strip() for x in r):
        yield [thing(x) for x in r]

# ---- 1. Columns ----
class Num:
  """Summarize stream of numbers."""
  def __init__(i, txt: str="", a: int=0):
    i.txt, i.at, i.n   = txt, a, 0
    i.mu, i.m2, i.sd   = 0.0, 0.0, 0.0
    i.heaven           = txt[-1:] != "-"

class Sym:
  """Summarize stream of symbols."""
  def __init__(i, txt: str="", a: int=0):
    i.txt, i.at, i.n, i.has = txt, a, 0, {}

Col_t = Num | Sym

def Col(txt: str="", a: int=0) -> Col_t:
  """Num if txt starts uppercase, else Sym."""
  return (Num if txt[:1].isupper() else Sym)(txt, a)

def add(i, v, w: int = 1):
  """Incremental update for Num/Sym/Data. w=-1 removes."""
  if isa(i, Data):
    for c in i.cols.all: add(c, v[c.at], w)
    (i.rows.append if w > 0 else i.rows.remove)(v)
    return v
  if v == "?": return v
  if isa(i, Sym): i.has[v] = w + i.has.get(v, 0)
  elif w < 0 and i.n <= 2: i.n = i.mu = i.m2 = i.sd = 0
  else:
    i.n += w
    d     = v - i.mu
    i.mu += w * d / i.n
    i.m2 += w * d * (v - i.mu)
    i.sd  = (max(0, i.m2) / (i.n - 1))**.5 if i.n > 1 else 0
  return v

def sub(i, v): return add(i, v, -1)

def norm(i: Num, v: Atom) -> float | str:
  """Logistic normalize a Num value."""
  if v == "?": return v
  z = (v - i.mu) / (i.sd + 1e-32)
  z = max(-3, min(3, z))
  return 1 / (1 + exp(-1.7 * z))

def adds(vs: Iterable[Atom], n: Num | None=None) -> Num:
  """Summarize numbers into n (fresh Num by default)."""
  n = n or Num()
  for v in vs: add(n, v)
  return n

# ---- 2. Data ----
class Cols:
  """Organize columns from header names."""
  def __init__(i, names: list[str]):
    i.names = names
    i.all   = [Col(t, j) for j, t in enumerate(names)]
    i.xs    = [c for c in i.all if c.txt[-1] not in "+-!X"]
    i.ys    = [c for c in i.all if c.txt[-1] in "+-!"]

class Data:
  """Rows plus summarized columns."""
  def __init__(i, src: Iterable[Row] | None=None):
    src    = iter(src or [])
    i.cols = Cols(next(src))  # type: ignore[arg-type]
    i.rows: Rows = []
    for r in src: add(i, r)

def clone(d: Data, rs: Rows | None=None) -> Data:
  """New Data, same cols, optional rows."""
  return Data([d.cols.names] + list(rs or []))

# ---- 3. Distance ----
def minkowski(xs: Iterable[float]) -> float:
  """Minkowski p-distance from iterable."""
  tot, n = 0.0, 1e-32
  for x in xs: tot, n = tot + x**the.p, n + 1
  return (tot / n) ** (1 / the.p)

def aha(c: Col_t, u: Atom, v: Atom) -> float:
  """Distance between two values on one column."""
  if u == v == "?": return 1
  if isa(c, Sym): return float(u != v)
  u, v = norm(c, u), norm(c, v)  # type: ignore[arg-type]
  u = u if u != "?" else (0 if v > 0.5 else 1)
  v = v if v != "?" else (0 if u > 0.5 else 1)
  return abs(u - v)

def distx(d: Data, a: Row, b: Row) -> float:
  """Distance between rows on X-cols."""
  return minkowski(aha(c, a[c.at], b[c.at]) for c in d.cols.xs)

def disty(d: Data, r: Row) -> float:
  """Distance to heaven on Y-cols."""
  return minkowski(abs(norm(c, r[c.at]) - c.heaven)  # type: ignore[operator]
                   for c in d.cols.ys)

# ---- 4. Dim ----
class Dim:
  """Fastmap axis: two poles, gap, cut points."""
  def __init__(i, east: Row, west: Row, gap: float, cuts: list[float] | None=None):
    i.east, i.west, i.gap, i.cuts = east, west, gap, cuts or []

def proj(i: Dim, d: Data, r: Row) -> float:
  """Project r onto east-west via fastmap."""
  de, dw = distx(d, r, i.east), distx(d, r, i.west)
  return (de*de + i.gap*i.gap - dw*dw) / (2*i.gap + 1e-32)

def index(i: Dim, d: Data, r: Row) -> int:
  """Bin index of r on this axis."""
  p = proj(i, d, r)
  for k, c in enumerate(i.cuts):
    if p <= c: return k
  return len(i.cuts)

def poles(d: Data, rs: Rows, prev: list[Dim]) -> tuple[Row, Row, float]:
  """Pick 2 poles maximizing residual dist after prev dims."""
  t = sample(rs, min(the.few, len(rs)))
  def r2(a: Row, b: Row) -> float:
    dx2 = distx(d, a, b) ** 2
    for p in prev:
      ax  = proj(p, d, a) - proj(p, d, b)
      dx2 = max(0, dx2 - ax*ax)
    return dx2
  g2, e, w = max(((r2(a, b), a, b) for a in t for b in t),
                 key=lambda z: z[0])
  return e, w, g2 ** .5

def newDim(d: Data, rs: Rows, epsx: float, prev: list[Dim]) -> Dim | None:
  """Build next Dim with equal-width bins; None if gap < epsx."""
  e, w, gap = poles(d, rs, prev)
  if gap < epsx: return None
  step = gap / the.bins
  cuts = [step * (i + 1) for i in range(int(the.bins) - 1)]
  return Dim(e, w, gap, cuts)

def dims(d: Data, rs: Rows) -> list[Dim]:
  """Grow Dims; stop on first no-split or cap."""
  t    = sample(rs, min(the.few, len(rs)))
  xSd  = adds(distx(d, a, b) for a in t for b in t if a is not b).sd
  epsx = the.eps * xSd
  out: list[Dim] = []
  while len(out) < the.dims:
    dim = newDim(d, rs, epsx, out)
    if dim is None: break
    out.append(dim)
  return out

def clusters(d: Data, rs: Rows, ds: list[Dim]) -> dict[tuple[int, ...], Rows]:
  """Group rows by joint Dim-index tuple."""
  out: dict[tuple[int, ...], Rows] = {}
  for r in rs:
    k = tuple(index(x, d, r) for x in ds)
    out.setdefault(k, []).append(r)
  return out

def labelled(d: Data) -> Data:
  """Clone with the.budget random rows."""
  rs = d.rows[:]; shuffle(rs)
  return clone(d, rs[:the.budget])

# ---- 5. Tests ----
def test__the():
  """Print config."""
  print(o(the.__dict__))

def test__data():
  """Load file, show shape."""
  d = Data(csv(the.file))
  print("rows", len(d.rows),
        "x", len(d.cols.xs), "y", len(d.cols.ys))

def test__dim():
  """Build dims on labelled subset; report clusters."""
  d   = Data(csv(the.file))
  lab = labelled(d)
  ds  = dims(lab, lab.rows)
  cs  = clusters(lab, lab.rows, ds)
  print("dims", len(ds), "clusters", len(cs),
        "bins", [len(x.cuts) + 1 for x in ds])

def test__4d():
  """7D data, only 4 dims have spread. Should find ~4 dims, not 7."""
  import os
  path = "etc/spread4of7.csv"
  if not os.path.exists(path):
    print("skip: run `python3 etc/spread4of7.py > etc/spread4of7.csv` first"); return
  d   = Data(csv(path))
  lab = labelled(d)
  ds  = dims(lab, lab.rows)
  cs  = clusters(lab, lab.rows, ds)
  print("4d test: dims", len(ds), "clusters", len(cs),
        "bins", [len(x.cuts) + 1 for x in ds])
  assert len(ds) <= 5, f"expected ~4 dims, got {len(ds)}"

def test__simplex():
  """Noisy 5-simplex: 6 vertices x 20 rows. Rank-5; expect ~5 dims."""
  import os
  path = "etc/simplex5.csv"
  if not os.path.exists(path):
    print("skip: run `make etc/simplex5.csv` first"); return
  d   = Data(csv(path))
  lab = labelled(d)
  ds  = dims(lab, lab.rows)
  cs  = clusters(lab, lab.rows, ds)
  print("simplex dims", len(ds), "clusters", len(cs),
        "bins", [len(x.cuts) + 1 for x in ds])

def _tests():
  return [(k, v) for k, v in globals().items()
          if k.startswith("test__") and k != "test__all"]

def test__all():
  """Run every test__ function."""
  for k, v in _tests(): print(f"# {k}"); v()

# ---- Ready, set, go. ----
the = S()
for k, v in re.findall(r"([\w.]+)=(\S+)", __doc__ or ""):
  setattr(the, k, thing(v))

def parseArgs(argv: list[str]) -> list[str]:
  """Parse --k=v flags into `the`, collect bare --t as test names."""
  tests: list[str] = []
  for a in argv:
    if a.startswith("--"):
      if "=" in a:
        k, v = a[2:].split("=", 1); setattr(the, k, thing(v))
      else: tests.append(a[2:])
  return tests

if __name__ == "__main__":
  random.seed(the.seed)
  for t in parseArgs(sys.argv[1:]):
    if fn := globals().get(f"test__{t}"): fn()
