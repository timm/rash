#!/usr/bin/env python3 -B
# rash.py: recursive axis-split clustering
# (c) 2026 Tim Menzies, timm@ieee.org, MIT license

# Naming (data types):
#   c        Num|Sym       n        Num
#   cs       list[Num|Sym] r,rs     Row, Rows
#   d        Data          s        Sym
#   dim,ds   Dim, list[Dim]
#
# Naming ( vars):
#   at       col index     lab      labelled Data
#   cuts     bin edges     mu,sd,m2 moments
#   e,w,gap  poles         out      collector
#   fn       fn            prev     prior list
#   i,j      int           t        sample
#   it       self          txt      str
#   k        key           v        value
"""
Options:

    --seed=1          random seed
    --p=2             Minkowski p
    --eps=0.35        margin multiplier
    --budget=30       label budget
    --bins=6          bins per dim
    --few=32          pole sample size
    --dims=6          max dims
    --file=auto93.csv input csv
"""
import re, sys, random
from random import sample, shuffle
from math import exp
from types import SimpleNamespace as S

isa = isinstance

# ---- 0. Structs (simple to complex) ----
def Col(txt="", at=0):
  "Num if txt starts upper, else Sym."
  return (Num if txt[:1].isupper() else Sym)(txt, at)

class Num:
  "Summarize numbers."
  def __init__(it, txt="", at=0):
    it.txt, it.at, it.n = txt, at, 0
    it.mu, it.m2, it.sd = 0.0, 0.0, 0.0
    it.heaven           = txt[-1:] != "-"

class Sym:
  "Summarize symbols."
  def __init__(it, txt="", at=0):
    it.txt, it.at, it.n, it.has = txt, at, 0, {}

class Cols:
  "Organize columns from header names."
  def __init__(it, names):
    it.names = names
    it.all   = [Col(t, j) for j, t in enumerate(names)]
    it.xs    = [c for c in it.all if c.txt[-1] not in "+-!X"]
    it.ys    = [c for c in it.all if c.txt[-1]     in "+-!"]

class Data:
  "Rows plus summarized columns."
  def __init__(it, src=None):
    src = iter(src or [])
    it.rows, it.cols = [], Cols(next(src))
    for r in src: add(it, r)

class Dim:
  "Fastmap axis: poles, gap, cuts."
  def __init__(it, east, west, gap, cuts=None):
    it.east, it.west = east, west
    it.gap, it.cuts  = gap, cuts or []

# ---- 1. Basics ----
def sub(it, v): return add(it, v, -1)

def adds(vs, it=None):
  "Summarize numbers into it."
  it = it or Num()
  for v in vs: add(it, v)
  return it

def add(it, v, w=1):
  "Update Num/Sym/Data; w=-1 removes."
  if v != "?":
    match it:
      case Data():
        add(it.cols, v, w)
        (it.rows.append if w > 0 else it.rows.remove)(v)
      case Cols():
        for c in it.all: add(c, v[c.at], w)
      case Sym():
        it.has[v] = w + it.has.get(v, 0)
      case Num() if w < 0 and it.n <= 2:
        it.n = it.mu = it.m2 = it.sd = 0
      case Num():
        it.n  += w
        delta  = v - it.mu
        it.mu += w * delta / it.n
        it.m2 += w * delta * (v - it.mu)
        it.sd  = (max(0,it.m2)/(it.n-1))**.5 if it.n > 1 else 0
  return v

def clone(d, rs=None):
  "New Data, same cols, optional rows."
  return Data([d.cols.names] + list(rs or []))

def labelled(d):
  "Clone with the.budget random rows."
  rs = d.rows[:]; shuffle(rs)
  return clone(d, rs[:the.budget])

# ---- 2. Distance ----
def norm(n, v):
  "Logistic normalize Num value."
  if v == "?": return v
  z = (v - n.mu) / (n.sd + 1e-32)
  z = max(-3, min(3, z))
  return 1 / (1 + exp(-1.7 * z))

def minkowski(vs):
  "Minkowski p-distance from iterable."
  tot, i = 0.0, 1e-32
  for v in vs: tot, i = tot + v**the.p, i + 1
  return (tot / i) ** (1 / the.p)

def aha(c, u, v):
  "Distance between two values on one column."
  if u == v == "?": return 1
  match c:
    case Sym(): return float(u != v)
    case _:
      u, v = norm(c, u), norm(c, v)
      u = u if u != "?" else (0 if v > 0.5 else 1)
      v = v if v != "?" else (0 if u > 0.5 else 1)
      return abs(u - v)

def distx(d, r1, r2):
  "Distance between rows on X-cols."
  return minkowski(aha(c, r1[c.at], r2[c.at])
                   for c in d.cols.xs)

def disty(d, r):
  "Distance to heaven on Y-cols."
  return minkowski(abs(norm(c, r[c.at]) - c.heaven)
                   for c in d.cols.ys)

def proj(dim, d, r):
  "Project r onto east-west via fastmap (uses gap)."
  de = distx(d, r, dim.east)
  dw = distx(d, r, dim.west)
  return (de*de + dim.gap*dim.gap - dw*dw) / (2*dim.gap+1e-32)

# ---- 3. Fastmap ----
def index(dim, d, r):
  "Bin index of r on this axis."
  p = proj(dim, d, r)
  for k, cut in enumerate(dim.cuts):
    if p <= cut: return k
  return len(dim.cuts)

def poles(d, rs, prev):
  "Pick 2 poles maximizing residual dist."
  t = sample(rs, min(the.few, len(rs)))
  def gap2(r1, r2):
    dx2 = distx(d, r1, r2) ** 2
    for dim in prev:
      ax  = proj(dim, d, r1) - proj(dim, d, r2)
      dx2 = max(0, dx2 - ax*ax)
    return dx2
  g2,e,w = max(((gap2(r1,r2), r1, r2) for r1 in t for r2 in t),
               key=lambda z: z[0])
  return e, w, g2 ** .5

def newDim(d, rs, epsx, prev):
  "Build next Dim; None if gap < epsx."
  e, w, gap = poles(d, rs, prev)
  if gap < epsx: return None
  step = gap / the.bins
  cuts = [step*(i+1) for i in range(int(the.bins)-1)]
  return Dim(e, w, gap, cuts)

def dims(d, rs):
  "Grow Dims; stop on no-split or cap."
  t    = sample(rs, min(the.few, len(rs)))
  xSd  = adds(distx(d, r1, r2) for r1 in t for r2 in t
              if r1 is not r2).sd
  epsx = the.eps * xSd
  out  = []
  while len(out) < the.dims:
    dim = newDim(d, rs, epsx, out)
    if dim is None: break
    out.append(dim)
  return out

def clusters(d, rs, ds):
  "Group rows by joint Dim-index tuple."
  out = {}
  for r in rs:
    k = tuple(index(dim, d, r) for dim in ds)
    out.setdefault(k, []).append(r)
  return out

# ---- 4. Utilities ----
def o(x):
  "Format recursive, 2dp floats."
  if isa(x, float): return f"{x:.2f}"
  if isa(x, dict):
    return "{" + ", ".join(f"{k}={o(v)}"
           for k, v in sorted(x.items())) + "}"
  if isa(x, list): return "{" + ", ".join(map(o, x)) + "}"
  if isa(x, S):    return "S" + o(x.__dict__)
  if hasattr(x, "__dict__"):
    return x.__class__.__name__ + o(x.__dict__)
  return str(x)

def thing(txt):
  "Coerce str to int/float/bool/str."
  txt = txt.strip()
  for fn in [int, float]:
    try: return fn(txt)
    except ValueError: pass
  return {"true":True, "false":False}.get(txt.lower(), txt)

def csv(f):
  "Yield typed rows from file."
  with open(f, encoding="utf-8") as file:
    for line in file:
      r = line.partition("#")[0].split(",")
      if any(x.strip() for x in r):
        yield [thing(x) for x in r]

# ---- 5. Tests ----
def test__the():
  "Print config."
  print(o(the.__dict__))

def test__data():
  "Load file, show shape."
  d = Data(csv(the.file))
  print("rows", len(d.rows),
        "x", len(d.cols.xs), "y", len(d.cols.ys))

def test__dim():
  "Build dims on labelled subset."
  d   = Data(csv(the.file))
  lab = labelled(d)
  ds  = dims(lab, lab.rows)
  groups = clusters(lab, lab.rows, ds)
  print("dims", len(ds), "clusters", len(groups),
        "bins", [len(dim.cuts)+1 for dim in ds])

def test__4d():
  "7D data, 4 dims have spread."
  import os
  path = "etc/spread4of7.csv"
  if not os.path.exists(path):
    print("skip: gen etc/spread4of7.csv first"); return
  d   = Data(csv(path))
  lab = labelled(d)
  ds  = dims(lab, lab.rows)
  groups = clusters(lab, lab.rows, ds)
  print("4d test: dims", len(ds), "clusters", len(groups),
        "bins", [len(dim.cuts)+1 for dim in ds])
  assert len(ds) <= 5, f"want ~4 dims, got {len(ds)}"

def test__simplex():
  "Noisy 5-simplex: expect ~5 dims."
  import os
  path = "etc/simplex5.csv"
  if not os.path.exists(path):
    print("skip: make etc/simplex5.csv first"); return
  d   = Data(csv(path))
  lab = labelled(d)
  ds  = dims(lab, lab.rows)
  groups = clusters(lab, lab.rows, ds)
  print("simplex dims", len(ds), "clusters", len(groups),
        "bins", [len(dim.cuts)+1 for dim in ds])

def test__all():
  "Run every test__ function; reseed each."
  for k, v in globals().items():
    if k.startswith("test__") and k != "test__all":
      random.seed(the.seed); print(f"# {k}"); v()

# ---- 6. Main ----
the = S(**{k: thing(v) for k, v in
           re.findall(r"(\w+)=(\S+)", __doc__ or "")})

def cli(argv):
  "Parse --k=v into `the`; bare --t runs test."
  for a in argv:
    if a.startswith("--"):
      if "=" in a:
        k, v = a[2:].split("=", 1)
        setattr(the, k, thing(v))
      elif fn := globals().get(f"test__{a[2:]}"):
        random.seed(the.seed)
        fn()

if __name__ == "__main__": cli(sys.argv[1:])
