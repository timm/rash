#!/usr/bin/env python3 -B
# rash.py: recursive axis-split clustering
# (c) 2026 Tim Menzies, timm@ieee.org, MIT license

# Naming (data types):
#   col      Num|Sym       num      Num
#   cols     list[col]     row,rows Row, list[Row]
#   data     Data          sym      Sym
#   dim,dims Dim, list[Dim]
#
# Naming (vars):
#   at       col index     out      collector
#   lst      sample        fn       fn
#   it       self          txt      str
#   k,v      key,value     i        int
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
from math import exp,sqrt,atan2
from types import SimpleNamespace as S

isa = isinstance

### 0. Structs (simple to complex) ---------------------------------------
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

class Roles:
  "Organize columns from header names."
  def __init__(it, names):
    it.names = names
    it.all = [Col(t, j) for j, t in enumerate(names)]
    it.xs  = [col for col in it.all if col.txt[-1] not in "+-!X"]
    it.ys  = [col for col in it.all if col.txt[-1]     in "+-!"]

class Data:
  "Rows plus summarized columns."
  def __init__(it, src=None):
    src = iter(src or [])
    it.rows, it.roles = [], Roles(next(src))
    for row in src: add(it, row)

class Dim:
  "Fastmap axis: poles, gap, cuts."
  def __init__(it, east, west, gap, cuts=None):
    it.east, it.west = east, west
    it.gap, it.cuts  = gap, cuts or []

### 1. Basics ---------------------------------------
def sub(it, v): return add(it, v, -1)

def adds(vs, it=None):
  "Summarize numbers into it."
  it = it or Num()
  for v in vs: add(it, v)
  return it

def add(o, v, w=1):
  "Update Num/Sym/Data; w=-1 removes."
  if v == "?": return v
  match o:
    case Data():
      add(o.roles, v, w)
      (o.rows.append if w > 0 else o.rows.remove)(v)
    case Roles():
      for col in o.all: add(col, v[col.at], w)
    case Sym():
      o.has[v] = w + o.has.get(v, 0)
    case Num() if w < 0 and o.n <= 2:
      o.n = o.mu = o.m2 = o.sd = 0
    case Num():
      o.n  += w
      delta  = v - o.mu
      o.mu += w * delta / o.n
      o.m2 += w * delta * (v - o.mu)
      o.sd  = 0 if o.n<2 else sqrt(max(0,o.m2)/(o.n-1)) if o.n>1
  return v

def clone(data, rows=None):
  "New Data, same cols, optional rows."
  return Data([data.roles.names] + list(rows or []))

def labelled(data):
  "Clone with the.budget random rows."
  rows = data.rows[:]; shuffle(rows)
  return clone(data, rows[:the.budget])

### 2. Distance ---------------------------------------
def norm(num, v):
  "Logistic normalize Num value."
  if v == "?": return v
  z = (v - num.mu) / (num.sd + 1e-32)
  z = max(-3, min(3, z))
  return 1 / (1 + exp(-1.7 * z))

def minkowski(vs):
  "Minkowski p-distance from iterable."
  tot, i = 0.0, 1e-32
  for v in vs: tot, i = tot + v**the.p, i + 1
  return (tot / i) ** (1 / the.p)

def aha(col, u, v):
  "Distance between two values on one column."
  if u == v == "?": return 1
  match col:
    case Sym(): return float(u != v)
    case _:
      u, v = norm(col, u), norm(col, v)
      u = u if u != "?" else (0 if v > .5 else 1)
      v = v if v != "?" else (0 if u > .5 else 1)
      return abs(u - v)

def distx(data, row1, row2):
  "Distance between rows on X-cols."
  return minkowski(aha(col, row1[col.at], row2[col.at])
                   for col in data.roles.xs)

def disty(data, row):
  "Distance to heaven on Y-cols."
  return minkowski(abs(norm(col, row[col.at]) - col.heaven)
                   for col in data.roles.ys)

#   Two different jobs:
#
#   East (1 row): offaxis(row) uses cosine(...)[1] = y. Score = "how off-axis is this single row vs prev
#   dims". Pick row with max y → it points in a new direction.
#
#   West (1 row, paired with east): residual2(east, other) = pairwise distance between east and
#   candidate, MINUS distance already explained by prev dims. Picks row farthest from east in the new
#   direction.
#
#   Why different formulas:
#   - y/theta from cosine measures ONE row's deviation from ONE existing dim.
#   - For west you need a pair distance (east ↔ candidate) to define an axis with a real gap.
#
#   If you used y*theta for west: you'd get the second-most off-axis row. Risk: it sits next to east in
#   the new direction → tiny gap → useless axis.
#
#   residual2 ensures west is far from east specifically in the unexplained subspace.
#
#   If you want symmetry/simplicity:
#
#   # Pick top-2 by offaxis. Cheap but risky (poles may cluster).
#   east, west = sorted(rows[:n], key=offaxis, reverse=True)[:2]
#
#   Use only if your data spreads evenly off-axis. Otherwise keep residual2.
#
# ✻ Brewed for 46s
#
# ───────────────────────────────────────────────────────────────────────────────────────────────────────
# ❯ now i understand resutaul2.
#
#
def proj(dim, data, row):
  "Project row onto east-west via fastmap (uses gap)."
  de = distx(data, row, dim.east)
  dw = distx(data, row, dim.west)
  return (de*de + dim.gap*dim.gap - dw*dw) / (2*dim.gap+1e-32)

def cosine(data, dim, row):
  "x=along axis, y=perp dist, theta=angle from east."
  a = distx(data, row, dim.east)
  b = distx(data, row, dim.west)
  c = dim.gap
  x = (a*a + c*c - b*b) / (2*c)
  y = max(0, a*a - x*x) ** .5
  return x, y, atan2(y, x)

### 3. Fastmap ---------------------------------------
def index(dim, data, row):
  "Bin index of row on this axis."
  p = proj(dim, data, row)
  for k, cut in enumerate(dim.cuts):
    if p <= cut: return k
  return len(dim.cuts)

def nextDim(data, rows, prev):
  "Pick 2 poles maximizing residual dist."
  n = min(the.few, len(rows))
  def fn(row, dim):
    x, y, theta = cosine(data,dim,row)
    return y*theta, x, y 
  _, east, west = max((fn(rows[i], dim)
                      for i in range(n) for dim in prev)
                      key=lambda z: z[0])
  return east, west, g2 ** .5

def newDim(data, rows, epsx, prev):
  "Build next Dim; None if gap < epsx."
  e, west, gap = poles(data, rows, prev)
  if gap < epsx: return None
  step = gap / the.bins
  cuts = [step*(i+1) for i in range(int(the.bins)-1)]
  return Dim(e, west, gap, cuts)

def newDims(data, rows):
  "Grow Dims; stop on no-split or cap."
  lst  = sample(rows, min(the.few, len(rows)))
  xSd  = adds(distx(data, row1, row2) for row1 in lst for row2 in lst
              if row1 is not row2).sd
  epsx = the.eps * xSd
  out  = []
  while len(out) < the.dims:
    dim = newDim(data, rows, epsx, out)
    if dim is None: break
    out.append(dim)
  return out

def clusters(data, rows, dims):
  "Group rows by joint Dim-index tuple."
  out = {}
  for row in rows:
    k = tuple(index(dim, data, row) for dim in dims)
    out.setdefault(k, []).append(row)
  return out

### 4. Utilities ---------------------------------------
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
      row = line.partition("#")[0].split(",")
      if any(x.strip() for x in row):
        yield [thing(x) for x in row]

### 5. Tests ---------------------------------------
def test__the():
  "Print config."
  print(o(the.__dict__))

def test__data():
  "Load file, show shape."
  data = Data(csv(the.file))
  print("rows", len(data.rows),
        "x", len(data.roles.xs), "y", len(data.roles.ys))

def test__dim():
  "Build dims on labelled subset."
  data = Data(csv(the.file))
  lab  = labelled(data)
  dims = newDims(lab, lab.rows)
  groups = clusters(lab, lab.rows, dims)
  print("dims", len(dims), "clusters", len(groups),
        "bins", [len(dim.cuts)+1 for dim in dims])

def test__4d():
  "7D data, 4 dims have spread."
  import os
  path = "etc/spread4of7.csv"
  if not os.path.exists(path):
    print("skip: gen etc/spread4of7.csv first"); return
  data = Data(csv(path))
  lab  = labelled(data)
  dims = newDims(lab, lab.rows)
  groups = clusters(lab, lab.rows, dims)
  print("4d test: dims", len(dims), "clusters", len(groups),
        "bins", [len(dim.cuts)+1 for dim in dims])
  assert len(dims) <= 5, f"want ~4 dims, got {len(dims)}"

def test__simplex():
  "Noisy 5-simplex: expect ~5 dims."
  import os
  path = "etc/simplex5.csv"
  if not os.path.exists(path):
    print("skip: make etc/simplex5.csv first"); return
  data = Data(csv(path))
  lab  = labelled(data)
  dims = newDims(lab, lab.rows)
  groups = clusters(lab, lab.rows, dims)
  print("simplex dims", len(dims), "clusters", len(groups),
        "bins", [len(dim.cuts)+1 for dim in dims])

def test__all():
  "Run every test__ function; reseed each."
  for k, v in globals().items():
    if k.startswith("test__") and k != "test__all":
      random.seed(the.seed); print(f"# {k}"); v()

### 6. Main ---------------------------------------
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
