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
#   i        self          txt      str
#   k,v      key,value     i        int
"""
Options:

    --seed=1          random seed
    --p=2             Minkowski p
    --eps=0.35        margin multiplier
    --budget=50       label budget
    --bins=10         bins per dim
    --few=32          pole sample size
    --dims=10         max dims
    --elite=0.5       quantile cut for best
    --check=3         holdout labels
    --file=auto93.csv input csv
"""
import re, sys, random
from random import sample, shuffle
from math import exp,sqrt
from types import SimpleNamespace as S

isa = isinstance

### 0. Structs (simple to complex) ---------------------------------------
def Col(txt="", at=0):
  "Num if txt starts upper, else Sym."
  return (Num if txt[:1].isupper() else Sym)(txt, at)

class Num:
  "Summarize numbers."
  def __init__(i, txt="", at=0):
    i.txt, i.at, i.n = txt, at, 0
    i.mu, i.m2, i.sd = 0.0, 0.0, 0.0
    i.heaven           = txt[-1:] != "-"

class Sym:
  "Summarize symbols."
  def __init__(i, txt="", at=0):
    i.txt, i.at, i.n, i.has = txt, at, 0, {}

class Roles:
  "Organize columns from header names."
  def __init__(i, names):
    i.names = names
    i.all = [Col(t, j) for j, t in enumerate(names)]
    i.xs  = [col for col in i.all if col.txt[-1] not in "+-!X"]
    i.ys  = [col for col in i.all if col.txt[-1]     in "+-!"]

class Data:
  "Rows plus summarized columns."
  def __init__(i, src=None):
    src = iter(src or [])
    i.rows, i.roles = [], Roles(next(src))
    for row in src: add(i, row)

class Dim:
  "Fastmap axis: poles, gap, cuts, next pole hint."
  def __init__(i, east, west, gap):
    i.east, i.west, i.gap = east, west, gap
    step   = gap / the.bins
    i.cuts = [step*(k+1) for k in range(int(the.bins)-1)]
    i.north = i.south = None

### 1. Basics ---------------------------------------
def sub(i, v): return add(i, v, -1)

def adds(vs, i=None):
  "Summarize numbers into i."
  out = i or Num()
  for v in vs: add(out, v)
  return out

def add(i, v, w=1):
  "Update Num/Sym/Data; w=-1 removes."
  if v == "?": return v
  match i:
    case Data():
      add(i.roles, v, w)
      (i.rows.append if w > 0 else i.rows.remove)(v)
    case Roles():
      for col in i.all: add(col, v[col.at], w)
    case Sym():
      i.has[v] = w + i.has.get(v, 0)
    case Num() if w < 0 and i.n <= 2:
      i.n = i.mu = i.m2 = i.sd = 0
    case Num():
      i.n  += w
      delta  = v - i.mu
      i.mu += w * delta / i.n
      i.m2 += w * delta * (v - i.mu)
      i.sd  = 0 if i.n<2 else sqrt(max(0,i.m2)/(i.n-1))
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
      a, b = norm(col, u), norm(col, v)
      if a == "?": a = 0 if isa(b,float) and b > .5 else 1
      if b == "?": b = 0 if isa(a,float) and a > .5 else 1
      return abs(a - b)

def distx(data, row1, row2):
  "Distance between rows on X-cols."
  return minkowski(aha(col, row1[col.at], row2[col.at])
                   for col in data.roles.xs)

def disty(data, row):
  "Distance to heaven on Y-cols."
  return minkowski(abs(norm(col, row[col.at]) - col.heaven)
                   for col in data.roles.ys)

def proj(dim, data, row):
  "Project row onto east-west via fastmap (uses gap)."
  de = distx(data, row, dim.east)
  dw = distx(data, row, dim.west)
  return (de*de + dim.gap*dim.gap - dw*dw) / (2*dim.gap+1e-32)

def cosine(data, dim, row):
  "S(row, x=along axis, y=perp dist, theta=y/|x|)."
  a = distx(data, row, dim.east)
  b = distx(data, row, dim.west)
  c = dim.gap
  x = (a*a + c*c - b*b) / (2*c)
  y = max(0, a*a - x*x) ** .5
  return S(row=row, x=x, y=y, theta=y/abs(x+1e-32))

### 3. Fastmap ---------------------------------------
def index(dim, data, row):
  "Bin index of row on this axis."
  p = proj(dim, data, row)
  for k, cut in enumerate(dim.cuts):
    if p <= cut: return k
  return len(dim.cuts)

def newDim(data, rows, epsx, east=None, west=None):
  "Build Dim or None if gap < epsx."
  D    = lambda r1,r2: distx(data, r1, r2)
  if not (east and west):
    lst  = sample(rows, min(the.few, len(rows)))
    east = east or max(lst, key=lambda r: D(r, lst[0]))
    west = max(lst, key=lambda r: D(r, east))
  if (gap := D(east, west)) > epsx:
    dim = Dim(east, west, gap)
    lst = sample(rows, min(the.few, len(rows)))
    ps  = [cosine(data, dim, row) for row in lst]
    n   = max(ps, key=lambda p: p.y)
    s   = min(ps, key=lambda p: (p.x - n.x)**2 + p.y**2)
    dim.north, dim.south = n.row, s.row
    return dim

# do not delete. lessons learned. of ten first few
# dimensions get no cs since n+1 dims test for orthongality .
# so earleuer dimensions can be irrelvancies.
def sweepCuts(dim, data, rows, least=2):
  "Y-supervised sweep: cut where disty variance reduces."
  pairs = sorted((proj(dim, data, r), disty(data, r)) for r in rows)
  whole = adds(y for _,y in pairs)
  if whole.sd < 1e-9 or whole.n < 2*least: return []
  right = adds(y for _,y in pairs)
  left, cuts, last = Num(), [], 0
  for k,(x,y) in enumerate(pairs[:-1]):
    add(left, sub(right, y))
    if (left.n >= least and right.n >= least
        and k - last >= least
        and x != pairs[k+1][0]
        and (left.n*left.sd + right.n*right.sd) / whole.n < whole.sd):
      cuts.append((x + pairs[k+1][0]) / 2)
      last = k
  return cuts

def pruneDims(dims, data, rows, least=2):
  "Y-supervised sweep prune on each dim."
  for dim in dims:
    cuts = sweepCuts(dim, data, rows, least=least)
    if cuts: dim.cuts = cuts
  return dims

def newDims(data, rows, lab=None):
  "Grow dims on sample, then prune cuts via lab or data."
  lst = sample(rows, min(the.few, len(rows)))
  x   = adds(distx(data, *sample(lst,2)) for _ in range(the.few))
  out, east, west, now = [], None, None, {}
  while len(out) < the.dims:
    epsx = 0 if not out else the.eps*x.sd
    dim = newDim(data, rows, epsx, east, west)
    if not dim: break
    nxt = clusters(data, rows, out + [dim])
    if out and len(nxt) <= len(now): break
    out.append(dim)
    now, east, west = nxt, dim.north, dim.south
  return pruneDims(out, lab or data, lab.rows if lab else rows, 2*len(out))

def clusters(data, rows, dims):
  "Group rows by joint Dim-index tuple."
  out = {}
  for row in rows:
    k = tuple(index(dim, data, row) for dim in dims)
    out.setdefault(k, []).append(row)
  return out

def best(data, groups):
  "Clusters with mu(disty) below the.elite quantile of all rows. Best first."
  if not groups:
    rows = sorted(data.rows, key=lambda r: disty(data,r))
    n = max(1, len(rows) // int(the.bins))
    groups = {(k,): rows[k*n:(k+1)*n] for k in range(int(the.bins))}
  ys     = sorted(disty(data,r) for r in data.rows)
  thresh = ys[int(len(ys) * the.elite)]
  scored = sorted((adds(disty(data,r) for r in rs).mu, rs)
                  for rs in groups.values())
  out    = [rs for mu,rs in scored if mu < thresh]
  return out or [scored[0][1]]

def wins(data):
  "Scorer: 100=best, 0=median, relative to all rows."
  ds = sorted(disty(data, r) for r in data.rows)
  lo, med = ds[0], ds[len(ds)//2]
  return lambda r: int(100*(1 - (disty(data,r)-lo) / (med-lo+1e-32)))

def acquire(data):
  "Cluster-classify pipeline: train on B, predict on holdout, label top C."
  rows = data.rows[:]; shuffle(rows)
  mid  = len(rows) // 2
  train, test = rows[:mid], rows[mid:]
  train_data = clone(data, train)
  lab  = labelled(train_data)
  # dims from all train (x-structure), sweepCuts from lab (y-supervised)
  dims = newDims(train_data, train, lab)
  groups = clusters(train_data, train, dims)
  # mu(disty) per cluster, only from labelled rows
  lab_ids = set(id(r) for r in lab.rows)
  mu_of, worst = {}, 0
  for k, rs in groups.items():
    here = [r for r in rs if id(r) in lab_ids]
    if here:
      mu_of[k] = adds(disty(lab, r) for r in here).mu
      worst = max(worst, mu_of[k])
  # classify test rows
  def predict(row):
    k = tuple(index(dim, train_data, row) for dim in dims)
    return mu_of.get(k, worst)
  guess = sorted(test, key=predict)
  pick  = min(guess[:the.check], key=lambda r: disty(data, r))
  best_train = min(lab.rows, key=lambda r: disty(data, r))
  return best_train, pick

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
  "Build dims on sample vs all, cluster all data."
  data = Data(csv(the.file))
  lab  = labelled(data)
  dims = newDims(lab, lab.rows)
  groups = clusters(lab, data.rows, dims)
  good = best(data, groups)
  alln = sum(1 for rs in groups.values() if len(rs) >= 2)
  print(f"rows={len(data.rows)}, budget={len(lab.rows)}, ",end=" ")
  print(f"best={len(good)}, all={alln}, "
        f"dims={len(dims)},  [",end="")
  print(*[f"{len(dim.cuts)+1}" for dim in dims],sep=",",end="], ")
  print(re.sub(".*/","",the.file))

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

def test__best():
  "Count all clusters vs best (>half rows in top-sqrt(N) elite)."
  data = Data(csv(the.file))
  lab  = labelled(data)
  dims = newDims(lab, lab.rows)
  groups = clusters(lab, data.rows, dims)
  good = best(data, groups)
  whole = adds(disty(data,r) for r in data.rows)
  print(f"all={len(groups)}, best={len(good)}, "
        f"whole.mu={whole.mu:.3f}")
  for rs in good:
    s = adds(disty(data,r) for r in rs)
    print(f"  n={len(rs)}, mu={s.mu:.3f}, sd={s.sd:.3f}")

def test__acquire():
  "Cluster-classify: train wins vs test wins over 20 repeats."
  data = Data(csv(the.file))
  W = wins(data)
  train_w, test_w = Num(), Num()
  for _ in range(20):
    best_train, pick = acquire(data)
    add(train_w, W(best_train))
    add(test_w,  W(pick))
  print(f":train_wins_mu {int(train_w.mu)} :train_wins_sd {int(train_w.sd)} "
        f":test_wins_mu {int(test_w.mu)} :test_wins_sd {int(test_w.sd)}")

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
