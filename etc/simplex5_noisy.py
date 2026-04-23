#!/usr/bin/env python3
# Generate noisy samples around vertices of regular 5-simplex (pure Python).
# 6 clusters x `per` points, in R^5, side length 1, Gaussian noise.
import math, random

def simplex(d, side=1.0):
  """Return d+1 vertices of regular d-simplex in R^d, all pairwise dist=side."""
  if d == 1: return [[0.0], [side]]
  sub = simplex(d-1, side)
  c   = [sum(p[j] for p in sub)/len(sub) for j in range(d-1)]
  sub = [[p[j]-c[j] for j in range(d-1)] for p in sub]
  r   = math.sqrt(sum(sub[0][j]**2 for j in range(d-1)))
  h   = math.sqrt(side*side - r*r)
  return [p + [0.0] for p in sub] + [[0.0]*(d-1) + [h]]

random.seed(1)
d, per = 5, 20
verts  = simplex(d)

print(",".join(f"X{i+1}" for i in range(d)) + ",Goal-")
for v in verts:
  for _ in range(per):
    r = [v[j] + random.gauss(0, 0.02) for j in range(d)]
    print(",".join(f"{x:.4f}" for x in r) + ",0")
