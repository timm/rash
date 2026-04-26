#!/usr/bin/env python3
# Generate noisy samples around vertices of regular k-simplex embedded in R^N.
# Signal in first k cols; remaining N-k cols are near-zero-variance filler.
# Expect rash to discover ~k dims.
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
K, N, PER = 3, 8, 20       # signal rank, total cols, rows per vertex
SIG       = 0.02           # noise std on signal cols; filler cols are constant 0
verts     = simplex(K)

print(",".join(f"X{i+1}" for i in range(N)) + ",Goal-")
for i, v in enumerate(verts):
  for _ in range(PER):
    sig  = [v[j] + random.gauss(0, SIG) for j in range(K)]
    fill = [0.0] * (N - K)
    print(",".join(f"{x:.4f}" for x in sig + fill) + f",{i}")
