#!/usr/bin/env python3
# Generate 7D data: 4 dims with real spread, 3 dims near zero.
# Clusters arranged so fastmap axes align with live dims.
import random

random.seed(42)
N = 200

# 5 clusters along a gradient in 4D — each step moves in live dims
centers = [
    [0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0, 1.0],
]

print(",".join(f"X{i+1}" for i in range(7)) + ",Goal-")
for c in centers:
    goal_base = sum(c)  # 0, 1, 2, 3, 4
    for _ in range(N // len(centers)):
        live = [c[j] + random.gauss(0, 0.08) for j in range(4)]
        dead = [random.gauss(0, 0.001) for _ in range(3)]
        vals = live + dead
        goal = goal_base + random.gauss(0, 0.2)
        print(",".join(f"{v:.5f}" for v in vals) + f",{goal:.4f}")
