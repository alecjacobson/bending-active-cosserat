#!/usr/bin/env python3
"""Plot the lateral-buckling bifurcation curve from results/buckling.csv.

max out-of-plane deflection max|z| vs dimensionless load gamma*, one curve per
strip width. Below a critical gamma* the strip is planar (max|z|=0); above it,
it buckles out of plane. Reproduces the qualitative result of "Better Bending"
Section 10.3 with the BAC model.
"""
import csv
import sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "results/buckling.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "docs/buckling_curve.png"

data = defaultdict(list)
with open(src) as f:
    for row in csv.DictReader(f):
        data[float(row["width"])].append(
            (float(row["gamma"]), float(row["max_lateral"]), int(row["converged"]))
        )

plt.figure(figsize=(7, 5))
colors = plt.cm.viridis
widths = sorted(data)
for i, W in enumerate(widths):
    pts = sorted(data[W])
    g = [p[0] for p in pts]
    z = [p[1] for p in pts]
    c = colors(i / max(1, len(widths) - 1))
    plt.plot(g, z, "-o", color=c, ms=4, label=f"W/L = {W:.2f}")
    # mark any non-converged points
    for gg, zz, conv in pts:
        if not conv:
            plt.plot([gg], [zz], "x", color="red", ms=8, mew=2)

plt.xlabel(r"dimensionless load  $\gamma^*$")
plt.ylabel(r"max out-of-plane deflection  $\max|z|$  (m)")
plt.title("BAC lateral buckling of a cantilevered strip (Better Bending §10.3)")
plt.gca().invert_xaxis()  # sweep goes high -> low load
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(out, dpi=130)
print("wrote", out)
