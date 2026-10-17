#!/usr/bin/env python3
"""Octree-level histogram for a HybridOctree_Hex VTK hex mesh.

Every model in this repo is rescaled into a fixed 100-unit cube before
octree construction, so a leaf hex's edge length is exactly its octree
cell size, and its level is log2(100/edge). Established in analysis.md
(2026-08-17/18 sessions) that the METRIC used to turn a hex's 12 edges
into one length matters: the mean of all 12 edges is the one that
reproduces this repo's own previously-published per-model level counts
exactly (Bunny L4=204, Dragon Stand2 L5=735, Gargoyle L5=368, Bottle1's
Finding-11 reference column to ~1%) -- a single-edge proxy (first edge,
min edge, or max edge) distorts badly, especially at the deepest levels
where projection-stage smoothing squashes hexes unevenly. Use this tool,
not a hand-rolled single-edge version, for any future level-histogram
measurement in this analysis.

Usage: level_histogram.py <mesh.vtk> [--parents]

  --parents   also print refined-parent counts per tier (N/8 recursed
              upward from each level's leaf count) -- the sharper
              per-tier diagnostic introduced in the 2026-08-18 session,
              directly comparable across configs tier-by-tier rather
              than as a single mesh-total.
"""
import sys
import math
import collections

# The 12 edges of a VTK_HEXAHEDRON, as point-index pairs into the
# standard 8-node VTK hex connectivity order.
HEX_EDGES = [
    (0, 1), (1, 2), (2, 3), (3, 0),  # bottom face
    (4, 5), (5, 6), (6, 7), (7, 4),  # top face
    (0, 4), (1, 5), (2, 6), (3, 7),  # verticals
]


def load_vtk_hex(path):
    with open(path) as f:
        lines = f.readlines()
    i = 0
    pts = []
    cells = []
    while i < len(lines):
        if lines[i].startswith("POINTS"):
            npts = int(lines[i].split()[1])
            i += 1
            vals = []
            while len(vals) < npts * 3:
                vals.extend(float(x) for x in lines[i].split())
                i += 1
            for k in range(npts):
                pts.append((vals[3 * k], vals[3 * k + 1], vals[3 * k + 2]))
            continue
        if lines[i].startswith("CELLS"):
            ncells = int(lines[i].split()[1])
            i += 1
            for _ in range(ncells):
                vals = [int(x) for x in lines[i].split()]
                i += 1
                cells.append(vals[1:1 + vals[0]])
            continue
        i += 1
    return pts, cells


def dist(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} <mesh.vtk> [--parents]", file=sys.stderr)
        sys.exit(1)
    path = sys.argv[1]
    want_parents = "--parents" in sys.argv[2:]

    pts, cells = load_vtk_hex(path)
    levels = collections.Counter()
    skipped = 0
    for c in cells:
        if len(c) != 8:
            skipped += 1
            continue
        edge_len = sum(dist(pts[c[a]], pts[c[b]]) for a, b in HEX_EDGES) / 12.0
        if edge_len <= 0:
            skipped += 1
            continue
        lvl = round(math.log2(100.0 / edge_len))
        levels[lvl] += 1

    total = sum(levels.values())
    print(f"{path}: {total} hex cells ({skipped} skipped)")
    for lvl in sorted(levels):
        pct = 100.0 * levels[lvl] / total
        print(f"  level {lvl:2d}: {levels[lvl]:7d}  ({pct:5.2f}%)")

    if want_parents:
        print("refined-parent counts per tier (N/8 recursed upward):")
        max_lvl = max(levels)
        min_lvl = min(levels)
        # running[l] = number of level-l cells that exist, either as
        # measured leaves or as inferred parents of deeper leaves.
        running = dict(levels)
        for lvl in range(max_lvl, min_lvl, -1):
            parent_count = running.get(lvl, 0) / 8.0
            running[lvl - 1] = running.get(lvl - 1, 0) + parent_count
            print(f"  tier feeding level {lvl:2d} <- level {lvl - 1:2d}: "
                  f"{parent_count:8.1f} parents")


if __name__ == "__main__":
    main()
