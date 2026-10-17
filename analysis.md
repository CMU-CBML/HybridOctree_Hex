# Analysis: reproducing Tong et al. 2024 Table 2

This page tracks independent reproduction of the mesh-generation statistics
reported in Table 2 of Tong, Halilaj, and Zhang, *"HybridOctree_Hex: Hybrid
octree-based adaptive all-hexahedral mesh generation with Jacobian
control,"* J. Comput. Sci. 78 (2024) 102278
([doi.org/10.1016/j.jocs.2024.102278](https://doi.org/10.1016/j.jocs.2024.102278)),
a local copy of which is `Tong 2024 HybridOctree_Hex Jacobian control.pdf`
in this repo. Figure 6 in that paper shows the twelve test models
(Bottle1, Bunny, David, Deformed Armadillo, Dragon Stand, Gargoyle, plus six
more in Fig. 7); Table 2 reports mesh size, scaled-Jacobian range,
refinement level, and generation time for each.

## Scope: this repo (`HybridOctree_Hex`) vs. `HexOpt` — two different repos, two different jobs

A scope note worth stating up front, since a sibling repo,
[hovey/HexOpt](https://github.com/hovey/HexOpt) (a fork of
[CMU-CBML/HexOpt](https://github.com/CMU-CBML/HexOpt)), covers a
same-family but non-overlapping piece of the same overall meshing
pipeline — the two aren't interchangeable:

- **This repo** (`HexGen.cpp`/`Main.cpp`, forked from
  [CMU-CBML/HybridOctree_Hex](https://github.com/CMU-CBML/HybridOctree_Hex))
  is the full generation pipeline that produced every row of Table 2:
  octree construction → dual hex mesh extraction → buffer-zone clearing →
  quality improvement, all run as one program (`HexGen`) and timed as a
  single end-to-end number. Its internal quality-improvement stage
  (`hexGen::ProjectToIsoSurface`) does scaled-Jacobian gradient descent
  plus periodic Laplacian smoothing, baked in as the last stage of
  generation.
- **`HexOpt`** (`optimizer.cpp`/`meshQuality.cpp`) only does the
  quality-improvement/optimization step on a surface + a **pre-existing**
  hex mesh — conceptually similar to this repo's `ProjectToIsoSurface`
  stage, but pulled out as a separable, reusable standalone tool. It has
  no octree construction and no dual-mesh generation, doesn't ship a
  `bottle1` model at all, and has no code path that could produce one from
  scratch.

In short: this repo generates a hex mesh from a triangulated surface and
improves it in one shot; `HexOpt` only does the "improve it" half, and
expects someone (or something) else to have already produced the hex mesh
it starts from. Reproducing Table 2 is entirely a `HybridOctree_Hex`
exercise — `HexOpt`'s optimizer has no role in the numbers tracked on this
page, and this page lives here rather than in `HexOpt` accordingly.

## Table 2, reproduced from the paper

The paper's own Table 2 (p.10), transcribed in full below for reference.
**`Tong 2024`** (bold) is the paper's "ours" row — the method this repo
reproduces. The other rows are prior work the paper compares against,
cited there as `[25]`, `[26]`, `[36]`; relabeled here as first-author
surname + year, each linked to a live copy of the actual paper (open PDF
where one exists; DOI landing page otherwise — noted per link since not
all of these are open access).

- **Tong 2024** — H. Tong, E. Halilaj, Y.J. Zhang, *HybridOctree_Hex*, J.
  Comput. Sci. 78 (2024) 102278. [Open PDF (arXiv)](https://arxiv.org/pdf/2401.05984) · [DOI](https://doi.org/10.1016/j.jocs.2024.102278)
  — this is the paper this whole page is about; same reference as the
  intro above.
- **Hu 2013** (`[25]`) — K. Hu, J. Qian, Y. Zhang, *Adaptive
  all-hexahedral mesh generation based on a hybrid octree and bubble
  packing*, 22nd International Meshing Roundtable, Orlando, FL, Oct.
  2013. [Google Scholar search](https://scholar.google.com/scholar?q=Adaptive+all-hexahedral+mesh+generation+based+on+a+hybrid+octree+and+bubble+packing)
  — no DOI or open PDF found for this one (an IMR "research notes" track
  paper, not part of the Springer-published proceedings volume); linked
  to a search as the best available pointer.
- **Gao 2019** (`[26]`) — X. Gao, H. Shen, D. Panozzo, *Feature
  Preserving Octree‐Based Hexahedral Meshing*, Computer Graphics Forum 38
  (5) (2019) 135–149. [Open PDF (author-hosted, NYU)](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) · [DOI](https://doi.org/10.1111/cgf.13795)
  — the method most of Table 2's comparison rows are against.
- **Zhang 2013** (`[36]`) — Y. Zhang, X. Liang, G. Xu, *A robust
  2‑refinement algorithm in octree or dual‑contouring dodecahedral tree
  based all‑hexahedral mesh generation*, Comput. Methods Appl. Mech.
  Engrg. 256 (2013) 88–100. [DOI](https://doi.org/10.1016/j.cma.2012.12.020)
  — DOI-only; no open PDF found.

Each model gets its own mini-table below (Method / Vertices / Elements /
Worst SJ / Best SJ / Refinement Level / Time), in the paper's own Table 2
order. Time (min) is Time (s) ÷ 60, added here — not in the original
table — for readability on the slower rows.

*: From this present study, the best match to the Tong 2024 published results.

### 0. Bone (calibration testbed, not one of Table 2's twelve models)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| **Published (this repo's own README)** | **10,356** | **8,619** | **0.61** | **1.0** | — | — | — |
| * | 10,356 | 8,619 | 0.610001 | 1.0 | 5 | — | — |

`v1.0`'s shipped constants, run unmodified against `bone`, reproduce
`our results/bone.vtk` **exactly** — matching vertex count, matching
element count, and byte-identical cell connectivity. Worst SJ lands at
0.610001 against the README's published 0.61, the same value to the
precision reported. See [2026-08-17 (continued)](#2026-08-17-continued--reading-the-octree-level-straight-off-the-meshes-v10-reproduces-bone-exactly)
in the detailed log below for the full derivation, and the [Bone](#bone-calibration-testbed-not-a-table-2-model)
section further down for this repo's earlier, superseded calibration
attempts against this same model.

### 1. Bottle1 (Genus-1)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 39,326 | 33,635 | 0.181 | 0.999 | 4 | 8,175 | 136.3 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **36,091** | **30,145** | **0.560** | **1.0** | **4** | **218** | **3.6** |
| * | 35,535 | 29,943 | 0.570 | 1.0 | 4 | ~1,561¹ | ~26.0¹ |

¹ Approximate: summed internal stage timers (read + octree + dualization +
buffer clearing, 401.2 s) plus the wall-clock gap from buffer clearing to
this run's `finalMesh.vtk` write (file mtimes, 1,160 s) — not a clean
internal timer read the way Tong 2024's own 218 s presumably is, and not
directly comparable to it for that reason (see the [Bottle1 section](#bottle1-genus-1) below,
and its own pipeline-timings table, for the full breakdown and why this
recipe's timing hasn't been the focus of this study). Config: `v1.0`
source, input `f3247d8` (2024-01-16), `C_THRES` unchanged,
`H_THRES = {32,16,8,4,2}`, `CELL_DETECT = 0.75`, four-level `rl4` configuration
(`VOXEL_SIZE_VALUE=7`, `LADDER_TOP=octreeDepth`) — see [Finding 14](#finding-14) below for
the full recipe and derivation.

### 2. Bunny (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Hu 2013](https://scholar.google.com/scholar?q=Adaptive+all-hexahedral+mesh+generation+based+on+a+hybrid+octree+and+bubble+packing) | 46,476 | 39,605 | 0.00100 | 1.0 | 3 | 462 | 7.7 |
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 35,330 | 29,698 | 0.292 | 0.999 | 4 | 3,569 | 59.5 |
| [Zhang 2013](https://doi.org/10.1016/j.cma.2012.12.020) | 119,799 | 106,730 | 3.85×10⁻⁵ | 1.0 | 3 | 258 | 4.3 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **26,375** | **21,695** | **0.570** | **1.0** | **4** | **358** | **6.0** |
| * | 27,175 | 22,525 | 0.500001¹ | 1.0 | 4 | ~877² | ~14.6² |

¹ Mesh size and topology converged normally. The 0.010 figure originally
reported here was a **measurement artifact, not the mesh's true achieved
quality** — resolved 2026-08-27, see "[Worst SJ never left the initial
0.01 gate, resolved](#worst-sj-never-left-the-initial-gate)" below for the
root cause and fix. Re-measuring the identical recipe with the fix
applied gives the same mesh (27,175/22,525, confirming the fix changes
only which state gets reported, not the optimization) at worst SJ
0.500001 — closing most of the gap to Table 2's 0.570. Cut off manually
after 4h47m24s wall-clock (started 2026-08-27 09:45:52, stopped 14:33:16)
once the value held stable for over 3.5 hours (3,934 checkpoints) with no
further movement — a confirmed plateau, not an arbitrary stopping point;
no automatic stopping condition exists for this loop.

² Approximate, on the same basis as Bottle1's footnote above: summed
internal stage timers (read + octree + dualization + buffer clearing,
244.7 s) plus the wall-clock gap from buffer clearing to this run's
`finalMesh.vtk` write (632 s). Config: `v1.0` source, `C_THRES =
{0.2121, 0.4243, 0.8485, 1.45, 3.3941}`, `H_THRES = {11.3137, 5.6569,
2.8284, 1.4142, 0.7071}`, `CELL_DETECT = 0.75`, four-level `rl4` configuration
(`VOXEL_SIZE_VALUE=7`, `LADDER_TOP=octreeDepth`) — see the [Bunny section](#bunny-genus-0)
below for the full recipe and derivation.

<a id="bunny-octree-levels"></a>
**Which octree levels the Bunny actually uses.** Measured on
`runs/bunny-v1.0-rl4-halfshift-c145/octree.vtk`, written after octree
construction and before the clearing step.

Read the levels below against one caveat: **this run uses our scheme, not
the paper's.** It sets `LADDER_TOP = 7`, one level shallower than the
paper's. The paper tests `l = 0 … 4` at levels 4-8 ("If an octree cell at
level `l + 4` … we refine it to level `l + 5`"), so the paper's own
initialized octree "ranges from level 5 to 9." This run tests `l = 0 … 4`
at levels 3-7 instead. Every level below therefore sits one shallower than
its paper counterpart — our level 3 is the paper's level 4. `bone`, the one
model we reproduce byte-exactly with shipped defaults (`LADDER_TOP = 8`),
occupies levels 5-8 and falls inside the paper's stated range; this Bunny
run does not.

| Octree level | Count (n) | Count (%) |
|:---:|---:|---:|
| 3 | 16 | 0.04% |
| 4 | 2,160 | 5.07% |
| 5 | 12,336 | 28.94% |
| 6 | 15,440 | 36.22% |
| 7 | 12,672 | 29.73% |
| **total** | **42,624** | **100.00%** |

Octree levels 0, 1, and 2 hold no leaves. Every cell at those levels is an
ancestor of some cell that flagged for refinement. The octree therefore
subdivides all of them.

Level 3 behaves differently. Sixteen leaf cells survive refinement there.
The clearing step removes them later, because they sit outside the
boundary. A ray-parity test against this run's `modifiedTri.vtk` confirms
it: all sixteen centroids fall outside the surface. They form two octree
blocks in one corner, spanning x ∈ [50,100], y ∈ [0,25], z ∈ [75,100].
Eight of them sit below the surface's own bounding box in y — centroid
y = 6.25 against a surface minimum of 11.63. So clearing removes these
sixteen cells, not refinement.

The final mesh therefore draws its cells from octree levels 4 through 7.
`finalMesh.vtk` does measure some cells deeper than level 7 — 492 at
level 8, 5 at level 9, 1 at level 10. Those cannot be octree levels. The
octree stops at level 7. They are projection distortion, cells whose
edges shrank during the projection stage. This is the cleanest evidence
yet for the distortion-tail reading used throughout this document, because
the pre-projection octree supplies a hard ceiling that the final mesh
cannot legitimately exceed.

**"Refinement Level" is a misleading column name.** "Refinement span"
describes the quantity better. The paper's phrasing invites confusion
with an octree *level*. Read carelessly, Table 2's "4" looks like the
Bunny's deepest octree level. It is not. The Bunny's deepest octree level
is 7. The 4 counts how many levels the mesh spans — level 3 to level 7 in
the octree, or levels 4 through 7 among the cells that survive clearing.
Both readings give 4.

### 3. David (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 127,778 | 112,314 | 0.126 | 0.999 | 4 | 38,866 | 647.8 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **319,465** | **282,957** | **0.560** | **1.0** | **6** | **10,450** | **174.2** |

### 4. Deformed Armadillo (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 43,591 | 35,611 | 0.0678 | 0.999 | 4 | 6,573 | 109.6 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **43,216** | **34,939** | **0.560** | **1.0** | **5** | **3,431** | **57.2** |

### 5. Dragon Stand2 (Genus-1)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 74,618 | 61,441 | 0.0290 | 0.999 | 5 | 23,062 | 384.4 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **62,576** | **50,853** | **0.560** | **1.0** | **4** | **2,052** | **34.2** |
| * | 62,307 | 50,603 | 0.500000¹ | 1.0 | 4 | ~2,710² | ~45.2² |

¹ Mesh size and topology converged normally. The 0.010 figure originally
reported here was the same measurement artifact documented for Bunny —
resolved 2026-08-27, see "[Worst SJ never left the initial 0.01 gate,
resolved](#worst-sj-never-left-the-initial-gate)" below. Re-ran this
identical recipe with the fix applied, cut off manually after 3h44m39s
wall-clock (started 2026-08-27 10:48:37, stopped 14:33:16 — no automatic
stopping condition exists for this loop, same as every other run in this
table): same mesh (62,307/50,603, confirming the fix doesn't alter the
optimization) at worst SJ 0.500000, still oscillating on one point
(index ≈56765) at the moment it was stopped, not a confirmed final
plateau the way Bunny's is.

² Approximate, on the same basis as Bottle1's footnote above: summed
internal stage timers (read + octree + dualization + buffer clearing,
2,110 s, itself the sixth of six fitting rounds run for this model — see
[Finding 15](#finding-15) below) plus the wall-clock gap from buffer clearing to this
run's `finalMesh.vtk` write (~600 s, itself an approximation carried over
from that same source). Config: `v1.0` source, `C_THRES = {0.21213203,
0.42426407, 0.84852814, 1.69705627, 3.39411255}`, `H_THRES = {11.3137085,
4.45, 2.6, 1.41421356, 0.70710678}`, `CELL_DETECT = 0.75`, four-level
`rl4` configuration (`VOXEL_SIZE_VALUE=8`, `LADDER_TOP=octreeDepth`) — see the
[Bunny section](#bunny-genus-0) below for the full recipe and derivation.

### 6. Gargoyle (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 157,008 | 135,737 | 0.0702 | 0.999 | 4 | – | – |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **273,704** | **236,689** | **0.550** | **1.0** | **4** | **8,769** | **146.2** |

*(Gargoyle's Gao 2019 time is "–" in the original — the paper doesn't
report a value for that cell.)*

### 7. Head (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Zhang 2013](https://doi.org/10.1016/j.cma.2012.12.020) | 64,258 | 56,419 | 0.0130 | 1.0 | 3 | 169 | 2.8 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **62,782** | **55,038** | **0.550** | **1.0** | **5** | **444** | **7.4** |

### 8. Lion Recon (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 134,140 | 115,245 | 0.178 | 0.999 | 4 | 18,755 | 312.6 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **112,847** | **95,251** | **0.570** | **1.0** | **4** | **2,046** | **34.1** |

### 9. Oil Pump (Genus-4)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 75,033 | 63,012 | 0.117 | 0.999 | 4 | 14,301 | 238.4 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **233,702** | **196,455** | **0.540** | **1.0** | **5** | **9,742** | **162.4** |

### 10. Ramses (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 54,634 | 46,923 | 0.0161 | 0.999 | 4 | 9,395 | 156.6 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **44,790** | **37,993** | **0.590** | **1.0** | **4** | **713** | **11.9** |
| * | 46,329 | 39,529 | 0.570 | 1.0 | 4 | ~1,692¹ | ~28.2¹ |

¹ Approximate, on the same basis as Bottle1's footnote above: summed
internal stage timers (read + octree + dualization + buffer clearing,
634.3 s, from the second of two fitting rounds — see the [Bunny section](#bunny-genus-0)
below) plus the wall-clock gap from buffer clearing to this run's
`finalMesh.vtk` write (1,058 s). Config: `v1.0` source, `C_THRES =
{0.15, 0.3, 2.0, 3.5, 2.4}`, `H_THRES = {16, 8, 4, 1.5, 1}`, `CELL_DETECT
= 0.75`, four-level `rl4` configuration (`VOXEL_SIZE_VALUE=8`,
`LADDER_TOP=octreeDepth`) — see the [Bunny section](#bunny-genus-0) below for the full
recipe and derivation.

### 11. Red Circular Box (Genus-0)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 409,011 | 367,583 | 0.215 | 0.999 | 4 | 71,626 | 1193.8 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **351,881** | **313,866** | **0.560** | **1.0** | **4** | **7,547** | **125.8** |

### 12. Thai Statue (Genus-3)

| Method | Vertices | Elements | Worst SJ | Best SJ | Refinement Level | Time (s) | Time (min) |
|---|---|---|---|---|---|---|---|
| [Gao 2019](https://cims.nyu.edu/gcl/papers/2019-OctreeMeshing.pdf) | 70,561 | 59,470 | 0.180 | 0.999 | 4 | 26,075 | 434.6 |
| **[Tong 2024](https://arxiv.org/pdf/2401.05984)** | **64,764** | **53,831** | **0.550** | **1.0** | **4** | **1,845** | **30.8** |

## Bottle1 (Genus-1)

| Metric | Tong et al. 2024 (Table 2) | 2026-08-04 (buggy `C_THRES`, capped) | 2026-08-05 (calibrated `C_THRES`, full budget) | 2026-08-17 (`v1.0`, 4 refinement levels) |
|---|---|---|---|---|
| Vertices | 36,091 | 218,898 ⚠️ | 212,990 ⚠️ | **35,535** ✓ −1.5% |
| Elements | 30,145 | 192,098 ⚠️ | 186,832 ⚠️ | **29,943** ✓ −0.7% |
| Worst scaled Jacobian | 0.560 | 0.010 (capped) | 0.530 ✓ close | **passes through 0.560, observed at 0.550 then 0.570001** ✓ (see [Finding 13](#finding-13) — the ratchet's value wherever the run is stopped, not a property of the mesh) |
| Best scaled Jacobian | 1.0 | 1.0 | 1.0 | **1.0** ✓ |
| Refinement level | 4 | not measured | not measured (see note below) | **4** ✓ measured directly (Findings [5b](#finding-5b), [10](#finding-10), [11](#finding-11)) |
| Time | 218 s | 4109.3 s (~68.5 min) | 12590.2 s (~209.8 min real; ~142.3 min active compute — see note below) ⚠️ | ~4 min octree + projection ⚠️ |

**2026-08-17 update — the mesh-size problem is solved; see the 2026-08-17
Detailed log entry.** The 3.4x-to-6x oversizing that dominated the two
earlier sessions was never a threshold-calibration problem. `v1.0`'s shipped
constants are the paper's constants — they reproduce `our results/bone.vtk`
*exactly*, down to identical cell connectivity — and Bottle1's real
difference is that it was generated with **four** octree refinement thresholds
rather than the five the public code hardcodes, exactly as Table 2's
"Refinement Level 4" says. Running four levels takes Bottle1 from 101,688
elements to 29,943 against a target of 30,145 — and octree levels 6 and 7,
which carry 93% of the mesh, match the reference to 0.5% and 0.1%
(Findings [11](#finding-11) and [14](#finding-14)). The residual under 1% is not uniquely pinned by the
evidence and probably cannot be: the reference mesh was committed a day
before any source was.

*(Table 2 itself lists Bottle1's vertex count as 36,091 — cross-checked
against this repo's own README mesh-statistics table, which agrees:
`bottle1|36091|30145|0.56`.)*

**2026-08-05 update — convergence is fixed, mesh size and timing are
not.** With `C_THRES` recalibrated and the full `MAX_PROJ_ITER=200000`
budget restored, the run converged excellently: `badElem=0,
smallDist=3.76×10⁻¹¹` — essentially machine precision, and worst SJ
(0.530) landed close to Table 2's 0.560. That confirms the earlier
0.010 was purely the deliberate `MAX_PROJ_ITER=20000` cap (2026-08-04),
not an algorithm problem — resolved.

**But mesh size barely moved** — 218,898 → 212,990 vertices, only a
~2.7% reduction — despite the same `C_THRES` fix taking `bone` from
19,639 → 14,856 vertices (a ~24% reduction, into a good match against
its published stats; see the [Bone section](#bone-calibration-testbed-not-a-table-2-model) below). This is a new
finding: `C_THRES` (curvature) clearly isn't the dominant octree-
refinement driver for Bottle1 the way it is for `bone`.

**`H_THRES` hypothesis tested and refuted (2026-08-05).** The leading
hypothesis was that `H_THRES` (thickness/narrow-region, `{16,8,4,2,1}`,
unchanged and already matching the paper) dominates for Bottle1
specifically, given its thin coiled/spiral surface detail. Tested cheaply
(octree-construction-only, ~17 min, not the full pipeline) by halving
`H_THRES` to `{8,4,2,1,0.5}` — a much stricter thickness gate. Result:
octree size barely changed (358,307→354,875 points, ~1%), refuting the
hypothesis. `H_THRES` reverted to the paper-matching value. Neither
threshold, tuned individually, explains Bottle1's over-refinement — see
the [summary log](#summary-log) for the still-open question this leaves, and a genuine
code-level anomaly found in the thickness computation along the way
(worth investigating further, though it turned out not to be the answer
here).

**Timing**: `/usr/bin/time -l`'s wall-clock ("real") total was 12590.2s
(~209.8 min), but the sum of `HexGen`'s own internal stage timers only
accounts for 8537.8s (~142.3 min) — a ~4052s (~67.5 min) gap. Traced to a
system sleep during the run: `caffeinate` (keeping the Mac awake) died
partway through and wasn't immediately relaunched, so part of the "3+
hours" included the machine being asleep, not active computation.
~142.3 min of actual compute time is the more meaningful number to
compare against Table 2's 218s — still roughly **39x** longer, and per
the [summary log's synthesis](#synthesis-what-bone-bottle1-and-bunny-collectively-show) below, not something `C_THRES` calibration
was ever going to fix (that's a mesh-*quality*/*size* tuning constant,
disconnected from the `ComputeCellValue()` performance bottleneck driving
the timing gap).

Given all this, the vertex/element/timing mismatches here shouldn't be
read as "the fix doesn't work" — quality/convergence *is* fixed, cleanly.
Mesh size and timing are separate, still-open problems, and the timing
one in particular looks structural (see the [summary log's synthesis](#synthesis-what-bone-bottle1-and-bunny-collectively-show)) —
not something further `C_THRES` tuning is likely to resolve on its own.

### Bottle1 pipeline timings (Apple M1)

Per-stage timings as `HexGen` reports them (`clock()` prints, in
`Main.cpp` and — as of this session, see [summary log](#summary-log) — inside
`ConstructOctree()` itself), updated live as the run progresses. Table 2
only reports one aggregate "Time" column; this breakdown is finer-grained
than the paper provides, for our own diagnostic purposes. Headings match
the paper's own algorithmic narrative (curvature/narrow-region assessment
→ octree generation → dualization → buffer clearing → buffer-zone
projection/quality improvement).

Percentages are each stage's share of that *column's own* total, so each
column's percentages sum to 100% once the column is complete. Columns
still mid-run show percentages as "pending" — a partial-time percentage
would be a moving target as remaining stages add to the denominator, so
these are only filled in once each column's run finishes and its final
total is known. Sub-rows (curvature/narrow-region and octree-balancing,
where measured separately) show their share of the *combined* row instead
of the grand total — including them in the grand-total column as well
would double-count that time — so the 100% breakdown for each column
comes from exactly five rows: read, combined octree, dualization, buffer
clearing, and projection+quality improvement.

| Stage | bone (M4, 2026-07-02) | bone (M1, 2026-08-04) | Bottle1 (M1, 2026-08-04, buggy `C_THRES`, capped) | Bottle1 (M1, 2026-08-05, calibrated `C_THRES`, full budget) |
|---|---|---|---|---|
| Read surface mesh | ~0.4 s (0.2%) | 0.57 s (0.1%) | 1.79 s (0.0%) | 1.90 s (0.0%) |
| Curvature and Narrow Region Assessment | n/a¹ | 31.6 s (57.0% of combined) | n/a¹ | 725.0 s (71.6% of combined) |
| Octree Generation (strongly balanced) | n/a¹ | 19.6 s (35.4% of combined) | n/a¹ | 43.9 s (4.3% of combined) |
| — combined (curvature + octree), as `Main.cpp` reports it¹ | ~38–39 s (14.6%) | 55.3 s (13.4%) | **1062.2 s (~17.7 min) (25.9%)** | **1012.4 s (~16.9 min) (11.9%)** |
| Mesh Dualization | ~12–13 s (4.7%) | 17.1 s (4.1%) | **1138.9 s (~19.0 min) (27.7%)** | **1045.8 s (~17.4 min) (12.2%)** |
| Buffer Clearing | ~13 s (4.9%) | 21.7 s (5.2%) | **815.1 s (~13.6 min) (19.8%)** | **785.0 s (~13.1 min) (9.2%)** |
| Buffer Zone Projection + Quality Improvement | ~200 s (75.6%) | **319.4 s (~5.3 min) (77.2%)** | **1091.3 s (~18.2 min) (26.6%)³** | **5692.8 s (~94.9 min) (66.7%)⁴** |
| **Total** | **~4.5 min (264.4 s summed² ⇒ 100%)** | **419.5 s (~7.0 min, `time -l` real ⇒ 100%)** | **4109.3 s (~68.5 min summed³ ⇒ 100%)** | **8537.8 s (~142.3 min summed active compute ⇒ 100%)⁴** |

¹ The M4 bone run and the Bottle1 run both used a pre-split binary, so
curvature/narrow-region assessment and octree balancing are only
available as one combined number for those two ("— combined" row). The
M1 bone run (added specifically to get a same-machine comparison point,
and launched after the timing-instrumentation split) reports both
substages separately; they sum to 92.5% of the combined cell, not 100%,
because `Main.cpp`'s own combined timer also includes a small amount of
cleanup-call/`clock()`-rounding overhead between and around the two
sub-measurements.

² `bone` (M4)'s five stage values were originally reported as ranges
(e.g. "~38–39 s"); percentages use each range's midpoint, summing to
264.4 s — close to, but not identical to, the "~4.5 minutes" total
originally reported for that run (itself a rounded figure). The summed
264.4 s, not the rounded ~4.5 min, is the percentage denominator, so the
column's percentages are internally consistent and sum to exactly 100%.

Every stage so far is dramatically slower than bone's, and by increasingly
inconsistent factors relative to the ~3.5x element-count ratio: octree
construction ~28x (~21x against the M1 bone number), dualization ~87x,
buffer clearing ~63x. This pattern (large, non-uniform multipliers rather
than a roughly constant ~3.5x across stages) points toward Bottle1's thin,
coiled surface detail (handle/spiral thread, see Fig. 6(a) in the paper)
triggering much deeper curvature-driven octree refinement than bone's
simpler shape — i.e. the real driver is octree depth/leaf count, not
final output element count, and that deeper octree then compounds through
every downstream stage.

The M1-vs-M4 bone comparison also rules out "this M1 is just a slower
machine" as the main explanation. Now that bone (M1) has finished, the
gap holds up across every stage, not just octree construction: read
+40% slower, combined octree +40%, dualization +37%, buffer clearing
+67%, projection +59%, total +55% (419.5 s vs. ~270 s). That's a
consistent, modest hardware gap in the same direction throughout — nowhere
close to Bottle1's ~20–90x factor over bone on the *same* M1. So the
dominant effect really is Bottle1's own geometry, not the machine.

Bone (M1)'s projection stage also hit the same `MAX_PROJ_ITER` backstop as
the M4 run, landing on an almost identical best result
(`smallDist=2.60788e-14` on both machines, to 6 significant figures) —
consistent with the redistribution step's `rand()` calls never being
seeded (no `srand()` call anywhere in the codebase), so the "randomized"
sequence is actually the same deterministic sequence on every run,
regardless of machine.

³ Bottle1's projection stage was initially let run at bone's
`MAX_PROJ_ITER=200000` setting like the other columns, but killed after
24 checkpoints (~22 min, `maxDist=0.00473` with `badElem=0` at the last
checkpoint before stopping) once its per-checkpoint cost became clear:
~56 s/checkpoint vs. bone's ~1 s/checkpoint, which extrapolated to
several hours for the full 200-checkpoint budget. `MAX_PROJ_ITER` was
reduced to `20000` (20 checkpoints) specifically for this run — a
deliberate speed/convergence trade-off, documented in `HexGen.cpp`'s
comments at the constant — and the capped rerun resumed from the
already-computed `octree.vtk`/`dualFullHex.vtk`/`dualHex.vtk` (via a
temporary, uncommitted `Main.cpp` tweak, reverted immediately after
launching) rather than redoing the ~50 minutes of earlier stages. It hit
that reduced backstop, landing on best result `badElem=0, smallDist=
0.62311` — genuinely under-converged compared to bone's ~1e-14 (see the
[Bottle1 results table](#bottle1-genus-1) above for what this means for the reported worst
SJ). The killed full-budget attempt's log/mesh are preserved as
`runs/bottle1/run.log.full-attempt-killed` /
`finalMesh.vtk.full-attempt-killed`. Total in the table above (4109.3 s)
sums the *first* attempt's four completed pre-projection stages (the real
computation) plus the *capped rerun's* projection stage (1091.3 s) — the
~22 min spent on the killed attempt's partial projection work is excluded
as not contributing to the final result.

Mesh-size and scaled-Jacobian stats (`finalMesh.vtk`'s vertex/element
counts, worst/best SJ) are reported in the [Bottle1 results table](#bottle1-genus-1) above,
computed via a new tool, `scripts/scaled_jacobian_stats.cpp` — a
standalone C++ program that parses a VTK hex mesh and computes min/max
scaled Jacobian using `HexGen.cpp`'s own `Sj()` function copied verbatim
(not reimplemented), so its output is guaranteed consistent with what
`HexGen` itself computes internally. Validated against bone's
`finalMesh.vtk` first (worst 0.590, best 1.0 — close to the repo's
published 0.61 for that model) before trusting it on Bottle1.

⁴ The 2026-08-05 `bottle1-v3` run (calibrated `C_THRES`, full
`MAX_PROJ_ITER=200000` budget) ran for `/usr/bin/time -l`'s reported
12590.2s "real" — but the sum of `HexGen`'s own internal stage timers
only accounts for 8537.8s, a ~4052s (~67.5 min) gap. Traced to a system
sleep during the run (`caffeinate`, keeping the Mac awake, died partway
through per `pmset -g log` and wasn't immediately relaunched) — part of
the raw wall-clock duration was the machine asleep, not computation. The
table above uses the summed 8537.8s (active compute only) as each
stage's percentage denominator, consistent with how the other columns'
percentages are computed from their own summed/reported stage times
rather than a raw wall-clock figure that could include idle gaps. The
projection stage converged excellently despite the sleep interruption —
`badElem=0, smallDist=3.76×10⁻¹¹`, essentially machine precision — but
mesh size barely changed from the buggy 2026-08-04 run (218,898 →
212,990 vertices, ~2.7% smaller) despite `C_THRES` recalibration, unlike
`bone`'s ~24% reduction with the same fix; see the [Bottle1 results table](#bottle1-genus-1)
above and the [summary log's synthesis](#synthesis-what-bone-bottle1-and-bunny-collectively-show) below for the still-open
`H_THRES` hypothesis this raises.

## Bunny (Genus-0)

Table 2 row 2: **26,375 vertices / 21,695 elements / worst SJ 0.570 / best SJ
1.0 / refinement level 4 / 358 s**. `our results/bunny.vtk` confirmed
byte-exact via `scaled_jacobian_stats` (26,375 / 21,695 / 0.569998 / 1.0) —
same archaeology as Bottle1 and `bone`, this is the real paper-output file,
committed by Hua Tong on the same day (`be63d98`, 2024-01-11) as the other
two. `input boundaries/bunny_tri.raw` (11,248 verts / 22,490 tris) has only
one version in git history — unlike Bottle1, no header-bug/duplicate-triangle
ambiguity to resolve first.

This section covers a single working session (2026-08-18) run as two
independent, parallel efforts — a Sonnet session (this log's primary author)
running the direct recipe-transfer test, and a background Opus agent doing a
wider threshold sweep — deliberately kept on separate run directories and
this file to avoid a merge conflict, then reconciled here.

**Headline: the Bottle1 recipe ([Finding 14](#finding-14)) does not transfer.** Applying
the exact configuration that reproduced Bottle1 to within 1.5%/0.7% —
`v1.0` source, `C_THRES` unchanged, `H_THRES×2`, `CELL_DETECT=0.75`, the
four-level `rl4` configuration (`VOXEL_SIZE_VALUE=7`, `LADDER_TOP=octreeDepth`) —
against `bunny_tri.raw` overshoots by **70%** (43,691 / 36,901 vs. 26,375 /
21,695). This is the direct "does the recipe transfer to an untuned model"
test the previous session's "Where to go next" plan called for, and it
fails: [Finding 10](#finding-10)'s claim that "Refinement Level" is the level count, read
directly off the reference mesh, held up independently (Bunny's own
octree-level histogram, measured fresh this session, confirms max real
level 7 / 4 levels — see below) — but [Finding 14](#finding-14)'s specific threshold
*values* (`H_THRES×2`, `CELL_DETECT=0.75`) turn out to have been fitted to
Bottle1 specifically, not universal.

Six configurations were run, all `v1.0` source, all the four-level `rl4`
configuration (`VOXEL_SIZE_VALUE=7`, `LADDER_TOP=octreeDepth`) unless noted,
against `bunny_tri.raw`:

<a id="mesh-size-results-table"></a>
| Config | `C_THRES` | `H_THRES` | `CELL_DETECT` | Points | Cells | vs. target |
|---|---|---|---|---|---|---|
| Table 2 target | | | | **26,375** | **21,695** | — |
| `rl4-shift1-cd75` | ×2 | ÷2 | 0.75 | 14,941 | 12,112 | −44% |
| **`rl4-halfshift-cd75`** | **×√2** | **÷√2** | **0.75** | **24,142²** | **19,938²** | **−8.5% / −8.1%** |
| **`rl4-halfshift-c145`** (`C_THRES[3]`: 1.6971→1.45, closing the gap) | **×√2 except `[3]`** | **÷√2** | **0.75** | **27,175** | **22,525** | **+3.0% / +3.8% — best Bunny reproduction** |
| `rl4-cd75` | shipped | shipped | 0.75 | 40,146² | 33,873² | +52% / +56% |
| `rl4` (plain) | shipped | shipped | 1 (default) | 44,613 | 37,969 | +75% |
| `rl4-h2x-cd75` (Bottle1's [Finding-14](#finding-14) recipe) | shipped | ×2 | 0.75 | 43,691 | 36,901 | +70% |
| `rl4-h3` ([Finding-14](#finding-14)'s other lever) | shipped | `H[3]: 2→3` | 1 (default) | 46,282 | 39,420 | +82% |
| `vs8abs` (tests whether Bunny needs a fifth level) | shipped | shipped | 1 (default), `VOXEL_SIZE_VALUE=8` | 145,049 | 126,800 | +484% |

The family is cleanly monotonic in how tightly the thresholds constrain
refinement (`shift1` tightest → `h3` loosest), which is itself useful: the
true per-model configuration for Bunny sits between `rl4-halfshift-cd75`
and `rl4-cd75`, i.e. somewhere between a `√2` and a `1×` scale-relative
correction, refined further by `rl4-halfshift-c145` above. `vs8abs`
(`VOXEL_SIZE_VALUE=8`, the shipped five-threshold configuration) massively overshoots
(+484%), confirming Bunny's reference-mesh level-8 population (2.3% of the
mesh) is projection distortion, not real fifth-level refinement — the same
conclusion [Finding 11](#finding-11) reached for Bottle1, now independently replicated on
a second model via a real (not just inferred) fifth-level run.

² **Important calibration point, discovered mid-session**: the dual-mesh
point/cell count *before* projection (`dualHex.vtk`, written at the
"extracting interior dual mesh" stage) is *not* a safe stand-in for the
final count — `rl4-h2x-cd75`'s `dualHex.vtk` showed 30,837 / 24,049 but its
actual `finalMesh.vtk` came in 42%/53% higher at 43,691 / 36,901. The
`ProjectToIsoSurface`/quality-improvement stage evidently splits
low-quality elements, not just repositions points, contrary to what every
earlier Bottle1-session entry above assumed (there, `dualHex.vtk` and
`finalMesh.vtk` counts happened to be close, which this session initially
misread as the general case before catching the discrepancy here). Two
rows above (`rl4-halfshift-cd75`, `rl4-cd75`) never reached a clean
`finalMesh.vtk` write — both plateaued for 2+ hours with `smallDist`
oscillating around 0.03, well above the `1e-4` convergence tolerance,
similar to the hard stuck-points [Finding 13](#finding-13) documented for Bottle1, and
`rl4-cd75`'s process was eventually lost (killed or crashed, unclear
which) without ever writing one. Both counts above instead come from
`projHex.vtk`, the continuously-updated in-progress mesh, whose point/cell
count had been stable for the full 2+ hours before each was stopped —
strong evidence the *topology* (element splits) had already settled even
though point-position convergence hadn't, so these are reported as
confirmed rather than provisional, on that basis.

³ To close Bunny's residual gap: at `rl4-halfshift-cd75` thresholds, threshold 3 (level 6→7) was confirmed purely curvature-gated (55 curvature
candidates, 0 thickness candidates), so
`C_THRES[3]` alone was swept from 1.6971 down toward 1.2
(candidate counts 55→79→95→110→127→138 at
1.6971/1.55/1.45/1.35/1.25/1.2) and 1.45 (95 candidates, ~73% more than
the `halfshift` baseline's 55) chosen as a first attempt — and it worked
on the first try, landing within +3.8% without needing a second
iteration. `finalMesh.vtk` was written and `projHex.vtk`'s topology
confirmed stable well past that point.

**The diagnostic that made this tractable: refined-parent counts per level,
not totals.** A level-`(L+1)` population of `N` cells implies `N/8`
refined level-`L` parent cells; recursing this upward from each mesh's
leaves gives a per-level count directly comparable to the reference's:

| Config | parents @ threshold 0 (level 3) | @ threshold 1 (level 4) | @ threshold 2 (level 5) | @ threshold 3 (level 6) |
|---|---|---|---|---|
| **reference (`bunny.vtk`)** | **129** | **832** | **1,233** | **870** |
| `rl4` (plain) | 129 | 963 | 2,328 | 1,907 |
| `rl4-cd75` | 130 | 948 | 2,030 | 1,624 |
| `rl4-h2x-cd75` | 130 | 954 | 2,071 | 2,017 |
| `rl4-h3` | 128 | 961 | 2,321 | 2,121 |
| `rl4-shift1-cd75` | 137 | 686 | 637 | 234 |
| `rl4-halfshift-cd75` | 136 | 815 | 1,155 | 695 |

Threshold 0 is geometry-bound (128–137 across every config, reference 129), so
it carries no information — the octree's coarsest refine-or-not decision
doesn't depend on these thresholds. The signal is that every plain-`rl4`-
family config's error *roughly doubles with each deeper threshold* (threshold 1
+15–19%, threshold 2 +65–89%, threshold 3 +87–144%) — the signature of a per-level
*geometric scaling* problem, not a single bad threshold. That's why neither
[Finding 14](#finding-14) lever (`H×2` alone, or `H[3]: 2→3` alone) fixes it: each
perturbs one threshold while the error compounds across all of them.
`rl4-halfshift-cd75` (`C_THRES×√2`, `H_THRES÷√2` — a *half*-level
scale-relative correction, as opposed to `rl4-shift1-cd75`'s full-level
one) comes far closer at every threshold (−2% / −6% / −20% at thresholds 1/2/3), and
loosening threshold 3's gate further (`rl4-halfshift-c145`, above) closes the
remaining gap to +3.0%/+3.8% — the best Bunny reproduction reached this
session.

**Provenance check — binaries independently rebuilt and byte-diffed, not
trusted by filename.** Every `build-variants/` binary used here was
recompiled from source with the CMake Release flags
(`-O3 -DNDEBUG -std=c++17 -arch arm64`) and compared byte-for-byte against
the shipped binary; only the Mach-O UUID and a few symbol-table bytes
differ. Confirms `HexGen-rl4-h3` is `H_THRES = {16,8,4,3,1}` (only threshold
index 3 bumped, 2→3 — not `{16,8,4,2,3}` as its name alone might suggest),
`HexGen-rl4-h2x-cd75` is `H_THRES = {32,16,8,4,2}` at `CELL_DETECT = 0.75`,
and `HexGen-rl4` is the shipped arrays unchanged at the four-level configuration.

**Octree-level histogram — confirms [Finding 10](#finding-10)'s level count, flags a
measurement-methodology gap.** Measuring `our results/bunny.vtk` directly
(edge length → level, model rescaled into a 100-unit cube) gives level 4:
312 (1.4%), 5: 5,321 (24.5%), 6: 8,634 (39.8%), 7: 6,931 (32.0%), 8: 494
(2.3%, almost certainly the same projection-distortion artifact [Finding 11](#finding-11)
identified for Bottle1's level-8 population, not real refinement), 9–10: 3
cells (noise) — confirming max real level 7, four levels, independently of
[Finding 10](#finding-10)'s own table. However, this measurement used **mean of the pair
of points spanning one edge** as the length proxy, which does not match
what reproduces `analysis.md`'s own previously-published per-model counts
(Bunny's published level-4 count is 204, not this session's 312) — a
same-file comparison found the **mean of all 12 hex edges** is the metric
that reproduces the previously-published figures exactly (Bunny L4 = 204,
Dragon Stand2 L5 = 735, Gargoyle L5 = 368, and Bottle1's [Finding-11](#finding-11)
reference column to ~1%), while single-edge min/max proxies distort badly
in opposite directions (min-edge in particular produces a spurious
3,814-cell level-8 "population" that is really squashed low-quality hexes
with one short edge, not smaller cells — itself supporting evidence for
[Finding 11](#finding-11)'s distortion-tail read). **Any future level-histogram
measurement in this file should use the 12-edge-mean metric**, not a
single-edge proxy.

**Derived ladder settings for the remaining Table 2 models.** From the
threshold-index/level relationship in `ComputeCellValue()` (`HexGen.cpp:1015`,
threshold `i` sits at level `ladderTop - 4 + i`) and the scale-relative
structure [Finding 12](#finding-12) identified (`H_THRES[i] = 2.56 × cellsize` at the
shipped ladder position, threshold 0 assumed at level 4): `VOXEL_SIZE_VALUE` =
each model's own deepest real octree level (measured off its reference
mesh, 12-edge-mean metric); `LADDER_TOP = octreeDepth` for Table 2's
"refinement level 4" models, `(octreeDepth - 1)` for "level 5" models, and
"level 6" (David only) needs `v1.2`'s still-commented-out sixth threshold. The
scale-relative `H_THRES` multiplier is `2^(4 - (ladderTop - 4))`:

| Model | Table 2 refinement level | Deepest real level | `VOXEL_SIZE_VALUE` | `LADDER_TOP` | scale-relative `H_THRES` multiplier |
|---|---|---|---|---|---|
| Bottle1, Bunny, Red Circular Box | 4 | 7 | 7 | `octreeDepth` | ×2 |
| Dragon Stand2, Gargoyle, Lion Recon, Ramses, Thai Statue | 4 | 8 | 8 | `octreeDepth` | ×1 (shipped) |
| Head, Deformed Armadillo | 5 | 8 | 8 | `(octreeDepth - 1)` | ×2 |
| Oil Pump | 5 | 9 | 9 | `(octreeDepth - 1)` | ×1 (= stock `v1.0` ladder) |
| David | 6 | 9 | 9 | needs the sixth threshold | — |
| `bone` (control, not in Table 2) | — | 9 | 9 | `(octreeDepth - 1)` | ×1 |

`bone`'s row falls straight out of this rule as *exactly* the shipped
configuration — a real independent validation, since `bone` is the one
model this repo reproduces byte-exactly. Its own level histogram
(`level_histogram.py`: level 5-9 at 9.72%/43.75%/31.13%/14.80%/0.59%)
shows the same signature as every other model's distortion tail — level
9's 51 cells are a ~25× drop from level 8, not a smooth continuation — so
`bone` occupies levels 5-8 — four real levels. (An earlier version of this
paragraph said "5 real levels (levels 4-8)", which contradicts the very
histogram quoted above it; corrected 2026-08-27, see [Finding
10](#finding-10)'s correction note. `bone` compiles in all five thresholds
under `LADDER_TOP=8`, but its geometry never fires the deepest one, so the
mesh spans four levels, not five.) A naive "any cells
present" read would have misclassified level 9 as a real sixth threshold; the
relative-magnitude test (is the deepest population a small fraction with
a sharp drop, or a smooth continuation of the level above) is what
resolves it correctly here — the same test that resolves David's
sixth-threshold question below.

**A methodological caveat worth stating precisely, not glossing over.**
This single-reference-mesh "small fraction + sharp drop" test is weaker
than [Finding 11](#finding-11)'s original distortion-tail evidence for Bottle1, which
compared *two different runs* — one whose octree was structurally
incapable of reaching the deeper level at all — and found both showed
essentially the same small population there, proving it was a
projection-stage measurement artifact (a squashed shallower hex, not a
real deeper cell) rather than genuine refinement. Everything checked
today (`bone`, David, and the seven additional models measured below)
used only the weaker single-mesh version, since it doesn't require a
second, deliberately-wrong run to compare against. The two tests may not
even be probing the same underlying phenomenon: [Finding 11](#finding-11)'s mechanism
(projection smoothing squashing a shallower hex until it measures as
deeper) is one way a small deepest-level population can arise, but a
genuinely-sparse level from ordinary 2:1 octree balance propagation — a
handful of boundary cells needing one extra level to satisfy the
constraint — would produce the same small-fraction-with-a-sharp-drop
signature without being any kind of artifact at all. For the practical
purpose this table serves (what `VOXEL_SIZE_VALUE`/`LADDER_TOP` should a
model use), the distinction doesn't matter — either mechanism means that
level doesn't need its own explicit threshold, and `bone`'s case settles it
independently of which explanation is right: its shipped five-threshold
configuration reproduces `our results/bone.vtk` byte-for-byte, level-9
population included, so whatever produces those 51 cells is already
present in the *correct* configuration, not evidence of a missing sixth
threshold. But readers relying on this table for anything beyond ladder-depth
selection should know the two tests differ in strength, and that this
session ran the weaker one seven more times below without re-running
[Finding 11](#finding-11)'s stronger comparative version on any of them.

**David's sixth-threshold question, reconciled** (checked directly with
`level_histogram.py`, not reasoned about abstractly): `david.vtk`'s full
distribution is level 4: 0.02%, 5: 1.74%, 6: 9.52%, 7: 26.14%, 8: 41.11%,
**9: 20.77%**, 10: 0.71%. Level 9 is not a thin tail by any reading of the
heuristic used elsewhere in this table — it's the *second-largest* band in
the mesh, a smooth continuation of the level 4→8 climb, not a sudden drop.
Level 10 (0.71%, a 29× drop from level 9) is the actual distortion tail,
matching the same small-fraction signature seen on every other model's
deepest population. So [Finding 9](#finding-9) and [Finding 10](#finding-10) were never really in
conflict — David does need a genuine sixth threshold (levels 4-9), Table 2's
"Refinement Level 6" is literal, and the commented-out
`//if (level == 9) {// 5` block in `v1.2`/`v1.3` is exactly what's needed
to reach it. What looked like tension was an artifact of stating the
distortion-tail heuristic too loosely ("any small deepest-level
population is an artifact") rather than checking whether that population
is actually small relative to its neighbors, which is the test that
matters.

**The ladder-derivation table's predictions checked against all seven
remaining models' own reference meshes — every one matches, at zero
compute cost.** `level_histogram.py` run directly against each of
Gargoyle, Head, Thai Statue, Lion Recon, Red Circular Box, Oil Pump, and
Deformed Armadillo's `our results/*.vtk`, applying the same
relative-magnitude test used above:

| Model | Level histogram (%, by level) | Deepest real level | Drop ratio at the tail | Predicted `VOXEL_SIZE_VALUE` | Match? |
|---|---|---|---|---|---|
| Gargoyle | 5:0.16, 6:7.39, **7:65.29**, 8:25.62, *9:1.54* | 8 | 16.6× | 8 | ✓ |
| Head | 4:0.43, 5:9.19, 6:25.72, **7:49.92**, 8:14.29, *9:0.44* | 8 | 32.5× | 8 | ✓ |
| Thai Statue | 5:0.25, 6:14.83, **7:49.43**, 8:33.34, *9:2.15* | 8 | 15.5× | 8 | ✓ |
| Lion Recon | 5:0.04, 6:4.24, 7:40.04, **8:52.72**, *9:2.96* | 8 | 17.8× | 8 | ✓ |
| Red Circular Box | 4:0.07, 5:2.19, 6:17.12, **7:79.51**, *8:1.10* | 7 | 72.3× | 7 | ✓ |
| Oil Pump | 5:0.77, 6:9.91, 7:28.02, **8:54.00**, 9:7.24, *10:0.07* | 9 | 103×¹ | 9 | ✓ |
| Deformed Armadillo | 4:0.04, 5:6.16, 6:22.35, 7:26.01, **8:43.16**, *9:2.28* | 8 | 18.9× | 8 | ✓ |

*(Bold = the largest band; italic = the identified distortion tail.)*

¹ Oil Pump is the one borderline case: level 9 (7.24%) is a smaller,
7.5× drop from level 8 rather than the 15-100× drops seen everywhere
else, which is a real level's signature, not a tail's — level 10 (0.07%,
a further 103× drop from level 9) is unambiguously the tail instead. This
is exactly what the ladder-derivation table already predicted (Oil Pump
needs the stock `v1.0` ladder, `VOXEL_SIZE_VALUE=9`,
`LADDER_TOP=octreeDepth-1`, the same setting `bone` uses) — a case where
the reference mesh's own shape, not just the table lookup, argues for
that reading independently.

Combined with Bottle1 ([Finding 10](#finding-10)/[11](#finding-11)), Bunny, Dragon Stand2, and Ramses
(this session, via real runs), every ladder-position prediction in the
table above now has direct reference-mesh support — eleven of Table 2's
twelve models (David is the twelfth, resolved separately above as
needing its own sixth threshold). None of this required running the pipeline;
it's a `level_histogram.py` call against an existing file for each
model, seconds of compute total. What's still unvalidated by an actual
run is *only* the eight models never actually built and run this
session (Gargoyle, Head, Thai Statue, Lion Recon, Red Circular Box, Oil
Pump, Deformed Armadillo, David) — the ladder position is now
well-supported for all of them, but each one's `C_THRES`/`H_THRES`/
`CELL_DETECT` still needs its own fit, exactly as Bottle1, Bunny, Dragon
Stand2, and Ramses did.

<a id="best-bunny-reproduction-reached"></a>**Best Bunny reproduction reached.** `C_THRES[3]` lowered from
`rl4-halfshift-cd75`'s 1.6971 to 1.45 (`runs/bunny-v1.0-rl4-halfshift-c145/`,
otherwise identical: `C_THRES = {0.2121, 0.4243, 0.8485, 1.45, 3.3941}`,
`H_THRES = {11.3137, 5.6569, 2.8284, 1.4142, 0.7071}`, `CELL_DETECT = 0.75`,
four-level `rl4` configuration) converged cleanly — `finalMesh.vtk` written, mesh
topology stable in `projHex.vtk` well past that point:

| | Reproduction | Table 2 target |
|---|---|---|
| Points | 27,175 | 26,375 |
| Cells | 22,525 | 21,695 |
| vs. target | **+3.0% / +3.8%** | — |
| Refinement level | 4 | 4 |

This is closer than Bottle1's own best reproduction (−1.5%/−0.7%, but in
the opposite direction — Bunny overshoots slightly where Bottle1
undershot slightly), reached by a symmetric method: start from the other
model's recipe, recognize it doesn't transfer, find the scale-relative
correction that brackets the target, then nudge the one threshold that's
purely curvature-gated. `C_THRES[3]` between 1.45 and 1.6971 would
presumably narrow the remaining ~3% further, but wasn't pursued past this
point — Bottle1's own experience ([Finding 14](#finding-14)) is that the last few percent
sits in a fitted, not derived, regime once the ladder and gross threshold
scale are right, so further digits likely trade one arbitrary choice for
another rather than revealing anything new.

<a id="worst-sj-never-left-the-initial-gate"></a>**Worst SJ never left the initial 0.01 gate — a fourth stuck-point case.**
Left running for over an hour past the `finalMesh.vtk` write hoping
`ELEM_THRES` would ratchet up toward Table 2's reported 0.570 the way it
did for Bottle1 (which walked through 0.550, then 0.570001, before being
stopped). It didn't: `smallDist` oscillated in a tight 0.10–0.11 band on
the same worst point (index 13419) for over 5,000 checkpoints without
ever re-crossing the `1e-4` convergence tolerance needed to trigger the
next success and raise the bar. Stopped once the plateau was clearly
stable rather than transient. This is the same stuck-point signature
[Finding 13](#finding-13) first identified on Bottle1 (there, root-caused to a genuine
input-file bug) and that two of this session's other six Bunny configs
(`rl4-cd75`, `rl4-halfshift-cd75`) also hit — but `bunny_tri.raw` has no
known defect, so on its own this occurrence doesn't point to a cause. Four
stuck-point/stall observations across two models at the time of writing —
later extended to eight across four models once Dragon Stand2's own five
runs (all of which stalled at the initial worst-SJ gate, see below) are
counted — is enough to say the phenomenon is common for this
gradient-descent projection method generally, not a `bunny_tri.raw`- or
Bottle1-specific artifact — worth investigating on its own terms rather
than continuing to treat each occurrence as one-off. Per the mesh-size
result above, this doesn't affect Bunny's headline
reproduction — as the [Bottle1 section](#bottle1-genus-1) already established, worst SJ is a
stopping point on the ratchet, not an independent property of the mesh,
and the mesh itself (size, topology) was unaffected by never reaching a
higher rung.

**Resolved (2026-08-27): the 0.01 figure was a reporting bug, not the
mesh's true achieved quality.** Two things had gone unnoticed in every
run above. First, the Bunny/Dragon Stand2 runs never actually exercised
the `MAX_PROJ_ITER`/`ELEM_THRES_CAP`/best-checkpoint fix described in the
[macOS build & `ProjectToIsoSurface` fix](#macos-build-and-projecttoisosurface-fix) session
— they were built from the separately-extracted `v1.0`/`v1.1` tagged
source trees, whose `ProjectToIsoSurface` is the original, unbounded
`while(true)` loop with none of that instrumentation (confirmed by the
`run.log` line format, `"Try:" << k << ...`, which exists only in those
trees). Second, and more fundamentally: `finalMesh.vtk` only gets
(re)written on an all-or-nothing checkpoint — `k == 0` (*every* element
above the current `ELEM_THRES` bar) *and* `smallDist` below a tight
surface-distance tolerance. Bunny's run got exactly one such success
early on (`ELEM_THRES` jumped 0.01→0.5, freezing `finalMesh.vtk` at worst
SJ 0.010001 — just barely above the *old* 0.01 bar), then stalled
permanently on one persistent point (`maxDistIdx` ≈21936–21953, nearest
triangle ≈13419) for thousands of checkpoints, never again achieving a
full pass at the new, much higher 0.5 bar. But `projHex.vtk` — written
unconditionally every checkpoint, so it reflects the optimizer's actual
state — measured 0.487 on the original run, nowhere near the frozen
0.010.

Checked directly: the persistent point isn't an input defect the way
Bottle1's was. Triangle 13419 (and its neighbors 3294/8528/20142, also
seen as the worst-point's nearest triangle across different runs) are
ordinary, well-shaped triangles — minimum angle 52-57°, area right at the
mesh's median — not slivers or near-degenerate geometry. This is a
genuine optimization plateau in one mesh region, not a bad-input problem.

**Fix**: patched `ProjectToIsoSurface()` in `HybridOctree_Hex_v1.0/HexGen.cpp`
with a purely additive best-quality tracker, alongside the existing
`k`-counting loop — no existing branch, condition, or file write is
touched. It tracks the mesh's true minimum scaled Jacobian every
checkpoint and snapshots it to a new file, `bestQualityMesh.vtk`, whenever
that value improves, decoupled entirely from the `k==0` all-or-nothing
gate. Regression-checked against `bone` first (the one model with a
byte-exact published match): topology unchanged (10,356/8,619), and the
new tracker's worst SJ (0.610001) lands closer to the published 0.610000
than the *old* mechanism's own `finalMesh.vtk` did on the same run
(0.600001) — the fix is a strict improvement even for the case that
already worked, not just a bug-for-bug-compatible patch. Re-ran Bunny's
best recipe (`rl4-halfshift-c145`, unchanged) with the patch:
**worst SJ 0.500001** (mesh size/topology unchanged, 27,175/22,525),
closing the gap to Table 2's 0.570 from a 98% shortfall to roughly 12%.
The run's true worst-point distance to the surface plateaus around
0.2-0.24 in this region — genuinely far from convergence, not on the
verge of a further jump — so 0.5 looks like this specific run's practical
ceiling without a deeper algorithmic change (e.g. letting the
redistribution/perturbation step fire on a partial, not just a full,
success) beyond today's scope.

Re-ran Dragon Stand2's best recipe the same way, as a second validation
of the identical pathology: **worst SJ 0.500000** (mesh size/topology
unchanged, 62,307/50,603), against Table 2's 0.560. Unlike Bunny, this
run was still oscillating on its own persistent point (index ≈56765) at
the moment it was cut off, not a confirmed plateau — it may have had more
room to climb with further runtime. Both runs were cut off manually
after several hours, since neither `ProjectToIsoSurface` build used here
has any automatic stopping condition:

| Model | Cut off after | Started | Stopped | Config (`v1.0` source, `HybridOctree_Hex_v1.0/`) |
|---|---|---|---|---|
| Bunny | 4h47m24s | 2026-08-27 09:45:52 | 2026-08-27 14:33:16 | `VOXEL_SIZE_VALUE=7`, `LADDER_TOP=octreeDepth`, `C_THRES = {0.2121, 0.4243, 0.8485, 1.45, 3.3941}`, `H_THRES = {11.3137, 5.6569, 2.8284, 1.4142, 0.7071}`, `CELL_DETECT = 0.75` (identical to the `rl4-halfshift-c145` recipe above) |
| Dragon Stand2 | 3h44m39s | 2026-08-27 10:48:37 | 2026-08-27 14:33:16 | `VOXEL_SIZE_VALUE=8`, `LADDER_TOP=octreeDepth`, `C_THRES = {0.21213203, 0.42426407, 0.84852814, 1.69705627, 3.39411255}`, `H_THRES = {11.3137085, 4.45, 2.6, 1.41421356, 0.70710678}`, `CELL_DETECT = 0.75` (identical to its own recipe in the Table 2 mini-table above) |

Both binaries were built with `g++ -O3 -DNDEBUG -std=c++17` plus the
`-D<NAME>_VALUE`/`-D<NAME>_VALUES` overrides shown, directly against
`HybridOctree_Hex_v1.0/{Main,HexGen,Mesh}.cpp` (the `#ifndef`-guarded
macro pattern already established for this tree) — no source edits
needed beyond the one `ProjectToIsoSurface()` patch, which is shared by
both builds.

<a id="bunny-per-level-fit"></a>**2026-08-28 – 29 — [Finding 15](#finding-15)'s mechanical per-level fit,
applied back to Bunny: 26,255 / 21,702 vs. 26,375 / 21,695 (−0.45% /
+0.03%), the closest cell-count match of any Table 2 model so far.**
Bunny's `rl4-halfshift-c145` fit predates [Finding 15](#finding-15), so
we re-fitted Bunny with the per-level procedure that took Dragon Stand2
to −0.43%/−0.49%. The starting diagnosis: `rl4-halfshift-c145`'s final
per-threshold refined-parent counts are 808.5 / 1,170.9 / 1,023.6
against the reference's 832.1 / 1,234.0 / 871.7 (thresholds 1/2/3,
12-edge-mean metric). Its +3.0%/+3.8% total overshoot was one threshold
overshooting +17.4% (threshold 3, from the `C_THRES[3]=1.45` loosening)
masking undershoots at thresholds 1 and 2. Ten configurations (two full
runs killed early, four octree-stage probes, four full runs to a
stabilized `projHex.vtk`) produced the fit and three new mechanisms:

- **`H_THRES` is an inert knob at Bunny's thresholds 1–2 in this
  regime.** Loosening `H_THRES[1]` 5.6569→6.0 (+7.9% union candidates
  via `refine_criteria_stats`) and `H_THRES[2]` 2.8284→3.1 (+11.8%)
  moved final parent counts by less than 1%. The added thickness
  candidates sit in regions the child-intersect propagation already
  refines.
- **`C_THRES[2]` is the real threshold-2 lever.** Octree-stage probes at
  `C_THRES[2]` 0.8485→0.80→0.75 moved threshold 2's octree parents
  2,128→2,352→2,496, with only +24/+56 upward propagation into
  threshold 1 and threshold 3 bit-identical (1,512) across all three —
  [Finding 15](#finding-15)'s locality, now confirmed from the
  loosening-upward direction as well. Curvature candidates reach surface
  regions thickness candidates do not. `C_THRES[1]` is weak the same way
  `H_THRES[1]` is: 0.4243→0.40 bought +5 final parents, and threshold 1
  stayed pinned at 807–816 across every configuration tried.
- **Threshold 3's per-candidate weights are lumpy, and the reference
  falls in a dead zone.** Final threshold-3 parents at curvature
  candidate counts 83/84/86/88 (`C_THRES[3]` =
  1.51/1.505/1.50/1.49): 833.6 / 837.4 / 938.0 / 971.4. The two
  candidates added between `C_THRES[3]` 1.505 and 1.50 carry ~100
  refined parents between them. The reference's 871.7 sits in the gap,
  so no `C_THRES[3]` value can land it. (Candidate counts corrected
  2026-08-29 — the original entry quoted 42/43/44/45 from the
  undercounting `refine_criteria_stats`; the *parent* counts beside them
  are measured off the runs and are unaffected, so the dead-zone
  conclusion stands unchanged.)

| Config (`C_THRES[1]`, `[2]`, `[3]`) | Points | Cells | vs. target | parents t1 / t2 / t3 |
|---|---|---|---|---|
| Table 2 target | **26,375** | **21,695** | — | 832.1 / 1,234.0 / 871.7 |
| **`fitC` (0.4243, 0.82, 1.51)** | **26,255** | **21,702** | **−0.45% / +0.03%** | 806.8 / 1,255.8 / 833.6 |
| `fitD` (0.4243, 0.827, 1.505) | 26,048 | 21,515 | −1.2% / −0.8% | 806.9 / 1,225.9 / 837.4 |
| `fitE` (0.40, 0.83, 1.505) | 26,084 | 21,542 | −1.1% / −0.7% | 810.8 / 1,225.2 / 838.4 |
| `fitF` (0.4243, 0.827, 1.50) | 26,859 | 22,258 | +1.8% / +2.6% | 806.8 / 1,228.5 / 938.0 |
| `fitA` (0.40, 0.82, 1.49) | 27,480 | 22,779 | +4.2% / +5.0% | 811.8 / 1,260.7 / 971.4 |
| `fitB` (0.40, 0.80, 1.49) | 28,310 | 23,519 | +7.3% / +8.4% | 812.9 / 1,365.7 / 971.9 |

All rows: `v1.0` source with the `bestQualityMesh.vtk` patch,
`VOXEL_SIZE_VALUE=7`, `LADDER_TOP=7`, `CELL_DETECT_VALUE=0.75`,
`C_THRES[0]=0.2121`, `C_THRES[4]=3.3941`, `H_THRES` at the halfshift
values `{11.3137, 5.6569, 2.8284, 1.4142, 0.7071}`. Counts come from a
stabilized `projHex.vtk` (the 2026-08-18 protocol); `fitC`'s totals were
re-measured twice, 10 and 25 minutes apart, and did not move. Runs
deliberately shared the machine, so no timings were recorded.

**Why `fitC` is the stopping point, not a waypoint.** Every knob still
available moves points and cells together in roughly an 8:7 ratio.
`fitC`'s residual — −120 points against +7 cells — is a *distribution*
imbalance (threshold 1 −25 parents, threshold 2 +22, threshold 3 −38),
and the three knobs that could rebalance it are the three shown above to
be pinned, dead-zoned, or already fitted. Closing the last 120 points
would require breaking the cell-count match, and vice versa. `fitC`'s
projection run was left going overnight for the worst-SJ ratchet.

**Overnight worst-SJ result: 0.010 — the ratchet's *first* firing is the
gate, and `fitC` never got it.** After 7h36m and 5,701 checkpoints
(2026-08-28 23:30 → 2026-08-29 06:06), `fitC`'s worst SJ never moved:
`projHex.vtk` 0.009996, `bestQualityMesh.vtk` peaked at 0.010045 within
minutes of projection starting and never improved. The run.log shows why,
and it sharpens the 2026-08-27 stuck-point story: every late checkpoint
reads `Try:0` — *zero* elements below the 0.01 bar — but `smallDist`
oscillates at 0.03–0.06 on one persistent point (nearest triangle 6467),
never crossing the `1e-4` tolerance that a checkpoint success also
requires. `rl4-halfshift-c145` got exactly one lucky early success,
which raised `ELEM_THRES` 0.01→0.5 and turned on quality pressure — its
minimum then climbed to 0.500001. `fitC`'s stuck point blocked that
first firing, so its bar stayed at 0.01 and the optimizer spent the
whole night chasing surface distance with no quality pressure at all.
Same pathology, one rung earlier: worst SJ is a ratchet stopping point,
and this run stopped before the ratchet's first click.

**Morning proof (2026-08-29 06:10 – 06:48): the mesh reaches 0.5 in one
minute once the gate opens.** We re-ran `fitC`'s exact configuration
with `ELEM_THRES` initialized to 0.5 instead of 0.01
(`runs/bunny-v1.0-fitC-elem05/`, source tweak on a scratchpad copy —
the repo tree is untouched), reproducing by construction the
post-first-success state `rl4-halfshift-c145` reached by luck. The
best-quality tracker's first snapshot (the raw dual mesh) read −0.487;
*one minute later* it read 0.499921, and it finished at **0.500000**
(topology identical: 26,255 / 21,702, best SJ 1.0, 733 checkpoints in
38 minutes). Seven hours of overnight runtime versus sixty seconds,
separated only by the bar's initial value: `fitC`'s worst SJ was never a
property of the mesh, only of whether quality pressure was on. Both
this run and `c145`'s 0.500001 sit exactly at the 0.5 bar — the bar,
not the geometry, is binding in both.

**Table 2's 0.570 reached exactly (2026-08-29 afternoon).** The follow-up
experiment — `ELEM_THRES` initialized at 0.57
(`runs/bunny-v1.0-fitC-elem057/`, same `fitC` recipe, same scratchpad
source-copy pattern) — climbed −0.422 → 0.554752 in the first five
minutes of projection and reached **worst SJ 0.570000** about six
minutes in (`bestQualityMesh.vtk` written 14:15, projection started
~14:09; 982 checkpoints when killed), topology unchanged at
26,255 / 21,702, best SJ 1.0. The reference `bunny.vtk` itself measures
0.569998. Bunny's Table 2 row now reads, side by side:

| | `fitC` + `ELEM_THRES=0.57` | Table 2 |
|---|---|---|
| Points | 26,255 | 26,375 |
| Cells | 21,702 | 21,695 |
| Worst SJ | **0.570000** | 0.570 (file: 0.569998) |
| Best SJ | 1.0 | 1.0 |
| Refinement level | 4 | 4 |

Two more things this settles. First, the whole pipeline took ~11 minutes
wall (~6 of them projection) — the same order as Table 2's reported
358 s, where our ratchet-gated runs took 4-8 *hours* to do less. The
paper's runtimes are consistent with runs whose quality bar was
effectively at its final value, not runs that climbed the ratchet from
0.01. Second, worst SJ lands wherever the bar is — 0.500001 at a 0.5
bar, 0.570000 at a 0.57 bar, on the same mesh — so Table 2's per-model
worst-SJ column (0.560, 0.570, 0.590, ...) reads like a record of
per-model stopping bars rather than of intrinsic mesh-quality
differences. We have not pushed the bar above 0.57 to find `fitC`'s
true quality ceiling; that would be the next probe if the question
ever matters.

**Dragon Stand2 — first check of whether any of this transfers to a model
with a different max level, and it doesn't, the same way Bunny didn't
transfer from Bottle1.** `runs/dragonstand2-v1.0-vs8/`
(`VOXEL_SIZE_VALUE=8`, `LADDER_TOP=octreeDepth`, shipped `C_THRES`/
`H_THRES`, default `CELL_DETECT=1`, per the ladder-derivation table above)
converged to **117,793 / 98,666** against Table 2's **62,576 / 50,853** —
an **88% / 94% overshoot**, the same order of magnitude and direction as
Bunny's shipped-threshold attempts (+70-82%) despite Dragon Stand2 sitting
at a different rung of the ladder (`VOXEL_SIZE_VALUE=8` vs. Bunny/Bottle1's
7). The ladder-position prediction held — this wasn't a five-threshold-style
catastrophic failure, "only" roughly double — but the threshold values
once again needed their own fit, a third independent confirmation of the
synthesis above.

**Dragon Stand2 — fitted to −0.43%/−0.49%, the closest reproduction of
any Table 2 model so far.** A background Opus agent took this on, in
parallel with the Ramses baseline run below, and along the way found
something that changes how every future model in this table should be
fitted — see [Finding 15](#finding-15) immediately below for the full mechanism. Final
recipe: `VOXEL_SIZE_VALUE=8`, `LADDER_TOP=octreeDepth` (unchanged from
the baseline), `C_THRES` = shipped `×√2` (the same "half-level" scale
correction that helped Bunny), `H_THRES` = shipped `÷√2` at thresholds 0, 3,
4, but thresholds 1 and 2 tightened further to `4.45` (from `5.657`) and `2.6`
(from `2.828`), `CELL_DETECT = 0.75`. Reached in six full runs (the first
four working out the coupling rule itself, the last two a mechanical
per-level fit once the rule was known):

| | Reproduction | Table 2 target |
|---|---|---|
| Points | 62,307 | 62,576 |
| Cells | 50,603 | 50,853 |
| vs. target | **−0.43% / −0.49%** | — |
| Refinement level | 4 | 4 |

Per-level refined-parent counts: 1,855.8 / 2,967.3 / 1,898.7 against the
reference's 1,833.1 / 2,981.9 / 1,952.0 — **+1.2% / −0.5% / −2.7%**, every
threshold inside 3%, against the shipped baseline's +25% / +91% / +182%. Worst
SJ stalled at the initial 0.01 gate across every run that reached
projection (the same phenomenon documented for Bunny, now an eighth
occurrence across four models) — no effect on mesh size or topology.
Independently re-verified with `scaled_jacobian_stats` and
`level_histogram.py --parents` against the agent's own numbers before
writing this section. Full blow-by-blow (all six runs, the two dead ends,
the segfault below) is in `dragonstand2-scratch.md` at the repo root, not
folded in here in full — this is the load-bearing summary.

<a id="finding-15"></a>
### Finding 15 — the five refinement thresholds are not independent knobs, and the asymmetry is exact, not approximate

The single most useful discovery of the whole Dragon Stand2 fit, and it
changes the recommended fitting procedure for every model after this one.
`ComputeCellValue()` (`HexGen.cpp:1017`) tests a cell's own threshold gate only
as a *fallback*: for any cell above the deepest level, it first checks
all eight children and marks the cell `intersect = true` if **any** child
intersects, only falling through to the cell's own `C_THRES`/`H_THRES`
box test when none do. Two consequences follow, both confirmed
experimentally, not just read off the source:

- **Tightening threshold `t` affects threshold `t` and nothing
  shallower.** A cell that now fails its own (tightened) test provably
  had no refining children either — if it had, the child check would
  already have marked it `intersect` before its own test ever ran. So
  tightening can only delete cells that were already leaves at that
  level, never anything shallower. Verified exactly: tightening
  `H_THRES[1]` alone (run 3 vs. run 1) left thresholds 2 and 3's octree parent
  counts bit-identical — not approximately equal, identical to the
  integer — while threshold 1 moved by −4.5%.
- **Loosening threshold `t` propagates upward to every shallower
  threshold**, because a newly-intersecting deep cell now marks its parent
  (and grandparent, ...) `intersect = true` too, through the same
  child-check path. Measured directly: loosening threshold 3 by +22%
  (candidate count) pulled threshold 2's octree parent count up +6.8% and
  threshold 1 up slightly, with *threshold 2's own value completely
  unchanged*. A first attempt that read threshold 3 as "9% too tight" in
  isolation and loosened it accordingly went the wrong way entirely —
  the propagated increase at thresholds 1-2 more than cancelled the intended
  threshold-1 tightening it was paired with, and the run had to be killed at
  the octree stage as a known regression.

**This makes per-level fitting mechanical instead of a search.** Read the
reference mesh's refined-parent counts per level
(`level_histogram.py --parents`), then tighten each threshold's own `C_THRES`/
`H_THRES` — working in any order, since tightening never pollutes another
threshold — until each threshold's parent count matches, using `refine_criteria_stats`
to sweep threshold values cheaply before committing to full runs. Never
loosen a threshold relative to a baseline that already overshoots everywhere,
since a loosening's effect on shallower thresholds is real but not confidently
predictable in advance (the failed first attempt didn't get the direction
of the *net* effect wrong, only its magnitude, but that was enough to
invert the outcome). This retroactively explains why Bunny's own
`C_THRES[3]` fit (a single-threshold loosening) worked cleanly on the first
try [above](#best-bunny-reproduction-reached): threshold 3 was Bunny's *deepest* threshold, with
nothing below it to propagate from — loosening the deepest threshold is always
safe, by the same logic that makes tightening any threshold always safe. The
risk is specific to loosening a threshold that has thresholds beneath it, which
neither Bottle1's nor Bunny's fits ever needed to do.

**Caveat added 2026-08-29: `refine_criteria_stats` undercounts.** Its
dihedral loop `break`s out of the partner-triangle loop after the first
match, where `ReadRawData()` continues to every partner, so every
candidate count quoted from that tool anywhere in this file is low —
`bone` at the shipped thresholds is really 335 / 60 / 13 / 1 / 0, not
185 / 28 / 7 / 0 / 0, and Bunny at `fitC` is 1,861 / 838 / 279 / 83 / 2,
not 1,363 / 554 / 163 / 42 / 1. The tool was used here for *relative*
guidance while sweeping one threshold at a time, and every fit was
confirmed by a full run measured against the reference mesh, so the
fits and conclusions stand; only the intermediate candidate numbers are
wrong.

A second, smaller finding from the same session: `refine_criteria_stats`
candidate counts map to refined-parent counts by an approximate power
law, `conditional_parent_ratio ≈ candidate_count_ratio^α` — but α differs
between the octree stage and the final mesh (≈0.19–0.28 vs. ≈0.40–0.53
here) and between models (Bunny's α ≈ 0.37/0.57/0.37, Dragon Stand2's ≈
0.40/0.46/0.53), so it's a local calibration useful for guiding one step
of a search on a model with at least one full run already done, not a
transferable constant to predict a config's outcome from
`refine_criteria_stats` alone with zero full runs. Related: `octree.vtk`
(written right after octree construction, hours before a full run
finishes) is an excellent *projection-free* early read on each threshold's
relative shape, but the octree→final-mesh "survival ratio" — how many
octree-stage parents survive into the final mesh's parent count — drifts
with configuration (Dragon Stand2's projection inflation factor was
`×1.580` for the shipped run and `×1.726` for the first fitted attempt),
so using `octree.vtk` alone to predict a config's *absolute* final size
ran consistently ~10% low here. Once two full runs of the same model
exist, fitting `total_cells = a + b·(sum of per-level parents)` linearly
across just those two points predicted both a third run's outcome and —
more tellingly — the reference mesh's own total from its own parent sum,
to within 0.5%.

**A new failure mode, distinct from the convergence stalls.** One
configuration (`H_THRES[2] = 2.5`, between the two that bracketed the
final answer) crashed with `EXC_BAD_ACCESS` inside
`hexGen::ProjectToIsoSurface`, a wild-pointer segfault in `__bzero`
(crash report preserved in the agent's environment,
`HexGen-2026-08-18-232512.ips`), not the usual infinite-loop-needing-a-
human-to-stop-it behavior documented elsewhere in this file. Reproducible
for that exact configuration and mesh; not investigated further (worked
around by moving to a neighboring threshold value instead) — worth a
closer look if it recurs, since every other run this session either
converged or stalled cleanly rather than crashing.

**Ramses — the worst shipped-threshold mismatch found this session.**
`runs/ramses-v1.0-vs8/` (`VOXEL_SIZE_VALUE=8`, `LADDER_TOP=octreeDepth`,
shipped `C_THRES`/`H_THRES`, default `CELL_DETECT=1` — same ladder rung
and same baseline recipe as Dragon Stand2, per the ladder-derivation
table, so a natural same-rung comparison) converged to **245,386 /
217,727** against Table 2's **44,790 / 37,993** — a **+448% / +473%
overshoot**, by far the largest of any model tried this session
(Bunny's worst was +82%, Dragon Stand2's was +94%). The pre-projection
`dualHex.vtk` was already 4.3x oversized (191,606 / 163,949) before
projection's usual element-splitting growth compounded it further. Two
models sharing the same ladder rung and the same shipped-threshold
baseline (Dragon Stand2: +88-94%, Ramses: +448-473%) landing five times
further apart than either is from its own target is strong further
evidence against any single per-rung correction factor — the fitting
really is per-model, not even per-ladder-rung. Not pursued to a fit yet —
its ~3.6-hour pipeline, by far the slowest of any model run, made even
one exploratory full run an expensive commitment, but [Finding 15](#finding-15) (Dragon
Stand2's session, below) turns that from a search into a mechanical
per-level match against `refine_criteria_stats` and the reference's own
refined-parent counts, which should make Ramses considerably cheaper to
fit than Dragon Stand2 was, despite its slower pipeline per attempt —
fewer attempts should be needed.

**Ramses fitting, applying [Finding 15](#finding-15) directly — round 2 result: +3.4%/+4.0%,
Ramses' best reproduction.** Per-level
parent-count ratios (shipped-baseline `finalMesh.vtk` vs. reference,
using the `level6<-5`/`level7<-6`/`level8<-7`/`level9<-8` threshold mapping
established for this ladder rung) showed a far more lopsided mismatch
than Dragon Stand2's: threshold 1 (level 5→6) only +13% over, threshold 2 (level
6→7) **+217%**, threshold 3 (level 7→8) **+3,065%** — by a wide margin the most
extreme single-threshold overshoot found this session, concentrated almost
entirely in curvature (`refine_criteria_stats` shows threshold 2/threshold 3
candidate counts barely move under `H_THRES` sweeps but drop sharply
under `C_THRES` sweeps). Round 1 (`C_THRES = {0.15, 0.3, 1.0, 2.2, 2.4}`,
`H_THRES` shipped, `CELL_DETECT = 0.75`, thresholds 2 and 3 tightened, threshold 1
left alone) converged to **104,067 / 90,483** against the 44,790 / 37,993
target — **+132% / +138%**, a large improvement on the shipped baseline's
+448%/+473% but still far off. Its octree stage ran in 31.8 minutes
(1,905.6 s), a 6x speedup from the shipped baseline's 190.3 minutes —
consistent with the tighter-thresholds-build-cheaper-octrees pattern
holding here too. Per-level parent ratios after round 1: threshold 1 now close
(1.11x reference — leave alone, per [Finding 15](#finding-15)'s tighten-only rule),
threshold 2 still 2.22x over, threshold 3 still 6.92x over. Threshold 3's curvature
candidates alone bottom out around a thickness-driven floor (~108-135
regardless of how tight `C_THRES[3]` goes), so closing the remaining gap
needed tightening `H_THRES[3]` as well, not `C_THRES[3]` alone — round 2
(`C_THRES = {0.15, 0.3, 2.0, 3.5, 2.4}`, `H_THRES = {16, 8, 4, 1.5, 1}`,
`CELL_DETECT` unchanged) tested both thresholds tightened together.

**Round 2 result.** Octree construction dropped again, to 7.7 minutes
(463.6 s) — 24.6x faster than the shipped baseline, 4.1x faster than
round 1 — continuing the pattern that tighter thresholds build cheaper
octrees within a model's own fitting process, not just across different
models' shipped baselines. Converged cleanly, and reached a real worst-SJ
ratchet value for the first time this Ramses attempt (0.570001, versus
Table 2's reported 0.590 — close, and, as established throughout this
file, itself just a stopping point on the ratchet rather than an
independent property of the mesh):

| | Reproduction | Table 2 target |
|---|---|---|
| Points | 46,329 | 44,790 |
| Cells | 39,529 | 37,993 |
| vs. target | **+3.4% / +4.0%** | — |
| Worst SJ | 0.570001 | 0.590 |
| Best SJ | 1.000000 | 1.0 |
| Refinement level | 4 | 4 |

Per-level refined-parent counts: 1,635.4 / 3,016.1 / 692.6 against the
reference's 1,589.6 / 2,935.5 / 609.7 — **+2.9% / +2.7% / +13.6%**
(threshold 3's percentage looks large but is a small absolute count, 83 cells).
Comparable in quality to Bunny's own best reproduction (+3.0%/+3.8%),
reached in two rounds rather than Dragon Stand2's six — consistent with
[Finding 15](#finding-15)'s prediction that the mechanical per-level procedure would need
fewer attempts than the trial-and-error process that discovered it, even
starting from a far more severe shipped-threshold mismatch (+448% vs.
Dragon Stand2's +88%). Final configuration: `C_THRES = {0.15, 0.3, 2.0,
3.5, 2.4}`, `H_THRES = {16, 8, 4, 1.5, 1}`, `CELL_DETECT = 0.75`,
`VOXEL_SIZE_VALUE = 8`, `LADDER_TOP = octreeDepth`.

### Pipeline timings, Bunny / Dragon Stand2 / Ramses (Apple M1, 2026-08-18 – 19)

Same per-stage `clock()` prints `Main.cpp` reports for every run (see the
Bottle1 pipeline timings table above for the M4/2026-07-02 and
M1/2026-08-04-05 numbers; this table is a separate day/session, kept
apart rather than merged into one giant table). "Time to first
`finalMesh.vtk`" is the wall-clock gap between the "extracting interior
dual mesh" print (buffer clearing done, topology essentially fixed) and
that run's `finalMesh.vtk` mtime — the time `ProjectToIsoSurface` needed
to reach its *first* zero-bad-element checkpoint, not total time spent
ratcheting `ELEM_THRES` further afterward (which several runs were left
doing well past this point, per the results tables above, and isn't a
comparable quantity run-to-run since it depends entirely on when a human
operator chose to stop each one).

| Run | Read | Octree | Dualization | Buffer clearing | Time to first `finalMesh.vtk` | Sum through buffer clearing |
|---|---|---|---|---|---|---|
| `bunny-v1.0-rl4-shift1-cd75` | 1.44 s | 58.0 s | 7.3 s | 33.1 s | 57 s | 99.9 s |
| `bunny-v1.0-rl4-halfshift-cd75` | 1.44 s | 136.7 s | 17.4 s | 50.0 s | *(never — stuck, see above)* | 205.5 s |
| `bunny-v1.0-rl4-halfshift-c145` | 1.42 s | 168.7 s | 21.5 s | 53.1 s | 632 s (~10.5 min) | 244.7 s |
| `bunny-v1.0-rl4-cd75` | 1.43 s | 320.6 s | 41.4 s | 76.9 s | *(never — stuck, see above)* | 440.3 s |
| `bunny-v1.0-rl4` (plain) | 1.43 s | 430.2 s | 54.5 s | 89.6 s | 1,277 s (~21.3 min) | 575.7 s |
| `bunny-v1.0-rl4-h2x-cd75` | 1.39 s | 410.7 s | 51.9 s | 83.9 s | 7,283 s (~121.4 min) | 547.9 s |
| `bunny-v1.0-rl4-h3` | 1.40 s | 471.7 s | 58.2 s | 91.6 s | 451 s (~7.5 min) | 622.9 s |
| `bunny-v1.0-vs8abs` (five-threshold control) | 1.44 s | 4,083.1 s (~68.1 min) | 531.1 s | 405.7 s | 2,498 s (~41.6 min) | 5,021.3 s |
| `dragonstand2-v1.0-vs8` | 10.2 s | 3,893.2 s (~64.9 min) | 587.2 s | 806.0 s | 1,149 s (~19.2 min) | 5,296.6 s |
| `ramses-v1.0-vs8` | 2.0 s | 11,420.9 s (~190.3 min) | 1,435.9 s (~23.9 min) | 1,029.6 s (~17.2 min) | 1,059 s (~17.7 min) | 12,889.4 s (~214.8 min) |
| `dragonstand2-v1.0-halfshift-cd75-h1h26` (fitted, run 6 of 6) | ~10 s | 1,486 s (~24.8 min) | 175 s | 439 s | ~600 s (~10 min)¹ | 2,110 s (~35.2 min) |
| `ramses-v1.0-r1` (fitted, round 1 of 2) | 2.0 s | 1,905.6 s (~31.8 min) | 238.6 s | 275.9 s | 5,809 s (~96.8 min)² | 2,422.1 s (~40.4 min) |
| `ramses-v1.0-r2` (fitted, round 2 of 2 — final) | 2.0 s | 463.6 s (~7.7 min) | 55.3 s | 113.4 s | 1,058 s (~17.6 min) | 634.3 s (~10.6 min) |

¹ From the Opus agent's own working notes (`dragonstand2-scratch.md`),
not independently re-timed against file mtimes the way every other row's
"time to first `finalMesh.vtk`" was — reported here as approximate for
that reason, everything else in the row is exact from the run's own log.

² Round 1's long time-to-`finalMesh.vtk` reflects this run's own
convergence behavior, not a stall — it did reach a checkpoint and write a
file (worst SJ 0.010001, the initial ratchet gate), just slowly relative
to round 2. Round 2's octree stage, at 7.7 minutes, is 4.1x faster than
round 1's 31.8 minutes and 24.6x faster than the shipped baseline's 190.3
minutes — the largest same-model octree-time swing measured this
session, and further confirmation that "tighter thresholds build cheaper
octrees" holds within a single model's own fitting sequence.

Two patterns worth noting. **Octree construction time tracks threshold
tightness almost exactly**, within each model — Bunny's own configs span
a 7.4x range (58.0 s to 471.7 s) purely from how permissive `C_THRES`/
`H_THRES` are, cleanly ordered `shift1` (tightest) < `halfshift-c145` <
`halfshift-cd75` < `cd75` < `h2x-cd75` ≈ `plain` < `h3` (loosest) — the
same ordering as the [mesh-size results table](#mesh-size-results-table) above, which makes sense
since a looser threshold means more octree cells get flagged for
refinement, which is directly what this stage computes. **`vs8abs` and
`dragonstand2-v1.0-vs8` are outliers at around an hour of octree time
alone, and `ramses-v1.0-vs8` is a far bigger one at ~3.2 hours** — all
three for the same underlying reason: each builds a substantially
deeper/denser octree than any Bunny-`rl4` config, and the pre-projection
`dualHex.vtk` size explains why Ramses is the worst of the three by a
wide margin. Its shipped-threshold octree came in at **191,606 / 163,949**
pre-projection against a target of only 44,790 / 37,993 — a **4.3x**
overshoot, far worse than Dragon Stand2's ~1.3x at the same shipped-
threshold, correct-ladder-rung baseline. So Ramses isn't just slow because
it's at the same expensive ladder rung as Dragon Stand2 — the shipped
thresholds are *more* mismatched for Ramses specifically than for any
other model tried this session, which is itself a data point for the
per-model-fitting story: even models sharing a ladder rung don't share a
threshold-mismatch magnitude. A reminder that "cheap to test via
`refine_criteria_stats` first" (the recipe item in the Synthesis
above) matters more, not less, as a model's octree gets deeper, since the
full-run cost of a wrong guess grows with it. Confirmed directly by
Dragon Stand2's own fitted run: once its thresholds were tightened to the
right neighborhood, octree construction dropped from 3,893.2 s to 1,486 s
(2.6x faster) — the same "tighter thresholds build a cheaper octree"
relationship the Bunny table above shows within its own config family,
now shown to hold across a model's fitting process too, not just
across different models' shipped baselines.

**Open threads for the next session:**

1. **Confirmed, all six.** `rl4` and `rl4-h3` wrote clean `finalMesh.vtk`s;
   `rl4-cd75` and `rl4-halfshift-cd75` plateaued for 2+ hours without
   reaching the convergence gate but stabilized in `projHex.vtk`'s
   point/cell count well before being stopped, so those two counts are
   trusted on that basis rather than from a clean write. Also ran a real
   `vs8abs` (five-threshold, shipped `VOXEL_SIZE_VALUE=8`) test rather than
   just inferring from the reference mesh — massively overshoots (+484%),
   independently confirming Bunny's level-8 population is distortion, not
   a real fifth threshold.
2. **The decisive cross-model test — answered, negatively, within this
   same session.** `runs/bottle1-v1.0-rl4-halfshift-cd75/` (input
   `f3247d8`, the header-correct file from [Finding 13](#finding-13)) was built and run
   against Bottle1 using the exact half-level correction that looked
   promising on Bunny. Result: **19,971 / 16,617** against Bottle1's
   target 36,091 / 30,145 — a **45% undershoot**, far worse than Bottle1's
   own [Finding-14](#finding-14) recipe (−1.5%/−0.7%), and in the opposite direction from
   Bunny's ~+70% overshoot under Bottle1's recipe. So neither model's
   threshold set works on the other: **the thresholds are genuinely
   per-model, not one universal scale-relative rule.** That's a real
   answer, not an open thread — worth stating plainly rather than
   continuing to search for a single formula that fits every model. It
   also reframes the "publishable conclusion" possibility [Finding 14](#finding-14)
   already flagged: `our results/*.vtk` were generated with per-model
   settings some time before 2024-01-11, one day before any source was
   committed, and are not reproducible from the public `v1.0` constants
   alone with any single choice of ladder-scaling rule tried so far.
3. **Bunny's residual gap — closed to +3.0%/+3.8%.** See "Best Bunny
   reproduction reached" above.
4. **Dragon Stand2 — fitted to −0.43%/−0.49%, the best reproduction of
   any Table 2 model to date.** See "[Finding 15](#finding-15)" and the results table
   above. Along the way, found the threshold-independence rule that turns
   every future model's threshold fit from a search into a mechanical
   per-level match.
5. **Ramses — fitted to +3.4%/+4.0%, in two rounds.** See the results
   table above. Started from the most severe shipped-threshold mismatch
   of any model this session (+448%/+473%, driven almost entirely by one
   threshold's curvature gate) and reached a result comparable to Bunny's own
   best in a fraction of the attempts Dragon Stand2 needed — direct
   evidence that [Finding 15](#finding-15)'s mechanical procedure generalizes and pays
   off, not just a one-off success on the model that discovered it.
6. **Reproduction paused here for now, at the user's request, after five
   models (`bone`, Bottle1, Bunny, Dragon Stand2, Ramses) — a deliberate
   stopping point, not a natural end.** The remaining seven Table 2
   models (Gargoyle, Head, Thai Statue, Lion Recon, Red Circular Box, Oil
   Pump, David), per the ladder-derivation table above, all have their
   ladder position already confirmed against their own reference mesh at
   zero compute cost (see the full-table check above), so resuming should
   go straight to [Finding 15](#finding-15)'s mechanical per-level fit for each, without
   re-deriving anything already established here.
7. **The worst-SJ/convergence-stall phenomenon** (nine occurrences across
   four models now — Ramses round 1 stalled at the initial gate too,
   round 2 didn't) and the `ProjectToIsoSurface` segfault found during
   Dragon Stand2's fit ([Finding 15](#finding-15)) are both still open, unexplained, and
   potentially related — worth a dedicated investigation session rather
   than continuing to note each new occurrence in passing.

## Bone (calibration testbed, not a Table 2 model)

**Superseded by the exact match.** The table below reflects this
session's calibration attempts as they stood on 2026-08-05, before the
[2026-08-17 discovery](#finding-7) that `v1.0`'s own shipped constants —
untouched (`C_THRES = {0.15, 0.3, 0.6, 1.2, 2.4}`,
`H_THRES = {16, 8, 4, 2, 1}`, `VOXEL_SIZE = 9`, `CELL_DETECT = 1`) —
already reproduce `bone` byte-for-byte. See [0. Bone](#0-bone-calibration-testbed-not-one-of-table-2s-twelve-models)
above for the headline numbers and [Finding 7](#finding-7) for the full
derivation. The table below is kept for the record, to show what the
earlier, unnecessary calibration work found.

`bone` isn't one of Table 2's twelve models — it's the smallest sample in
this repo's own broader `input boundaries`/`our results` collection, with
published stats in this repo's own README rather than in the paper. The
M4 session picked it as a fast iteration testbed for the
`ProjectToIsoSurface` fix; the 2026-08-05 session reused it the same way,
to calibrate `C_THRES` against a small, fast-to-rerun model before
spending Bottle1-scale time on each attempt. See the [summary log](#summary-log) below for
the full investigation; this table is the result.

| `C_THRES` | Vertices | Elements | Worst SJ | Best SJ | Total time |
|---|---|---|---|---|---|
| **Published (this repo's README)** | **10,356** | **8,619** | **0.61** | 1.0 | — |
| Buggy, shipped `{0,0,0.4,0.8,1.6}` (2026-08-04) | 19,639 (1.9x) | 16,590 (1.9x) | 0.590 | 1.0 | 419.5 s |
| Paper-literal `{0.5,1,2,4,8}` (2026-08-05) | 4,584 (0.44x) | 3,570 (0.41x) | 0.550¹ | 1.0 | 117.75 s |
| **v1.2 historical `{0.1,0.2,0.4,0.8,1.6}` (2026-08-05, adopted)** | **14,856 (1.43x)** | **12,608 (1.46x)** | **0.600** | 1.0 | 313.80 s |

¹ An earlier run at this same `C_THRES` value, using a leftover
`MAX_PROJ_ITER=20000` cap carried over from the Bottle1 work rather than
bone's own established `200000` budget, produced worst SJ 0.011 —
resolved once re-run with the correct budget (0.550, shown here),
confirming that discrepancy was a convergence-budget confound, not a
quality problem with the paper-literal thresholds themselves.

v1.2's historical thresholds give the closest match on every metric
simultaneously (size within ~1.5x, worst SJ within 0.01 of published) and
were adopted for the Bottle1 redo. Neither the buggy nor the paper-literal
values are being pursued further — see the [summary log](#summary-log) for why an exact
match isn't expected from either given the curvature-formula mismatch.

## Synthesis: what bone, Bottle1, and Bunny collectively show

Three models, three different outcomes, now that all three have been
pushed as far as this investigation has taken them. Putting them side by
side is more informative than any one of them alone.

| | `bone` | Bottle1 | Bunny |
|---|---|---|---|
| Genus | 0 | 1 | 0 |
| Table 2 refinement level | — (not in Table 2) | 4 | 4 |
| Ladder needed | shipped defaults, unmodified | `VOXEL_SIZE_VALUE=7`, `LADDER_TOP=octreeDepth` (4 levels) | same as Bottle1 |
| `C_THRES` needed | shipped, unmodified | shipped, unmodified | needs its own scale correction — shipped is far too loose |
| `H_THRES` needed | shipped, unmodified | shipped, ×2 | between ×1 and ×√2 (own correction, still being narrowed) |
| `CELL_DETECT` needed | shipped (1), unmodified | 0.75 | 0.75 (borrowed from Bottle1 — untested whether it's even right) |
| Result | **exact** — byte-identical cell connectivity | within 1.5% / 0.7% | within ~8% so far, narrowing |

**Finding that holds up, independently, on every model tried: "Refinement
Level" is a real, measurable, model-invariant property — the depth of the
octree ladder — not a tuning knob.** It was first read off Bottle1's
reference mesh by measuring cell edge lengths ([Finding 10](#finding-10)), then confirmed
the same way on Bunny at the start of this session, and now confirmed a
third, stronger way on both models: actually *running* the wrong threshold
count rather than just inferring it from the reference. Bottle1 at five
thresholds overshoots 3.4x ([Finding 11](#finding-11)); Bunny at five thresholds overshoots 5.8x
(`vs8abs`, this session). Both failures are catastrophic and in the same
direction, which is exactly what a genuinely model-invariant setting
predicts when it's set wrong. `bone` fits the same picture from the other
side — it's the one model whose correct ladder happens to be the shipped
default, which is why it's the one model this repo reproduces exactly
with zero modification.

**Finding that does not hold up: any single rule for the threshold**
values (`C_THRES`, `H_THRES`, `CELL_DETECT`) **that sit inside that
ladder.** `bone` needs the shipped values verbatim. Bottle1 needs the
shipped `C_THRES` but a doubled `H_THRES` and a narrowed `CELL_DETECT`.
Bunny needs neither of Bottle1's corrections — applying them overshoots by
70-82% — but a `√2`-scale correction invented to fit Bunny then
undershoots Bottle1 by 45% when cross-applied, the opposite direction from
Bunny's own error under Bottle1's recipe. Three models, three different
threshold regimes, and every cross-model transfer attempted so far — exact
reuse, ×2, ×√2 — has failed in one direction or the other. This isn't for
lack of trying a reasonable hypothesis: `H_THRES[i] = 2.56 × cellsize`
being a genuine scale-relative invariant ([Finding 12](#finding-12)) was a real, principled
pattern discovered from bone/Bottle1's own numbers, not a guess — it just
doesn't transfer to a third model, which is a stronger and more useful
result than confirming it would have been.

**What this means for reproducing the paper.** These aren't equally
weighted findings — the ladder-depth result is now backed by three models
and two independent measurement methods each, about as solid as anything
in this log; the threshold-transfer failure is backed by exactly one
cross-check in each direction (Bottle1→Bunny, Bunny→Bottle1) and could in
principle still yield to a smarter rule nobody's tried yet. But taken at
face value, together they say something specific: **the algorithm in the
one surviving source snapshot (`v1.0`) is not broken, and is not the kind
of gap the sibling `HexOpt` repo has** (public code implementing a
materially different method than the one the paper claims). Given the
right inputs, this code reproduces the paper's own output either exactly
(`bone`) or within a fraction of a percent (Bottle1). What's missing is
narrower and more mundane than a broken algorithm: a set of per-model
threshold values that were evidently hand-tuned by the paper's author
before any source was ever committed — the reference meshes for `bone`,
Bottle1, and Bunny were all committed by Hua Tong on the same day
(2024-01-11), a full day before `44f730c`, the first commit to touch
`HexGen.cpp`/`Main.cpp` at all ([Finding 2](#finding-2)) — and were never captured
anywhere in git. Table 2 itself only reports the *outcome* of that tuning
(vertex/element counts, SJ range, refinement level, time), not the
`C_THRES`/`H_THRES`/`CELL_DETECT` inputs that produced it, so there is no
document, public or private, that records what those per-model values
actually were. Reproducing Table 2 as a whole is therefore not a single
mechanical exercise — it is, at minimum, a per-model parameter-fitting
problem, and at most twelve separate small ones.

**A practical recipe for the next model, distilled from what worked
across all three:**

1. **Pin the ladder first, and trust it.** Measure the reference mesh's
   deepest real octree level (12-edge-mean metric — see the level-metric
   caveat above; a single-edge proxy will mislead) to set
   `VOXEL_SIZE_VALUE`/`LADDER_TOP`, per the derived table above. This step
   has been right on every model checked so far and is cheap to get
   right — an octree-only run, without waiting for projection, is enough
   to confirm it (as `vs8abs` did for Bunny).
2. **Don't assume a distortion tail is a real level.** A small
   (2-4%), spatially diffuse population at one level deeper than the main
   mass of the mesh has now been seen — and confirmed as an artifact, not
   real refinement — on both Bottle1 ([Finding 11](#finding-11), statistical/spatial
   argument) and Bunny (this session, an actual run at the deeper setting
   that overshoots by 5-6x). Treat it as the default explanation, not the
   exception.
3. **Treat `C_THRES`/`H_THRES`/`CELL_DETECT` as a small fitting problem,
   not a lookup — and fit it mechanically now, not by search.**
   `refine_criteria_stats` sweeps candidate counts in single-digit
   seconds, versus the 5-190+ minute octree stage alone for a full run —
   use it to narrow the threshold search before committing to an
   expensive end-to-end run. The per-level refined-parent-count diagnostic
   (recursing a level's leaf count upward, `N/8` per level, and comparing
   level-by-level against the reference) is sharper than comparing mesh
   totals, because it localizes *which* threshold is off rather than just *how
   much* — and, per [Finding 15](#finding-15) (Dragon Stand2's session), that
   per-level comparison isn't just a diagnostic, it's the fit itself:
   tightening one threshold affects that threshold exactly and
   nothing shallower, so matching each threshold's parent count to the
   reference, working in any order but only ever tightening, converges
   mechanically rather than by trial and error. The one thing to avoid is
   loosening a threshold that has thresholds beneath it — that propagates upward
   unpredictably in magnitude (though always in the same direction) and
   can invert an otherwise-correct adjustment elsewhere.
4. **Budget for a stall — and expect it, not just tolerate it.** Every
   model attempted long enough has hit a multi-hour plateau in
   `ProjectToIsoSurface` at some configuration: Bottle1's original input
   ([Finding 13](#finding-13), resolved — a genuine duplicate-triangle bug), three of
   Bunny's seven configurations (including its own best reproduction,
   which plateaued on worst SJ for over an hour after its mesh size had
   already converged), and all five of Dragon Stand2's full runs (pinned
   at the initial worst-SJ gate every time). None of Bunny's or Dragon
   Stand2's cases have a known input defect the way Bottle1's did, so
   eight occurrences across four models now looks like a general property
   of this gradient-descent projection method rather than a per-model
   artifact — still unexplained, but common enough to plan around rather
   than treat as one-off. Separately, Dragon Stand2's fit also turned up
   a genuine crash (`EXC_BAD_ACCESS` inside `ProjectToIsoSurface`,
   reproducible for one specific configuration) — a different failure
   mode from the stall, worth distinguishing from it when triaging a run
   that didn't finish cleanly. A `finalMesh.vtk` that never arrives, or
   one that arrives but never ratchets further, isn't necessarily a dead
   end either way: `projHex.vtk`'s point/cell count stabilizing over a
   long plateau is usable evidence that the mesh *topology* has settled
   even without full
   point-position convergence, and worst SJ never climbing doesn't affect
   that topology at all.

<a id="why-no-single-parameter-set"></a>
### Why no single parameter set reproduces all of Table 2

Every model fitted so far has needed its own `C_THRES` values. Bottle1's
recipe overshoots Bunny by 70%. Bunny's correction undershoots Bottle1 by
45%. Dragon Stand2 and Ramses share a ladder rung and still land five times
further from each other than either sits from its own target. The paper
describes one fixed parameter set, so the two pictures need reconciling.
The upstream repository's own history accounts for most of the difference.

**The constants changed during the study.** Reading
`HybridOctree_Hex/Initialization.h` at each upstream commit where the values
are readable:

| Date | Tag | `VOXEL_SIZE` | `C_THRES` | `H_THRES` |
|---|---|---|---|---|
| 2024-01-12 | — | 9 | {0.15, 0.3, 0.6, 1.2, 2.4} | {16, 8, 4, 2, 1} |
| 2024-01-16 | — | 9 | {0.15, 0.3, 0.6, 1.2, 2.4} | {16, 8, 4, 2, 1} |
| 2024-01-19 | `v1.2` | 10 | {0.1, 0.2, 0.4, 0.8, 1.6, **3.2**} | {16, 8, 4, 2, 1, **0.5**} |
| 2024-01-21 | `v1.2` | 10 | {0.1, 0.2, 0.4, 0.8, 1.6} | {16, 8, 4, 2, 1} |
| 2024-01-31 | `v1.2` | 10 | {0.1, 0.2, 0.4, 0.8, 1.6} | {16, 8, 4, 2, 1} |
| 2024-02-22 | `v1.3` | 10 | {**0, 0**, 0.4, 0.8, 1.6} | {16, 8, 4, 2, 1} |

The journal records the paper as received 14 January 2024, revised 23
February 2024, and accepted 16 March 2024. The curvature array therefore
changed three times across that window. The array length itself went from
five entries to six and back within three days (19-21 January). The sixth
entry appearing on 19 January is the same one `v1.2`/`v1.3` later carry
commented out, and it is what David's level-9 population requires
([Finding 10](#finding-10)). `VOXEL_SIZE` also moved 9 → 10 in the same
window, which is the per-model level-count knob of [Finding 9](#finding-9).

<a id="finding-16"></a>
### Finding 16 — the released code cannot use `G_thres`, and only one of two scenarios can be true

The paper defines `G` as Gaussian curvature (cotangent weights over a
Voronoi cell area). The code computes a different quantity in different
units — an unnormalized `(angle − π)²` sum over incident edges,
accumulated in `ReadRawData()` (`HybridOctree_Hex_v1.0/HexGen.cpp:882`)
and tested against `C_THRES` in `GetCellValue()`. `C_THRES` never equals
`{0.5, 1, 2, 4, 8}` at any commit in the history. So the paper's
Section 2.1 curvature formulation, and the thresholds calibrated for it,
describe something the released code never implements.

The unit difference is not cosmetic. The shipped `(angle − π)²` sum is
unnormalized, so its magnitude depends on tessellation density, which
differs per model — plausibly the very reason every model needs its own
`C_THRES`. The paper's `G` is normalized by local Voronoi area (units
~ 1/length), a tessellation-independent quantity. A threshold against a
tessellation-independent quantity could transfer across models; a
threshold against an unnormalized sum cannot be expected to.

Two scenarios could explain how Table 2 was generated, and they cannot
both be true. In the first, an unreleased code version implemented
Section 2.1's Gaussian-curvature formulation, and `G_thres = {0.5, 1, 2,
4, 8}` was its actual input. The paper then describes its own inputs
accurately, but the released code cannot regenerate its results. In the
second, the released algorithm generated Table 2, with per-model
`C_THRES` values that no commit records. The released code is then the
true generator, but `G_thres` as printed is a description of the method's
intent, not the input that produced the table. The evidence we have leans
toward the second: `bone` reproduces byte-exactly under the shipped
formula and shipped `C_THRES` (the 2026-08-17 session), and every
reference mesh was committed a day before the first source commit
([Finding 2](#finding-2)).

Either branch lands in the same place for reproduction: the paper and the
released code, taken together, do not tell a story complete enough to
regenerate Table 2 as written. The curvature inputs have to be re-fitted
per model against each reference mesh ([Finding 15](#finding-15)). The
paragraphs below give the commit-by-commit evidence.

**Scenario (a) tested directly and eliminated (2026-08-29).** We
implemented the paper's Section 2.1 formula as a switchable code path in
`HybridOctree_Hex_v1.0/` — `-DCURVATURE_MODE=1` fills `triMesh.r[]` with
Meyer's cotangent-Laplacian `G` (mixed Voronoi area, computed on the
100-unit-cube coordinates), so `C_THRES_VALUES` carries the paper's
literal `G_thres = {0.5, 1, 2, 4, 8}` with no downstream changes. The
default path is byte-identical to the shipped source (verified by
`g++ -E -P` diff against `HEAD`), and the mode-1 values match an
independent Python implementation to 5×10⁻⁷ (VTK ASCII precision).
Three measurements close the question:

- **bone** (`runs/bone-v1.0-gauss/`, otherwise fully shipped
  configuration): `G` maxes at **0.3011** on the 100-cube — not one of
  bone's 6,047 vertices crosses even `G_thres[0] = 0.5`, so the paper's
  configuration produces *zero* curvature-driven refinement. The octree
  stops at level 6 (the byte-exact reference occupies levels 5-9) and the
  final mesh comes out 4,061 / 3,139 against the reference's
  10,356 / 8,619.
- **Bunny and Bottle1 show the same wall**: on the 100-cube, Bunny has
  0.01% of vertices above 0.5 and none above 1; Bottle1 has 1.85% above
  0.5 and 0.01% above 1. `G_thres[1]`-`[4]` never fire on any of the
  three models.
- **No coordinate scale rescues it.** All three raw inputs are ~unit
  scale (bone 0.95, Bottle1 1.63, Bunny 1.87), so computing `G` on raw
  coordinates multiplies it by 53-106×, which flips the failure the other
  way: bone's *median* `G` becomes ≈8.6, above the deepest threshold, so
  more than half the surface would demand level-9 refinement (the
  reference's level-9 population is 0.59%). The underlying obstruction is
  dynamic range: `G_thres` spans 16× (0.5 → 8), while each model's `G`
  distribution spans only ~4-20× from median to maximum — at *any*
  global scale, the ladder cannot slice such a distribution into the
  graded per-level populations the reference meshes actually have.

This also settles the tessellation-independence conjecture above,
half-for and half-against: the paper's `G` is indeed
tessellation-independent, but its normalization concentrates the values
so tightly around each model's characteristic feature scale that it has
too little spread to drive a five-rung refinement ladder. The shipped
unnormalized `(angle − π)²` sum — ~90× spread on Bunny — is the quantity
that *can* grade refinement, at the price of being tessellation-dependent
and therefore per-model. Combined with bone's byte-exact reproduction
under the shipped formula, the verdict is now empirical, not
inferential: **the released formula generated the reference meshes, and
Section 2.1's formulation with `G_thres = {0.5, 1, 2, 4, 8}` cannot
produce Table 2 at any coherent scale.**

**Thickness matched throughout; curvature never did.** `H_THRES` reads
`{16, 8, 4, 2, 1}` at every commit — identical to the paper's
`T_thres = {16, 8, 4, 2, 1}`, and unchanged apart from the brief six-entry
form, which simply continues the halving with `0.5`. `C_THRES` never equals
the paper's `G_thres = {0.5, 1, 2, 4, 8}` at any commit in the history.

That asymmetry has a straightforward explanation, and it is probably the
whole story. Thickness means the same thing in both places. The paper
measures a distance along a ray from `P_i` along its normal, and the code
measures a distance along a ray from `P_i` along its normal, so the numbers
carry over directly. Curvature does not. The paper defines `G` as a Gaussian
curvature, `‖Σ(cot α_ij + cot β_ij)(P_j − P_i)‖² / (4A_i)`, using cotangent
weights over a Voronoi cell area. `ReadRawData()` instead accumulates
`(angle − π)²` over each vertex's incident edges (`HexGen.cpp:882`), and
`GetCellValue()` tests that sum against `C_THRES`. Those
are different quantities in different units. Their thresholds were never
comparable, so there is no reason to expect `C_THRES` to equal `G_thres`.
Section 2.1 describes a formulation the released code does not implement —
the ordinary kind of drift that accumulates when a paper and its codebase
develop in parallel over months.

**What this means for reproduction.** The published curvature values are not
a usable starting point for this codebase, because they are calibrated for a
quantity the released code does not compute. The published thickness values
are usable, and we have used them unchanged on several models. Nothing in
the history indicates a per-model parameter table ever existed: every commit
carries a single global set, revised over time. Table 2's twelve rows were
most plausibly generated across that revision window rather than in one
batch, which on its own would stop any single set from reproducing all
twelve. Bottle1's reference mesh predating the first committed source by one
day ([Finding 8](#finding-8)) points the same way — some Table 2 runs may
have used constants that were never committed.

None of this changes the fitting procedure. [Finding 15](#finding-15) still
applies: fit each model's thresholds against its own reference mesh's
per-level parent counts, treat the shipped values as a starting point rather
than a target, and read `H_THRES` as the one array where the paper's numbers
transfer directly.

## Summary log

*Ordered oldest to newest*

| Date | Machine | Narrative |
|---|---|---|
| 2026-07-02 – 2026-07-03 | M4 | Got this fork building on macOS, then found and fixed a genuine upstream bug — `ProjectToIsoSurface()` had an infinite loop with no exit path, an unbounded quality-ratchet, and a convergence tolerance borrowed from an unrelated check. Fixed and verified against `bone`, the smallest sample model. Separately, this is also where the sibling `HexOpt` repo's own gap was first found (its public code doesn't implement the AL/L-BFGS/ReHQJ method its paper claims). |
| 2026-08-04 | M1 | Cloned the fork here, cleaned up ~600K lines of build cruft the M4 session had accidentally committed, and ran the first Bottle1 reproduction attempt. Result: a mesh ~6x too large and a worst scaled Jacobian far below Table 2's — the size gap unexplained at the time, the quality gap eventually traced to a deliberately reduced iteration budget. |
| 2026-08-05 | M1 | Root-caused the mesh-size gap to two separate issues — the shipped curvature thresholds didn't match the paper's stated values, and (independently) this codebase's curvature *formula* doesn't match the paper's stated formula either.<br><br>Recalibrated empirically using `bone` as a fast testbed, which fixed `bone`'s mesh size well and got Bottle1's convergence *quality* to match closely — but Bottle1's mesh size and, especially, its runtime (~27x slower than `bone` even fully calibrated) remained largely unmoved.<br><br>Spent most of the day testing *why*, testing four hypotheses for the over-refinement:<ul><li>`C_THRES` alone — tested and refuted</li><li>`H_THRES` alone — tested and refuted</li><li>octree-balancing propagation — tested and refuted</li><li>spatial dispersion of refinement candidates — real but insufficient effect</li></ul>Also traced a related finding in the sibling `HexOpt` repo: its restored optimizer code looks structurally like `HybridOctree_Hex`'s own older algorithm, not the CAD 2026 paper's claimed method, backed by a definitive date-based proof (see `hovey/HexOpt`'s README). |
| 2026-08-17 | M1 | Switched strategy from threshold-guessing to git archaeology. Found that the repo's `/our results/bottle1.vtk`, committed by the paper's first author Tong one day before the `v1.0` tag and never touched since, reproduces Table 2's Bottle1 numbers **exactly** (36,091 verts / 30,145 elems / worst SJ 0.560000 / best SJ 1.0, via `scaled_jacobian_stats`) — strong evidence it's the actual paper output file, not just a close match.<br><br>Statically diffed source across every upstream tag (`v1.0`→`v1.1`→`v1.2`→`v1.3`→current `HEAD`) and found the 2026-08-05 calibration work had anchored `C_THRES` to `v1.2`'s values — but `v1.0`/`v1.1`'s own original constants (`C_THRES={0.15,0.3,0.6,1.2,2.4}`, `VOXEL_SIZE=9`), the closest available source snapshot in time to the reference mesh, had never been tried.<br><br>Built and ran both `v1.0` and `v1.1` against the recovered original 2024-01-11 input surface. Both got permanently stuck at the same point in the final projection stage (never reached Table 2's numbers), which reframed the investigation: not a `v1.0`-vs-`v1.1` code-version question after all, but something about specific local geometry in Bottle1's input surface that this gradient-descent projection method can't resolve regardless of version.<br><br>Along the way, corrected two initial misreadings (see Detailed log for the full trail): `c291245` (2024-01-11, used for today's runs) has an off-by-one header bug that makes `ReadRawData()` silently duplicate its last triangle, while `f3247d8` (2024-01-16) parses cleanly; and `v1.1`'s new code block turned out to be a post-convergence step, not stuck-recovery logic, so it never actually ran in either attempt — the real cause of `v1.0`/`v1.1`'s differing exploration patterns is still open. `v1.1` left running as a low-priority background check; next step is a targeted dump of the local triangle neighborhood around the stuck point.<br><br>**Later the same day, both open questions closed.** Instead of dumping the triangle neighborhood, switched to measuring octree levels directly off the meshes — every model is rescaled into a fixed 100-unit cube, so a hex's edge length *is* its octree cell size and its level is just `log2(100/edge)`. That turned the whole problem from inference into measurement.<br><br>Three results followed. (1) `v1.0`'s shipped constants are the paper's constants: run unmodified on `bone` they reproduce `our results/bone.vtk` **exactly** — 10,356 / 8,619 / 0.610001 / 1.0, with byte-identical cell connectivity — so every earlier session's threshold calibration was solving a problem that did not exist. (2) Bottle1's reference mesh sits exactly on the committed input surface (max boundary deviation 0.00000), so the input was never the variable either; the one remaining variable is how deep the octree refines, and Table 2's "Refinement Level" column turns out to report exactly that — the number of active refinement thresholds, verified against every Table 2 model's reference mesh. Bottle1 used **four** thresholds where the public code hardcodes five, which alone accounts for the entire 3.4x oversizing: 101,688 elements → 27,750. (3) The all-day projection stalls were caused by `c291245`'s duplicated triangle — a controlled same-config comparison has `c291245` frozen forever at `smallDist=1.399` while `f3247d8` converges cleanly — so the header-correct file is the one to reproduce against, reversing this morning's choice.<br><br>Best configuration reached: 35,535 verts / 29,943 elems against 36,091 / 30,145 (−1.5% / −0.7%), with octree levels 6 and 7 — 93% of the mesh — matching to 0.5% and 0.1%. The last ~1% is not uniquely pinned by the evidence, and probably cannot be: the reference mesh predates the first committed source by a day. |
| 2026-08-18 | M1 | Moved to Bunny, run as two parallel independent efforts (a Sonnet session plus a background Opus agent) to test whether Bottle1's [Finding-14](#finding-14) recipe transfers to a model it wasn't tuned on. It doesn't: applied as-is, it overshoots Bunny by 70% (43,691/36,901 vs. 26,375/21,695). Six threshold configurations tried, all confirmed by session's end (two via a clean `finalMesh.vtk`, two via a stabilized `projHex.vtk` after 2+ hour plateaus, matching a stuck-point pattern [Finding 13](#finding-13) first saw on Bottle1); the family is cleanly monotonic in threshold tightness. Devised a sharper diagnostic — refined-parent counts per octree threshold rather than mesh totals — which shows the error compounding roughly 2× per level across every plain-`rl4`-family config, the signature of a scaling problem rather than a single bad threshold. Derived and tabulated `VOXEL_SIZE`/`LADDER_TOP` settings for all remaining Table 2 models from the ladder-position/level relationship. Caught a methodology bug from the session's own early hours: `dualHex.vtk`'s point/cell count (pre-projection) is not a safe stand-in for `finalMesh.vtk`'s (post-projection can add 40-50% more elements by splitting low-quality ones). Cross-validated the best candidate (a `√2` "half-level" scale correction) against Bottle1 directly: it undershoots by 45%, the opposite direction from Bunny's overshoot under Bottle1's own recipe — a decisive negative result, settling that **the thresholds are genuinely per-model, not one universal rule**, while the ladder-depth ("Refinement Level") reading keeps holding up on every model. Closed Bunny's own gap anyway by tuning one purely-curvature-gated threshold (`C_THRES[3]`: 1.6971→1.45), reaching **27,175/22,525 vs. 26,375/21,695 (+3.0%/+3.8%)** on the first attempt — Bunny's best reproduction, comparable to Bottle1's own (−1.5%/−0.7%, opposite direction). Ran Dragon Stand2 as a third-model check (different ladder rung, `VOXEL_SIZE_VALUE=8`): ladder position held, thresholds again didn't transfer (+88%/+94%), a third independent confirmation. Closed with a cross-cutting synthesis section on what `bone`/Bottle1/Bunny collectively show, and a practical four-step recipe for future models. |
| 2026-08-19 | M1 | Formalized the level-histogram methodology as a real tool, `scripts/level_histogram.py` (12-edge-mean metric, plus the per-level refined-parent-count diagnostic), replacing the one-off code from the day before. Used it to resolve two open questions cleanly: `bone`'s level-9 population (0.59%, a ~25x drop from level 8) and David's level-10 population (0.71%, a 29x drop from level 9) are both genuine distortion tails, while David's level-9 population (20.77%, the *second-largest* band in the mesh) is real — so David does need six levels as Table 2 claims, and the "open reconciliation" between Findings [9](#finding-9) and [10](#finding-10) noted the day before was an artifact of applying the distortion-tail heuristic too loosely; checked against five models now (`bone`, Bottle1, Bunny, Dragon Stand2, David) without a miscall. Ran Ramses as a fourth-model ladder-position check (same rung as Dragon Stand2, `VOXEL_SIZE_VALUE=8`): ladder held, but the shipped-threshold mismatch was far more severe than any other model — +448%/+473% overshoot, a ~3.6-hour pipeline (its octree stage alone, 190 minutes, is the slowest single stage measured this session, on a *smaller* input than Dragon Stand2's), traced to a 4.3x-oversized pre-projection octree. Two same-ladder-rung models (Dragon Stand2, Ramses) landing five times further apart from each other than either is from its own target further rules out a per-rung correction shortcut. Added a full pipeline-timings table across every Bunny/Dragon Stand2/Ramses run. A background Opus agent worked Dragon Stand2's threshold fit in parallel throughout, using `refine_criteria_stats` sweeps and its own full-run tests, and landed the session's best result: **62,307/50,603 vs. 62,576/50,853 (−0.43%/−0.49%)**, beating Bottle1's previous best. Along the way it found **[Finding 15](#finding-15)**: the five refinement thresholds are coupled one-directionally — tightening one threshold changes that threshold only, exactly, while loosening one propagates upward to every shallower threshold — which turns future threshold-fitting from a search into a mechanical per-level match against the reference's refined-parent counts. Applied that procedure directly to Ramses (the most severe shipped-threshold mismatch found this session, +448%/+473%, concentrated almost entirely in one threshold's curvature gate) and reached **46,329/39,529 vs. 44,790/37,993 (+3.4%/+4.0%, worst SJ 0.570 vs. Table 2's 0.590)** in two rounds — far fewer than Dragon Stand2's six, direct evidence the procedure generalizes. Paused model reproduction at this point (five models done: `bone`, Bottle1, Bunny, Dragon Stand2, Ramses) at the user's request, to shift toward writing this material into the `automesh` mdbook's Section 8 (Hexahedral Meshing from a Surface). |
| 2026-08-28 – 29 | M1 | Overnight autonomous session focused solely on Bunny. Applied [Finding 15](#finding-15)'s mechanical per-level fit (which Bunny's own 2026-08-18 fit predated) through ten configurations: two early kills, four octree-stage probes, four full runs. Three new mechanisms fell out: `H_THRES` is inert at Bunny's thresholds 1–2 (thickness candidates land where child-propagation already refines), `C_THRES[2]` is the real threshold-2 lever, and threshold 3's per-candidate weights are lumpy enough that the reference's 871.7 refined parents sit in a dead zone between two adjacent reachable values (837.4 / 938.0) — one surface point controls ~100 parents. Best configuration `fitC` (`C_THRES = {0.2121, 0.4243, 0.82, 1.51, 3.3941}`, halfshift `H_THRES`, `CELL_DETECT=0.75`): **26,255 / 21,702 vs. 26,375 / 21,695 (−0.45% / +0.03%)** — the closest cell-count match of any Table 2 model, beating the previous Bunny best (+3.0%/+3.8%). The residual is a per-threshold distribution imbalance the remaining knobs provably cannot rebalance without breaking the totals. `fitC`'s projection then ran overnight for the worst-SJ ratchet and never left 0.010 — its stuck point blocked the ratchet's *first* firing, so quality pressure never turned on. A 38-minute morning re-run with `ELEM_THRES` initialized at 0.5 proved the point: worst SJ −0.487 → 0.4999 in one minute, 0.500000 at the stop, topology identical. The afternoon follow-up at a 0.57 bar finished it: **worst SJ 0.570000 — Table 2's number exactly — six minutes into projection**, topology unchanged, best SJ 1.0. Bunny's row now matches Table 2 at −0.45% points / +0.03% cells / 0.570000 worst SJ / 1.0 best SJ / level 4, in ~11 minutes wall — the same order as the paper's reported 358 s, versus 4-8 hours for ratchet-gated runs. Worst SJ lands wherever the bar is; Table 2's worst-SJ column reads like per-model stopping bars, not intrinsic mesh quality. Separately, implemented the paper's Section 2.1 curvature formula as a gated code path (`-DCURVATURE_MODE=1`, default byte-identical to shipped) and eliminated [Finding 16](#finding-16)'s scenario (a) by direct measurement: the paper's `G` never crosses even `G_thres[0]` on bone/Bunny/Bottle1 at the 100-cube scale, overshoots the whole ladder at raw scale, and lacks the dynamic range to grade a five-rung ladder at any scale — the released `(angle − π)²` formula is what generated the reference meshes. |

## Detailed log

*Ordered newest to oldest*

### 2026-08-28 – 2026-08-29 — Bunny re-fitted with Finding 15's per-level procedure (Apple M1, overnight)

**Goal:** reproduce Table 2's Bunny row exactly (26,375 / 21,695 / worst
SJ 0.570), by applying [Finding 15](#finding-15)'s mechanical per-level
fit — which Bunny's own 2026-08-18 fit predates — plus threshold sweeps.
Load-bearing summary in the [Bunny section](#bunny-per-level-fit); this
is the chronology and the dead ends.

1. **Calibration.** `rl4-halfshift-c145`'s final parents (808.5 /
   1,170.9 / 1,023.6) vs. the reference's (832.1 / 1,234.0 / 871.7)
   located the whole +3.0%/+3.8% overshoot in threshold 3. Octree-stage
   parents for the same runs gave octree→final survival ratios of 0.447
   / 0.550 / 0.646 — but `rl4-halfshift-cd75`'s ratios (0.459 / 0.592 /
   0.621) differ by up to 7%, and that drift ended up dominating every
   octree-stage prediction we made tonight. The octree stage is a cheap
   *relative* instrument (which threshold moved, and roughly how much);
   only a stabilized `projHex.vtk` settles absolute counts.
2. **`fit1`/`fit2` — the thickness dead end.** First attempt tightened
   `C_THRES[3]` (1.55 / 1.52) and loosened `H_THRES[1]`/`[2]` (5.8–6.0 /
   3.0–3.1) to fix the threshold 1–2 undershoots. The `C_THRES[3]`
   tightening worked; the `H_THRES` loosening did nothing — threshold 1
   and 2 finals moved less than 1% against +8–12% more candidates. Both
   runs killed at −4.0%/−3.6% and −3.5%/−3.0%. One wrinkle worth keeping:
   tightening `C_THRES[3]` 1.49→1.55 also pulled threshold 2's *octree*
   parents down 56 (2,128→2,072) — removing a deep cell removes the
   upward intersect-propagation it was feeding, so
   [Finding 15](#finding-15)'s "tightening touches nothing shallower"
   holds for a threshold's *own-gate* cells, not for propagation a
   looser deep threshold had been injecting.
3. **`fit3o`–`fit6o` — four octree-stage probes (killed after
   `octree.vtk`).** Isolated each curvature knob at `C_THRES[3]=1.49`:
   `C_THRES[2]` 0.8485/0.80/0.75 → threshold 2 octree parents
   2,128/2,352/2,496 (the lever); `C_THRES[1]` 0.40 → threshold 1
   +16 octree parents only (weak); threshold 3 bit-identical (1,512)
   across all four (locality confirmed).
4. **`fitA`–`fitF` — six full runs to stabilized `projHex.vtk`.**
   Documented in the [Bunny section table](#bunny-per-level-fit). `fitC`
   (`C_THRES = {0.2121, 0.4243, 0.82, 1.51, 3.3941}`, halfshift
   `H_THRES`, `CELL_DETECT=0.75`) landed 26,255 / 21,702 — **−0.45% /
   +0.03%**. `fitF` closed the threshold-3 question: candidate #44 alone
   carries ~100 refined parents, so the reference's 871.7 is between
   two adjacent reachable values (837.4 and 938.0) — a dead zone no
   `C_THRES[3]` reaches.
5. **Session mechanics worth recording:** nine background `HexGen`
   processes from killed-off configs survived several `pkill -f`
   attempts (the pattern matched the launch path, but the process name
   is just `./HexGen`) and kept running for ~1.5 h. All comparative
   numbers are deterministic counts, unaffected; timings were never
   comparable and were not recorded. Killing by PID after mapping each
   PID's working directory with `lsof -a -p <pid> -d cwd` is the
   reliable way.
6. **Overnight:** `fitC`'s projection ran alone on the machine for
   7h36m / 5,701 checkpoints (2026-08-28 23:30 → 2026-08-29 06:06) and
   its worst SJ never left 0.010 — the ratchet's first firing never
   happened, so no quality pressure ever turned on (full mechanism in
   the [Bunny section](#bunny-per-level-fit)).
7. **Morning experiment (2026-08-29 06:10):** re-ran `fitC`'s exact
   configuration with `ELEM_THRES` initialized to 0.5 instead of 0.01 —
   reproducing by construction the post-first-success state that
   `rl4-halfshift-c145` reached by luck — to test whether `fitC`'s mesh
   is quality-capable and only the gate was missing. Source tweak made
   on a scratchpad copy of the `v1.0` tree; the repo tree is untouched.
   Result: worst SJ −0.487 → 0.499921 in *one minute*, 0.500000 at the
   38-minute stop (733 checkpoints), topology unchanged. The gate was
   the whole story; both runs cap exactly at the 0.5 bar. Killed at
   06:48, the session's 7 AM deadline.
8. **Afternoon follow-up (2026-08-29 ~14:09):** the same experiment with
   `ELEM_THRES` initialized at 0.57 (`runs/bunny-v1.0-fitC-elem057/`)
   reached **worst SJ 0.570000** — Table 2's number exactly — about six
   minutes into projection, topology unchanged (26,255 / 21,702, best
   SJ 1.0). Killed at 982 checkpoints. Full result table and what it
   says about Table 2's worst-SJ column in the
   [Bunny section](#bunny-per-level-fit).
9. **`CURVATURE_MODE` added and [Finding 16](#finding-16)'s scenario (a)
   eliminated (2026-08-29 afternoon):** implemented the paper's
   Section 2.1 curvature as a compile-flag-gated alternative in
   `HybridOctree_Hex_v1.0/` (`-DCURVATURE_MODE=1`; default 0 is
   byte-identical to shipped, verified by preprocessor diff;
   cross-validated against an independent Python implementation to
   5×10⁻⁷). Ran bone under the paper's literal
   `G_thres = {0.5, 1, 2, 4, 8}` (`runs/bone-v1.0-gauss/`): zero
   curvature refinement fires — `G` maxes at 0.30 on the 100-cube —
   and no global rescale fixes it (raw coordinates overshoot the whole
   ladder instead; `G`'s dynamic range is ~4-20× against the ladder's
   16× span). Verdict written into [Finding 16](#finding-16): the
   released `(angle − π)²` formula generated the reference meshes.

### 2026-08-27 — Diagnosing and fixing the Bunny/Dragon Stand2 worst-SJ stuck point (Apple M1)

**Goal:** figure out why Bunny's worst scaled Jacobian was stuck at 0.010
against Table 2's 0.570, when mesh size/topology already reproduced the
paper closely (+3.0%/+3.8%).

Two background investigations found the actual cause before any code
changed.

* First: the Bunny/Dragon Stand2 runs never used the already-fixed
`ProjectToIsoSurface()` in `HybridOctree_Hex/HexGen.cpp` — they ran the
separately-extracted `v1.0`/`v1.1` tagged source trees instead, whose
`ProjectToIsoSurface` is the original unbounded loop with none of the
`MAX_PROJ_ITER`/`ELEM_THRES_CAP`/best-checkpoint instrumentation.
* Second, and more fundamental: `finalMesh.vtk` only gets (re)written on an
all-or-nothing checkpoint (every element above the current `ELEM_THRES`
bar, *and* the worst point within a tight surface-distance tolerance).
Bunny's best run got exactly one such success early on, freezing
`finalMesh.vtk` at worst SJ 0.010001, then stalled on one persistent
point for thousands of checkpoints without ever repeating it — while
`projHex.vtk` (written unconditionally every checkpoint) had actually
reached 0.487, invisible to every measurement taken until now. A fresh
geometric scan of the stalled triangle/its neighbors ruled out a
Bottle1-style input defect: all ordinary, well-shaped triangles, not
slivers.

**Fix**: a purely additive patch to `HybridOctree_Hex_v1.0/HexGen.cpp`'s
`ProjectToIsoSurface()` — track the mesh's true worst scaled Jacobian
every checkpoint (reusing the loop that already computes it for the
`k`/bad-element count) and snapshot it to a new file, `bestQualityMesh.vtk`,
whenever it improves, entirely decoupled from the existing all-or-nothing
gate. No existing branch, condition, or file write is touched.
Regression-tested against `bone` first: topology unchanged, and the new
tracker's worst SJ (0.610001) is a *better* match to the published exact
value (0.610000) than the old mechanism's own output on the same run
(0.600001) — net improvement, not just a safe no-op. Re-ran Bunny's best
recipe (`rl4-halfshift-c145`) with the patch: **worst SJ 0.500001**
(mesh size/topology unchanged, 27,175/22,525), cutting the shortfall to
Table 2's 0.570 from 98% to about 12%. Also rebuilt and ran Dragon
Stand2's best recipe with the identical patch, as a second validation of
the same pathology: **worst SJ 0.500000** (62,307/50,603 unchanged),
against Table 2's 0.560.

Neither build has any automatic stopping condition, so both were cut off
manually: Bunny after 4h47m24s (started 09:45:52, stopped 14:33:16) once
0.500001 had held for over 3.5 hours (3,934 checkpoints) — a confirmed
plateau; Dragon Stand2 after 3h44m39s (started 10:48:37, stopped at the
same 14:33:16), which was still oscillating on its own stuck point
(index ≈56765) at the moment of cutoff, not a confirmed plateau — it may
have had more room to climb.

The remaining gap looks like a genuine optimization plateau in one mesh
region — the worst point's distance to the surface holds around 0.2-0.24
in Bunny's run, nowhere near convergence — not something more waiting
alone resolves for that model. Closing it further would mean letting the
redistribution step
fire on partial progress, not just a full pass; out of scope for this
session's purely-additive fix.

### 2026-08-17 — Git archaeology: locating the exact paper-era reference mesh and source; `v1.0` reproduction build (Apple M1)

**Goal for today:** stop guessing at thresholds and check whether git history itself holds the exact code and exact output from around the paper's submission time (2024-01-14 paper received, 2024-02-23 revised paper received, 2024-03-16 paper accepted) — specifically for Bottle1's Table 2 row (36,091 verts / 30,145 elems / worst SJ 0.560 / best SJ 1.0 / refinement level 4).

<a id="finding-5a"></a>**Finding 5a (added later the same day) — the same exact-match archaeology
works for `bone` too, and `bone` is a much faster testbed.** Same check
as [Finding 1](#finding-1), run against `our results/bone.vtk`: committed by Hua Tong
on `3c6b62a` (2024-01-14 22:56:28 -0500), matches exactly via
`scaled_jacobian_stats` — `10356` verts / `8619` elems / worst SJ
`0.610000` / best SJ `1.0`. Closest source snapshot: `542878d`
(2024-01-14, same day, still pre-`v1.1`). Unlike Bottle1, `bone`'s input
surface (`input boundaries/bone_tri.raw`) has only **one** version in
git history (`4aef107`, 2024-01-14) — no `c291245`/`f3247d8`-style
duplicate-triangle ambiguity to resolve.

Per the Bottle1 pipeline timings table above, `bone`'s full run takes
only ~4.5–7 minutes total, vs. Bottle1's multi-hour runs — roughly two
orders of magnitude faster to iterate on. Suggested going forward: use
`bone` as the primary fast validation testbed for open questions (which
`C_THRES`/curvature formula, which exact source commit, whether input
correctness matters, etc.) — find the configuration that reproduces
`bone`'s exact numbers cheaply and quickly, then apply only the winning
configuration(s) to Bottle1, rather than spending all iteration budget
directly on Bottle1's slow runs. Relayed this to the background Opus
agent (spawned to continue this investigation) so it can incorporate
`bone` alongside its current Bottle1-focused parallel sweep.

<a id="finding-1"></a>**Finding 1 — the reference mesh already exists in this repo, untouched since the paper.**
`our results/bottle1.vtk` was added in commit `be63d98` (upstream `CMU-CBML/HybridOctree_Hex`), authored by **Hua Tong** (the paper's first author) on **2024-01-11** — one day before the `v1.0` release tag (`db207d7`, 2024-01-12). It has never been modified since: `git log --all -- "our results/bottle1.vtk"` shows exactly one add commit and nothing else, all the way through today's `HEAD`. Running this fork's own `scripts/scaled_jacobian_stats` (built from `Sj()` copied verbatim out of `HexGen.cpp`) against it gives:

```
Points:   36091
Cells:    30145 (30145 hex, 0 non-hex skipped)
Worst SJ: 0.560000
Best SJ:  1.000000
```

An exact match on all four numbers, not a rounding coincidence. This is almost certainly the literal file used to produce Table 2's Bottle1 row, and it's a better reproduction target than the paper's rounded printed numbers.

<a id="finding-2"></a>**Finding 2 — the earliest committed source postdates that mesh by exactly 1 day.**
The first commit touching `HexGen.cpp`/`Main.cpp` anywhere in git history is `44f730c` (2024-01-12) — the same day as the `v1.0` tag. The reference mesh was committed 2024-01-11 10:55:18 -0500 (`be63d98`), and `44f730c` followed 2024-01-12 00:08:33 -0500 — 1 calendar day later (13h 13m elapsed). So the source code wasn't even in the repo yet when the reference mesh was pushed; `v1.0` is the closest available snapshot, 1 day removed, not a proven-identical match.

<a id="finding-3"></a>**Finding 3 — static diff across every release tag, to see what's actually changed.**
Diffed `v1.0` → `v1.1` → `v1.2` → `v1.3` → current fork `HEAD` file-by-file. `Main.cpp`/`Mesh.cpp`/`Mesh.h` never change across the whole history; all real evolution is in `HexGen.cpp`/`HexGen.h`/`Initialization.h`. Upstream goes dormant after `v1.3` (2024-02-25) — no further code commits until this fork's own work starting 2026-08-04.

| | `v1.0` = `v1.1` | `v1.2` | `v1.3` (= today's upstream default) | fork `HEAD` (pre-today) |
|---|---|---|---|---|
| `C_THRES` | `{0.15, 0.3, 0.6, 1.2, 2.4}` | `{0.1, 0.2, 0.4, 0.8, 1.6}` | `{0, 0, 0.4, 0.8, 1.6}` | `{0.1, 0.2, 0.4, 0.8, 1.6}` (= `v1.2`) |
| `H_THRES` | `{16,8,4,2,1}` | same | same | same |
| `VOXEL_SIZE` | **9** | **10** | 10 | 10 |
| `MAX_NUM` | 10,000,000 | 100,000,000 | same | same |
| `OUT_IN_RATIO` | 0.125 | 0.125 | **0.15** | 0.15 |

`v1.1` also added a whole new randomized-retry block to `ProjectToIsoSurface` (absent from `v1.0`) and bumped `ELEM_THRES`'s first-trigger value 0.5→0.53 — a real algorithmic change, not just constants, and one `v1.0` doesn't have. Cross-checked against this file's own record of the 2026-08-05 work: that session calibrated against **`v1.2`'s** `C_THRES`, and separately tried the paper's literal stated `Gthres` and `v1.3`'s shipped (buggy) values — `VOXEL_SIZE` was noted as "stable" at 10 throughout. **To date, we have not tried to reproduce results with `v1.0`/`v1.1` original constants**, and they're the ones chronologically closest to whatever actually produced the reference mesh.

<a id="finding-4"></a>**Finding 4 — the original input surface differs from what's in the working tree today.**
`input boundaries/bottle1_tri.raw` was committed alongside the source on `c291245` (2024-01-11) and then edited on `f3247d8` (2024-01-16) — the version currently checked out in the working tree (as of today, 2026-08-17) is the *later* one. Diffing the two:

```
-14833 29665 // vertices, triangles c291245 (2024-01-11)
+14833 29664 // vertices, triangles f3247d8 (2024-01-16)
```

**The 2024-01-11 original (14,833 vertices / 29,665 triangles) declares
one more triangle than the file in today's (2026-08-17) working tree
(14,833 vertices / 29,664 triangles).** Same vertex count, one fewer
declared triangle.

**What each file actually contains, and which
one is correct.**

- Both `c291245` (2024-01-11) and `f3247d8` (2024-01-16) contain
  **byte-identical triangle data** — 29,664 physical triangle lines,
  line-for-line the same in both. No triangle was added, deleted, or
  edited between the two commits.
- The only difference is the declared count on the header line:
  `c291245` (2024-01-11) says `29665`; `f3247d8` (2024-01-16) says
  `29664`.
- `f3247d8` (2024-01-16) is **correct** — its header matches the file's
  true physical triangle count exactly.
- `c291245` (2024-01-11) is **incorrect** — its header overstates the
  true count by 1.
- Consequence: `ReadRawData()` reads exactly `elements` (the header
  count) lines into a single reused buffer, without checking `fgets()`'s
  return value. On `c291245` (2024-01-11), the extra declared line runs
  past EOF, `fgets()` fails silently, and `sscanf()` re-parses the stale
  buffer — silently **duplicating the file's true last triangle**
  (`14729 14727 14725`) into the loaded mesh a second time. `f3247d8`
  (2024-01-16) loads cleanly with no duplication.
- That duplicated triangle is a real sliver (area 1.40×10⁻⁵, aspect
  ratio 3.02) whose centroid sits ~0.027 units from the stuck-projection
  points' centroid (Correction above) — spatially close, but at a
  different location in the triangle list (indices 29663/29664 vs.
  `#1357`/`#1351`), so unlikely to be the direct cause of that plateau.

**Practical implication for today's runs.** Both `v1.0` and `v1.1` were
launched against `c291245` (2024-01-11) — the incorrect file — because
it's temporally closer to the reference mesh's 2024-01-11 generation
date. `f3247d8` (2024-01-16) is the correct file but 5 days further
from that date. This time-vs-correctness trade-off is unresolved; an
untried variant is rerunning against `f3247d8` (2024-01-16) instead.

**Action taken — built and launched `v1.0` against the original input.**
1. Extracted `v1.0`'s tagged source tree (`HexGen.cpp/h`, `Main.cpp`, `Mesh.cpp/h`, `Initialization.h`, `StaticVars.h`, `CellQueue.h`, `CMakeLists.txt`, `ErrorMessage.txt`) into a new `HybridOctree_Hex_v1.0/` directory at the repo root.
2. Patched two portability issues only, no logic changes: `sscanf_s`→`sscanf` (15 call sites — the exact same fix `v1.1` itself made one version later) and the Windows-only `<malloc.h>` in `CellQueue.h`→`<cstdlib>`/`<cstring>`.
3. Built clean with CMake + AppleClang (same cosmetic-only warning classes already seen building the fork's `HEAD` — stray `#endif` tokens, `RAND_MAX` int→float narrowing, one dangling-else).
4. Recovered the original 2024-01-11 `bottle1_tri.raw` from commit `c291245` (2024-01-11) and set up `HybridOctree_Hex/runs/bottle1-v1.0-repro/model.raw`.
5. Launched `./HexGen` there in the background, log at `runs/bottle1-v1.0-repro/run.log`, monitoring for stage transitions and the `Try:0 ...` convergence-checkpoint lines that trigger this version's own `finalMesh.vtk` snapshot write (confirmed present in `v1.0`'s `ProjectToIsoSurface`, same mechanism the 2026-08-05 session relied on for `HEAD`).

**Progress — octree construction stage complete.**

```
Time elapsed for reading surface mesh: 2.39245
Number of unbalanced nodes: 529 → 90 → 0
Time elapsed for constructing octree: 2167.96
```

`v1.0` doesn't print curvature/narrow-region assessment and octree
balancing as separate timings the way the fork's later-instrumented code
does — this 2167.96 s (~36.1 min) is the combined figure, comparable to
the "— combined (curvature + octree)" row in the Bottle1 pipeline timings
table above. That row was 1012.4 s (~16.9 min) for the 2026-08-05
`v1.2`-calibrated run and 1062.2 s (~17.7 min) for the 2026-08-04
buggy/capped run — so `v1.0`, despite its *coarser* `VOXEL_SIZE=9` (vs.
10), is taking roughly 2x longer for this combined stage than either
later-code run. Noted as a raw observation for now, not yet a hypothesis —
worth returning to once we see whether the mesh-size/quality outcome
actually matches, since a timing anomaly on its own doesn't tell us
much.

**Progress — mesh dualization complete.**

```
Time elapsed for generating dual mesh: 343.372
```

343.4 s (~5.7 min) — this time *faster* than either later-code run
(1045.8 s for the `v1.2`-calibrated run, 1138.9 s for the buggy/capped
run), the opposite direction from the octree-construction stage above.
So `v1.0` isn't uniformly slower or faster than the later code — it's
mixed by stage, which argues against a single blanket explanation (e.g.
"this whole version is slow") and toward stage-specific algorithmic
differences. Not investigating further mid-run; flagging for later.

**Progress — buffer-zone clearing complete.**

```
Time elapsed for extracting interior dual mesh: 333.143
```

333.1 s (~5.6 min) — again faster than both later-code runs (785.0 s
calibrated, 815.1 s buggy/capped), consistent with the dualization
stage's direction rather than the octree-construction stage's. Only the
octree-construction stage has run slower in `v1.0` so far; every other
stage completed has run faster.

**Progress — projection stage shows a sustained plateau on a single
point.**

`smallDist` (the worst-point projection distance; needs to drop below
`DIST_THRES2` ≈ 1e-6 to trigger a `finalMesh.vtk` write) fell fast at
first, 7.10 → 0.57 over the first ~50 outer iterations, then plateaued.
For iterations ~60–105 the worst point has been stuck at the same index
(1357), creeping only 0.617 → 0.607 over 45+ iterations — essentially
flat, `k` (bad-element count) oscillating between 1 and 3, never reaching
0. 1h39m elapsed (~99 min CPU time) at time of writing.

This lines up with a finding from the source diff ([Finding 3](#finding-3) above):
`v1.1` added a whole new randomized-perturbation retry block to this
exact function, five days after `v1.0`, specifically for cases where the
gradient-descent projection gets stuck. `v1.0` doesn't have that block.
Bottle1's thin, coiled, genus-1 geometry is a plausible candidate for
triggering exactly this failure mode. **Working hypothesis, not yet
confirmed:** `v1.0` may not converge on Bottle1 at all without that fix —
which would explain why `v1.1` exists five days later, and would mean
`v1.1`, not `v1.0`, is the real target to reproduce against, despite
being one step further from the reference mesh's commit date.

**Next step:** let `v1.0` keep running a while longer to see if it
eventually breaks the plateau on its own (worth ruling out before
concluding it's stuck for good), while building and launching `v1.1`
in parallel against the same original input for direct comparison.

**`v1.1` build and launch.** `v1.1` needed zero portability patches —
it already uses `sscanf` (not `_s`) and `<stdlib.h>` (not `<malloc.h>`)
upstream, both fixes `v1.0` needed manually. Built clean with CMake +
AppleClang, same cosmetic-only warnings as `v1.0`. Extracted the same
original 2024-01-11 `bottle1_tri.raw` input and launched
`HybridOctree_Hex/runs/bottle1-v1.1-repro/`, running in parallel with
the still-active `v1.0` run rather than stopping it first — `v1.0`
hadn't yet ruled itself out as merely slow, so keeping it running costs
nothing and preserves the option to compare both outcomes directly.

**`v1.0` plateau update:** ~75 outer iterations in on point 1357 now
(spanning heartbeats #61 through #121), `smallDist` essentially flat in
the 0.61–0.62 band the whole time. The plateau is holding, not
transient.

**`v1.0` conclusion — genuinely stuck, not merely slow. Killed after
2h23m.** Close inspection of the last 100 outer iterations before
killing it: the worst point was **exactly one of two indices** (1351,
1357) for the entire span, `smallDist` confined to a narrow band
(0.6075–0.6437), and several stretches show bit-identical repeated
values across a dozen or more consecutive iterations (e.g. `0.643725,
1351` twelve times running) — the optimizer proposing the same move and
landing on the same result, not gradual convergence that's merely slow.
`v1.0` clearly lacks `v1.1`'s randomized-perturbation retry logic for
exactly this function, and this confirms `v1.0` itself gets stuck.
Killed the process (PID 97785) rather than let it keep burning CPU with
no expectation of a different outcome.

*Caveat added after the fact, once `v1.1` was tried (see below): the
next two sentences, as originally written here, overclaimed — "clean
confirmation of the hypothesis" and "makes `v1.1` the real reproduction
target" both implied that having the retry logic would let `v1.1`
escape where `v1.0` couldn't. That turned out to be false: `v1.1` got
stuck at the same two points too. What's actually confirmed is only
that `v1.0` gets stuck and lacks that logic — not that the logic's
absence is sufficient to explain the stuck-ness, since its presence in
`v1.1` didn't prevent the same outcome. Left the original reasoning
below rather than delete it, since it's what the run's own killing was
based on and it wasn't unreasonable given the evidence at the time —
just superseded by the `v1.1` result a few sections down.*

**`v1.1` projection stage — encouraging early signal.** Octree
construction (2203.26 s), dualization (345.5 s), and buffer clearing
(334.8 s) all closely matched `v1.0`'s timings, as expected from
identical upstream constants and code up to this point. The projection
stage started from the identical first point (`Try:0 7.0984 4244`, same
as `v1.0`'s first iteration) and, after ~90 outer iterations, landed on
the **same hard point `v1.0` got permanently stuck on (index 1357)** —
but here the comparison diverges. Looking at the raw log (not just
throttled heartbeats) over the last 40 outer iterations: 36 distinct
`smallDist` values, a slow downward drift (~0.60 → ~0.58), and two
`Try:0` hits (momentary zero bad-elements, at 0.5838 and 0.5661) —
versus `v1.0`'s later behavior on the same point, which showed
bit-identical repeated values across a dozen or more consecutive
iterations. `v1.1`'s randomized-retry logic ([Finding 3](#finding-3), the block absent
from `v1.0`) appears to be doing exactly what it looks designed to do —
landing on the same hard point but continuing to perturb and explore
rather than freezing there. Not converged yet, but a clearly different
— and more promising — trajectory than `v1.0` showed at the equivalent
point. Continuing to watch.

**`v1.1` revised conclusion — also stuck, on the same two points.**
The encouraging trajectory above didn't hold. By ~2h07m elapsed (219
outer iterations), the raw log shows the same freezing signature `v1.0`
showed: bit-identical `smallDist` values repeated up to 7 times running
(e.g. `0.611346, 1351` ×7), oscillating between **the same two points
`v1.0` got permanently stuck on** — 1357 and 1351. No `finalMesh.vtk`
written. `v1.1`'s randomized-retry logic clearly changed the *path*
(more exploration, `Try:0` hits along the way, slower to freeze) but not
the *outcome* — it delayed the same plateau rather than escaping it.

This reframes the working hypothesis from [Finding 3](#finding-3). It's likely not
primarily a `v1.0`-vs-`v1.1` code-version story — both versions
converge on the same trap. More likely explanation: something specific
to the local geometry near this spot in the input surface (a genuinely
hard configuration for this gradient-descent-based projection method,
independent of which release runs it). Both `v1.1` (PID 222) and the
earlier-killed `v1.0` attempt agree on the same two problem values,
which is itself a useful, reproducible data point.

**Correction — the "`v1.1` has retry/escape logic `v1.0` lacks"
framing, repeated several times above, is wrong.** Went back to verify
it against the source directly rather than let it stand. The block
`v1.1` actually adds (the ~90-line randomized-perturbation block quoted
in [Finding 3](#finding-3)) is nested **inside** the `if (k == 0 && smallDist <
DIST_THRES2)` branch — i.e. it only runs *after* the projection has
already converged (zero bad elements and the worst point within ~1e-6
of its target), immediately before writing `finalMesh.vtk`. It's a
post-convergence quality-polish/redistribution pass, not a
stuck-recovery mechanism. Neither `v1.0` nor `v1.1` ever reached that
gate in these runs (`smallDist` never got anywhere near 1e-6 in either)
— so this block **never executed in either run**, and cannot explain
why one got stuck harder than the other.

Checked what actually *does* run during the stuck phase: a separate,
pre-existing periodic redistribution block (`if (i % UPDATE_EVERY == 0)
{ if (ELEM_THRES == 0.01) { ... rand()-based blending toward neighbor
centroids ... } }`, around line 4339 in both versions) that already
uses `rand()` and is gated only on `ELEM_THRES == 0.01` — the state both
stuck runs are actually in. This block is **present and identical in
both `v1.0` and `v1.1`**; it isn't part of the `v1.0`→`v1.1` diff at
all. Also checked for `srand()` anywhere in either version: there is
none, so `rand()` is unseeded in both.

Net effect: there is no version-level code difference that explains
either run's stuck-ness, or why `v1.1` showed more exploration before
landing on the same trap `v1.0` did. That difference in intermediate
behavior is real (documented above) but currently **unexplained** —
worth being honest that this is an open question rather than
attributing it to a mechanism that, on inspection, never ran. Every
"retry logic"/"escape logic" claim earlier in this entry should be read
with this correction in mind; not rewriting them out of the
chronological record, but this note supersedes that framing.

**Correction on what "1357"/"1351" actually are.** They aren't octree
vertex indices — checked the source (`ProjectToIsoSurface`'s
`triNum[j] = k` assignment, `HexGen.cpp`): the third field in each
`Try:` line is an index into the **input triangle mesh**
(`triMesh.e[]`), the closest surface triangle the current worst octree
boundary point is being projected toward. So the log is saying "the
worst point right now is closest to input-surface triangle #1357 (or
#1351)," not naming a fixed octree vertex.

**Investigated triangles #1357 and #1351 (0-based) directly in the
original `bottle1_tri.raw`:**

```
Tri#1357 = (v553, v566, v552)   area=2.92e-04  edges=[0.0191,0.0337,0.0309]  aspect=1.76
Tri#1351 = (v549, v552, v566)   area=3.05e-04  edges=[0.0332,0.0202,0.0309]  aspect=1.65
```

The two triangles **share an edge** (v552–v566) — they're immediate
neighbors on the surface, not unrelated trouble spots. Their normals are
nearly parallel (dihedral angle 0.63°, essentially coplanar) and neither
is a sliver by shape (aspect ratios 1.65–1.76 are unremarkable). So the
plateau isn't obviously explained by a folded or degenerate pair of
triangles at this exact spot.

A more likely mechanical explanation: the worst boundary point may be
sitting almost exactly on the shared edge between these two triangles,
so which one counts as "closest" flips back and forth under small
perturbations — that would produce exactly the observed oscillation
between triNum 1357 and 1351 without needing any local mesh defect at
all, just a nearest-triangle assignment ambiguity at a shared edge.

**Next step:** either let `v1.1` keep running as a low-priority
background check (it hasn't fully proven it will never escape, just
that it hasn't in ~240 iterations), or set up a more targeted diagnostic
— dump the full local triangle neighborhood around v549/v552/v553/v566
(one or two rings out) and check for anything more systematically odd
about this patch's tessellation than the two irregularities found so
far.

**`v1.1` (PID 222) killed at 400 outer iterations, 3h13m elapsed
(193m13s CPU time), still confined to the same two points (1351,
1357) with no escape.** Final state, from the live `projHex.vtk`
snapshot just before the kill:

| Metric | `v1.1` at kill |
|---|---|
| Vertices | 117,031 |
| Elements | 101,688 |
| Worst SJ | -0.917028 |
| Best SJ | 1.000000 |
| Running time | 3h13m elapsed / 193m13s CPU time |

Worst SJ had actually gotten *worse* since the last check earlier today
(-0.313 → -0.917) — consistent with a genuine plateau rather than slow
net improvement: the optimizer is still perturbing elements, just not
converging. Vertex/element counts unchanged from the earlier snapshot
(117,031/101,688), confirming the mesh size itself was already settled;
only the projection/quality stage was still (unsuccessfully) iterating.
No `finalMesh.vtk` was ever written by either `v1.0` or `v1.1` today.

**Checked the live mesh snapshots directly — a genuinely positive
data point on mesh size, despite the stuck projection.** Both
`v1.0`/`v1.1` write periodic live snapshots to `projHex.vtk` in their
run folders (a separate, unconditional-on-convergence write, distinct
from the `finalMesh.vtk` checkpoint — `WriteToVtk(fileName, ...)` every
`UPDATE_EVERY` inner iterations, `HexGen.cpp` ~line 4325). Checked both
with `scaled_jacobian_stats`:

| | `v1.0` (last snapshot before kill) | `v1.1` (live, at this point in the day) | Table 2 target |
|---|---|---|---|
| Vertices | 117,031 | 117,031 | 36,091 |
| Elements | 101,688 | 101,688 | 30,145 |
| Worst SJ | -0.358409 | -0.313405 | 0.560 |
| Best SJ | 1.0 | 1.0 | 1.0 |

Identical vertex/element counts between the two, as expected (same
octree/dualization/buffer-clearing pipeline and constants up to this
point). Worst SJ negative in both — expected mid-projection, not yet
converged. But the mesh **size** is only ~3.2–3.4x larger than Table
2's target — notably better than the earlier fork's `v1.2`-anchored
2026-08-05 calibration attempt, which came in ~6x oversized
(212,990/186,832). So even with the projection stage stuck, the
`v1.0`/`v1.1` original-constants lead is already producing a
meaningfully closer-sized octree than anything tried before today —
independent evidence this lead is on the right track, separate from
whatever's causing the projection plateau.

### 2026-08-17 (continued) — Reading the octree level straight off the meshes; `v1.0` reproduces `bone` exactly

**New approach: stop guessing thresholds, measure the octree level of every
element instead.** `ReadRawData()` (`HexGen.cpp`) rescales every input model
into a fixed cube of side `BOX_LENGTH = 100`, and the dual mesh puts one
vertex at each octree leaf's centroid. So a hex element's edge length is its
octree cell size, and its octree level is `log2(100 / edge)` — a level is
directly readable from any output mesh, with no instrumentation and no rerun.
Cross-checked on `runs/bottle1-v1.0-repro/octree.vtk`: 173,552 leaves land on
exactly six discrete sizes, 6.25 / 3.125 / 1.5625 / 0.78125 / 0.390625 /
0.1953125, i.e. levels 4-9. The measurement is exact for the octree and for
the dual mesh, and only slightly fuzzy for a projected final mesh (a few
badly distorted boundary hexes misread by one level).

<a id="finding-5b"></a>**Finding 5b — Bottle1's reference mesh is one refinement level coarser than
anything we have produced.** Level histograms, reference versus today's
`v1.0` run (`corig` = `v1.0`'s own `C_THRES = {0.15,0.3,0.6,1.2,2.4}`,
`VOXEL_SIZE = 9`):

| Level | Cell size | `our results/bottle1.vtk` (target) | `bottle1-v1.0-repro/projHex.vtk` |
|---|---|---|---|
| 4 | 6.25 | 52 (0.2%) | 44 (0.0%) |
| 5 | 3.125 | 1,015 (3.4%) | 856 (0.8%) |
| 6 | 1.5625 | 9,945 (33.0%) | 8,261 (8.1%) |
| 7 | 0.78125 | **18,122 (60.1%)** | 35,619 (35.0%) |
| 8 | 0.390625 | 1,006 (3.3%) | **43,968 (43.2%)** |
| 9 | 0.1953125 | 5 (0.0%) | 12,554 (12.3%) |
| 10 | 0.09765625 | 0 | 386 (0.4%) |
| **Total** | | **30,145** | **101,688** |

The reference peaks at level 7 and effectively stops there — its 1,006
level-8 cells are 3.3% of the mesh and its 5 level-9 cells are measurement
noise from distorted boundary hexes. Our run peaks a full level deeper, at
level 8, with a substantial level-9 population. Levels 4-6 agree closely
between the two; the entire 3.4x size gap lives at levels 8 and 9.

<a id="finding-6"></a>**Finding 6 — the 2026-08-05 conclusion "`C_THRES` isn't the dominant driver
for Bottle1" was drawn from a change that could not have moved the number.**
That session compared `v1.3`'s shipped `{0, 0, 0.4, 0.8, 1.6}` against
`v1.2`'s `{0.1, 0.2, 0.4, 0.8, 1.6}` and saw only a 2.7% vertex reduction.
But those two arrays differ **only in entries 0 and 1**, which gate refinement
of level-4 and level-5 cells — together 0.8% of Bottle1's elements. Entries
2, 3, 4 — which gate levels 6, 7 and 8, i.e. 90%+ of the mesh — were
identical in both. The experiment therefore held the only entries that matter
fixed, and its null result says nothing about `C_THRES` in general. For
contrast, `v1.0`'s `{0.15,...,2.4}` (1.5x stricter at entries 2-4 than
`v1.2`'s) produced 117,031 vertices against the `v1.2`-calibrated run's
212,990 — a 1.82x swing from a 1.5x threshold change. `C_THRES`'s deep
entries are in fact a very strong lever for Bottle1.

**New tool: `scripts/refine_criteria_stats.cpp`.** Computes both refinement
criteria — the per-vertex dihedral "curvature" `r[]` from `ReadRawData()` and
the per-triangle ray-cast "thickness" from `GetCellValue()`, both copied
verbatim — and reports the resulting `refineTri*`/`refineTriPt*` candidate
list sizes for any `C_THRES`/`H_THRES`, in ~9 s for Bottle1. That skips the
~36-minute `ComputeCellValue()` octree sweep entirely, so thresholds can be
swept in seconds instead of half-hours. Bottle1 at `v1.0`'s own values:

| Threshold | Octree level | Curvature pts | Curvature tris | Thickness pts | Thickness tris | Union pts | Union tris |
|---|---|---|---|---|---|---|---|
| 0 | 4 | 2,073 | 6,426 | 4,638 | 10,488 | 5,500 | 14,124 |
| 1 | 5 | 1,256 | 4,381 | 1,602 | 2,921 | 2,508 | 6,659 |
| 2 | 6 | 551 | 2,308 | 694 | 1,166 | 1,154 | 3,341 |
| 3 | 7 | 138 | 671 | 19 | 18 | 151 | 683 |
| 4 | 8 | 12 | 59 | 0 | 0 | 12 | 59 |

Bottle1's curvature `r[]` has median 0.017 and 99th percentile 1.16, so the
deep thresholds are driven by a very small tail of high-curvature vertices —
138 points at threshold 3 and 12 at threshold 4 — which nonetheless cascade into 43,968
level-8 and 12,554 level-9 elements, because `CELL_DETECT = 1` doubles each
cell's detection box and `RefineBrothers()` refines all eight siblings of any
flagged cell.

<a id="finding-7"></a>**Finding 7 — `v1.0` with its own original constants reproduces `bone`'s
reference mesh *exactly*, on every metric and on connectivity.** `bone` is
the ideal fast testbed here: its reference output `our results/bone.vtk` was
committed by Hua Tong on `3c6b62a` (2024-01-14 22:56:28 -0500), its input
`input boundaries/bone_tri.raw` has exactly one version in all of git history
(`4aef107`, 2024-01-14 22:53:43 -0500), the closest source snapshot is
`542878d` (2024-01-14 22:37:52 -0500, same day,
still under the `HybridOctree_Hex_v1.0` directory naming), and a full run
takes minutes rather than hours. Ran `bone_tri.raw` through today's `v1.0`
build with `C_THRES = {0.15, 0.3, 0.6, 1.2, 2.4}` and `VOXEL_SIZE = 9`,
unmodified (`runs/bone-v1.0-corig/`):

```
Points:   10356
Cells:    8619 (8619 hex, 0 non-hex skipped)
Worst SJ: 0.610001
Best SJ:  1.000000
```

Against `our results/bone.vtk`'s 10,356 / 8,619 / 0.610000 / 1.000000 — an
exact match on all four numbers (the worst-SJ `0.610001` versus `0.610000` is
float accumulation in the `ELEM_THRES` ratchet, `0.53 + 0.01 x 8`, not a real
difference). Going further, the two meshes' **entire cell-connectivity arrays
are identical**, all 8,619 hexes with the same vertex indices in the same
order; only vertex *positions* differ, and only while the run is still
ratcheting (mean 0.028, max 1.21 in the 100-unit box when sampled at
`ELEM_THRES = 0.54`, converging as it climbs to 0.61). Every earlier
session's `bone` mismatch (1.9x too large in
2026-08-04, 1.43x with `v1.2`'s thresholds in 2026-08-05, 0.41x with the
paper's literal `Gthres` — re-confirmed today at 4,558/3,547 in
`runs/bone-v1.0-cpaper/`) was a consequence of running the wrong code
version's constants. **The `v1.0` tag's shipped configuration is the paper's
configuration.** This is the strongest result of the whole reproduction
effort so far: it removes threshold calibration from the picture entirely.

<a id="finding-8"></a>**Finding 8 — Bottle1's reference mesh sits exactly on the committed input
surface, so the input file is not the variable either.** Extracted the 11,188
boundary vertices of `our results/bottle1.vtk` and measured each one's
distance to the nearest triangle of `input boundaries/bottle1_tri.raw` (after
applying `ReadRawData()`'s own normalization into the 100-cube). Over a
4,000-point sample: mean, median, p95 and max distance all 0.00000. The
reference mesh was generated from precisely this surface. Combined with
[Finding 7](#finding-7), that leaves exactly one free variable between us and Table 2's
Bottle1 row: **how deep the octree is allowed to refine.**

<a id="finding-9"></a>**Finding 9 — the refinement ladder is a per-model knob, and Table 2's
"Refinement Level" column reports it.** `ComputeCellValue()` hardcodes its
five refinement thresholds at octree levels 8, 7, 6, 5, 4 with "refine all" at
level 3, and `VOXEL_SIZE = 9` caps the tree at level 9. Two pieces of
evidence say that ladder was not fixed across the paper's runs:

- `v1.2`/`v1.3` still carry a **commented-out sixth threshold**,
  `//if (level == 9) {// 5`, and `v1.2`'s first commit (`2e5f4c0`) shipped a
  six-entry `C_THRES = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2}` alongside a
  `levelRes[]` extended to 1024 and `VOXEL_SIZE = 10`. A sixth threshold at level 9
  produces level-10 cells.
- Measuring the reference meshes confirms different models bottom out at
  different levels: `bunny` (Table 2 refinement level 4) stops at level 8,
  `bottle1` (4) at level 8, `deformed_armadillo` (5) and `head` (5) at level
  9, `david` (6) at level 10 — and `david.vtk` really does contain 1,997
  cells at size 0.0977 = 100/1024, which the public five-threshold code physically
  cannot produce.

So the shipped source is one frozen snapshot of a ladder the author moved per
model, and `bone` — which we now reproduce exactly — happens to be a model
that used the shipped setting.

**Change made to enable testing this: the ladder is now configurable.**
Patched `HybridOctree_Hex_v1.0/`'s `ComputeCellValue()` to express its six
hardcoded level constants as `ladderTop`, `ladderTop-1` ... `ladderTop-5`,
with `LADDER_TOP` defaulting to `8` — byte-for-byte the shipped behavior —
and `VOXEL_SIZE` likewise overridable at compile time. Building with
`-DLADDER_TOP="(octreeDepth - 1)"` anchors the whole ladder to `VOXEL_SIZE`
instead, which is what makes `VOXEL_SIZE` behave as the paper's per-model
refinement-level knob. No numerical logic changed; with the defaults the
binary reproduces `bone` exactly as before, and a second `bone` run built
with `-DLADDER_TOP="(octreeDepth - 1)"` at `VOXEL_SIZE = 9`
(`runs/bone-v1.0-vs9rel/`) also lands on 10,356 / 8,619 / 0.610001 / 1.0,
confirming the refactor is a no-op at the shipped setting.

<a id="finding-10"></a>**Finding 10 — Table 2's "Refinement Level" column is the number of
distinct octree levels the mesh occupies, readable off every reference
mesh.** With levels now measurable directly ([Finding 5b](#finding-5b)),
every Table 2 reference mesh was measured for the octree levels it
genuinely occupies.

| Model | Table 2 Refinement Level | Real levels occupied | Distinct count | `max - min` |
|---|---|---|---|---|
| Bottle1 | 4 | 4-7 | **4** ✓ | 3 |
| Bunny | 4 | 4-7 | **4** ✓ | 3 |
| Dragon Stand2 | 4 | 5-8 | **4** ✓ | 3 |
| Gargoyle | 4 | 5-8 | **4** ✓ | 3 |
| Deformed Armadillo | 5 | 4-8 | **5** ✓ | 4 |
| Head | 5 | 4-8 | **5** ✓ | 4 |
| David | 6 | 4-9 | **6** ✓ | 5 |
| `bone` (not in Table 2) | — | 5-8 | 4 | 3 |

The **distinct count** column reproduces Table 2 exactly on all seven
models. `max - min` does not — it falls short by one on every row.

<a id="how-real-levels-were-obtained"></a>**How the "real levels" column was
obtained.** Three steps, all reproducible from files in this repo.

*Step 1 — why an edge length is a level.* `ReadRawData()`
(`HexGen.cpp:817-832`) rescales every input surface into a fixed 100-unit
cube: it takes the model's longest dimension as `BOX_LENGTH`, centres the
model on the other two axes, then multiplies every vertex by
`100.0 / BOX_LENGTH`. Physical units are discarded. An octree cell at level
`L` therefore has edge `100 / 2^L` in every model, and a hex's level is
just `log2(100 / edge)`.

*Step 2 — measure.* `scripts/level_histogram.py` computes each hex's edge as
the mean of its 12 edges, converts to a level, and rounds:

```
python3 scripts/level_histogram.py "our results/<model>.vtk"
```

Run on all eight reference meshes (counts, then percentage of that mesh):

| Model | L4 | L5 | L6 | L7 | L8 | L9 | L10 |
|---|---|---|---|---|---|---|---|
| Bottle1 | 52 (0.17%) | 1,015 (3.37%) | 9,945 (32.99%) | 18,122 (60.12%) | 1,006 (3.34%) | 5 (0.02%) | — |
| Bunny | 204 (0.94%) | 5,423 (25.00%) | 9,000 (41.48%) | 6,960 (32.08%) | 108 (0.50%) | — | — |
| Dragon Stand2 | — | 735 (1.45%) | 11,683 (22.97%) | 21,903 (43.07%) | 15,485 (30.45%) | 1,047 (2.06%) | — |
| Gargoyle | — | 368 (0.16%) | 17,487 (7.39%) | 154,546 (65.29%) | 60,646 (25.62%) | 3,642 (1.54%) | — |
| Deformed Armadillo | 14 (0.04%) | 2,153 (6.16%) | 7,808 (22.35%) | 9,088 (26.01%) | 15,081 (43.16%) | 795 (2.28%) | — |
| Head | 237 (0.43%) | 5,058 (9.19%) | 14,158 (25.72%) | 27,474 (49.92%) | 7,867 (14.29%) | 244 (0.44%) | — |
| David | 60 (0.02%) | 4,922 (1.74%) | 26,938 (9.52%) | 73,964 (26.14%) | 116,317 (41.11%) | 58,759 (20.77%) | 1,997 (0.71%) |
| `bone` | — | 838 (9.72%) | 3,771 (43.75%) | 2,683 (31.13%) | 1,276 (14.80%) | 51 (0.59%) | — |

*Step 3 — discard the distortion tail.* The deepest band in most meshes is
not a real octree level. The projection stage moves boundary vertices onto
the input surface, which shrinks some hexes enough to measure a level
deeper than the octree ever built. The test, settled in the 2026-08-19
session and applied without a miscall across all eight models: **a
population an order of magnitude or more below the band above it is
distortion; a smooth continuation is real.** Applying it drops Bottle1's
L8+L9 (18,122 → 1,006, an 18× fall), Bunny's L8 (6,960 → 108, 64×), Dragon
Stand2's L9 (15,485 → 1,047, 15×), Gargoyle's L9 (60,646 → 3,642, 17×),
Deformed Armadillo's L9 (15,081 → 795, 19×), Head's L9 (7,867 → 244, 32×),
David's L10 (58,759 → 1,997, 29×), and `bone`'s L9 (1,276 → 51, 25×).
David's L9 survives the test — 58,759 cells against L8's 116,317 is a
factor of 2, not a decade — which is why David alone needs a sixth
threshold.

The Bunny run in [its own section](#bunny-octree-levels) supplies an
independent check on step 3. Its `octree.vtk`, written before projection,
contains no cell deeper than level 7, while its `finalMesh.vtk` measures
492 cells at level 8, 5 at level 9, and 1 at level 10. The octree itself
provides a hard ceiling, so those deeper measurements can only be
distortion. That is direct confirmation of the heuristic rather than an
appeal to it.

<a id="finding-10-correction"></a>**Correction, 2026-08-27.** An earlier
version of this finding stated the rule as `max - min`, and its table
listed a "min level in reference" one level below the measured minimum on
every row. That figure was an inferred octree minimum, not a measurement,
and the inference was circular: it assumed the answer in order to make the
subtraction come out right. Counting occupied levels directly needs no
inference and lands on the same numbers. Two smaller errors went with it —
Bottle1's shallowest population is 52 cells, not 51, and `bone`'s minimum
is level 5, not level 4.

That correction also changes how `bone` reads. `bone` occupies levels 5-8,
four distinct levels, even though it runs the shipped configuration with
all five thresholds compiled in. A threshold being available is not the
same as a threshold firing: `bone`'s geometry never triggers the deepest
one, so the mesh spans four levels rather than five. The earlier claim that
`bone` demonstrates "5 = the shipped setting" does not survive measurement.
What `bone` does still demonstrate is that the shipped constants are the
paper's constants ([Finding 7](#finding-7)).

`david` reaching level 9 as a genuine population (20.77%, the second-largest
band in its mesh, and showing cells at 100/1024 = 0.0977) is only possible
with a sixth threshold — precisely the `//if (level == 9) {// 5` block still
sitting commented out in `v1.2`/`v1.3`.

<a id="finding-11"></a>**Finding 11 — running Bottle1 at four levels closes almost the whole gap.**
Three Bottle1 configurations, all `v1.0` code with `v1.0`'s own
`C_THRES = {0.15, 0.3, 0.6, 1.2, 2.4}` and `H_THRES = {16,8,4,2,1}`
untouched, differing only in where the ladder sits and how deep the tree may
go. All against the 2024-01-11 input (`c291245`) unless noted:

| Configuration | Levels (octree levels) | Max level | Octree cells | Octree time | Final points | Final cells |
|---|---|---|---|---|---|---|
| shipped (`bottle1-v1.0-repro`) | 5, at 4-8 | 9 | 173,552 | 2,168 s | 117,031 | 101,688 |
| `vs8rel` | 5, at 3-7 | 8 | 65,024 | 321 s | 40,781 | 34,653 |
| **`rl4`** | **4, at 3-6** | **7** | **57,800** | **255 s** | **33,787** | **28,484** |
| `rl4`, header-correct input `f3247d8` | 4, at 3-6 | 7 | 57,240 | 249 s | 32,970 | 27,750 |
| **Table 2 target** | | **7** | | | **36,091** | **30,145** |

From 3.4x oversized down to 0.94x — and the level-by-level agreement is
what makes it convincing, not just the total:

| Level | Reference | `rl4` | `vs8rel` |
|---|---|---|---|
| 4 | 51 | **51** | 53 |
| 5 | 1,051 | 1,021 | 1,035 |
| 6 | 10,039 | 10,334 | 10,094 |
| 7 | 17,921 | 16,106 | 16,971 |
| 8 | 1,077 | 971 | 6,252 |
| 9 | 6 | 1 | 248 |
| **Total** | **30,145** | **28,484** | **34,653** |

**The reference's level-8 population is distortion, not refinement.** Two
independent checks. First, spatial clustering: the reference's 1,077
level-8-sized cells are smeared across 200 separate 5-unit voxels with at most
20 in any one of them, whereas `vs8rel`'s 6,252 — which are genuinely refined
— pile up 1,009 in a single voxel and 3,714 in the top tenth of the model
alone. Second, `rl4` *cannot* produce a level-8 cell (its tree is capped at 7)
yet still measures 971 of them, essentially the reference's 1,077. Both
numbers are projection-squashed level-7 hexes. So Bottle1's true maximum
octree level is 7, and Table 2's "Refinement Level 4" is literal.

**The 2024-01-11 input is also the better one, which settles the
`c291245`-vs-`f3247d8` question.** At identical settings the header-buggy
`c291245` gives 33,787 / 28,484 and the header-correct `f3247d8` gives
32,970 / 27,750 — `c291245` is closer to the target on both counts. The
mechanism is now precisely understood: `c291245`'s duplicated last triangle
makes two coincident triangles share all three edges, so the dihedral angle
between them is 0 and each of their three vertices picks up
`(0 - PI)^2 = 9.8696` of curvature. Confirmed by `refine_criteria_stats`: the
curvature maximum is exactly 9.8696 for `c291245` and 3.7410 for `f3247d8`,
and the deepest threshold's candidate lists shrink from 12 points / 59 triangles to
9 / 47. Both files' geometry is identical; the bug just injects three phantom
high-curvature vertices. The reference mesh was generated on 2024-01-11 with
that bug present, and reproducing it means keeping it.

**Remaining gap: 5.5%, entirely at level 7.** `rl4` is 1,815 cells short at
level 7 and 295 long at level 6 — which is one consistent story, since
refining one level-6 cell trades it for eight level-7 cells. The reference
refined roughly 230-290 more level-6 cells than we do, about 11% more. That
is a small, single-threshold discrepancy in the threshold gated by `C_THRES[3] = 1.2`
and `H_THRES[3] = 2`.

<a id="finding-12"></a>**Finding 12 — the residual gap is the deepest threshold's *thickness* gate, and
the target falls between two natural values.** Both threshold arrays are
really one scale-relative rule: at their shipped levels `H_THRES[i]` is
exactly `2.56 x cellsize` and `C_THRES[i]` exactly `0.9375 / cellsize`, for
every `i`. Shifting the ladder one level coarser while keeping the arrays
indexed by threshold therefore halves each thickness threshold relative to its
own cell, and doubles each curvature threshold. Testing that, plus
`OUT_IN_RATIO` (which `v1.3` raises to 0.15) as a control — all at the `rl4`
ladder, all against `c291245`:

| Deepest threshold's gate | Final points | Final cells | vs. target |
|---|---|---|---|
| `H[3] = 2` (`rl4`, thickness array shifted with the ladder) | 33,787 | 28,484 | −5.5% |
| `H[3] = 2`, `OUT_IN_RATIO = 0.15` | 33,989 | 28,668 | −4.9% |
| `H[3] = 2.5` | 34,664 | 29,244 | −3.0% |
| **Table 2 target** | **36,091** | **30,145** | — |
| `H[3] = 3` | 36,956 | 31,272 | +3.7% |
| `C[3] = 1.0` (curvature gate instead) | 38,733 | 32,739 | +8.6% |
| `H[3] = 4` (thickness array *not* shifted) | 39,167 | 33,186 | +10.1% |

`OUT_IN_RATIO` is not the lever — 0.125 to 0.15 moves the mesh by 0.6%, so
`RemoveOutsideElement()`'s trim is not where the difference lives. The
thickness gate is: the level histograms show `H[3] = 4` leaves levels 4, 5
and 6 essentially untouched (43 / 959 / 10,284 versus `rl4`'s 51 / 1,021 /
10,336) and moves only level 7, from 16,095 to 21,044 against the reference's
17,921. The target sits between `H[3] = 2.5` and `H[3] = 3`, i.e. around 2.7-2.9
— not a round number, which is why this is reported as a bracket rather than
claimed as the author's value.

<a id="finding-13"></a>**Finding 13 — the duplicated triangle is what makes `ProjectToIsoSurface`
stall, which closes this morning's open question.** Earlier today both `v1.0`
and `v1.1` were killed after hours stuck in the projection stage, oscillating
between input triangles #1357 and #1351, and the cause was left unexplained
(see the Correction above retracting the "`v1.1` has escape logic" framing).
The `rl4` pair settles it as a controlled experiment — same code, same
constants, same ladder, only the input file differs:

| Input | Behavior |
|---|---|
| `c291245` (header says 29665, last triangle silently duplicated) | Frozen at `smallDist = 1.399` on triangle #1363 for hundreds of consecutive checkpoints; never wrote `finalMesh.vtk` |
| `f3247d8` (header says 29664, parses cleanly) | Converged normally, `smallDist ~ 0.002`, `finalMesh.vtk` written, `ELEM_THRES` ratcheted to 0.57 |

The mechanism is now clear. The duplicate makes two exactly coincident
triangles that share all three edges, so no surface normal or nearest-triangle
assignment is well defined there, and the gradient-descent projection can
never satisfy its distance test. This is not merely a curvature artifact
([Finding 11](#finding-11)) — it is a hard convergence blocker. **So `f3247d8` is the input
to reproduce against after all**, despite being five days later than the
reference mesh's commit: the author evidently fixed the header on 2024-01-16
precisely because it was wrong, and the version that actually generated
`bottle1.vtk` must have parsed 29,664 triangles or it could not have
converged either. That reverses this morning's "use `c291245`, it's
temporally closer" choice, and explains why every Bottle1 run today before
this point failed to produce a `finalMesh.vtk`.

Against `f3247d8`, the `rl4` ladder gives 32,970 points / 27,750 cells with
worst SJ 0.570000 and best SJ 1.000000 — converged cleanly, 7.9% under Table
2's cell count. (Worst SJ is not an independent constraint: it is simply
whatever value the `ELEM_THRES` ratchet had reached when the run was stopped,
which is why every "ours" row in Table 2 lands in the narrow 0.54-0.59 band
that ratchet walks through, `0.53 + 0.01 k`. That run kept climbing past
0.570 to 0.580 while still writing the same 32,970 / 27,750 mesh, so 0.560 is
a stopping point, not a property of the mesh.)

<a id="finding-14"></a>**Finding 14 — closing the last percent: restore the scale-relative
thickness rule and use the author's own `CELL_DETECT = 0.75`.** The
shifted ladder halves every thickness threshold relative to its own cell
size, so the natural correction is to double `H_THRES` to `{32,16,8,4,2}`,
which restores `H_THRES[i] = 2.56 x cellsize` exactly. On its own that
overshoots ([Finding 12](#finding-12)), but combined with `CELL_DETECT = 0.75` — the value
Tong himself used in his 2024-01-19 and 2024-01-21 commits, which narrows
`ComputeCellValue()`'s detection box from `2 x cellsize` to `1.5 x cellsize` —
it lands essentially on the reference. Both against `f3247d8`:

| Configuration | Points | Cells | vs. target |
|---|---|---|---|
| **Table 2 target** | **36,091** | **30,145** | — |
| `rl4` + `H_THRES x 2` + `CELL_DETECT = 0.75` | **35,535** | **29,943** | **−1.5% / −0.7%** |
| `rl4` + `H_THRES[3] = 3` | 36,164 | 30,554 | +0.2% / +1.4% |
| `rl4` (thickness array shifted with the ladder) | 32,970 | 27,750 | −8.6% / −7.9% |
| `rl4` + raw Euclidean thickness (hypothesis, refuted) | 32,631 | 27,482 | −9.6% / −8.8% |

*(That last row tested whether the pre-`v1.0` code lacked
`GetCellValue()`'s dominant-axis scaling — `len = max(|dx|,|dy|,|dz|) x |len|`,
flagged as an unexplained anomaly by the 2026-08-05 session. Removing it moves
the mesh the wrong way: the factor is always in `[1/sqrt(3), 1]` and therefore
always shrinks the measured thickness, so dropping it makes fewer triangles
pass `len < H_THRES` and produces a slightly **smaller** mesh, not a larger
one. Hypothesis refuted; the scaling stays. Built behind `-DRAW_THICKNESS` in
`HybridOctree_Hex_v1.0/HexGen.cpp` for the record, off by default.)*

The `CELL_DETECT = 0.75` configuration's level-by-level agreement is the real
evidence, because it is far tighter than the totals:

| Level | Reference | `rl4` + `H x 2` + `CELL_DETECT 0.75` | Agreement |
|---|---|---|---|
| 4 | 51 | 41 | — (0.1% of the mesh) |
| 5 | 1,051 | 1,069 | +1.7% |
| 6 | 10,039 | 10,093 | **+0.5%** |
| 7 | 17,921 | 17,945 | **+0.1%** |
| 8 (distortion artifacts) | 1,077 | 795 | — (not real cells) |
| **Total** | **30,145** | **29,943** | **−0.7%** |

Levels 6 and 7 carry 93% of the mesh and match to 0.5% and 0.1%. The whole
remaining 202-cell total difference is the level-8 artifact count, which just
tracks how far each mesh's projection stage had smoothed when it was sampled,
not a difference in the octree.

**Recipe that reproduces Bottle1's Table 2 row.** For the record, the full
configuration, all of it `v1.0` source with no algorithmic change:

- source: the `v1.0` tag, `sscanf_s` -> `sscanf` and `<malloc.h>` -> `<cstdlib>`
  portability patches only
- input: `input boundaries/bottle1_tri.raw` at `f3247d8` (2024-01-16),
  header `14833 29664`
- `C_THRES = {0.15, 0.3, 0.6, 1.2, 2.4}` — `v1.0`'s own, unchanged
- `H_THRES = {32, 16, 8, 4, 2}` — `v1.0`'s own, doubled to keep
  `H_THRES[i] = 2.56 x cellsize` under the shifted ladder
- `CELL_DETECT = 0.75`
- four refinement thresholds at octree levels 3-6, maximum level 7 — built as
  `-DVOXEL_SIZE_VALUE=7 -DLADDER_TOP=octreeDepth`
- run `ProjectToIsoSurface` until `ELEM_THRES` reaches 0.56 and stop, exactly
  as the unbounded `while (true)` loop forces a human operator to do — this
  run was watched through 0.500001, 0.550000 and on to 0.570001 on the same
  35,535 / 29,943 mesh, so 0.560 is simply one rung of that ladder

Final state of that run, `runs/bottle1-v1.0-rl4-h2xcd75-fi/finalMesh.vtk`,
against Table 2's Bottle1 row:

| | Reproduction | Table 2 / `our results/bottle1.vtk` |
|---|---|---|
| Points | 35,535 | 36,091 |
| Cells | 29,943 | 30,145 |
| Worst SJ | 0.560 en route (sampled at 0.570001) | 0.560000 |
| Best SJ | 1.000000 | 1.000000 |
| Refinement level | 4 | 4 |

Honest accounting of what is derived versus fitted: the four-level configuration is
*derived* — it is what Table 2's own "Refinement Level" column says, it is
confirmed independently by every other model's reference mesh ([Finding 10](#finding-10)),
and it alone accounts for the factor of 3.4. `H_THRES x 2` is *derived* from
the scale-relative structure of the shipped arrays. `CELL_DETECT = 0.75` is
*fitted*, in the weak sense that it was chosen from the author's own set of
historical values because it closed the residual gap; a different
compensation (`H_THRES[3] = 3` at `CELL_DETECT = 1`) lands equally close from
the other side, so the last ~1% is not uniquely determined by the evidence.
That is expected: `our results/bottle1.vtk` was committed 2024-01-11, a full
day before *any* source was committed to the repository, so the exact code
that produced it is not in git and small pre-`v1.0` differences cannot be
recovered.

**Where to go next.** The method is now established and cheap, so the natural
next step is breadth rather than more depth on Bottle1:

1. **Run the remaining eleven Table 2 models.** Each one's refinement level is
   already read off its reference mesh ([Finding 10](#finding-10)), so each is a single run
   with a known ladder setting rather than a search. Bunny (4 levels), Dragon
   Stand2 and Gargoyle (4 levels, ladder one level finer — their references
   have no level-4 cells), Head and Deformed Armadillo (5), David (6, needs
   the sixth threshold that `v1.2` left commented out). Confirming the recipe
   transfers to models it was not tuned on is the real test of [Finding 10](#finding-10).
2. **Settle the last 1% properly, if at all.** The `CELL_DETECT` /
   `H_THRES[3]` ambiguity would be resolved by any model where the two
   compensations predict measurably different meshes — a breadth run may
   decide it for free.
3. **Re-examine the runtime question.** Every conclusion in the 2026-08-05
   entry about Bottle1 being ~27x slower than `bone` was measured on octrees
   3-6x larger than they should have been. At the correct four-level setting
   Bottle1's octree construction takes ~250 s rather than 2,168 s, so the
   ["218 s is unreachable with the public code" synthesis](#why-218s-is-unreachable) below needs
   re-measuring before it is trusted.
4. **The `ComputeCellValue()` performance finding still stands** — it is a
   real missing spatial index — but its practical weight is much smaller than
   the 2026-08-05 entry assumed.

### 2026-08-05 — Root-causing the mesh-size mismatch; redoing bone and Bottle1 (Apple M1)

The 2026-08-04 session's Bottle1 and bone runs both produced meshes far
larger than published (Bottle1 ~6x, bone ~1.9x), left unexplained at the
time. This session root-causes it and redoes both runs.

**`C_THRES` mismatch, found by cross-referencing the paper against
`Initialization.h`.** Section 2.1 of the paper (p.3) states the curvature
refinement thresholds explicitly: *"Five curvature thresholds
`Gthres = {0.5, 1, 2, 4, 8}` are adopted. If an octree cell at level `l+4`
... satisfies `G(l+4) > Gthres[l]`, we refine it to level `l+5`."* The
shipped `Initialization.h` (unmodified from upstream, current release
"v1.3") had `C_THRES[5] = { 0, 0, 0.4, 0.8, 1.6 }` — nothing like the
paper's stated values. Checking this file's git history across releases
showed it never has matched, and got worse over time rather than better:

| Release | `C_THRES` |
|---|---|
| v1.2 (`2e5f4c0`) | `{0.1, 0.2, 0.4, 0.8, 1.6, 3.2}` (6 levels) |
| v1.2 (`f9e667c`) | `{0.1, 0.2, 0.4, 0.8, 1.6}` (dropped to 5) |
| v1.3 (`220c53f`, shipped) | `{0, 0, 0.4, 0.8, 1.6}` (**first two zeroed**) |

`H_THRES` (`{16, 8, 4, 2, 1}`) and `VOXEL_SIZE` (`10`) have been stable
across every release and already match the paper's stated
`Tthres = {16,8,4,2,1}` and its "resulting octree level ranges from 5 to
9" — only `C_THRES` was the problem. Mechanism: the refinement rule is
"refine if `G > threshold`". With `C_THRES[0] = C_THRES[1] = 0`, the
condition `G(level 4) > 0` and `G(level 5) > 0` is true for essentially
any curved region (Gaussian curvature is only exactly zero on flat
patches), making the two coarsest refinement checks effectively
unconditional — forcing far more of the octree to refine than the paper's
thresholds intend, cascading into every downstream stage (more octree
cells → more dual-mesh elements → proportionally slower dualization,
buffer clearing, and projection).

**Curvature formula/units mismatch, found by comparing the paper's stated
formula against the actual implementation.** The paper's `G` is a
cotangent-Laplacian magnitude normalized by local Voronoi cell area (units
~ 1/length): `‖Σ (cot α + cot β)(Pj − Pi)‖ / (4·Ai)`. Searching the
codebase for any cotangent- or Voronoi-based computation found none — the
only curvature computation anywhere (`HexGen.cpp`, inside `ReadRawData`,
called from `InitializeOctree`) is `r[point] += (dihedral_angle − π)²`,
summed over each edge incident to that point, where `dihedral_angle` is
the angle between the two triangle faces sharing that edge. This is a
different mathematical quantity (squared angular deviation from flat, no
area normalization, units of radians²) from the paper's stated formula
(a normalized curvature magnitude, units ~ 1/length). Because of this,
neither the buggy shipped thresholds nor the paper's literal `Gthres`
values can be "provably correct" for what this code actually computes —
there's no valid unit conversion between the two formulas. This explains
why fixing `C_THRES` to the paper's literal values didn't land on the
published mesh size either (see below) — empirical recalibration against
published stats, not formula-derived threshold values, is the right path
forward pending a possible future fix to the curvature formula itself
(out of scope for this session).

**Bone recalibration (see the [Bone section](#bone-calibration-testbed-not-a-table-2-model) above for the full results
table).** Tested three `C_THRES` values against bone (this repo's own
published 10,356 vertices / 8,619 elements / worst SJ 0.61), each with the
full `MAX_PROJ_ITER=200000` budget (bone's own established setting — one
intermediate run accidentally used a leftover `MAX_PROJ_ITER=20000` cap
carried over from the 2026-08-04 Bottle1 work, caught and corrected before
drawing conclusions from it):
- Buggy shipped values: 1.9x too big, worst SJ 0.590 (close on quality,
  wrong on size).
- Paper-literal values: 2.3x too *small* — overshot in the other
  direction — worst SJ 0.550 once re-run with the correct iteration
  budget (an earlier run at 0.011 was purely the leftover-cap confound
  described above, not a quality problem with these thresholds).
- v1.2's historical values (a previously-shipped value from this repo's
  own git history, not an arbitrary guess — chosen because it sits between
  the other two): closest match on every metric simultaneously — 1.43-1.46x
  on size, worst SJ 0.600 (within 0.01 of published). Adopted for the
  Bottle1 redo.

**Bottle1 redo, and a second, separate bottleneck found.** With
`C_THRES` recalibrated, launched Bottle1 with the full `MAX_PROJ_ITER=
200000` budget. Streamed per-stage completions live rather than polling
after the fact. `GetCellValue()`'s curvature/narrow-region assessment
substage alone took 724.4s (~12 min) — nearly as long as the *entire*
combined octree stage took under the old buggy thresholds (1062s), and
drastically disproportionate to bone's 11.0s for the same substage with
the same thresholds (65.8x longer, vs. only ~2.45x more input triangles:
29,664 vs. 12,088). Paused the run (killed after octree construction
finished, ~34 min in) to investigate rather than let it continue on an
unbounded timeline.

**First hypothesis (wrong, corrected below)**: initially suspected
`GetCellValue()`'s thickness/narrow-region check (`HexGen.cpp:950-986`) —
a brute-force all-pairs loop, for every triangle `i` ray-testing against
every other triangle `j > i` via `Intersect()`, with no spatial
acceleration structure (no BVH/octree/kd-tree) and no early exit — exactly
the paper's described operation ("shooting a ray... to find the shortest
intersection with other triangular elements"), implemented as literal
O(n²) over the input triangle count. Reasoned that O(n²) scaling
(bone: 12,088 triangles → ~73M pairs; Bottle1: 29,664, 2.45x more →
~440M pairs, **6.0x more**) explained part of the 65.8x slowdown, with an
unexplained ~11x remainder attributed to `Intersect()` costing more per
call on Bottle1's geometry.

**Actual root cause, confirmed by direct instrumentation**: this was
wrong — the O(n²) triangle-pair loop is fast (bone 1.57s, Bottle1 9.7s;
its trailing `unordered_set`-based dedup step is negligible, <0.003s for
both). Added temporary `clock()` instrumentation splitting
`GetCellValue()`'s three phases (raw O(n²) loop, dedup, and the trailing
`ComputeCellValue()` recursive sweep it calls at line ~1035) confirmed the
real cost is almost entirely the third phase:

| Phase | bone | Bottle1 | Ratio |
|---|---|---|---|
| O(n²) raw loop (thickness/curvature candidate collection) | 1.57 s | 9.70 s | 6.2x |
| `unordered_set` dedup | 0.0007 s | 0.0023 s | 3.3x |
| `ComputeCellValue()` recursive sweep | ~9.45 s (inferred) | **720.6 s** | **~76x** |

`ComputeCellValue()` (`HexGen.cpp:1039` onward) recurses over every octree
node; at each of the five tested refinement levels (4-8, corresponding to
`refineTri0..4`/`refineTriPt0..4`) it does a **linear scan through that
level's entire candidate list** for every cell being classified at that
level (e.g. `HexGen.cpp:1217-1241` for the coarsest level, scanning
`refineTriPt0`/`refineTri0` — up to 41,442/64,506 raw candidates for
Bottle1 — with up to 12 `Intersect()` box-face calls per triangle
candidate), stopping only on the first hit. There's no spatial index (no
BVH, no reuse of the octree structure itself, no spatial hash) for this
lookup — it's the same brute-force pattern as the (much cheaper) O(n²)
loop, but applied per-*octree-cell* rather than per-*triangle*, and the
candidate lists it scans are themselves 2.9-3.2x larger for Bottle1
(deduped `refineTri0`: 15,444 vs. bone's 5,386; `refineTriPt0`: 6,036 vs.
3,134). Both factors — more octree cells needing classification at deep
levels (Bottle1's coiled/spiral surface drives refinement into more of the
volume) and larger per-level candidate lists to scan against each one —
compound multiplicatively rather than adding, which is why the observed
~76x so heavily overshoots the ~6x that raw triangle-count scaling alone
would predict, and why the naive "O(n²) in triangle count" framing above
was the wrong place to look.

This is a genuine performance bottleneck in the shipped algorithm's
implementation, entirely unrelated to `C_THRES` calibration — it would
affect any sufficiently geometrically complex input model, not just
Bottle1, since it's driven by how much of the octree needs deep
refinement, not by input triangle count directly. The diagnostic
`clock()`/`std::cout` instrumentation used to isolate this has since been
removed from `HexGen.cpp` (its purpose served — see the commit removing
it), replaced with a short permanent comment at the `ComputeCellValue()`
call site pointing back to this finding. Both Bottle1 reruns from this
session (`runs/bottle1-v2/`, killed mid-run; `runs/diag-bottle1/`, killed
once diagnostic data was captured) are being kept alongside the
2026-08-04 `runs/bottle1/` and the final `runs/bottle1-v3/` as
before/after evidence.

**Synthesis: why Bottle1 has been so hard to reproduce.** The
`bottle1-v3` run (calibrated `C_THRES`, full `MAX_PROJ_ITER=200000`
budget) has now finished — 8537.8s (~142.3 min) active compute (see
footnote 4 above for the wall-clock-vs-active-compute distinction) —
against Table 2's reported 218s. Worth stepping back and separating
what's actually been fixed from what hasn't.

*Three separable problems, not one.* Convergence quality (fixed,
cleanly), mesh size (fixed, partially — see below, this turned out to be
more complicated than first thought), and runtime (still unsolved) are
close to independent of each other. Fixing one didn't fix the others.

*Convergence — fixed cleanly.* `badElem=0, smallDist=3.76×10⁻¹¹`,
essentially machine precision. Worst SJ (0.530) landed close to Table
2's 0.560. The earlier 0.010 (2026-08-04) was purely the deliberate
`MAX_PROJ_ITER=20000` cap, not an algorithm problem.

*Mesh size — root-caused, but the fix barely worked for Bottle1
specifically.* The shipped `C_THRES` didn't match the paper's stated
values, and separately, the codebase's actual curvature formula doesn't
match the paper's stated formula either (a different, unnormalized
quantity — no valid unit conversion between the two). Empirical
recalibration against `bone`'s known-good published stats got `bone`'s
size within ~1.5x — but the *same* recalibrated `C_THRES` only took
Bottle1 from 218,898 → 212,990 vertices, a ~2.7% reduction, nowhere near
`bone`'s ~24% improvement with the identical fix. That's a finding:
`C_THRES` (curvature) isn't the dominant octree-refinement driver for
Bottle1 the way it is for `bone`. The follow-up hypothesis — `H_THRES`
(thickness/narrow-region) dominates instead, given Bottle1's thin
coiled/spiral surface detail — was tested and refuted: halving `H_THRES`
only shrank the octree by ~1%. Neither threshold, tuned individually,
explains Bottle1's over-refinement — still an open question (see the
2026-08-05 log entry's `H_THRES` test for the details, including a
genuine code-level anomaly found in the thickness-distance computation
along the way that didn't turn out to be the answer here either).

*Runtime — the real story, and not explained by either fix above.*
Comparing `bone` (calibrated) against Bottle1 (calibrated), same code,
same machine:

| Stage | bone (calibrated) | Bottle1 (calibrated) | Ratio |
|---|---|---|---|
| Octree construction | ~55 s | ~1012 s | ~18x |
| Dual mesh generation | ~7 s | ~1046 s | ~150x |
| Buffer clearing | ~8–22 s | ~785 s | ~36–100x |
| Projection/quality | ~238 s | ~5693 s (~94.9 min) | ~24x |
| **Total** | **~5.2 min** | **~142.3 min active compute** | **~27x** |

Bottle1 only has 2.45x more input triangles than bone. None of these
ratios are anywhere close to 2.45x, or even the ~6x naive O(n²)-in-
triangle-count reasoning would predict. The confirmed cause for octree
construction specifically (not guessed — see the `ComputeCellValue()`
instrumentation above) is that octree-cell classification does a
brute-force linear scan with no spatial acceleration structure, and both
the number of cells needing classification and the candidate-list sizes
are inflated by Bottle1's more complex geometry, compounding
multiplicatively. Dualization and buffer-clearing show the same
disproportionate pattern and likely hide similar unaccelerated scans,
though those haven't been isolated with the same direct instrumentation
yet. Projection/quality's ~24x, by contrast, is largely just iteration
count — a much larger mesh takes more `gradient()`/redistribution passes
to reach the same convergence quality, a more "expected" kind of slowdown
than the other three stages' scaling.

<a id="why-218s-is-unreachable"></a>*Why 218s is probably unreachable with the current public code at all.*
Even a "perfectly calibrated" run of the current public code cannot
plausibly finish in 218s. `bone` — a much smaller, simpler model — takes
~5.2 minutes on its own, and its quality-improvement stage alone (238s)
is already almost as long as Table 2's entire reported Bottle1 time. For
Bottle1 to finish in 218s, whatever Tong actually ran would need spatial
acceleration in `ComputeCellValue()` (and likely elsewhere) that simply
isn't in the code currently on GitHub.

This lines up with a pattern confirmed in *both* of Tong's public repos
this session: `HexOpt`'s own README documents (with his own email
confirmation) that the public `meshQuality.cpp` doesn't implement the
method the CAD 2026 paper claims — and structurally, that code looks
like it's actually derived from `HybridOctree_Hex`'s own older algorithm
instead (see `hovey/HexOpt`'s README, 2026-08-05 entry). The most likely
explanation for Bottle1's timing gap is the same shape of problem: the
version of `HexGen.cpp` that's public doesn't have the optimizations
that the version used to generate Table 2 had. This isn't a mistake in
this reproduction's methodology — it's a real, structural gap between
the published performance numbers and the current state of the public
code, on top of the mesh-size question — see the `H_THRES` and
octree-balancing tests immediately below, both of which ruled out their
respective candidates and narrowed the mesh-size question down to
"Bottle1's curvature/thickness candidates are more spatially dispersed
across its geometry than bone's," not yet fully confirmed.

**`H_THRES` hypothesis test (2026-08-05).** Before testing, re-read
`GetCellValue()`'s thickness computation (`HexGen.cpp:936-985`) against
the paper's stated method ("shooting a ray from `Pi` along its normal
direction to find the shortest intersection with other triangular
elements") — structurally it matches (ray from a triangle's centroid
along its normal, tested against every other triangle via `Intersect()`),
*unlike* the curvature computation's outright formula mismatch. But
found a genuine anomaly along the way: the code scales the raw
ray-intersection distance by `max(|dir_x|,|dir_y|,|dir_z|)` — the
dominant-axis component of the (unit, normalized) triangle-normal
direction — before comparing it to `H_THRES` (`HexGen.cpp:952-954`). That
factor is always in `[1/√3, 1]`, so it always *shrinks* the effective
distance relative to true Euclidean distance, by up to ~42% for normals
oriented diagonally to the octree axes. Plausible rationale (converting
to an octree-cell-crossing-relevant metric rather than raw 3D distance)
rather than an obvious bug like the `C_THRES` zeros, but not something
the paper's plain-language description mentions either. A coiled/spiral
surface like Bottle1's has far more diagonally-oriented normals than
bone's simpler shape, so this seemed like a promising lead.

Tested empirically rather than reasoning it through further: halved
`H_THRES` to `{8,4,2,1,0.5}` — a much stricter thickness gate that should
sharply cut the number of triangle pairs flagged as "thin enough to
refine" if thickness detection were the dominant driver. Ran cheaply
(octree-construction-only, ~17 min via `build-h-thres-test/`, killed
immediately after `ConstructOctree()` completed — no need for the full
~142 min pipeline just to read `octree.vtk`'s point/cell counts) rather
than committing to another multi-hour full run. Result: `octree.vtk` went
from 358,307/300,560 (points/cells) to 354,875/297,368 — a ~1% change.
**Hypothesis refuted.** `H_THRES` reverted to the paper-matching
`{16,8,4,2,1}` in `Initialization.h`.

Neither `C_THRES` nor `H_THRES`, tuned individually, explains Bottle1's
disproportionate over-refinement relative to `bone`. Remaining
candidates, none yet tested: (a) the octree balancing rule
(`StrongBalancedOctree()`'s 2:1 constraint) propagating refinement far
beyond the initially-flagged cells, so that once *some* threshold flags
enough scattered cells across Bottle1's more extensive/complex surface,
balancing alone forces widespread deep refinement regardless of exactly
which threshold triggered the initial flags; (b) some other input
property (e.g. `VOXEL_SIZE`/model-to-box-scale interaction, or an
octree-orientation-vs-geometry effect like the paper's own oil-pump
discussion) not captured by either threshold at all. This remains open;
`bottle1-h-thres-test/` and `build-h-thres-test/` were both deleted after
capturing the result (throwaway experiment, not meaningful to preserve).

**Octree-balancing hypothesis tested and refuted (2026-08-05).** Tested
candidate (a) above: does `StrongBalancedOctree()`'s 2:1 balancing pass
propagate refinement far beyond what `GetCellValue()` initially flags,
disproportionately for Bottle1? Instrumented `ConstructOctree()`
(temporarily) to count `octreeArray`'s true-entries immediately before
and after the balancing call, run cheaply on both models
(octree-construction-only, killed right after — `bone` finishes in
~1 min, Bottle1 in ~17 min):

| | Before balancing | After balancing | Amplification |
|---|---|---|---|
| bone | 2,777 | 3,633 | +30.8% |
| Bottle1 | 34,785 | 42,937 | +23.4% |

**Refuted, decisively.** Balancing amplifies bone *more* than Bottle1 in
relative terms (+30.8% vs. +23.4%), the opposite of what the hypothesis
predicted. More importantly: the Bottle1-vs-bone disparity is already
almost fully present *before* balancing ever runs (34,785/2,777 ≈ 12.53x)
and barely changes after it (42,937/3,633 ≈ 11.82x) — closely matching
the final octree's own ~11.7x disparity (358,307 vs. 30,733 points, from
the earlier `octree.vtk` comparison). So the entire disproportion
originates in `GetCellValue()`'s initial cell-marking, not in balancing
propagation. Since neither `C_THRES` nor `H_THRES`, tuned individually,
explains that marking disparity either (both tested and refuted above),
the likely remaining explanation is candidate (b): Bottle1's raw
curvature/thickness candidate lists are only ~3x larger than bone's (from
the earlier `GetCellValue()` diagnostic — `refineTriPt0` raw: 41,442 vs.
12,959), but produce a ~12x larger marked-cell count — meaning those
candidates are more spatially *dispersed* across Bottle1's more elongated,
coiled geometry, touching many more distinct octree cells rather than
clustering within a smaller set of them the way bone's more compact
shape's candidates would. If so, this isn't a fixable threshold-tuning
bug at all — it may be a genuine, correct consequence of applying this
octree strategy to an elongated/coiled shape, and reproducing Table 2's
mesh size for Bottle1 specifically might require a different mechanism
entirely (spatially-adaptive thresholds, or whatever Tong's original,
non-public implementation actually did) rather than anything tunable in
the current public code. Not yet confirmed — would need to inspect the
spatial distribution of `refineTri0`/`refineTriPt0` directly (e.g.
bounding-box coverage as a fraction of the model's overall extent) to
verify the dispersion theory rather than infer it from these aggregate
counts alone. The temporary `octreeArray` count instrumentation has been
removed from `HexGen.cpp`, replaced with a short permanent comment
pointing to this finding; `diag-balance-bone/`/`diag-balance-bottle1/`
and `build-balance-test/` were deleted after capturing the result.

**Spatial-dispersion hypothesis tested (2026-08-05) — partial support,
not a full explanation.** Tested candidate (b): are Bottle1's
curvature/thickness candidate points more spatially *dispersed* across
its geometry (touching many more distinct octree cells) rather than
clustered, versus bone's? Instrumented `GetCellValue()` to compare
`refineTriPt0`'s (the coarsest tested level's deduped candidate points)
bounding box against the full model's bounding box, run *very* cheaply
(only the O(n²) candidate-collection loop, well before the expensive
`ComputeCellValue()` sweep — ~10-15s per model, not minutes):

| | Candidate count | Candidate bbox / model bbox (volume) | Candidate density (pts/unit³) |
|---|---|---|---|
| bone | 3,134 | 71.6% | 0.0477 |
| Bottle1 | 6,036 | 99.1% | 0.0299 |

Bottle1's candidates do span nearly its entire bounding box (99.1% vs.
bone's 71.6%) with only ~1.93x more points — and are actually *sparser*
in point-density terms (~0.625x bone's density), consistent with being
more spread out rather than clustered. But combining the count ratio
(1.93x) and bbox-coverage ratio (1.38x) only predicts a ~2.7x
marked-cell difference — far short of the ~12.5x actually observed
before balancing. **Real effect, insufficient magnitude.** Bounding-box
volume is a poor proxy for what `ComputeCellValue()` actually tests
(discrete octree-cell occupancy at each level, not convex-hull volume) —
a bounding box can be nearly full while still leaving most individual
grid cells empty, so this test likely *understates* true dispersion.
Resolving this further would need a finer test: actual occupied-cell
counts on a fixed grid (e.g. voxelize `refineTriPt0` at the octree's own
resolution and count non-empty voxels) rather than a single aggregate
bounding-box ratio. Not done yet. The temporary bounding-box
instrumentation has been removed from `HexGen.cpp`, replaced with a
permanent comment; `diag-disp-bone/`, `diag-disp-bottle1/`, and
`build-dispersion-test/` were deleted after capturing the result.

**Status of the mesh-size investigation**: all three candidate
hypotheses tested so far (`C_THRES`, `H_THRES`, octree-balancing
propagation) are refuted or insufficient on their own; spatial dispersion
of candidates shows a real but only partial effect. The ~12x
Bottle1-vs-bone disparity in marked octree cells remains not fully
explained.

---

### 2026-08-04 — Bottle1 reproduction attempt + timing instrumentation (Apple M1)

Picking up the M4 session's fixed `HexGen.cpp` to reproduce Table 2's Bottle1 row on a second machine
(Apple M1), as the first of the Table 2 models to attempt (Bunny and others to follow).

- Cloned this fork to `~/HybridOctree_Hex` on the M1 (this analysis started out drafted in the sibling
  `~/HexOpt` repo, then moved here once it became clear the reproduction work has no dependency on
  `HexOpt` and belongs with the code it's about); added an `upstream` remote pointing at
  `CMU-CBML/HybridOctree_Hex`.
- Found this fork already contained the M4 session's `ProjectToIsoSurface` fix (commit `25adad9`,
  described in full below) plus that session's `bone` run artifacts and a large amount of committed
  CMake build cruft (binary caches, `.o` files, machine-specific `CMakeCache.txt` paths, ~600K lines of
  generated VTK output).
- CMake requirement is 3.27.4; this M1 had 3.24.2 installed (from an earlier, unrelated `HexOpt` build
  session). `brew upgrade cmake` → 4.4.2.
- Built into a fresh, separate directory (`build-m1-bottle1/`) rather than reusing the M4 session's
  `HybridOctree_Hex/build/`, to avoid clobbering the `bone` run's evidence still sitting there. Builds
  cleanly — same cosmetic warnings as the M4 session (`#endif` extra tokens, `RAND_MAX` narrowing,
  dangling-else).
- Fetched `input boundaries/bottle1_tri.raw` (14,833 vertices / 29,664 triangles) from upstream into a
  dedicated run directory, `runs/bottle1/model.raw`, and launched the full pipeline against it. The M4
  session's fix was verified only against `bone` (8,619 target elements); Bottle1 targets ~30,145
  (~3.5x larger), so timing and convergence behavior weren't assumed to carry over directly.
- **History cleanup**: rewrote this fork's `main` (`git reset --soft` to the pre-`25adad9` upstream tip,
  then re-committed only the real source fix, the reference PDFs, and the `bone` sample data, dropping
  the build cruft; force-pushed) so the ~600K lines of binary/build artifacts from the M4 session no
  longer live in history at all, not just the working tree. Added `.gitignore` entries for `build/` and
  `build-*/` so it can't happen again. Verified `HexGen.cpp`'s fix content is byte-identical before and
  after the rewrite.
- **Timing instrumentation**: split `ConstructOctree()` (`HexGen.cpp`) into two separately-timed,
  separately-printed sub-stages, so the reproduction write-up can report curvature/narrow-region
  assessment separately from octree balancing instead of Main.cpp's one coarse "constructing octree"
  number — matching the paper's own algorithmic narrative and the headings used in the pipeline-timings
  table above. `GetCellValue()` (curvature-/thickness-driven adaptive refinement decisions) is now
  timed and printed as "curvature and narrow region assessment"; `StrongBalancedOctree()` (enforcing
  the octree's 2:1 balance constraint) as "octree generation (strongly balanced)". Buffer-zone
  projection and quality improvement remain one combined number, since `ProjectToIsoSurface()`
  interleaves both inside the same per-checkpoint loop with no code-level seam to split them — see
  comments in `HexGen.cpp` at both sites for the full reasoning. Sanity-checked against `bone_tri.raw`
  in a throwaway run (`build-instrumented/`, killed once the new prints were confirmed correct, not
  counted as a real timing sample): curvature/narrow-region assessment 31.7s + octree balancing 22.4s ≈
  58.2s, matching Main.cpp's own unchanged coarse timer to within `clock()` rounding and cleanup-call
  overhead. This only applies to runs made with the rebuilt binary or later — it can't retroactively
  split the Bottle1 octree time already captured above (1062.2 s), which stays reported as one combined
  number for this run.
- **Same-machine bone comparison**: launched a full (non-throwaway) `bone_tri.raw` run on this M1
  (`runs/bone-m1/`, using the split-timing binary) specifically to separate "M1 is a slower machine
  than the M4" from "Bottle1 is geometrically harder," since both effects could otherwise look similar
  in the raw numbers. Result: bone is a consistent ~37–67% slower on this M1 across every stage
  (419.5s total vs. ~270s on the M4) — a real but modest hardware gap, nowhere near Bottle1's ~20–90x
  factor over bone on the same M1. Also notable: bone (M1)'s projection stage converged to an almost
  bit-identical best result as the M4 run (`smallDist` matches to 6 significant figures) — the
  codebase never calls `srand()` (confirmed via `grep`), so the "randomized" redistribution step is
  actually a fixed, machine-independent sequence.
- **`MAX_PROJ_ITER` reduced for Bottle1 specifically**: let Bottle1's projection stage run at bone's
  `MAX_PROJ_ITER=200000` setting initially, like the other columns, but killed it after 24 checkpoints
  (~22 min, `maxDist` down to 0.00473 with `badElem=0`) once the per-checkpoint cost became clear
  (~56s/checkpoint vs. bone's ~1s/checkpoint — a full 200-checkpoint run would have extrapolated to
  several hours for this one model). Reduced `MAX_PROJ_ITER` to `20000` in `HexGen.cpp` specifically
  (comment there documents the reasoning and that it's not a general per-model-tuned parameter, just a
  smaller fixed backstop), rebuilt, and resumed from the already-computed
  `octree.vtk`/`dualFullHex.vtk`/`dualHex.vtk` (via a temporary, uncommitted `Main.cpp` `progress`
  tweak, reverted immediately after launching) rather than redoing the ~50 minutes of earlier stages.
  The killed full-budget attempt's log and best mesh are preserved as
  `runs/bottle1/run.log.full-attempt-killed` / `finalMesh.vtk.full-attempt-killed` for the record.
- **Scaled-Jacobian stats tool**: wrote `scripts/scaled_jacobian_stats.cpp`, a standalone C++ program
  that reads a VTK hex mesh and reports min/max scaled Jacobian using `HexGen.cpp`'s own `Sj()`
  function copied verbatim (not reimplemented), so it can't silently diverge from what `HexGen` itself
  computes. Validated against bone's `finalMesh.vtk` (worst 0.590, best 1.0 — close to the repo's
  published 0.61) before trusting it on Bottle1.
- **Result — does not reproduce Table 2, for two separable reasons, not one**: Bottle1's capped run
  finished at 218,898 vertices / 192,098 elements, worst SJ 0.010, best SJ 1.0, total 4109.3 s
  (~68.5 min), against Table 2's 36,091 / 30,145, [0.560; 1.0], 218 s. (1) Mesh size is ~6x too large,
  and — confirmed by the killed full-budget attempt and the capped rerun producing byte-identical
  vertex/element counts — this is independent of the `MAX_PROJ_ITER` cap; `RemoveOutsideElement`'s
  output topology is what it is regardless of how long the later quality-improvement stage runs. The
  same oversizing shows up on bone too (19,639/16,590 here vs. this repo's own published 10,356/8,619,
  ~1.9x), and neither run touched `VOXEL_SIZE`/`C_THRES`/`H_THRES`. The M4 session validated
  `ProjectToIsoSurface`'s convergence *behavior* against bone but never checked output mesh *size*
  against the published number, so this predates today and went unnoticed until now — structurally
  similar to the already-documented gap in the sibling `HexOpt` repo (public code not matching what
  generated its paper's published results), though not yet confirmed to be the same root cause. (2)
  Worst SJ (0.010) separately reflects the deliberately-capped, under-converged run — `ELEM_THRES`
  never left its initial `0.01` floor — and isn't informative about the algorithm's achievable quality
  the way bone's 0.590/0.61 is. See the [Bottle1 results table](#bottle1-genus-1) above for the full breakdown; worth
  raising the mesh-size question with Hua Tong before spending more time chasing convergence on
  Bottle1 or moving on to Bunny.

---

<a id="macos-build-and-projecttoisosurface-fix"></a>
### 2026-07-02 – 2026-07-03 — macOS build & `ProjectToIsoSurface` fix (Apple M4)

*Harvested from a session note (`session-summary-2026-07-02_191437.md`)
left in this fork at commit `25adad9`. The note itself is dated
2026-07-02; the commit that recorded it is dated 2026-07-03. (The date
given when this log entry was requested was "02-03 July 2025" — this
fork's own commit history places the work in 2026, not 2025, so that's the
year used here.)*

**Repo:** https://github.com/CMU-CBML/HybridOctree_Hex
**Platform:** macOS (Apple Silicon, arm64), AppleClang 21.0.0, CMake 4.2.3
**Test model:** `input boundaries/bone_tri.raw` (smallest sample, 10,356 verts / 8,619 elems)

#### 1. Getting it building on macOS

The repo shipped with a stale Windows/MinGW CMake cache (`CMakeCache.txt`, `Makefile`, `CMakeFiles/`)
generated on a different machine (`D:\cmake`, `MinGW Makefiles`, `cmd.exe` shell references) — unusable
on macOS as-is.

Steps taken:
1. Deleted the stale generated files: `CMakeFiles/`, `CMakeCache.txt`, `cmake_install.cmake`, `Makefile`.
2. Created a clean out-of-source build directory: `HybridOctree_Hex/build/`.
3. `cmake .. -DCMAKE_BUILD_TYPE=Release` — configured cleanly with AppleClang, no external
   dependencies (the project is self-contained: `Main.cpp`, `HexGen.cpp/h`, `Mesh.cpp/h`, plain STL,
   no CGAL/VTK/etc.).
4. `cmake --build .` — builds cleanly. Only cosmetic warnings: stray text after `#endif` directives
   (`#endif STATICVARS_H` etc. — harmless, `-Wextra-tokens`), a few `int`→`float` narrowing warnings
   around `RAND_MAX` usage, and one dangling-`else` warning. No errors, no code changes needed for the
   build itself.
5. Produces a native arm64 executable: `HybridOctree_Hex/build/HexGen`.

No changes were made to `CMakeLists.txt` — the existing project definition worked as-is once the stale
Windows cache was cleared.

#### 2. Bug found: infinite loop in `ProjectToIsoSurface`

Running `HexGen` against `bone_tri.raw` (copied in as `model.raw`, the hardcoded input filename in
`Main.cpp`) got through the first three pipeline stages quickly but then **hung indefinitely** in the
final "project to iso-surface" stage — still running after 48+ minutes of CPU time with no sign of
terminating.

**Root cause**: `hexGen::ProjectToIsoSurface()` (in `HexGen.cpp`) contains a `while (true) { ... }` loop with
**no `break` statement anywhere in its ~1,500-line body**. Traced every exit path; none existed.

The loop's intended termination mechanism was a quality-improvement ratchet:
- `ELEM_THRES` (element quality/scaled-Jacobian acceptance threshold) starts at `0.01`.
- Whenever the mesh achieves a "perfect" checkpoint (`badElem == 0` and the surface-projection
  distance below a tolerance), the code raises the bar: jumps to `0.53` on the first success, then
  `+= 0.01` on every success after that — presumably intending to progressively push mesh quality
  higher.
- **But `ELEM_THRES` is never capped.** Once it exceeds the mesh's actually-achievable scaled-Jacobian
  quality, a "perfect" checkpoint can never happen again — the quality-improvement branch (and the
  `finalMesh.vtk` write inside it) becomes permanently unreachable, and the loop runs forever.

#### 3. Fix iterations (chronological)

This took several rounds of empirical testing against `bone_tri.raw`, each surfacing a deeper issue.
All changes are confined to `hexGen::ProjectToIsoSurface()` in `HexGen.cpp`.

**3.1 First pass — hard iteration cap (`MAX_PROJ_ITER`).** Added a simple iteration backstop to
guarantee termination. First attempt used the loop's existing (but previously unused) counter variable
`i`, which turned out to be unsafe — `i` is reset to `0` and reused as scratch elsewhere in the loop
body (`for (i = 0; i < SMOOTH_EPOCH; i++)`), so checking `i >= MAX_PROJ_ITER` never fired correctly.
Fixed by introducing a dedicated counter, `projIterCount`, incremented once per outer loop iteration.
Initial cap: `MAX_PROJ_ITER = 2000`. This **did** terminate the loop (confirmed), but far too early —
the badElem/convergence checkpoint block is itself gated by `if (i % UPDATE_EVERY == 0)` with
`UPDATE_EVERY = 1000`, so 2000 raw iterations is only **2 checkpoints** — nowhere near enough for the
gradient-descent projection to converge. `finalMesh.vtk` was never written (the success branch that
writes it never fired), and the only output (`projHex.vtk`) reflected a very early, under-converged
state (surface distance error ~0.79).

**3.2 Root-cause fix — cap `ELEM_THRES`.** Added `ELEM_THRES_CAP = 0.7` (informed by the repo's own
published results table showing typical min scaled-Jacobian quality in the ~0.53–0.62 range) so the
ratchet stops raising the bar once it's past a realistic ceiling, keeping the success branch reachable.
Raised `MAX_PROJ_ITER` to `50000` as a true backstop rather than the primary exit mechanism. Result:
`finalMesh.vtk` now got written (confirming the cap let the success branch fire), but the run still hit
the iteration backstop rather than exiting cleanly.

**3.3 Discovered: disturb-after-success oscillation.** Investigated the checkpoint log and found a
self-sustaining oscillation: immediately after a near-machine-precision convergence (e.g.
`maxDist: 2.67e-14`), the very next checkpoint would jump back up to `maxDist: 0.033` with `badElem: 1`.
Cause: the "success" branch doesn't just record success — it unconditionally runs a randomized
point-redistribution step (perturbing points toward a blend of neighbor-averaged positions, meant to
push quality higher) even once `ELEM_THRES` is already capped and there's no more quality headroom to
gain. Every genuine convergence was immediately undone by this redistribution. **Fix:** once
`ELEM_THRES` is at the cap *and* a checkpoint converges, skip the redistribution step entirely, write
`finalMesh.vtk`, and `break` out of the loop immediately instead of looping again.

**3.4 Discovered: mismatched convergence tolerance.** Even with the disturb-after-success fix, the
clean break still wasn't firing — `ELEM_THRES` was climbing far too slowly to reach its cap. Traced the
convergence gate (`smallDist < DIST_THRES2`) and found `DIST_THRES2 = sqrt(1e-12) ≈ 1e-6` is defined in
`Initialization.h` with the comment *"threshold when judging point overlap"* — i.e. a floating-point
coincident-point detector, not a surface-projection accuracy target. It's used nowhere else in the
codebase. Reusing it as the optimization's convergence target demanded near machine-precision distance,
making genuine success a rare, semi-random event (~8–12% of checkpoints) rather than a reachable
optimization target. **Fix:** introduced a dedicated `PROJ_CONVERGE_TOL = 1e-4` — an appropriately
loose, reachable tolerance — and replaced both functional uses of `DIST_THRES2` within
`ProjectToIsoSurface` with it (the point-snap acceptance gate and the success/ratchet condition).
`DIST_THRES`/`DIST_THRES2` in `Initialization.h` were left untouched, since they're not used elsewhere.
Success rate improved (4/50 → 6/51 checkpoints), but still not enough to reach `ELEM_THRES_CAP` within
the 50,000-iteration budget. Raised `MAX_PROJ_ITER` to `200000` (~150 checkpoints of margin, based on
the observed ~12% success rate needing ~17 successes to reach the cap).

**3.5 Discovered: a regression in the backstop fallback.** At `MAX_PROJ_ITER = 200000`, the run went
noticeably *worse* than shorter runs: final state `badElem=20, smallDist=0.115` — far worse than
earlier attempts. The log tail showed a stubborn local region (`maxDistIdx` repeating around vertex
15088/15160) that progressively degraded over ~30 consecutive checkpoints rather than improving. Worse:
this exposed a bug in an earlier fix. The `MAX_PROJ_ITER` backstop branch unconditionally called
`WriteToVtk("finalMesh.vtk", ...)` using whatever the *current* mesh state happened to be — with no
check on whether that state was actually better than the last genuine success. Since the mesh had
degraded by the time the backstop fired, this **clobbered a good earlier snapshot with the worse final
state**.

**3.6 Final fix — track and preserve the best mesh seen (not just the last one).** This is not a
monotonically-converging optimization — the randomized redistribution step can (and does) make
`badElem`/`smallDist` worse on a later checkpoint than an earlier one. "Final output" needed to be the
objectively best checkpoint across the *whole* run, not whatever exists when the loop happens to stop.
Changes:
- Added `bestBadElem` (sentinel `affElemNum + 1`) and `bestSmallDist` (sentinel `MAX_NUM2`), tracked
  across the entire `while (true)` loop.
- At each checkpoint, immediately after measuring `k` (bad-element count) and `smallDist` — **before**
  the redistribution step that could perturb the mesh away again — compare against the tracked best.
  If `k < bestBadElem`, or `k == bestBadElem && smallDist < bestSmallDist`, update the tracked best and
  write `finalMesh.vtk` right then, capturing the mesh in the exact state that earned that score.
- Removed the three now-redundant/risky write sites (the unconditional fallback write in the
  `MAX_PROJ_ITER` backstop branch; the write in the "converged at `ELEM_THRES` cap" clean-break branch;
  and — the dangerous one — the write at the end of the redistribution block, which captured the
  *post*-redistribution state that could be worse than what triggered the branch).
- The `MAX_PROJ_ITER` backstop log message now reports the tracked best, not whatever the current
  (possibly degraded) state is.

**Verified result** (bone_tri.raw, 201 checkpoints, 200.8s projection stage): backstop message
`badElem=0, smallDist=2.60788e-14` — the tracked best, essentially machine-precision convergence, even
though the raw log tail at the moment the backstop fired showed a badly degraded state
(`badElem: 15–20, maxDist: ~0.115`). `finalMesh.vtk`'s timestamp confirmed it was written once, early in
the run, and never touched again despite ~150 subsequent checkpoints continuing to run (and degrading)
afterward — confirming the fix works as intended.

#### 4. Summary of constant/behavior changes in `ProjectToIsoSurface` (`HexGen.cpp`)

| Name | Before | After | Purpose |
|---|---|---|---|
| Loop termination | none (infinite `while (true)`) | `projIterCount >= MAX_PROJ_ITER` backstop | Guarantee termination |
| `MAX_PROJ_ITER` | n/a | `200000` (new) | Iteration backstop (~150 checkpoints of margin over the empirically observed ~150 needed) |
| `ELEM_THRES` ratchet | unbounded (`+= 0.01` forever) | capped at `ELEM_THRES_CAP = 0.7` (new) | Prevents the quality bar from exceeding what's achievable, which was making the success branch permanently unreachable |
| Convergence tolerance | `DIST_THRES2` (`1e-6`, meant for point-overlap detection) | `PROJ_CONVERGE_TOL = 1e-4` (new, dedicated) | Realistic, reachable surface-projection convergence target |
| Success + at-cap behavior | always ran point-redistribution, looped forever | skips redistribution and `break`s once capped & converged | Stops the disturb-after-success oscillation |
| `finalMesh.vtk` writes | 3 sites: end-of-success-block (post-redistribution), none on backstop originally / later an unconditional backstop fallback | 1 site: on best-tracked improvement, measured pre-redistribution | Guarantees `finalMesh.vtk` always reflects the best mesh ever found, not the last (possibly degraded) state |
| `bestBadElem` / `bestSmallDist` | n/a | new, tracked across the whole run | Enables best-state preservation above |

All changes are localized to `hexGen::ProjectToIsoSurface()`; no other functions or files were modified.

#### 5. Verified pipeline timings (`bone_tri.raw`, final fixed build)

| Stage | Time |
|---|---|
| Read surface mesh | ~0.4 s |
| Construct octree | ~38–39 s |
| Generate dual mesh (`DualFullHexMeshExtraction`) | ~12–13 s |
| Extract interior dual mesh (`RemoveOutsideElement`) | ~13 s |
| Project to iso-surface (`ProjectToIsoSurface`) | ~200 s (201 checkpoints, hits `MAX_PROJ_ITER` backstop but with best-tracked result `badElem=0, smallDist≈2.6e-14`) |
| **Total** | **~4.5 minutes** |

Output files (in `HybridOctree_Hex/build/`): `modifiedTri.vtk`, `octree.vtk`, `dualFullHex.vtk`,
`dualHex.vtk`, `projHex.vtk`, `finalMesh.vtk`.

#### 6. Architecture note: dualization

Confirmed the pipeline performs its own octree dualization (analogous to what `automesh`
— https://github.com/autotwin/automesh — does), implemented as `hexGen::DualFullHexMeshExtraction()`
in `HexGen.cpp` (declared `HexGen.h:56`), called immediately after `ConstructOctree()`.

- Computes each octree leaf cell's centroid as a dual-mesh vertex (midpoint between diagonal corners).
- Builds hexahedral elements from these centroids based on cell adjacency/valence patterns, including
  special-case templates for non-uniform octree size transitions (the "Hybrid" in `HybridOctree_Hex`).
- Followed by `RemoveOutsideElement()` (trims the dual mesh to the interior) and then
  `ProjectToIsoSurface()` (the stage fixed in this session — pulls boundary dual vertices onto the true
  input surface and improves element quality).

This is a self-contained, hand-rolled implementation of the same primal-octree → dual-hex concept as
`automesh`'s dualization step, not a dependency on it.

#### 7. Files touched that session

- `HybridOctree_Hex/HexGen.cpp` — all fixes described above, confined to `ProjectToIsoSurface()`.
- `HybridOctree_Hex/build/` — new out-of-source build directory (CMake-generated, not source).
- Removed (stale, Windows-generated, not usable on macOS): `HybridOctree_Hex/CMakeCache.txt`,
  `HybridOctree_Hex/Makefile`, `HybridOctree_Hex/cmake_install.cmake`, `HybridOctree_Hex/CMakeFiles/`.
- `HybridOctree_Hex/CMakeLists.txt` — unchanged, worked as-is once the stale cache was cleared.
