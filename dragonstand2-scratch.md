# Dragon Stand2 threshold fit — scratch notes (2026-08-18, Opus session)

Not for `analysis.md`. Plain working notes.

Target (Table 2 row 5): **62,576 verts / 50,853 elements / worst SJ 0.560 /
best SJ 1.0 / refinement level 4 / 2,052 s**.
Input: `input boundaries/dragonstand2_tri.raw` (30,562 verts / 61,122 tris).

## 0. Ladder position — confirmed, and the tool's level labels ARE directly usable here

`octreeDepth = VOXEL_SIZE = VOXEL_SIZE_VALUE` (`HexGen.cpp:776`, `Main.cpp:26`).
`ComputeCellValue()` places the five tiers at levels `ladderTop-4 … ladderTop`.
With `VOXEL_SIZE_VALUE=8`, `LADDER_TOP=octreeDepth=8`, tier index `t` sits at
level `4+t` — which is exactly `refine_criteria_stats`'s own hardcoded
`level = t+4` column. **So for Dragon Stand2 there is no off-by-one**: the
tool's printed level labels are literally correct, unlike Bottle1/Bunny
(`VOXEL_SIZE_VALUE=7` → tier `t` at level `3+t`, tool labels one too deep).
Checked, not assumed.

Corollary used throughout: a refined level-`L` cell is gated by tier `L-4`.

## 1. Reference vs. shipped-threshold overshoot, per tier

`level_histogram.py --parents`:

| | ref `our results/dragonstand2.vtk` | shipped `runs/dragonstand2-v1.0-vs8/finalMesh.vtk` |
|---|---|---|
| total cells | 50,853 | 98,666 |
| L5 leaves | 735 | 195 |
| L6 leaves | 11,683 | 12,668 |
| L7 leaves | 21,903 | 39,974 |
| L8 leaves | 15,485 | 43,838 |
| L9 leaves (distortion tail) | 1,047 | 1,989 |
| refined L5 parents (tier 1) | 1,833.1 | 2,294.2 |
| refined L6 parents (tier 2) | 2,981.9 | 5,685.6 |
| refined L7 parents (tier 3) | 1,952.0 | 5,510.8 |

Total level-5 population (leaves + parents) is 2,568 (ref) vs 2,489 (shipped) —
within 3%, so tier 0 / level-4 refinement is geometry-bound and carries no
information, same as Bunny's tier 0. All the leverage is in tiers 1-3.

**Conditional refine fraction** `parents_t / (8 × parents_{t-1})` — the number
that actually needs fixing, because parent counts compound:

| tier | ref | shipped | ratio needed |
|---|---|---|---|
| 1 (L5→L6) | 0.7138 | 0.9217 | **×0.774** |
| 2 (L6→L7) | 0.2033 | 0.3098 | **×0.656** |
| 3 (L7→L8) | 0.0818 | 0.1212 | **×0.675** |

Roughly uniform ~×0.7 tightening per tier. Not a single bad tier.

## 2. `refine_criteria_stats` candidate profile

Shipped `C={0.15,0.3,0.6,1.2,2.4}`, `H={16,8,4,2,1}` (unionTri per tier):
`t0=44109  t1=23665  t2=10413  t3=2264  t4=425`

**Dragon Stand2 is thickness-dominated, unlike Bunny.** At tier 3, curvTri=537
vs thickTri=1,947. Bunny's tier 3 was essentially pure curvature (thickTri=2),
which is why Bunny's fit reduced to a single `C_THRES[3]` nudge. That lever
does *not* exist here — tightening `C_THRES[3]` alone would move ~20% of the
tier-3 candidate set. Both `C` and `H` have to move together.

Also confirmed: because shipped `C` doubles and `H` halves per tier, a
"`C×2, H÷2`" (shift1) config is *exactly* a one-tier index shift of the same
nested candidate family (verified numerically: ds2-shift1 tier0 ≡ ds2-shipped
tier1, to the last digit). `C×√2, H÷√2` ("halfshift") is the half-step.

Candidate-count ratios vs. shipped, unionTri:

| | t1 | t2 | t3 |
|---|---|---|---|
| DS2 halfshift | 0.712 | 0.498 | 0.417 |
| DS2 shift1 | 0.440 | 0.217 | 0.188 |
| Bunny halfshift | 0.608 | 0.503 | 0.425 |
| Bunny shift1 | 0.342 | 0.212 | 0.150 |

Dragon's response to halfshift is very close to Bunny's at tiers 2/3, a bit
weaker at tier 1.

## 3. The transfer function: candidate count → refined-parent count

New idea this session, and the thing that made a single full run enough.
Bunny is the one model where we have both candidate counts (cheap tool) and
refined-parent counts (four full runs, from `analysis.md`'s Bunny table), so
it can calibrate the *elasticity* `α` in
`conditional_frac_ratio ≈ candidate_ratio^α`:

Bunny halfshift-cd75 vs rl4-cd75 → conditional ratios 0.822 / 0.662 / 0.752
against candidate ratios 0.608 / 0.503 / 0.425 → α = 0.393 / 0.601 / 0.334.
Bunny shift1-cd75 vs rl4-cd75 → conditional ratios 0.687 / 0.434 / 0.459
against candidate ratios 0.342 / 0.212 / 0.150 → α = 0.350 / 0.538 / 0.411.

Two independent config pairs, consistent. Blended **α = (0.37, 0.57, 0.37)**.
Sub-linear, as expected: refinement saturates — a cell refines if *any*
candidate lands in its dilated box, so halving the candidate set does not
halve the refined-cell count.

`CELL_DETECT` 1 → 0.75 measured separately off Bunny's `rl4` vs `rl4-cd75`
parents: conditional ratios **0.977 / 0.886 / 0.976**. Small, and concentrated
at tier 2. (Worth noting: `CELL_DETECT` is a box half-width in *cellsize*
units, so at the default 1 the detection box is 2 cells wide = 8× the cell
volume. 0.5 would be the exact cell.)

## 4. Prediction for Dragon Stand2

Apply α to DS2's halfshift candidate ratios, times the `CELL_DETECT=0.75`
factors, on top of the shipped run's measured conditional fractions:

| tier | shipped frac | × halfshift | × cd75 | predicted frac | predicted parents | reference |
|---|---|---|---|---|---|---|
| 1 | 0.9217 | 0.882 | 0.977 | 0.794 | 1,975 | 1,833 |
| 2 | 0.3098 | 0.672 | 0.886 | 0.184 | 2,914 | 2,982 |
| 3 | 0.1212 | 0.724 | 0.976 | 0.086 | 1,995 | 1,952 |

Predicted leaves: L5 514, L6 12,886, L7 21,317, L8 15,960 → **~50,700 cells vs.
the target 50,853 (−0.3%)**, every tier within ±8%.

For contrast, same method predicts halfshift at `CELL_DETECT=1` → ~56,800
(+12%), and the sensitivity to α is real: α=1 would predict ~30,700 (−40%),
α=0 predicts no change (+94%). So the whole prediction rests on the Bunny-
calibrated α. Worth stating plainly as the load-bearing assumption.

## 5. Run launched

`runs/dragonstand2-v1.0-halfshift-cd75/`, binary
`HybridOctree_Hex_v1.0/build-dragonstand2/HexGen-ds2-halfshift-cd75`:

```
-DVOXEL_SIZE_VALUE=8 -DLADDER_TOP=octreeDepth
-DC_THRES_VALUES="{0.21213203,0.42426407,0.84852814,1.69705627,3.39411255}"
-DH_THRES_VALUES="{11.3137085,5.65685425,2.82842712,1.41421356,0.70710678}"
-DCELL_DETECT_VALUE=0.75
```

i.e. numerically Bunny's own `rl4-halfshift-cd75` threshold set, but at Dragon
Stand2's own ladder rung (`VOXEL_SIZE_VALUE=8`, one deeper than Bunny's 7), so
the *physical* cell sizes these thresholds gate are half Bunny's. The
coincidence of numbers is not the argument — the argument is section 4.

Launched 20:20.

### Robustness check on the level-9 distortion tail

Both reference (1,047 cells, 2.06%) and the shipped overshoot (1,989, 2.02%)
carry a level-9 tail at the same proportion. `level_histogram.py --parents`
recurses those as if they were real children of refined level-8 cells. Redid
every fraction above under the opposite assumption — that the tail is
misclassified level-8 cells, so L8 population = L8+L9 leaves — and the needed
ratios come out **0.7746 / 0.6567 / 0.6879** vs. the 0.774 / 0.656 / 0.675
above. Immaterial. Under that variant the predicted total is 51,336
(**+0.95%** vs. target) rather than 50,700 (−0.3%). Either way the config is
predicted to land inside ±1-2%.

### An early, cheap validation channel: `octree.vtk`

`octree.vtk` is written right after `ConstructOctree()` — hours before the
final mesh — and is the full balanced octree's leaves, so `level_histogram.py
--parents` reads its per-tier refined counts with **zero projection
distortion**. Shipped run's:

```
level 3: 8   4: 1696   5: 11336   6: 46344   7: 90128   8: 77184  (226,696 leaves)
refined parents: L3 504, L4 2336, L5 7352, L6 12472, L7 9648
conditional fractions: t0 0.579, t1 0.393, t2 0.212, t3 0.0967
```

Octree parents run 1.7-2.2× the final-mesh-derived parents (the interior-dual
extraction keeps only a boundary band), with survival ratios
final/octree = 0.312 / 0.456 / 0.571 at tiers 1/2/3. If those survival ratios
are roughly config-invariant, the new run's `octree.vtk` alone predicts its
final size. Predicted new-run octree conditional fractions under this config:
t0 ~0.519, t1 ~0.338, t2 ~0.126, t3 ~0.068.

### Survival ratio (final-mesh parents / octree parents) is stable across configs

Checked on Bunny's four surviving `octree.vtk`s against `analysis.md`'s
final-mesh parent table (Bunny ladder: tier `t` gates refined level `3+t`):

| config | octree t0/t1/t2/t3 | final t0/t1/t2/t3 | survival |
|---|---|---|---|
| `rl4` | 504 / 2304 / 3824 / 2992 | 129 / 963 / 2328 / 1907 | 0.256 / 0.418 / 0.609 / 0.637 |
| `rl4-cd75` | 504 / 2128 / 3248 / 2376 | 130 / 948 / 2030 / 1624 | 0.258 / 0.446 / 0.625 / 0.683 |
| `rl4-halfshift-cd75` | 496 / 1776 / 1952 / 1120 | 136 / 815 / 1155 / 695 | 0.274 / 0.459 / 0.592 / 0.620 |

Stable to ~10% across configs spanning a 2.7× range in octree size. So
`octree.vtk` is a legitimate early predictor, not just a relative signal.

**Therefore the octree-stage acceptance test for this run.** Dragon Stand2's
shipped survival ratios are 0.312 / 0.456 / 0.571 at tiers 1/2/3, so the
reference's final parents (1,833 / 2,982 / 1,952) imply target octree parents
of **≈ 5,875 / 6,540 / 3,418** (vs. the shipped run's 7,352 / 12,472 / 9,648).
The α-model predicts this config's octree at 6,330 / 6,386 / 3,493 — i.e.
+8% / −2% / +2% against those targets.

## Results

### Octree stage (t = 1,641 s, vs. 3,893 s for the shipped-threshold run)

`octree.vtk`: 126,960 leaves (shipped: 226,696).

```
level 3: 8   4: 1864   5: 11832   6: 37328   7: 51032   8: 24896
refined parents: L3 504, L4 2168, L5 5512, L6 6768, L7 3112
```

| tier | new run | derived target | vs. target | α-model prediction |
|---|---|---|---|---|
| 1 (L5) | 5,512 | 5,875 | −6.2% | 6,330 |
| 2 (L6) | 6,768 | 6,540 | +3.5% | 6,386 |
| 3 (L7) | 3,112 | 3,418 | −9.0% | 3,493 |

The α-model was right in shape and slightly loose in magnitude — the run came
in tighter than predicted at tiers 1 and 3, looser at tier 2. Every tier is
inside ±10% of target, against a shipped-threshold baseline that was
+25% / +91% / +182%.

Propagating through the shipped survival ratios gives an expected final mesh
of ≈ **48,400 cells (−4.8% vs. the 50,853 target)**, with per-tier final
parents ≈ 1,720 / 3,086 / 1,777 against the reference's 1,833 / 2,982 / 1,952.

### Dual-mesh stage (193 s, vs. 587 s shipped)

`dualHex.vtk`: 43,172 points / 30,922 cells. Shipped run's was 81,587 /
62,460, and its `finalMesh.vtk` came in at 117,793 / 98,666 — i.e. projection
inflated it ×1.444 in points and ×1.580 in cells (the split-on-low-quality
effect `analysis.md` flags). Applying the same inflation here predicts
**62,340 points / 48,857 cells**, against the 62,576 / 50,853 target: **−0.4%
points, −3.9% cells**.

Two independent estimators — octree survival ratios (≈48,400 cells) and
dual-mesh inflation (≈48,900 cells) — agree.

### Final mesh — `runs/dragonstand2-v1.0-halfshift-cd75/finalMesh.vtk`

| | reproduction | Table 2 target | error |
|---|---|---|---|
| Points | **65,497** | 62,576 | **+4.7%** |
| Cells | **53,364** | 50,853 | **+4.9%** |
| Worst SJ | 0.010003 | 0.560 | stuck on the initial gate |
| Best SJ | 1.000000 | 1.0 | ✓ |
| Refinement level | 4 | 4 | ✓ |

Stage times: read 10 s, octree 1,641 s, dual mesh 193 s, interior dual 456 s,
projection ~600 s to the `finalMesh.vtk` write. ~48 min wall to a written
final mesh, against Table 2's reported 2,052 s.

Level histogram is a close shape match:
`5: 1.03% / 6: 24.08% / 7: 44.44% / 8: 28.34% / 9: 2.11%` vs. the reference's
`1.45 / 22.97 / 43.07 / 30.45 / 2.06`.

Per-tier refined parents:

| tier | reproduction | reference | error |
|---|---|---|---|
| 1 (L5) | 2,006.4 | 1,833.1 | +9.5% |
| 2 (L6) | 3,203.1 | 2,981.9 | +7.4% |
| 3 (L7) | 1,907.7 | 1,952.0 | −2.3% |

vs. the shipped-threshold run's +25% / +91% / +182%.

**Both size predictors under-called by ~9%.** Actual 53,364 vs. predicted
48,400-48,900. Cause: projection inflated this run ×1.726 in cells, where the
shipped run inflated ×1.580 — the tighter the octree, the larger the
projection-stage split factor. So the octree/dual early-signal channel is good
for *shape* and for ruling a config in or out, but its absolute calibration
drifts with the config. Worth knowing before anyone leans on it as a
stopping rule.

**Worst SJ never left the initial 0.01 gate** — the same stall
`analysis.md` documents on Bunny and Bottle1 (five occurrences now, three
models). `smallDist` bounced around 0.56-0.59 without re-crossing `1e-4`.
Mesh size and topology are unaffected.

## Second iteration: recalibrating α on Dragon Stand2's own two runs

With two full DS2 runs in hand the elasticity no longer has to be borrowed
from Bunny. Threshold-only conditional ratios (dividing out Bunny's measured
`CELL_DETECT` factors 0.977 / 0.886 / 0.976) against the candidate ratios
0.712 / 0.498 / 0.417 give **α = 0.40 / 0.46 / 0.53** for Dragon Stand2,
versus Bunny's 0.37 / 0.57 / 0.37. Tier 1 agrees closely; tier 3 does not.

Reading the residual per-tier error in *conditional* rather than cumulative
terms localizes it cleanly:

- tier 1 fraction 0.7858 vs. reference 0.7138 → needs ×0.908 (too loose)
- tier 2 fraction 0.1996 vs. reference 0.2033 → ×1.019, **already correct**;
  its +7.4% parent excess is entirely inherited from tier 1
- tier 3 fraction 0.0745 vs. reference 0.0818 → ×1.098 (too *tight*)

So the correction is not "tighten everything" — it is tighten tier 1, leave
tier 2 alone, loosen tier 3. Solving through α and interpolating the
`refine_criteria_stats` sweep (tier 1 unionTri: 16,850 at `H[1]`=5.657 →
13,639 at 4.5 → 13,203 at 4.35 → 12,152 at 4.0; tier 3 unionTri: 944 at
`H[3]`=1.414 → 1,121 at 1.52 → 1,444 at 1.7) gives **`H_THRES[1] = 4.45`,
`H_THRES[3] = 1.56`**, everything else unchanged. Both tiers are
thickness-dominated here, so `H` is the right lever and `C` stays put.

Predicted: refined parents 1,836 / 2,931 / 1,977 against the reference's
1,833 / 2,982 / 1,952, i.e. essentially dead-on, with ±3% model slop.

`runs/dragonstand2-v1.0-halfshift-cd75-tuned/`, binary
`build-dragonstand2/HexGen-ds2-tuned`. Launched 21:18.

### It went the wrong way — and why. **The tiers are not independent.**

Tuned octree (1,804 s): 132,336 leaves, refined parents 504 / 2,160 / 5,304 /
6,952 / 3,912, versus the previous run's 504 / 2,168 / 5,512 / 6,768 / 3,112.
Bigger, not smaller. Octree conditional fractions:

| tier | halfshift-cd75 | tuned | ratio | intended |
|---|---|---|---|---|
| 1 | 0.3178 | 0.3069 | 0.966 | 0.908 |
| 2 | 0.15347 | 0.16383 | **1.068** | **1.000 — thresholds untouched** |
| 3 | 0.05748 | 0.07033 | 1.224 | 1.098 |

**Tier 2 moved by +6.8% with its own `C_THRES[2]`/`H_THRES[2]` unchanged.**
The cause is in `ComputeCellValue()` (`HexGen.cpp:1017`): for any cell with
`level < octreeDepth - 1`, the function first loops over the eight children
and returns `intersect = true` if *any* child intersects — it only falls
through to its own tier's box test when no child does. So `intersect`
propagates *upward*: loosening a deep tier drags every shallower tier along
with it. My model had treated the five tiers as independent knobs. They are
not — deeper tiers are strictly coupled into shallower ones.

That inverts the sign of the tier-3 correction: tier 3 read as 9% "too tight"
in isolation, but loosening it added refinement at tiers 1 and 2 that more
than cancelled the tier-1 tightening. Projected final size for this config,
via the previous run's measured survival ratios, is ≈56,900 cells (+12%) —
worse than the run it was trying to improve. Killed at the octree stage rather
than spending another 40 min confirming a known regression; the octree
measurement is the informative part and it is recorded above.

Silver lining: this is a clean, quantified measurement of the coupling
(one level of propagation ≈ +6.8% on the tier above, for a tier-3 loosening of
+22%), which is new and is not in `analysis.md`.

## Third iteration: tier 1 only

Corrected reading of run 1's residual: leave tiers 2 and 3 alone entirely
(their conditional fractions were 1.9% and 9% off, and any deliberate change
to tier 3 pollutes tiers 1-2), and tighten **only** `H_THRES[1]`, which is
safe because tightening a shallow tier cannot propagate downward.

`H_THRES[1] = 4.45`, everything else exactly the run-1 halfshift-cd75 values.
Bracketing the unknown propagation factor `p` (1.00-1.03) that run 2's tier-1
number was contaminated by, the tier-1-only conditional ratio should be
0.938-0.966, giving a final mesh of **50,200-51,600 cells, i.e. −1.2% to
+1.6%** of the 50,853 target.

`runs/dragonstand2-v1.0-halfshift-cd75-h1/`, binary
`build-dragonstand2/HexGen-ds2-h1only`. Launched 21:50.

### Octree result — and the structural finding that follows from it

`octree.vtk` (1,595 s): 125,168 leaves, refined parents
504 / 2,160 / **5,264** / **6,768** / **3,112**, against the halfshift-cd75
run's 504 / 2,168 / **5,512** / **6,768** / **3,112**.

**Tiers 2 and 3 are bit-identical.** Tightening `H_THRES[1]` moved tier 1 by
−4.5% and left the two tiers below it untouched — to the cell.

That is not a coincidence, and reading `ComputeCellValue()` again explains it
exactly. A level-`L` cell (for `L < octreeDepth-1`) intersects if **any child
intersects**, and only falls through to its own tier test when none does. So a
level-5 cell that fails its own tier-1 test necessarily had no refining
children either — otherwise the child check would have marked it. Therefore
tightening tier 1 can only ever delete level-6 cells that were *already
leaves*. Each deletion is −7 cells (8 children replaced by 1 leaf) and nothing
deeper is affected.

**The rule, both directions:**

- **Tightening tier `t` affects tier `t` only.** Confirmed exactly here.
- **Loosening tier `t` propagates upward** to every shallower tier. Measured
  in run 2: a +22% tier-3 loosening pulled tier 2 up 6.8% and tier 1 up ~0.8%.

So the five tiers *are* separately fittable as long as you only ever tighten,
and each tier's refined-parent count should be matched against the reference
independently. That is a much better fitting procedure than the
`analysis.md` sessions had — Bunny's `C_THRES[3]` nudge worked precisely
because it was a single-tier loosening at the *deepest* tier, where there is
nothing below to propagate from.

It also corrects my run-1 error analysis: tier 2's +7.4% parent excess is
**not** inherited from tier 1 (as I wrote above — wrong), it is tier 2's own,
and tightening tier 1 will not fix it. Revised expectation for this run:
P1 2,006 → ~1,916, P2 and P3 unchanged → ≈ **52,700 cells, +3.65%**.
An improvement on +4.9%, but a small one.

## Fourth iteration: tier 1 and tier 2, both tightened

With the independence rule, the fit becomes a straightforward per-tier match
against the reference's refined-parent counts (1,833 / 2,982 / 1,952):

| tier | run-1 | after `H[1]`=4.45 | needed |
|---|---|---|---|
| 1 | 2,006 | 1,916 | 1,833 |
| 2 | 3,203 | 3,203 | 2,982 (ratio 0.931) |
| 3 | 1,908 | 1,908 | 1,952 — leave alone, and *must* leave alone, since loosening it would pollute tiers 1-2 |

Octree-level elasticity, now measurable directly from the uncontaminated
tier-1 pair (candidate ratio 13,482/16,850 = 0.800 → conditional ratio 0.9585):
**α_octree ≈ 0.19**, much flatter than the final-mesh-level 0.40-0.53. Scaling
that to tier 2 (α₂ ≈ 0.19-0.22) and solving for a 0.931 conditional ratio puts
the tier-2 candidate target at 3,550-3,750 unionTri, from 5,186.

Sweep: `H_THRES[2]` 2.828 → 5,186; **2.4 → 3,708**; 2.1 → 2,969. Take **2.4**.

`H_THRES = {11.3137085, 4.45, 2.4, 1.41421356, 0.70710678}`, `C_THRES` at the
halfshift values, `CELL_DETECT` 0.75. Predicted P1/P2/P3 ≈ 1,916 / 2,976 /
1,908 → **≈ 51,100 cells (+0.5%) and ≈ 62,700 points (+0.2%)**.

`runs/dragonstand2-v1.0-halfshift-cd75-h1h2/`, binary
`build-dragonstand2/HexGen-ds2-h1h2`. Launched 22:26. Run 3 left running in
parallel as an end-to-end check of the tier-1-only prediction.

### Run 3 final mesh — independence confirmed end-to-end

**64,519 points / 52,479 cells** (+3.1% / +3.2%), worst SJ 0.010001 (stalled
again), best SJ 1.0. Predicted +3.65%; actual +3.2%.

| tier | run 1 | run 3 (`H[1]`=4.45) | reference |
|---|---|---|---|
| 1 | 2,006.4 | **1,866.5** | 1,833.1 |
| 2 | 3,203.1 | **3,210.3** | 2,981.9 |
| 3 | 1,907.7 | **1,907.6** | 1,952.0 |

Tiers 2 and 3 moved by +0.2% and −0.01% — i.e. not at all — while tier 1 moved
−7.0%. The independence rule holds in the *final* mesh, not just the octree,
which is what makes it usable as a fitting procedure.

Fitting the two runs to `total = a + b·ΣP` gives `a` = 5,935, `b` = 6.664, and
that line predicts the reference's own total from its own ΣP to within 0.35%
(51,030 predicted vs. 50,853 actual). Run 4's predicted ΣP is 6,756 against
the reference's 6,767 → **≈ 50,960 cells (+0.2%) and ≈ 62,600 points
(+0.1%)**.

### Run 4 result

Octree (1,467 s): parents 504 / 2,144 / 5,136 / **6,008** / **3,112**. Tier 3
identical again (3,112, third run running). Tier 2 conditional went
0.16072 → 0.14621, ratio 0.9097 at candidate ratio 0.715 — so the true
tier-2 octree elasticity is **α₂ = 0.282**, not the 0.19-0.22 I assumed, and
`H[2]`=2.4 tightened ~7% harder than intended. Tier 1 also fell 2.4%, the
expected small downward drag from tightening a deeper tier.

Final mesh: **61,138 points / 49,620 cells** → **−2.30% / −2.42%**. Worst SJ
0.010000 (stalled), best SJ 1.0.

| tier | run 4 | reference | error |
|---|---|---|---|
| 1 | 1,854.6 | 1,833.1 | +1.2% |
| 2 | 2,822.4 | 2,981.9 | −5.3% |
| 3 | 1,895.5 | 1,952.0 | −2.9% |

## Fifth iteration — bracketing on `H_THRES[2]`

Run 3 (+3.2%) and run 4 (−2.4%) bracket the target, with `H[2]` = 2.828 and
2.4 the only difference. Refitting `total = a + b·ΣP` on runs 3 and 4 gives
`a` = 4,000, `b` = 6.941, which reproduces the reference's own total from its
own ΣP to +0.2% (50,969 predicted vs. 50,853). With the now-measured
α₂ = 0.282, `H[2]` = 2.5 predicts ΣP ≈ 6,654 → **≈ 50,200 cells, −1.3%**.

`runs/dragonstand2-v1.0-halfshift-cd75-h1h25/`, binary
`build-dragonstand2/HexGen-ds2-h1h25`. Launched 22:49.

### Run 5 crashed — a real, reproducible `v1.0` segfault

Octree came through fine (1,488 s; parents 504 / 2,144 / 5,136 / **6,152** /
3,112 — tier 1 and tier 3 bit-identical to run 4 again, tier 2 up 2.4% as
intended), dual mesh and interior dual both fine, then the process died
silently the moment `ProjectToIsoSurface` started.

`~/Library/Logs/DiagnosticReports/HexGen-2026-08-18-232512.ips`:

```
EXC_BAD_ACCESS (SIGSEGV), KERN_INVALID_ADDRESS at 0x0000000464854000
  __bzero + 88
  hexGen::ProjectToIsoSurface(char const*) + 6860
  main + 760
```

A wild pointer, not an allocation failure — so not a resource problem from
running three jobs at once, and deterministic for this mesh (nothing upstream
of it uses an unseeded `rand()`). Re-running the identical config would
reproduce it. Worth flagging on its own: `v1.0`'s projection stage can hard-
crash on a perfectly ordinary input, which is a different failure mode from
the convergence stalls `analysis.md` documents.

Worked around by moving to a neighbouring `H[2]` rather than debugging it.

## Sixth iteration — the final fit

Using the three measured `H[2]` points (2.4 → 6,008 octree tier-2 parents,
2.5 → 6,152, 2.828 → 6,768) and the measured octree→final survival ratio for
tier 2 (0.470-0.474), the reference's P2 = 2,982 implies an octree tier-2
target of ≈ 6,318, i.e. **`H_THRES[2]` ≈ 2.6**. At that setting P1 ≈ 1,866,
P2 ≈ 2,990, P3 ≈ 1,900 → ΣP ≈ 6,756 against the reference's 6,767 →
**≈ 50,900 cells (+0.1%)**.

`runs/dragonstand2-v1.0-halfshift-cd75-h1h26/`, binary
`build-dragonstand2/HexGen-ds2-h1h26`. Launched 23:50.

### Run 6 — best Dragon Stand2 reproduction

Octree tier-2 parents came in at **6,328** against the 6,318 target. Final
mesh:

| | reproduction | Table 2 target | error |
|---|---|---|---|
| Points | **62,307** | 62,576 | **−0.43%** |
| Cells | **50,603** | 50,853 | **−0.49%** |
| Worst SJ | 0.010004 | 0.560 | stalled (see below) |
| Best SJ | 1.000000 | 1.0 | ✓ |
| Refinement level | 4 | 4 | ✓ |

Octree-level histogram, reproduction vs. reference:

| level | reproduction | reference |
|---|---|---|
| 5 | 732 (1.45%) | 735 (1.45%) |
| 6 | 11,879 (23.47%) | 11,683 (22.97%) |
| 7 | 21,840 (43.16%) | 21,903 (43.07%) |
| 8 | 15,052 (29.75%) | 15,485 (30.45%) |
| 9 (distortion tail) | 1,099 (2.17%) | 1,047 (2.06%) |

Refined parents: **1,855.8 / 2,967.3 / 1,898.7** against the reference's
1,833.1 / 2,981.9 / 1,952.0 — **+1.2% / −0.5% / −2.7%**, every tier inside 3%.

Stage times: octree 1,486 s, dual mesh 175 s, interior dual 439 s, projection
to the `finalMesh.vtk` write ~600 s (~45 min wall, vs. Table 2's 2,052 s).

## Final configuration

```
/usr/bin/c++ -O3 -DNDEBUG -std=c++17 -arch arm64 \
  -DVOXEL_SIZE_VALUE=8 -DLADDER_TOP=octreeDepth \
  -DC_THRES_VALUES="{0.21213203,0.42426407,0.84852814,1.69705627,3.39411255}" \
  -DH_THRES_VALUES="{11.3137085,4.45,2.6,1.41421356,0.70710678}" \
  -DCELL_DETECT_VALUE=0.75 \
  Main.cpp HexGen.cpp Mesh.cpp -o HexGen
```

`C_THRES` = shipped × √2. `H_THRES` = shipped ÷ √2 at tiers 0, 3, 4, with
tiers 1 and 2 tightened further to 4.45 (from 5.657) and 2.6 (from 2.828).
`CELL_DETECT` = 0.75.

## Every run, in order

| run | `H_THRES` | points | cells | vs. target |
|---|---|---|---|---|
| target | | 62,576 | 50,853 | — |
| `vs8` (shipped, prior session) | {16,8,4,2,1}, `C` shipped, CD 1 | 117,793 | 98,666 | +88% / +94% |
| 1 `halfshift-cd75` | {11.314,5.657,2.828,1.414,0.707} | 65,497 | 53,364 | +4.7% / +4.9% |
| 2 `halfshift-cd75-tuned` | {11.314,4.45,2.828,**1.56**,0.707} | — | ≈56,900 est. | killed at octree, ≈+12% |
| 3 `halfshift-cd75-h1` | {11.314,**4.45**,2.828,1.414,0.707} | 64,519 | 52,479 | +3.1% / +3.2% |
| 4 `halfshift-cd75-h1h2` | {11.314,4.45,**2.4**,1.414,0.707} | 61,138 | 49,620 | −2.3% / −2.4% |
| 5 `halfshift-cd75-h1h25` | {11.314,4.45,**2.5**,1.414,0.707} | — | — | segfaulted in projection |
| **6 `halfshift-cd75-h1h26`** | {11.314,4.45,**2.6**,1.414,0.707} | **62,307** | **50,603** | **−0.43% / −0.49%** |

All six used `C_THRES` = {0.21213203, 0.42426407, 0.84852814, 1.69705627,
3.39411255} and `CELL_DETECT` = 0.75 at `VOXEL_SIZE_VALUE`=8,
`LADDER_TOP`=`octreeDepth`.

## What generalizes (candidates for `analysis.md`)

1. **Tightening a tier changes that tier and nothing below it — exactly.**
   Verified at the octree (bit-identical deeper-tier parent counts across
   runs 3/4/5/6) and end-to-end in the final mesh (run 3 moved tier 1 by −7.0%
   and tiers 2/3 by +0.2%/−0.01%). Follows from `ComputeCellValue()`'s
   child-first `intersect` check: a cell that fails its own tier test provably
   had no refining children.
2. **Loosening a tier propagates upward** to every shallower tier
   (+22% at tier 3 → +6.8% at tier 2, run 2). So a fit should only ever
   tighten, working shallow-to-deep, and should start from a configuration
   that overshoots at every tier.
3. **Therefore the fit is per-tier and mechanical**, not a global search:
   read the reference's refined-parent counts, then tighten each tier's own
   `C`/`H` until its count matches. Four full runs took Dragon Stand2 from
   +94% to −0.5%.
4. **`refine_criteria_stats` candidate counts map to refined-parent counts by
   a power law**, `conditional_ratio ≈ candidate_ratio^α`. Measured α differs
   between the octree and the final mesh (≈0.19-0.28 vs. ≈0.40-0.53) and
   between models (Bunny 0.37/0.57/0.37, Dragon Stand2 0.40/0.46/0.53), so it
   is a *local* calibration, good for one step of a search, not a transferable
   constant. Two runs of the same model are enough to fit it.
5. **`octree.vtk` is a good early signal for shape, a poor one for absolute
   size.** Its per-tier parent counts are exact and projection-free, but the
   octree→final survival ratio drifts with configuration (Dragon Stand2's
   projection inflation was ×1.580 for the shipped run and ×1.726 for run 1),
   so absolute-size predictions off the octree ran ~10% low. Once two runs of
   the same model exist, `total = a + b·ΣP` fitted on them predicted both the
   next run and the *reference's own* total to within 0.5%.
6. **The worst-SJ stall is universal here.** All five Dragon Stand2 runs that
   reached projection ended pinned at worst SJ ≈ 0.0100, never ratcheting.
   Eight occurrences now across four models. It has no effect on mesh size or
   topology.
7. **`v1.0`'s `ProjectToIsoSurface` can segfault** (wild pointer in a `bzero`,
   run 5). New failure mode, distinct from the stalls.

Whether the *specific* thresholds transfer to another model is untested — on
the `analysis.md` evidence they will not. But the *procedure* in items 1-3 is
model-independent and should make the remaining Table 2 models much cheaper
than Dragon Stand2's six runs, since most of those six went into discovering
the coupling rule rather than into the fit itself.
