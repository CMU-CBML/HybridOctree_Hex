[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jocs.2024.102278-blue)](https://doi.org/10.1016/j.jocs.2024.102278)
# HybridOctree_Hex
Please download the latest release of the software.

HybridOctree_Hex received the best technical poster award in [the 2024 International Meshing Roundtable](https://internationalmeshingroundtable.com/awards/).

Please also feel free to check out other works on grid-based hexahedral mesh generation:

1. https://github.com/CMU-CBML/MCHex
2. https://github.com/CMU-CBML/Element-Saving-Hexahedral-3-Refinement-Templates

## Build and reproduce legacy behavior

The default configuration preserves the meshing behavior of commit
`dfcb3db9b219bf8dbdf4e9130547d7c51f0d63d4`. Adaptive refinement, balancing decisions,
geometry predicates, random-call order and the projection quality ramp retain
that implementation. The earlier experimental BVH, stronger balancing and fixed
SJ 0.53 stopping policy have been removed from this PR.

```sh
cmake -S HybridOctree_Hex -B build/legacy -DCMAKE_BUILD_TYPE=Release -DHEXGEN_ABS_MODE=NATIVE
cmake --build build/legacy --config Release --parallel
# Run in a directory containing model.raw:
/path/to/build/legacy/HexGen --seed 1 --iterations 1000
```

`--iterations 1000` performs exactly 1,000 original outer-loop optimization
steps, writes the state at that step to `iterationMesh.vtk`, and returns. The
counter used for this limit is separate from the legacy loop variable; the
original smoothing/checkpoint schedule and quality ramp are unchanged. Reaching
the step limit does not imply convergence. `projHex.vtk` and `finalMesh.vtk` keep
their original checkpoint/milestone meanings. With no limit (or `--iterations 0`),
projection continues indefinitely, even after reaching a quality milestone, as
in the original implementation. Use a fresh run directory for each comparison.
The default input remains `model.raw`, depth 10, with the original curvature and
thickness thresholds. No parameter tuning is needed to enable compatibility.

Exact reproduction requires the same input bytes, constants, seed, compiler,
standard library, floating-point flags and optimization step as the legacy run.
Do not compare runs stopped after the same wall-clock duration: a faster program
will have performed more optimization steps. Historical `our results/*.vtk`
files lack their full build/run provenance, so they are not an exact regression
reference.

### Match the legacy numeric overload

The original source calls unqualified `abs(double)`. Its overload depends on the
build environment. On the tested GCC 11.4 build, `abs(0.5)` returns integer zero;
other legacy builds can select the floating overload. This PR does not silently
replace that behavior.

| CMake setting | Behavior |
| --- | --- |
| `HEXGEN_ABS_MODE=NATIVE` (default) | Keep the original compiler's unqualified `abs` resolution. Use for the same compiler/environment as the old build. |
| `HEXGEN_ABS_MODE=FLOATING` | Explicitly select the floating overload. Use only to match a legacy build using floating `abs(double)`. |

FLOATING selects an overload; it does not enable the removed balancing/BVH changes.
Matching this setting alone does not guarantee cross-compiler bitwise equality.

CMake 3.16+ and C++11 are required. Use an out-of-source build: the checked-in
cache/Makefiles refer to a Windows installation. Single-configuration generators
default to Release and respect an explicitly chosen build type.

## Result-preserving optimizations

- Index incident triangles for curvature, preserving triangle accumulation order.
- Store the dense level lookup as bytes, preserving every original value,
  including the zero-filled deepest level. This saves about 3.43 GiB at depth 10.
- Store only occupied balancing corner buckets. The original minimum-corner
  test, exactly-eight condition, sibling refinement, ordering, padded leaf arrays
  and recursive passes are unchanged. This removes the approximately 24 GiB
  dense corner-grid allocation per pass without strengthening balancing.
- Index integer octree vertices by exact lattice coordinates, preserving
  first-encounter numbering. On this integer lattice, the original squared
  distance `< 1e-12` test is equivalent to exact equality.
- Index incident cell/corner pairs during dual extraction, preserving their
  original cell-then-corner order.
- Skip the inactive level-9 refinement scan. Its original criterion is commented
  out; active levels retain their original descending traversal and predicates.

Compatibility deliberately retains legacy numerical and algorithmic behavior,
including its limitations. Approximately 2 GiB of initial octree/list storage and
expensive active cell/triangle scans remain. Bug fixes that alter geometry or
adaptive topology should be separate changes with separate acceptance criteria.

## Differential validation

The checked-in runner reads the actual baseline commit with `git show`, builds
both versions with the same flags, and compares exact output bytes. It does not
use published example meshes as a substitute for executing the reference.
Run from a clone containing the baseline commit, with Python 3, GCC/Clang and GNU
`time` on Linux:

```sh
ctest --test-dir build/legacy --output-on-failure
python3 tools/verify_legacy.py --suite fixtures --abs-mode FLOATING --output /tmp/hexgen-float-tests
python3 tools/verify_legacy.py --suite pipeline --abs-mode FLOATING --iterations 1000 --reference-release-grid \
  --output /tmp/hexgen-legacy-pipeline \
  --inputs "input boundaries/bone_tri.raw" "input boundaries/dtorus_tri.raw"
```

Fixtures compare three adaptive partitions, ordered leaf IDs/levels, octree and
dual-mesh coordinates/connectivity, and the 1,000th and 1,001st optimization
steps, including checkpoint state and the next random value. This checks both a
checkpoint boundary and a step between checkpoints. The pipeline suite compares
curvature, pre/post-balancing leaf lists, every mesh stage and projection
checkpoints, both as VTK and unrounded binary coordinate/connectivity data.
The reference records observations and returns after the requested number of
actual optimization steps.
FLOATING supplies the same floating overload to both reference and candidate.
No thresholds or refinement rules are patched.

The full reference requires substantial RAM (the runner caps address space at
36 GiB) and can be slow. `--reference-release-grid` frees its dead dense corner
table immediately before tail recursion, avoiding another approximately 24 GiB
of live allocation per recursive call. This changes allocation lifetime only;
all original balancing decisions, arrays and ordering remain intact. Fixture
tests also run without this reference adjustment. The report records its use.
A fixed-step comparison verifies that observation
boundary, not eventual convergence or byte identity with undocumented archived
meshes. Full run logs, compiler/source hashes, stage timings and output hashes
are written to the selected output directory.

### Recorded compatibility check (2026-10-05)

GCC 11.4, `-O3 -DNDEBUG -std=c++11`, FLOATING, default depth/thresholds:

| Case | Completed old pipeline | Completed new pipeline | Speedup | Stage outputs |
| --- | ---: | ---: | ---: | --- |
| bone | 136.11 s | 90.38 s | 1.51x | Byte-identical |
| dtorus | 209.48 s | 145.44 s | 1.44x | Byte-identical |

These already completed measurements include three original checkpoints (3,000
optimization steps), one run per version, with the reference grid-lifetime
adjustment described above. Subsequent short tests verify the actual-step limit
at 1,000 steps on both models, replaying the identical interior mesh, and at
1,000/1,001 steps on fixtures in both NATIVE and FLOATING modes. Coordinates and
connectivity match in unrounded binary form. Projection equations are unchanged;
the measured speedups come from preprocessing, indexing and octree bookkeeping.
See [machine-readable results](tools/legacy_validation.json) for hashes, stage
timings, matching leaf-level counts and the precise measurement scope.

To check a different projection step without rebuilding the octree, reuse the
binary interior mesh saved by a previous pipeline run:

```sh
python3 tools/verify_legacy.py --suite projection --abs-mode FLOATING --iterations 1000 \
  --state-dir /tmp/hexgen-legacy-pipeline/runs --output /tmp/hexgen-projection \
  --inputs "input boundaries/bone_tri.raw" "input boundaries/dtorus_tri.raw"
```

This is a stage replay: both versions load the same unrounded interior mesh and
reset the seed to 1 at the projection boundary. It verifies the fixed-step
projection behavior independently; its random state is not presented as a resume
of the full pipeline's random stream.

Preprocessing can be measured separately, with three alternating old/new runs
and six synthetic cases in addition to the supplied models:

```sh
baseline_dir="$(mktemp -d)"
git archive dfcb3db HybridOctree_Hex | tar -x -C "$baseline_dir"
python3 tools/benchmark_curvature.py "$baseline_dir/HybridOctree_Hex" \
  "input boundaries/bone_tri.raw" "input boundaries/bunny_tri.raw" \
  "input boundaries/fertility_tri.raw"
```

# Citation
```angular2html
@article{tong2024hybridoctree_hex,
  title={HybridOctree\_Hex: Hybrid octree-based adaptive all-hexahedral mesh generation with Jacobian control},
  author={Tong, Hua and Halilaj, Eni and Zhang, Yongjie Jessica},
  journal={Journal of Computational Science},
  volume={78},
  pages={102278},
  year={2024},
  publisher={Elsevier}
}
```

# Mesh Statistics:
|Model|#Vert|#Elem|Min SJ|
| :--: | :--: | :--: | :--: |
|airplane1|29899|24318|0.57|
|airplane2|29241|23646|0.61|
|ant2|57833|48120|0.60|
|armadillo|326447|282201|0.54|
|armchair|75577|65280|0.60|
|bimba|370664|330619|0.57|
|bird|25785|21067|0.58|
|blade|103933|88351|0.56|
|bone|10356|8619|0.61|
|botijo|161860|139327|0.57|
|bottle1|36091|30145|0.56|
|bottle2|89357|78366|0.58|
|bumpy_sphere|25466|21646|0.58|
|bumpy_torus|249852|219832|0.56|
|bunny|26375|21695|0.57|
|buste|71068|61441|0.56|
|camel|124573|105417|0.57|
|camille_hand|114453|100023|0.59|
|carter|190871|160122|0.55|
|chair|232460|192814|0.57|
|chair1|50650|41034|0.58|
|chinese_lion|351402|309525|0.55|
|cup|386677|342009|0.58|
|Cup1|89032|73266|0.58|
|dancer|58353|47994|0.56|
|dancer2|43774|35142|0.57|
|dancing_children|457139|396964|0.54|
|david|319465|282957|0.56|
|deformed_armadillo|43216|34939|0.56|
|dente|129398|114812|0.58|
|dilo|94803|79929|0.60|
|dino|113019|94303|0.55|
|dino2|75158|62562|0.56|
|dragonstand2|62576|50853|0.56|
|dragon_stand|240896|206695|0.55|
|dtorus|19515|15872|0.58|
|duck|65685|58538|0.58|
|duck2|218549|193450|0.58|
|eight|30157|25264|0.61|
|elephant|59595|49242|0.57|
|elk|355959|313725|0.56|
|eros|75907|65798|0.57|
|fertility|45514|38475|0.58|
|fish1|64760|53611|0.58|
|fish2|51367|42853|0.60|
|foot|64617|56468|0.58|
|gargoyle|273704|236689|0.55|
|genus3|54592|45176|0.57|
|glass1|42870|34034|0.57|
|glass2|14791|11602|0.62|
|grayloc|1102475|986401|0.56|
|greek_sculpture|187058|162674|0.57|
|hand|88887|75770|0.57|
|head|62782|55038|0.55|
|head1|153260|134611|0.56|
|head2|361054|327103|0.57|
|holes3|41611|36132|0.61|
|homer|166322|146056|0.57|
|horse|102832|87595|0.57|
|human1|148026|128154|0.60|
|human2|141982|119792|0.54|
|human3|187747|163852|0.55|
|igea|74420|64303|0.56|
|insect|108197|89154|0.54|
|isidore_horse|209974|182124|0.54|
|Kiss|116510|100375|0.55|
|kitten|38138|32258|0.55|
|lion_recon|112847|95251|0.57|
|master_cylinder|286569|247194|0.53|
|max|55222|47687|0.56|
|moai|145934|127950|0.54|
|mouse|99657|85890|0.53|
|oil_pump|233702|196455|0.54|
|oni|135607|117013|0.55|
|pear|19308|16122|0.58|
|pensatore|380202|340985|0.56|
|pierrot|170778|151683|0.55|
|pig|84664|72571|0.56|
|ramses|44790|37993|0.59|
|red_circular_box|351881|313866|0.56|
|retinal|19021|15950|0.58|
|rocker|93607|80266|0.54|
|rolling_stage|343540|297087|0.53|
|santa|93459|78597|0.54|
|screwdriver|147661|126098|0.53|
|sediapatch|56036|47445|0.55|
|sphereInShell|128626|108381|0.62|
|sphinx|206944|183216|0.55|
|teaport|20093|16726|0.57|
|thai_statue|64764|53831|0.55|
|toy1|20274|16660|0.58|
|toy2|15818|13268|0.58|
|uu-memento|34335|27228|0.58|
|venus|25794|21895|0.58|
|woodenfish|254517|214410|0.55|
|wrench|56712|45382|0.54|
