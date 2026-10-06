#!/usr/bin/env python3
"""Compare surface preprocessing against a baseline source directory (GCC/Linux).

Both versions use -O3, identical inputs and separate processes. The benchmark
accesses private members only in its temporary driver so it can time ReadRawData
without allocating the octree or entering the unbounded projection loop.
"""

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile


DRIVER = r'''
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <string>
#include <array>
#include <limits>
#include <numeric>
#include <stdexcept>
#define private public
#include "HexGen.cpp"
#undef private

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    hexGen generator(VOXEL_SIZE);
    const auto start = std::chrono::steady_clock::now();
    if (std::string(argv[1]) == "initialize")
        generator.InitializeOctree(argv[2], argv[3]);
    else
        generator.ReadRawData(argv[2], argv[3]);
    const auto end = std::chrono::steady_clock::now();
    std::cout.precision(10);
    std::cout << std::chrono::duration<double>(end - start).count() << '\n';

    if (std::string(argv[1]) == "initialize") {
#ifdef HEXGEN_DENSE_LEVEL_TABLE
        if (generator.getLevel.size() != static_cast<size_t>(levelId[VOXEL_SIZE + 1]))
            return 3;
        // Check every table value, outside the timed section. This also checks
        // the original zero-filled deepest level; changing that is out of scope.
        for (int level = 0; level <= VOXEL_SIZE; ++level)
            for (int i = levelId[level]; i < levelId[level + 1]; ++i)
                if (generator.getLevel[i] != (level < VOXEL_SIZE ? level : 0))
                    return 4;
#else
        for (int level = 0; level <= VOXEL_SIZE; ++level)
            if (generator.GetLevel(levelId[level]) != level ||
                generator.GetLevel(levelId[level + 1] - 1) != level)
                return 4;
#endif
    }

    FILE* output = fopen((std::string(argv[3]) + ".curvature").c_str(), "wb");
    if (!output) return 5;
    const size_t count = fwrite(generator.triMesh.r, sizeof(double),
                               generator.triMesh.vNum, output);
    fclose(output);
    return count == static_cast<size_t>(generator.triMesh.vNum) ? 0 : 6;
}
'''


def compile_driver(source, driver, executable, compiler):
    defines = (["-DHEXGEN_DENSE_LEVEL_TABLE"]
               if "getLevel" in (source / "HexGen.h").read_text() else [])
    result = subprocess.run(
        [compiler, "-std=c++11", "-O3", "-DNDEBUG"] + defines + ["-I", str(source),
         str(driver), str(source / "Mesh.cpp"), "-o", str(executable)],
        capture_output=True, text=True,
    )
    if result.returncode:
        raise RuntimeError(result.stderr)


def compare(executables, raw, mode, runs, work):
    samples = {name: [] for name in executables}
    peaks = {name: [] for name in executables}
    for run in range(runs):
        # Alternate order to reduce systematic warm-cache/order bias.
        order = list(executables) if run % 2 == 0 else list(reversed(executables))
        for name in order:
            output = work / (name + ".vtk")
            memory = work / (name + ".rss")
            result = subprocess.run(
                ["/usr/bin/time", "-f", "%M", "-o", str(memory),
                 str(executables[name]), mode, str(raw), str(output)],
                capture_output=True, text=True, check=True,
            )
            samples[name].append(float(result.stdout.strip()))
            peaks[name].append(int(memory.read_text().strip()))
        for suffix in (".vtk", ".vtk.curvature"):
            if ((work / ("baseline" + suffix)).read_bytes() !=
                    (work / ("candidate" + suffix)).read_bytes()):
                raise RuntimeError(f"Output mismatch: {raw.name}{suffix}, run {run}")
    medians = {name: statistics.median(values) for name, values in samples.items()}
    return {
        "input": raw.name, "mode": mode, "seconds": samples,
        "median_seconds": medians, "peak_rss_kib": peaks,
        "speedup": medians["baseline"] / medians["candidate"],
        "byte_identical": True,
    }


def fixtures(work):
    vertices = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    cases = {
        "boundary": [(0, 1, 2), (1, 4, 2)],
        "vertex_only": [(0, 1, 2), (0, 6, 3)],
        "nonmanifold_edge": [(0, 1, 2), (1, 0, 3), (0, 1, 5)],
        "duplicate_faces": [(0, 1, 2), (0, 1, 2), (2, 1, 0)],
        "repeated_vertex": [(0, 0, 1), (0, 1, 2), (1, 1, 0)],
        "mixed_order": [(7, 6, 5), (0, 1, 2), (2, 1, 4), (3, 0, 1), (0, 6, 3)],
    }
    for name, triangles in cases.items():
        raw = work / (name + ".raw")
        rows = [f"{len(vertices)} {len(triangles)}"]
        rows.extend(" ".join(map(str, values)) for values in vertices + triangles)
        raw.write_text("\n".join(rows) + "\n")
        yield raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path, help="Baseline C++ source directory")
    parser.add_argument("inputs", nargs="+", type=Path, help="Input .raw files")
    parser.add_argument("--candidate", type=Path,
                        default=Path(__file__).resolve().parents[1] / "HybridOctree_Hex")
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--mode", choices=("read", "initialize"), default="read")
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be positive")
    for raw in args.inputs:
        if not raw.is_file():
            parser.error(f"Missing input: {raw}")
    with tempfile.TemporaryDirectory(prefix="hexgen-curvature-") as temporary:
        work = Path(temporary)
        driver = work / "driver.cpp"
        driver.write_text(DRIVER)
        executables = {name: work / name for name in ("baseline", "candidate")}
        for name, source in (("baseline", args.baseline), ("candidate", args.candidate)):
            compile_driver(source.resolve(), driver, executables[name], args.compiler)
        for raw in fixtures(work):
            compare(executables, raw, "read", 1, work)
        print(json.dumps({"regression_fixtures_passed": 6}), flush=True)
        for raw in args.inputs:
            print(json.dumps(compare(executables, raw.resolve(), args.mode,
                                     args.runs, work)), flush=True)


if __name__ == "__main__":
    main()
