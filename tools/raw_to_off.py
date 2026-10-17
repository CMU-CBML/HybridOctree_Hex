#!/usr/bin/env python3
"""Converts HybridOctree_Hex "_tri.raw" triangle meshes to OFF, so tools
like MeshLab (which can't open .raw directly) can load them.

The .raw format (reverse-engineered against the companion bone_tri.obj,
2026-09-03):

    <num_vertices> <num_triangles>
    x y z          (repeated num_vertices times)
    i j k          (0-based triangle vertex indices, repeated num_triangles times)

OFF is the same data with a 2-line header (magic word, then counts
including an edge count OFF doesn't otherwise use) and each face prefixed
by its vertex count:

    OFF
    <num_vertices> <num_triangles> 0
    x y z          (repeated num_vertices times)
    3 i j k        (repeated num_triangles times)

Every file is validated before conversion. ReadRawData() (HexGen.cpp)
trusts the header's counts blindly, which is exactly what let
rocker_tri.raw's duplicated data block, and historically
bottle1_tri.raw's header/content mismatch (see analysis.md's Bottle1
section), go unnoticed -- a bad header doesn't error, it silently reads
too much or too little. Convert with --force to skip validation and
convert using only the header's declared nv + nt lines regardless of
what else the file contains; the result may be silently wrong for the
same reason the C++ reader can be, so treat --force as a debugging tool,
not a routine flag.

Usage
-----
    python tools/raw_to_off.py "input boundaries/bottle1_tri.raw"
    python tools/raw_to_off.py "input boundaries"          # every *_tri.raw
    python tools/raw_to_off.py "input boundaries" --check   # validate only
"""
import argparse
import sys
from pathlib import Path


def raw_validate(*, path):
    """Returns a list of problems with a .raw file; empty if it's clean.

    Checks the header's counts against what the file actually contains,
    then every vertex/triangle line's well-formedness, index range, and
    (for triangles) degeneracy and duplication.
    """
    problems = []
    lines = path.read_text().split("\n")
    if lines and lines[-1] == "":
        lines = lines[:-1]

    if not lines:
        return ["empty file"]

    header = lines[0].split()
    if len(header) != 2:
        return [f"header line has {len(header)} fields, expected 2: {lines[0]!r}"]
    try:
        nv, nt = int(header[0]), int(header[1])
    except ValueError:
        return [f"header fields aren't integers: {lines[0]!r}"]

    expected_lines = 1 + nv + nt
    actual_lines = len(lines)
    if actual_lines != expected_lines:
        problems.append(
            f"header/content mismatch: header says {nv} verts + {nt} tris "
            f"= {expected_lines} data lines, file actually has "
            f"{actual_lines - 1} data lines "
            f"({'short' if actual_lines < expected_lines else 'extra'} by "
            f"{abs(expected_lines - actual_lines)})"
        )
        return problems  # the split point below can't be trusted further

    vertex_lines = lines[1:1 + nv]
    triangle_lines = lines[1 + nv:1 + nv + nt]

    bad_vertices = [i for i, line in enumerate(vertex_lines)
                    if len(line.split()) != 3]
    if bad_vertices:
        problems.append(f"{len(bad_vertices)} malformed vertex line(s), "
                        f"first: line {bad_vertices[0]}")

    bad_triangles, out_of_range, degenerate = [], [], []
    seen, duplicates = {}, []
    for i, line in enumerate(triangle_lines):
        fields = line.split()
        if len(fields) != 3:
            bad_triangles.append(i)
            continue
        idx = [int(v) for v in fields]
        if any(v < 0 or v >= nv for v in idx):
            out_of_range.append((i, idx))
            continue
        if len(set(idx)) != 3:
            degenerate.append((i, idx))
            continue
        key = frozenset(idx)
        if key in seen:
            duplicates.append((i, seen[key]))
        else:
            seen[key] = i

    if bad_triangles:
        problems.append(f"{len(bad_triangles)} malformed triangle line(s), "
                        f"first: line {bad_triangles[0]}")
    if out_of_range:
        problems.append(f"{len(out_of_range)} triangle(s) with an "
                        f"out-of-range vertex index, first: triangle "
                        f"{out_of_range[0][0]} = {out_of_range[0][1]} "
                        f"(nv={nv})")
    if degenerate:
        problems.append(f"{len(degenerate)} degenerate triangle(s), first: "
                        f"triangle {degenerate[0][0]} = {degenerate[0][1]}")
    if duplicates:
        problems.append(f"{len(duplicates)} duplicate triangle(s), first: "
                        f"triangle {duplicates[0][0]} duplicates triangle "
                        f"{duplicates[0][1]}")
    return problems


def raw_to_off(*, raw_path, off_path, force=False):
    """Converts one .raw file to .off. Raises ValueError on a failed
    validation unless force=True, in which case it converts anyway using
    only the header's declared nv + nt lines."""
    problems = raw_validate(path=raw_path)
    if problems and not force:
        raise ValueError("; ".join(problems))

    lines = raw_path.read_text().split("\n")
    num_vertices, num_triangles = (int(v) for v in lines[0].split())
    vertices = lines[1:1 + num_vertices]
    faces = lines[1 + num_vertices:1 + num_vertices + num_triangles]
    with open(off_path, "w") as out:
        out.write("OFF\n")
        out.write(f"{num_vertices} {num_triangles} 0\n")
        for line in vertices:
            out.write(line.strip() + "\n")
        for line in faces:
            out.write("3 " + line.strip() + "\n")
    return num_vertices, num_triangles


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+",
                        help="One or more .raw files, or a directory to "
                             "convert every *.raw file within it")
    parser.add_argument("--check", action="store_true",
                        help="Validate only; don't write any .off files")
    parser.add_argument("--force", action="store_true",
                        help="Convert even if validation finds a problem "
                             "(uses only the header's declared nv + nt "
                             "lines; the result may be silently wrong)")
    args = parser.parse_args()

    targets = []
    for p in args.paths:
        p = Path(p).expanduser()
        targets.extend(sorted(p.glob("*.raw")) if p.is_dir() else [p])

    failures = []
    for raw_path in targets:
        problems = raw_validate(path=raw_path)
        if args.check:
            status = "clean" if not problems else "; ".join(problems)
            print(f"{raw_path.name}: {status}")
            if problems:
                failures.append(raw_path.name)
            continue
        try:
            off_path = raw_path.with_suffix(".off")
            nv, nt = raw_to_off(raw_path=raw_path, off_path=off_path,
                                force=args.force)
            print(f"{raw_path.name}: {nv} vertices, {nt} triangles -> "
                 f"{off_path.name}")
        except ValueError as e:
            print(f"SKIPPED {raw_path.name}: {e}", file=sys.stderr)
            failures.append(raw_path.name)

    if failures:
        verb = "failed validation" if args.check else \
              "failed validation and were not converted (use --force to " \
              "override, though the result may be silently wrong -- see " \
              "analysis.md's Bottle1 investigation for what a " \
              "trusted-but-wrong header can do downstream)"
        print(f"\n{len(failures)} file(s) {verb}: {', '.join(failures)}",
             file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
