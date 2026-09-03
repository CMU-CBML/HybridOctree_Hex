#!/usr/bin/env python3
"""Converts an OFF triangle mesh to OBJ.

OFF (0-based indices):

    OFF
    <num_vertices> <num_faces> <num_edges>
    x y z            (repeated num_vertices times)
    3 i j k          (repeated num_faces times; leading 3 = triangle)

OBJ (1-based indices, no header -- counts are implicit in how many v/f
lines there are):

    v x y z          (repeated num_vertices times)
    f i j k          (repeated num_faces times, each index OFF's + 1)

Validated against ground truth: this repo already carries bone_tri.obj
alongside bone_tri.raw (the file used to reverse-engineer the .raw format
in the first place). Converting bone_tri.off through this tool reproduces
bone_tri.obj's v/f lines exactly.

Every file is validated before conversion, the same policy as
tools/raw_to_off.py and for the same reason: a header ReadRawData()-style
code trusts blindly can silently read too much or too little instead of
erroring. Convert with --force to skip validation and convert anyway
using only the header's declared nv + nf lines; the result may be
silently wrong for the same reason, so treat --force as a debugging
tool, not a routine flag.

Usage
-----
    python tools/off_to_obj.py "input boundaries/bottle1_tri.off"
    python tools/off_to_obj.py "input boundaries"          # every *.off
    python tools/off_to_obj.py "input boundaries" --check   # validate only
"""
import argparse
import sys
from pathlib import Path


def off_validate(*, path):
    """Returns a list of problems with an .off file; empty if it's clean.

    Checks the magic word, the header's counts against what the file
    actually contains, then every vertex/face line's well-formedness,
    each face's declared vertex count against its actual index count,
    index range, and (for triangles) degeneracy and duplication.
    """
    problems = []
    lines = path.read_text().split("\n")
    if lines and lines[-1] == "":
        lines = lines[:-1]

    if len(lines) < 2:
        return ["file has fewer than 2 lines (missing magic word/header)"]

    if lines[0].strip() != "OFF":
        return [f"first line is {lines[0]!r}, expected exactly 'OFF' "
               f"(COFF/4OFF/nOFF variants aren't handled by this tool)"]

    header = lines[1].split()
    if len(header) != 3:
        return [f"header line has {len(header)} fields, expected 3 "
               f"(nv nf ne): {lines[1]!r}"]
    try:
        nv, nf, ne = (int(v) for v in header)
    except ValueError:
        return [f"header fields aren't integers: {lines[1]!r}"]

    expected_lines = 2 + nv + nf
    actual_lines = len(lines)
    if actual_lines != expected_lines:
        problems.append(
            f"header/content mismatch: header says {nv} verts + {nf} "
            f"faces = {expected_lines} total lines (+ magic word), file "
            f"actually has {actual_lines} lines "
            f"({'short' if actual_lines < expected_lines else 'extra'} by "
            f"{abs(expected_lines - actual_lines)})"
        )
        return problems  # the split point below can't be trusted further

    vertex_lines = lines[2:2 + nv]
    face_lines = lines[2 + nv:2 + nv + nf]

    bad_vertices = [i for i, line in enumerate(vertex_lines)
                    if len(line.split()) < 3]
    if bad_vertices:
        problems.append(f"{len(bad_vertices)} malformed vertex line(s), "
                        f"first: line {bad_vertices[0]}")

    bad_faces, wrong_count, out_of_range, degenerate = [], [], [], []
    seen, duplicates = {}, []
    for i, line in enumerate(face_lines):
        fields = line.split()
        if not fields:
            bad_faces.append(i)
            continue
        try:
            k = int(fields[0])
            idx = [int(v) for v in fields[1:1 + k]]
        except ValueError:
            bad_faces.append(i)
            continue
        if len(fields) != 1 + k:
            wrong_count.append((i, k, len(fields) - 1))
            continue
        if k != 3:
            continue  # not a triangle; degeneracy/dup checks below assume 3
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

    if bad_faces:
        problems.append(f"{len(bad_faces)} malformed face line(s), first: "
                        f"line {bad_faces[0]}")
    if wrong_count:
        problems.append(f"{len(wrong_count)} face(s) whose leading count "
                        f"doesn't match its index list, first: face "
                        f"{wrong_count[0][0]} says {wrong_count[0][1]} but "
                        f"has {wrong_count[0][2]} indices")
    if out_of_range:
        problems.append(f"{len(out_of_range)} triangle(s) with an "
                        f"out-of-range vertex index, first: face "
                        f"{out_of_range[0][0]} = {out_of_range[0][1]} "
                        f"(nv={nv})")
    if degenerate:
        problems.append(f"{len(degenerate)} degenerate triangle(s), first: "
                        f"face {degenerate[0][0]} = {degenerate[0][1]}")
    if duplicates:
        problems.append(f"{len(duplicates)} duplicate triangle(s), first: "
                        f"face {duplicates[0][0]} duplicates face "
                        f"{duplicates[0][1]}")
    return problems


def off_to_obj(*, off_path, obj_path, force=False):
    """Converts one .off file to .obj. Raises ValueError on a failed
    validation unless force=True, in which case it converts anyway using
    only the header's declared nv + nf lines."""
    problems = off_validate(path=off_path)
    if problems and not force:
        raise ValueError("; ".join(problems))

    lines = off_path.read_text().split("\n")
    num_vertices, num_faces, _ = (int(v) for v in lines[1].split())
    vertices = lines[2:2 + num_vertices]
    faces = lines[2 + num_vertices:2 + num_vertices + num_faces]
    with open(obj_path, "w") as out:
        for line in vertices:
            out.write("v " + line.strip() + "\n")
        for line in faces:
            fields = line.split()
            k = int(fields[0])
            idx = [str(int(v) + 1) for v in fields[1:1 + k]]
            out.write("f " + " ".join(idx) + "\n")
    return num_vertices, num_faces


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+",
                        help="One or more .off files, or a directory to "
                             "convert every *.off file within it")
    parser.add_argument("--check", action="store_true",
                        help="Validate only; don't write any .obj files")
    parser.add_argument("--force", action="store_true",
                        help="Convert even if validation finds a problem "
                             "(uses only the header's declared nv + nf "
                             "lines; the result may be silently wrong)")
    args = parser.parse_args()

    targets = []
    for p in args.paths:
        p = Path(p).expanduser()
        targets.extend(sorted(p.glob("*.off")) if p.is_dir() else [p])

    failures = []
    for off_path in targets:
        problems = off_validate(path=off_path)
        if args.check:
            status = "clean" if not problems else "; ".join(problems)
            print(f"{off_path.name}: {status}")
            if problems:
                failures.append(off_path.name)
            continue
        try:
            obj_path = off_path.with_suffix(".obj")
            nv, nf = off_to_obj(off_path=off_path, obj_path=obj_path,
                                force=args.force)
            print(f"{off_path.name}: {nv} vertices, {nf} faces -> "
                 f"{obj_path.name}")
        except ValueError as e:
            print(f"SKIPPED {off_path.name}: {e}", file=sys.stderr)
            failures.append(off_path.name)

    if failures:
        verb = "failed validation" if args.check else \
              "failed validation and were not converted (use --force to " \
              "override, though the result may be silently wrong)"
        print(f"\n{len(failures)} file(s) {verb}: {', '.join(failures)}",
             file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
