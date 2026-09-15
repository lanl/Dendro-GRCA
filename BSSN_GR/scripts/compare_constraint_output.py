#!/usr/bin/env python3
"""Compare BSSN constraint and psi4 output between two run directories.

Constraint differences are scaled to the peak of their own column. psi4 mode
differences are scaled to the peak over every mode and radius: the subdominant
modes of an equal-mass binary are numerically zero, so scaling those to
themselves measures the noise floor rather than the kernel.

Compare only builds made at the same CPU_ARCH. The derivative stencils are
plain loops the compiler auto-vectorizes, so a different -march changes the
summation order feeding the kernel and swamps what this is trying to measure.

Reference point for the tolerances: gcc 16.1 at -march=tigerlake, cascade vs
scalar, is 2.9e-15 on the constraints and 1.8e-15 on psi4.
"""
import argparse
import glob
import math
import os
import sys


def read_table(path):
    with open(path) as f:
        rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    header = [c.strip() for c in rows[0] if c.strip()]
    data = [[c.strip() for c in r if c.strip()] for r in rows[1:]]
    return header, data


def parse_cell(text):
    if text.startswith("("):
        re_s, im_s = text[1:-1].split(",")
        return complex(float(re_s), float(im_s))
    return complex(float(text), 0.0)


def column_values(data, col):
    return [parse_cell(r[col]) for r in data]


def compare_file(ref_path, new_path, tol, label, scale=None):
    h_ref, d_ref = read_table(ref_path)
    h_new, d_new = read_table(new_path)
    name = os.path.basename(ref_path)
    if h_ref != h_new:
        return [(name, "header", float("inf"), "column layout differs")]
    if len(d_ref) != len(d_new):
        return [(name, "rows", float("inf"),
                 f"{len(d_ref)} vs {len(d_new)} rows")]
    if not d_ref:
        return [(name, "rows", float("inf"), "no data rows")]

    worst = []
    ncol = len(h_ref)
    for col in range(ncol):
        if h_ref[col] in ("TimeStep", "t", "time"):
            # Step index and simulation time must match exactly; a drift here
            # means the runs diverged, not that the kernel is imprecise.
            for r, (ra, rb) in enumerate(zip(d_ref, d_new)):
                if ra[col] != rb[col]:
                    worst.append((name, h_ref[col], float("inf"),
                                  f"row {r}: {ra[col]} vs {rb[col]}"))
            continue
        vals_ref = column_values(d_ref, col)
        vals_new = column_values(d_new, col)
        peak = scale if scale else max(abs(v) for v in vals_ref)
        if peak == 0.0:
            peak = 1.0
        worst_rel, worst_row = 0.0, -1
        for r, (a, b) in enumerate(zip(vals_ref, vals_new)):
            rel = abs(a - b) / peak
            if math.isnan(rel) or rel > worst_rel:
                worst_rel, worst_row = rel, r
        worst.append((name, h_ref[col], worst_rel,
                      f"row {worst_row}, peak {peak:.3e}"))
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref_dir")
    ap.add_argument("new_dir")
    ap.add_argument("--prefix", default="bssn_profile")
    ap.add_argument("--tol-constraints", type=float, default=1e-12)
    ap.add_argument("--tol-psi4", type=float, default=1e-12)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    checks = []
    cons = f"{args.prefix}_Constraints.dat"
    if os.path.exists(os.path.join(args.ref_dir, cons)):
        checks.append((cons, args.tol_constraints, "constraints"))
    else:
        print(f"FAIL: {cons} missing from {args.ref_dir}")
        return 1
    for p in sorted(glob.glob(os.path.join(args.ref_dir,
                                           f"{args.prefix}_GW_l*.dat"))):
        checks.append((os.path.basename(p), args.tol_psi4, "psi4"))
    if len(checks) == 1:
        print("FAIL: no GW output found; psi4 would go unchecked")
        return 1

    # One scale for every psi4 mode and radius (see module docstring).
    psi4_scale = 0.0
    for fname, _, kind in checks:
        if kind != "psi4":
            continue
        header, data = read_table(os.path.join(args.ref_dir, fname))
        for col, name in enumerate(header):
            if name not in ("TimeStep", "t", "time"):
                psi4_scale = max(psi4_scale,
                                 max(abs(v) for v in column_values(data, col)))
    if not args.quiet:
        print(f"  psi4 scale   {psi4_scale:.3e}  (peak over all modes/radii)")

    failures, worst_by_kind = [], {}
    for fname, tol, kind in checks:
        ref_p = os.path.join(args.ref_dir, fname)
        new_p = os.path.join(args.new_dir, fname)
        if not os.path.exists(new_p):
            failures.append(f"{fname}: missing from {args.new_dir}")
            continue
        scale = psi4_scale if kind == "psi4" else None
        for name, col, rel, where in compare_file(ref_p, new_p, tol, kind,
                                                  scale):
            prev = worst_by_kind.get(kind, (0.0, ""))
            if rel > prev[0]:
                worst_by_kind[kind] = (rel, f"{name}:{col} ({where})")
            if math.isnan(rel) or rel > tol:
                failures.append(f"{name} col {col}: {rel:.3e} > {tol:.0e}"
                                f"  [{where}]")

    if not args.quiet:
        for kind, (rel, where) in sorted(worst_by_kind.items()):
            print(f"  worst {kind:12s} {rel:.3e}   {where}")
    for f in failures:
        print(f"  FAIL {f}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
