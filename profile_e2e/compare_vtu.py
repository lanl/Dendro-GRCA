#!/usr/bin/env python3
"""MAE gate: compare two DendroGR runs' fields point by point.

Companion to run_profile.sh, which measures speed but not correctness. Build two RHS
variants, run both on q1.verify.par.toml (BSSN_IO_OUTPUT_FREQ=1, so VTU lands every step),
then diff them here. Same parfile means identical meshes, so the comparison is
point-for-point with no interpolation.

  mkdir -p run_prod/dat run_casc/dat && cp q1.verify.par.toml run_prod/ && ...
  (cd run_prod && mpirun -np 4 ../build_production/BSSN_GR/bssnSolver q1.verify.par.toml)
  (cd run_casc && mpirun -np 4 ../build_cascade/BSSN_GR/bssnSolver    q1.verify.par.toml)
  ./compare_vtu.py run_prod run_casc

Expect a few ULP (~1e-16 on evolved fields, ~1e-15 on C_HAM, which amplifies as a derived
diagnostic) from the different FP operation ordering, flat across steps. Growth that does
not saturate is a real divergence, not roundoff. Requires python-vtk.
"""
import sys, os, glob, re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

GATE = 1e-10  # far above roundoff, far below anything structural


def read(path):
    r = (vtk.vtkXMLPUnstructuredGridReader() if path.endswith(".pvtu")
         else vtk.vtkXMLUnstructuredGridReader())
    r.SetFileName(path)
    r.Update()
    d = r.GetOutput()
    pts = vtk_to_numpy(d.GetPoints().GetData()) if d.GetPoints() else None
    pd = d.GetPointData()
    return {pd.GetArrayName(i): vtk_to_numpy(pd.GetArray(i)).astype(np.float64)
            for i in range(pd.GetNumberOfArrays())}, pts


def step_of(p):
    m = re.search(r"(\d+)\.p?vtu$", os.path.basename(p))
    return int(m.group(1)) if m else -1


def find(d):
    f = sorted(glob.glob(os.path.join(d, "**", "*.pvtu"), recursive=True), key=step_of)
    return f or sorted(glob.glob(os.path.join(d, "**", "*.vtu"), recursive=True), key=step_of)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    A, B = sys.argv[1], sys.argv[2]
    fa, fb = find(A), find(B)
    if not fa or not fb:
        sys.exit(f"no VTU output found (A={len(fa)}, B={len(fb)}); check BSSN_IO_OUTPUT_FREQ")

    sa = {step_of(p): p for p in fa}
    sb = {step_of(p): p for p in fb}
    common = sorted(set(sa) & set(sb))
    if not common:
        sys.exit(f"no common steps. A={sorted(sa)[:8]} B={sorted(sb)[:8]}")
    print(f"A: {len(fa)} files, B: {len(fb)} files, {len(common)} common steps: {common}\n")

    worst = (0.0, None, None)
    per_step = []
    for st in common:
        da, pa = read(sa[st])
        db, pb = read(sb[st])
        if pa is not None and pb is not None:
            if pa.shape != pb.shape:
                print(f"step {st}: GRID MISMATCH {pa.shape} vs {pb.shape} -- not comparable")
                continue
            gd = np.max(np.abs(pa - pb))
            if gd > 0:
                print(f"step {st}: WARNING grid coords differ by {gd:.3e}")
        names = sorted(set(da) & set(db))
        npts = len(next(iter(da.values())))
        print(f"--- step {st}  ({len(names)} arrays, {npts} points) ---")
        print(f"  {'variable':<16} {'MAE':>12} {'max|diff|':>12} {'max|A|':>12} {'rel MAE':>12}")
        su = sc = 0.0
        for n in names:
            a, b = da[n], db[n]
            if a.shape != b.shape:
                print(f"  {n:<16} shape mismatch {a.shape} vs {b.shape}")
                continue
            d = np.abs(a - b)
            mae, mx, s = d.mean(), d.max(), np.abs(a).max()
            print(f"  {n:<16} {mae:12.3e} {mx:12.3e} {s:12.3e} {mae/s if s else 0.0:12.3e}")
            if n.startswith("U_"):
                su = max(su, mx)
            elif n.startswith("C_"):
                sc = max(sc, mx)
            if mx > worst[0]:
                worst = (mx, n, st)
        per_step.append((st, su, sc))
        print()

    # A flat profile is roundoff; one that keeps climbing is a real divergence.
    print(f"{'step':>5} {'max|diff| evolved':>18} {'max|diff| constraints':>22}")
    for st, su, sc in per_step:
        print(f"{st:>5} {su:>18.3e} {sc:>22.3e}")
    print(f"\nWORST max|diff|: {worst[0]:.3e}  ({worst[1]}, step {worst[2]})")
    ok = worst[0] < GATE
    print("GATE:", "PASS" if ok else f"FAIL (above {GATE:g})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
