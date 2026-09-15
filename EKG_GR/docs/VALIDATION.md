# Validation record and required local checks

This record distinguishes equation/code-generation tests from a native Dendro
build. No native EKG evolution or actual wavelet AMR run was performed here.

## Environment limitations

The attached prior EKG package and private EMDA archive were available locally.
The requested public experimental branch could not be fetched in this execution
environment. The current implementation therefore consumes the user's local
BSSN_GR source at CMake configuration rather than pretending to include a
verified snapshot. DendroSym was not installed locally. Default generation uses
DendroSym on the user's machine; the tests below used the explicitly marked
SymPy reference emitter on the same mathematical expressions.

## Tests executed

### Symbolic identities (PASS)

`tests/test_symbolic.py` compared the conformal BSSN Ricci expression reconstructed
from the metric with the direct Christoffel Ricci expression on 12 general
positive-definite metric jets. It also checked the K source, all three Gamma-
driver source corrections, and the flat massive Klein-Gordon signs/operator.
These tests do not use mesh data or a Dendro time step.

### Generated C++ kernel algebra (PASS)

Reference-emitted C++ kernels were compiled with GCC and tested for:

* all 26 fused versus split RHS components on 30 general metric/scalar jets;
* the scalar propagation and mass-potential sign;
* an independent conformally flat Hamiltonian and all three momentum constraints;
* coupled massive/centered/puncture, coupled axion/upwind/puncture,
  passive massive/centered, and coupled massive/centered/synchronous variants.

For the default coupled massive build, the maximum absolute fused/split
discrepancy was 5.55112e-17. This tests algebraic consistency, NOT agreement with
the checked-in native vacuum kernel or a physical time evolution.

### Coupled radial initial data (PASS)

`tests/test_initial_data.cpp` exercised the nonlinear radial Hamiltonian shooting
solver, exterior matching, interpolation, and exact homogeneous FLRW formulas.
For its weak massive-scalar test at 4096 radial intervals:

```
ADM mass                 0.00348583
Robin residual          -1.93179e-14
radial finite-difference max |H|  1.04033e-10
```

The last value is a radial independent check, not the constraint residual on a
3D octree. Repeat radial and octree refinement independently on the target code.

### CMake adapter fixture (PASS)

`cmake -P tests/test_native_adapter.cmake` exercised minimal synthetic native
source fixtures. It checked zero-frequency event guards, field count, argument
forwarding, known literal
field-index array extension and preservation of the original source.
**This fixture is not a native Dendro checkout and does not establish ABI/API
compatibility.**

### TOML and Python syntax (PASS)

All five supplied parameter files parsed with Python's TOML parser; declared
refinement/output array counts and scalar indices were checked. Generator and
symbolic-test Python files passed byte compilation. This is not a run of the
native C++ parameter reader.

## Not executed / not established

* installed DendroSym low-level emission and its generated C++ compilation;
* configuration/compilation/linking of the actual experimental ekgSolver;
* runtime execution of the original native parameter reader with these inputs;
* native MPI ghost exchange and intergrid transfer for all 26 fields;
* native initialization WAMR convergence, dynamic remeshing/coarsening, ETS mesh
  synchronization, or CFL updates;
* a convergent full coupled time evolution or black-hole evolution;
* source equality with every optional gauge/modification in an arbitrary native
  vacuum build;
* checkpoint restore, LTS, GPU/SIMD/new derivative backends, horizon quantities
  or matter-corrected waveform extraction.

## Reproduce independent tests

With your DendroSym installed (preferred):

```sh
cmake -S EKG_GR -B ekg-algebra -DEKG_BUILD_SOLVER=OFF
cmake --build ekg-algebra -j 4
ctest --test-dir ekg-algebra --output-on-failure
```

Only to reproduce this environment's reference-emitter algebra test:

```sh
cmake -S EKG_GR -B ekg-reference-checks \
  -DEKG_BUILD_SOLVER=OFF -DEKG_REFERENCE_EMITTER=ON
cmake --build ekg-reference-checks -j 4
ctest --test-dir ekg-reference-checks --output-on-failure
```

The reference-emitter option is explicitly prohibited for ekgSolver. Passing
these tests is a prerequisite for, not a replacement for, the native run sequence
in README.md. Representative original test logs are in `docs/test_logs/`.
