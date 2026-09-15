# EKG_GR: native Dendro wavelet AMR with DendroSym-generated equations

This revision replaces the earlier custom EKG context/RK implementation. It is
an extension of the **local experimental BSSN_GR application**, not a new mesh
or time integrator. Add `add_subdirectory(EKG_GR)` after the parent's dependency
setup and existing BSSN application.

**Validation boundary:** the symbolic algebra, reference-emitted C++ kernels,
initial-data solver, and source-adapter fixture tests were run. DendroSym was not
installed in the execution environment, and a complete experimental Dendro
checkout could not be fetched. Therefore the DendroSym emitter, actual native
application build, MPI evolution, remeshing, and transfer have NOT been executed
here. The CMake adapter checks source hooks and stops on an incompatible layout.
See `docs/VALIDATION.md` for exactly what was and was not tested.

## 1. Application architecture

The sibling `BSSN_GR/` and its vacuum target are not edited or linked as a
24-field application. At configure time, `cmake/PrepareNative.cmake` copies your
local BSSN_GR into `build/EKG_GR/native/` and makes a small set of explicit
changes:

* append `U_SCALARPHI=24`, `U_SCALARPI=25`; compile the whole private application
  with `BSSN_NUM_VARS=26`;
* register new initial-data IDs and the `[EKG]` TOML table;
* dispatch the native per-block RHS and constraint entry points to generated EKG
  kernels;
* wrap native wavelet refinement tests with optional, fixed scalar normalization.

The adapter also guards zero-valued optional event frequencies (for example
GW extraction and checkpoint writing), preserving the native schedule for
nonzero values without modulo-by-zero.

The native `bssnCtx.cpp`, `bssngr_main.cpp`, DVector allocation, MPI exchange,
`is_remesh`, mesh construction, `grid_transfer`, remesh-and-transfer sequence,
and ETS mesh synchronization remain the basis of the application. Its native
RK implementation is used, not the bespoke RK4 loop from the previous package.
**All 26 fields go through the native stage storage and intergrid transfer.**

The source adapter does not parse arbitrary C++. It targets the experimental
layout described in `docs/PROVENANCE.md`; missing/ambiguous hooks cause a CMake
error. Compare those hooks with your checkout rather than forcing a mismatch.

The original entry-point bodies remain below an early return in the copied
`rhs.cpp` and `physcon.cpp`. They provide readable reference code and preserve
helper definitions; they are not executed a second time.

### Persistent, editable native application copy

Configuration recreates `build/EKG_GR/native`, so do not keep permanent edits
there. To take full ownership of the expanded source, copy it once:

```sh
cp -R build-ekg/EKG_GR/native /absolute/path/to/my-ekg-native
cmake -S . -B build-ekg \
  -DEKG_NATIVE_SOURCE_DIR=/absolute/path/to/my-ekg-native
```

This skips automatic source adaptation and compiles your edited 26-field native
source. `EKG_GR/src`, `include`, and `CodeGen` always remain ordinary editable
source. No Python installer or source-patching step is required.

## 2. Compile-time choices

| Option | Default | Meaning |
|---|---|---|
| `EKG_BACKREACTION` | ON | Generate Einstein equations including scalar stress-energy |
| `EKG_FREEZE_GEOMETRY` | OFF | Zero all 24 geometric RHS entries; requires backreaction OFF |
| `EKG_FUSED_RHS` | ON | Use one globally CSE-optimized 26-output GR+matter kernel |
| `EKG_POTENTIAL` | massive | `massive` or `axion`; mass and decay constant remain TOML parameters |
| `EKG_GAUGE` | puncture | `puncture` or `synchronous` |
| `EKG_ADVECTION` | centered | `centered` or native `upwind` derivative inputs |
| `EKG_REGENERATE` | ON | Run DendroSym during the build |
| `EKG_BUILD_TESTS` | ON | Build independent generated-kernel/initial-data tests |

`BACKREACTION=OFF, FREEZE_GEOMETRY=OFF` evolves vacuum BSSN plus a passive scalar.
`BACKREACTION=OFF, FREEZE_GEOMETRY=ON` evolves only the scalar on fixed geometry.
Disabling backreaction is therefore not the same as freezing the geometry.

`FUSED_RHS=OFF` is a debugging alternative using separately generated vacuum and
matter kernels plus explicit source addition. `B^i` receives the same matter
correction as the full Gamma-driver equation. The fused path must NOT add those
sources again. Scalar mass/self-interaction physics is not removed in test-field
mode; only its gravitational backreaction is removed.

## 3. Generate the equations

Use the Python interpreter in which you installed DendroSym:

```sh
python EKG_GR/CodeGen/generate_ekg.py \
  --out-dir EKG_GR/generated \
  --backreaction on \
  --potential massive \
  --gauge puncture \
  --advection centered
```

The two files to edit are:

* `CodeGen/ekg_equations.py`: canonical scalar, stress-energy, BSSN, gauge and
  constraint expressions;
* `CodeGen/generate_ekg.py`: DendroSym CSE/emission, array aliases, derivative
  dependency manifest, configuration checks and provenance.

The generator calls DendroSym's low-level `construct_cse_from_list` and
`generate_cpu_preextracted` API on SymPy expressions. This avoids depending on
NRConfig registration/global-metric initialization order. The input derivatives
are named explicitly and the finite-difference dependency list is generated from
those same symbols. See `CodeGen/README.md`.

The normal CMake build invokes that command automatically. To use your manual
output, configure with:

```sh
-DEKG_REGENERATE=OFF -DEKG_GENERATED_DIR=/absolute/path/to/EKG_GR/generated
```

Keep all generation options consistent with the CMake choices. Static assertions
reject mismatched coupling/gauge/potential/advection flags. Do not use the
`--reference-emitter` testing option for production: `ekgSolver` rejects it.

## 4. Build the fully coupled application

Keep the compiler and dependency settings of your working **CPU** vacuum build.
This first integration uses the native legacy explicit FD6 interface; it does
not port the GPU, SIMD, or new derivative backends.

```cmake
# Root CMakeLists.txt, after the existing dependency setup:
add_subdirectory(BSSN_GR)  # existing line
add_subdirectory(EKG_GR)   # new line
```

Example from the repository root:

```sh
cmake -S . -B build-ekg \
  -DEKG_BACKREACTION=ON \
  -DEKG_FREEZE_GEOMETRY=OFF \
  -DEKG_FUSED_RHS=ON \
  -DEKG_POTENTIAL=massive \
  -DEKG_GAUGE=puncture \
  -DEKG_ADVECTION=centered \
  -DDENDRO_USE_NEW_DERIVS=OFF \
  -DPython3_EXECUTABLE=/absolute/path/to/your/python

cmake --build build-ekg \
  --target bssnSolver ekgSolver ekgKernelTests ekgInitialDataTests -j 4

./build-ekg/EKG_GR/ekgKernelTests
./build-ekg/EKG_GR/ekgInitialDataTests
```

The example public vacuum target is `bssnSolver`; keep your local name if it
already differs. The new executable is `build-ekg/EKG_GR/ekgSolver`.

Do not link the original `bssn_common` or change the original BSSN field count.
The private native application must be compiled entirely with the EKG headers.
The generated `ekg_field_asserts.inc` checks all 26 enum positions. Also inspect
`native/EKG_FIXED_SIZE_AUDIT.txt` for literal 24-entry buffers in your checkout.
Known refinement/output index arrays are extended explicitly; unrelated arrays
and constants are not blindly rewritten.

## 5. The first run is fully coupled AND wavelet-adaptive

```sh
OMP_NUM_THREADS=1 mpirun -np 1 \
  ./build-ekg/EKG_GR/ekgSolver \
  EKG_GR/pars/coupled_massive_amr.toml 1
```

The final `1` selects the native uniform-step ETS, including its native remeshing
cycle. It does NOT select a uniform spatial grid. Local time stepping is not yet
validated by this integration and is rejected. The parameter file selects native
RK4; the application retains the native RK infrastructure rather than replacing
its methods.

Important settings in the supplied full file:

```toml
BSSN_ID_TYPE = 103
BSSN_ENABLE_BLOCK_ADAPTIVITY = 0
BSSN_REFINEMENT_MODE = 0
BSSN_REMESH_TEST_FREQ = 5
BSSN_REMESH_TEST_FREQ_AFTER_MERGER = 5
BSSN_NUM_REFINE_VARS = 5
BSSN_REFINE_VARIABLE_INDICES = [0,1,2,24,25]
BSSN_WAVELET_TOL = 2.0e-6
BSSN_ASYNC_COMM_K = 2
```

`0` for `ENABLE_BLOCK_ADAPTIVITY` follows the native wavelet/function-to-octree
path rather than the block-initialization path used by the earlier flat tests.
Both scalar fields are mandatory refinement inputs when remeshing is enabled.
The amplitude can be dynamically important while small compared with geometric
variables. Optional `[EKG] refine_scale_phi` and `refine_scale_pi` divide ONLY a
read-only view passed to native WAMR tests. They never rescale evolved fields.
Use fixed physical reference scales, not a time-dependent local maximum.
Defaults are both one, i.e. unchanged native absolute wavelet tolerance.

The initial function-to-octree seeding uses physical variables. The native
initialization convergence/remeshing tests and subsequent time-evolution WAMR
use the wrapper. With non-unit scales, verify initial mesh convergence as well
as subsequent remeshing. The defaults avoid this distinction.

`BSSN_ASYNC_COMM_K=2` divides 26. Do not retain a grouping of four from a 24-field
vacuum setup. The new C++ runtime checks this.

### ID 103: constraint-solved compact scalar pulse

This is NOT a nonzero scalar pasted onto Minkowski. It sets Pi=K=At=0,
gt=identity, chi=psi^(-4), and solves the nonlinear radial Hamiltonian equation:

```
psi'' + (2/r) psi' = -pi psi (phi')^2 - 2 pi psi^5 V(phi)
psi'(0) = 0
psi(R) + R psi'(R) = 1
```

The compact smooth scalar profile vanishes outside `support_radius`; `R` is
larger than that support. Beyond R, the conformal factor uses the asymptotically
flat exterior `1+M/(2r)`. The momentum constraint is satisfied by time symmetry.

This is a coupled, inhomogeneous code/AMR test with nonzero spatial stress. It is
not an initial black hole or an astrophysical dark-matter cloud model. The radial
solve and interpolation must be refined independently of the octree.

## 6. Other supplied parameter files

| File | ID | Build choices |
|---|---|---|
| `coupled_massive_amr.toml` | 103 | coupled, massive, puncture, dynamic WAMR |
| `coupled_axion_amr.toml` | 103 | coupled, axion, puncture, dynamic WAMR |
| `flat_gaussian_amr.toml` | 100 | backreaction OFF, freeze ON, massive with mu=0 |
| `flat_massive_plane.toml` | 101 | backreaction OFF, freeze ON, massive |
| `flrw_coupled.toml` | 102 | coupled, synchronous, massless, analytic boundary |

Use separate build directories for compile-time configurations:

```sh
cmake -S . -B build-ekg-testfield \
  -DEKG_BACKREACTION=OFF -DEKG_FREEZE_GEOMETRY=ON \
  -DEKG_POTENTIAL=massive -DDENDRO_USE_NEW_DERIVS=OFF \
  -DPython3_EXECUTABLE=/absolute/path/to/your/python
cmake --build build-ekg-testfield --target ekgSolver -j 4
mpirun -np 1 ./build-ekg-testfield/EKG_GR/ekgSolver \
  EKG_GR/pars/flat_gaussian_amr.toml 1
```

Use the full supplied TOML files: the inherited parameter reader still requires
its usual geometry/BH/TwoPunctures entries. Those placeholders do not turn ID103
into a black-hole initial-data solve. Geometry IDs below 100 remain native IDs
and initialize phi=Pi=0 unless you explicitly add a scalar-aware, constraint-
consistent initial-data path.

## 7. What must be checked on your machine

1. Run the generated C++ unit tests with YOUR DendroSym emission. They compare
   fused and split kernels over all 26 outputs, check scalar signs, the source
   terms and independent constraint specializations.
2. Run scalar amplitude zero and compare the generated vacuum-limit evolution
   against a native standard-gauge vacuum run. Match gauge, derivative, boundary
   and dissipation settings before interpreting differences.
3. Evolve ID103 on a fixed mesh first, then enable remeshing. Inspect initial H
   and M_i and their spatial/time convergence, not merely stability.
4. Confirm actual refine/coarsen events, all-26-field transfer, ETS mesh sync,
   and CFL changes in your native logs. Include a tolerance sweep such as
   8e-6, 4e-6, 2e-6, 1e-6 and increase MAXDEPTH to check saturation.
5. Repeat with two MPI ranks; compare constraints and observables, not only
   pictures or mesh counts. Compare fused and split builds as well.

A source-level hook test or pointwise algebra test is NOT an AMR evolution test.
No published-accuracy or production-readiness claim is made here.

## 8. Remaining scope boundaries

* CPU, explicit FD6, and native ETS are targeted. LTS/GPU/SIMD/new derivatives
  require their own source integration and tests.
* Native checkpoint writing paths are retained, but restore is deliberately
  rejected until EKG schema and compile-option compatibility checks are wired.
  Never restore an old 24-field checkpoint into this application.
* The generated constraint callback supplies H and the three covariant M_i.
  Native Psi4 slots are zeroed and must NOT be interpreted as waveforms;
  runtime validation disables GW extraction and excludes those slots from VTU.
  Matter-corrected waveform extraction and horizon analysis are subsequent work.
* Approximate Sommerfeld boundaries in the compact-pulse runs must be tested by
  increasing the domain. They are not exact nonreflecting massive-field BCs.
* ID103 is a weak-positive-branch radial elliptic solve. Failure to bracket a
  positive solution is reported, not disguised as constraint-satisfying data.
* The final source compatibility test is your local native build. The checked
  hooks are a guard against source drift, not a substitute for compilation.
