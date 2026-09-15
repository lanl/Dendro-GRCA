## Compile-time choices

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

## Generate the equations

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

## Build on LANL Darwin (current working configuration)

The Darwin login environment used here has a system Git/HTTPS OpenSSL--Kerberos 
library mismatch (`git-remote-https` can fail in `libk5crypto.so.3` while resolving
`EVP_KDF_ctrl`).  

The root project remains unchanged except for adding the EKG application after
the ordinary vacuum application:

```cmake
# Root CMakeLists.txt, after dependency setup:
add_subdirectory(BSSN_GR)  # existing vacuum application
add_subdirectory(EKG_GR)   # EKG application
```

4.1 Local dependency layout

A convenient Darwin layout is:

```text
/<your project space>/
  Dendro-GRCA/             # this repository
  Dendro-5.01/             # local Dendro source expected by Dendro-GRCA
  deps/
    spdlog/                # v1.14.1
    toml11/                # v4.4.0
    libxsmm/               # revision requested/validated with Dendro-5.01
```

The `spdlog` and `toml11` versions above are the versions requested by the
current Dendro-5.01 CMake configuration.  Dendro-5.01 currently requests the
`main` branch of LIBXSMM; record the exact commit used for every production
build.  If the local Dendro checkout is pinned by the parent project, use that
expected branch/commit rather than an unrelated Dendro-5.01 checkout.

Because HTTPS Git is broken in the affected Darwin environment, clone/update
these repositories using a working Git transport (for example GitHub SSH) or
another approved mechanism.  Once the source trees exist locally, the CMake
build itself does not need network access.

It is useful to record the dependency revisions before configuring:

```sh
git -C /<your project space>/Dendro-GRCA rev-parse HEAD
git -C /<your project space>/Dendro-5.01 rev-parse HEAD
git -C /<your project space>/deps/spdlog describe --tags --always
git -C /<your project space>/deps/toml11 describe --tags --always
git -C /<your project space>/deps/libxsmm rev-parse HEAD
```

### Generate equations before the HPC build

The recommended Darwin workflow is to generate the DendroSym kernels separately
and compile the generated files on Darwin.  This avoids making DendroSym a build
time dependency on the compute system.  From an environment with DendroSym
installed:

```sh
python EKG_GR/CodeGen/generate_ekg.py \
  --out-dir EKG_GR/generated \
  --backreaction on \
  --potential massive \
  --gauge puncture \
  --advection centered
```

Copy or commit the resulting `EKG_GR/generated/` directory into the Darwin
checkout.  The Darwin CMake configuration below therefore uses
`EKG_REGENERATE=OFF`.

### Configure the fully coupled EKG build on Darwin

From the repository root:

```sh
cd /<your project space>/Dendro-GRCA
rm -rf build-ekg

cmake -S . -B build-ekg \
  -DFETCHCONTENT_SOURCE_DIR_DENDROLIB=/<your project space>/Dendro-5.01 \
  -DFETCHCONTENT_SOURCE_DIR_SPDLOG=/<your project space>/deps/spdlog \
  -DFETCHCONTENT_SOURCE_DIR_TOML11=/<your project space>/deps/toml11 \
  -DUSE_LOCAL_XSMM=ON \
  -DLOCAL_XSMM_PATH=/<your project space>/deps/libxsmm \
  -DEKG_REGENERATE=OFF \
  -DEKG_GENERATED_DIR="$PWD/EKG_GR/generated" \
  -DEKG_BACKREACTION=ON \
  -DEKG_FREEZE_GEOMETRY=OFF \
  -DEKG_FUSED_RHS=ON \
  -DEKG_POTENTIAL=massive \
  -DEKG_GAUGE=puncture \
  -DEKG_ADVECTION=centered \
  -DDENDRO_USE_NEW_DERIVS=OFF
```

The local source overrides are important on this Darwin software stack.  Without
them CMake attempts HTTPS clones and can fail before either BSSN or EKG is
configured, e.g.

```text
/usr/libexec/git-core/git-remote-https: symbol lookup error:
/lib64/libk5crypto.so.3: undefined symbol: EVP_KDF_ctrl, version OPENSSL_1_1_1b
```

This error is a dependency-fetch/environment problem, not an EKG equation or
Dendro compilation error.

### Build vacuum GR and EKG independently

Build the vacuum target first, then the EKG target:

```sh
cmake --build build-ekg --target bssnSolver -j 8
cmake --build build-ekg --target ekgSolver  -j 8
```

Or build both in one command:

```sh
cmake --build build-ekg --target bssnSolver ekgSolver -j 8
```

If the optional EKG tests are enabled in the current checkout, also build/run:

```sh
cmake --build build-ekg --target ekgKernelTests ekgInitialDataTests -j 8
./build-ekg/EKG_GR/ekgKernelTests
./build-ekg/EKG_GR/ekgInitialDataTests
```

The expected executables are:

```text
build-ekg/BSSN_GR/bssnSolver
build-ekg/EKG_GR/ekgSolver
```

The vacuum application remains a 24-field application.  Do not change the
original BSSN field count or link its 24-field application objects into the EKG
solver.  The private EKG native application is compiled consistently with all 26
fields.  The generated `ekg_field_asserts.inc` checks their enum positions.
Also inspect `build-ekg/EKG_GR/native/EKG_FIXED_SIZE_AUDIT.txt` after
configuration for any checkout-specific literal 24-entry buffers.

### Reconfiguration rules

Use a fresh build directory whenever changing a compile-time physics choice,
especially:

```text
EKG_BACKREACTION
EKG_FREEZE_GEOMETRY
EKG_FUSED_RHS
EKG_POTENTIAL
EKG_GAUGE
EKG_ADVECTION
```

For example, do not reuse a massive/coupled build directory for an axion or
test-field build.  Also regenerate `EKG_GR/generated/` whenever the symbolic
physics/gauge/advection choices change, and keep those choices identical to the
CMake configuration.

## The first run is fully coupled AND wavelet-adaptive

Run from the Dendro-GRCA repository root so the relative TOML path is resolved
consistently.  Start with one MPI rank as a smoke/validation run; use the Darwin
batch scheduler and allocated compute resources for nontrivial production runs
rather than running them on a front-end/login node.

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

## Other supplied parameter files

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

## What must be checked on your machine

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

## Remaining scope boundaries

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

