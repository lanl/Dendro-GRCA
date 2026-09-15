# DendroSym generation interface

Edit `ekg_equations.py` for physics and `generate_ekg.py` for emission.
No evolution RHS is hand-coded in `src/ekg_rhs.cpp`; that file orchestrates native
finite differences, generated fragments, boundary conditions and KO dissipation.

## Invocation

```sh
python generate_ekg.py --out-dir ../generated \
  --backreaction on --potential massive --gauge puncture --advection centered
```

The default emitter imports your installed DendroSym and calls:

```python
from dendrosym import codegen
cse = codegen.construct_cse_from_list(expressions)
cpp = codegen.generate_cpu_preextracted(
    cse, names, '[pp]', original_operation_count,
    dtype='double', use_const=True,
)
```

The interface is the same low-level path discussed earlier. Its execution was
not available for testing in the delivery environment; the script checks that
the two functions exist and that the emitter returns a C++ string. An API mismatch
stops generation instead of silently substituting another emitter.

The generator uses plain SymPy symbols for fields and derivative arrays. This is
compatible with the low-level emitter without constructing NRConfig, registering
evolution fields, or initializing DendroSym's global metric state. The emitted
formulation uses the experimental BSSN variable order, not the private EMDA order.

## Outputs

| Output | Meaning |
|---|---|
| `ekg_full_rhs.inc` | All 26 RHS components passed together to CSE |
| `ekg_vacuum_rhs.inc` | All 24 vacuum geometry/gauge RHS components |
| `ekg_matter_rhs.inc` | Phi/Pi RHS, point-local Einstein sources and stress diagnostics |
| `ekg_apply_sources.inc` | Source additions used ONLY in split mode |
| `ekg_constraints.inc` | Geometrically evaluated H and three covariant M_i |
| `ekg_input_aliases.inc` | Native input-array aliases |
| `ekg_output_aliases.inc` | Native RHS-array aliases |
| `ekg_derivative_aliases.inc` | Derivative-workspace aliases |
| `ekg_derivative_specs.h` | Finite-difference dependencies, fields, directions and parents |
| `ekg_generated_config.h` | Coupling, potential, gauge and advection configuration tags |
| `ekg_field_asserts.inc` | Checks for the native 26-field enum and count |
| `ekg_potential.h` | Same symbolic potential for the radial elliptic initial-data solve |
| `manifest.json` | Emitter/version/configuration/source-hash provenance |

Include fragments inside the point loop, not as separately compiled sources.
The matter fragment uses local source variables, not permanent source arrays.
Fused and split temporary scopes are deliberately separate.

The full coupled expression is assembled BEFORE CSE:

```
K_full = K_vac + scalar_K_source
At_full = At_vac + scalar_At_source
Gt_full = Gt_vac + scalar_Gt_source
B_full uses Gt_full
full = [24 full geometry RHS, phi RHS, Pi RHS]
```

There is no second matter-source addition after `ekg_full_rhs.inc`.

## Testing fallback is not production generation

`--reference-emitter` emits the SAME symbolic expressions using direct SymPy CSE
and ccode. It exists to exercise algebra and C++ syntax without DendroSym and
records this distinction in generated headers and the manifest. The production
solver rejects such output. Local validation described in `docs/VALIDATION.md`
used this fallback, not the uninstalled DendroSym emitter.

The potential choice is compile time; `mu` and `fa` are runtime values. For axions,
V is written as `2 mu^2 fa^2 sin(phi/(2 fa))^2` to avoid small-field cancellation.
The derivative is the equivalent `mu^2 fa sin(phi/fa)`.
