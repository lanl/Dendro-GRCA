# Source basis and what was changed

## Requested public application

Target: dfvankomen/Dendro-GR, experimental branch, sibling BSSN_GR application.
Its live contents could not be retrieved in this execution environment. The
following paths/hook assumptions were available in the prior conversation's
source-review artifacts and are checked against the user's LOCAL files during
configuration:

```
BSSN_GR/CMakeLists.txt
BSSN_GR/include/grDef.h
BSSN_GR/include/parameters.h
BSSN_GR/src/bssngr_main.cpp
BSSN_GR/src/bssnCtx.cpp
BSSN_GR/src/grUtils.cpp
BSSN_GR/src/parameters.cpp
BSSN_GR/src/rhs.cpp
BSSN_GR/src/physcon.cpp
```

No public commit ID is fabricated here. CMake records the local commit and hashes
of the source files actually used in `ekg_build_manifest.txt` and the private
native copy's `EKG_SOURCE_HASHES.txt`.

## Attached archive inspected

The supplied private `emda-gr-main.zip` was inspected locally, particularly:

```
bssn_eqns_config.py
emdaSolver/src/emdaCtx.cpp
emdaSolver/include/dataUtils.h
emdaSolver/src/dataUtils.cpp
emdaSolver/src/parameters.cpp
emdaSolver/include/derivs.h
emdaSolver/CMakeLists.txt
```

These establish the attached application's native context, refinement wrappers,
transfer and symbolic BSSN pattern. They do not establish that every signature is
identical in the requested experimental fork. The private archive is not bundled
into this delivery and its 36-field ordering is NOT adopted.

## Prior EKG package

The previous delivered EKG_GR README explicitly described a separate hand-coded
RK4/context and a 24-field private geometry object library. That architecture is
removed here. The new code compiles the native application/context itself with
26 fields, and all evolution algebra is in the symbolic generator.

## Native source adaptation

`cmake/PrepareNative.cmake` makes an inspectable, configure-time private source
copy. It does not modify the original sibling application. It performs only:

1. field count and enum/name extension;
2. new ID and parameter registration;
3. native main validation and avoidance of BH-only actions for non-BH IDs;
4. read-only scalar-normalization hooks around native WAMR tests;
5. per-block generated RHS and constraint dispatch;
6. known literal field-index array extension and a remaining fixed-size audit;
7. zero-frequency event guards for native optional output/remeshing schedules.

Mesh ownership, initialization refinement iteration, actual remeshing, field
transfer, MPI-context rebuilding and native ETS/RK stay in the local native
source. Hooks require the expected functions and detect ambiguity. A successful
hook application is not a successful compile: local build/evolution validation
remains necessary.

## Files to review after configuration

```
build/EKG_GR/native/EKG_NATIVE_ADAPTER.txt
build/EKG_GR/native/EKG_SOURCE_HASHES.txt
build/EKG_GR/native/EKG_FIXED_SIZE_AUDIT.txt
build/EKG_GR/ekg_build_manifest.txt
build/EKG_GR/generated/manifest.json
```

Inspect remaining literal [24] occurrences. Some may belong to unreachable
vacuum generated code, while a live full-state buffer must be changed to 26 or
BSSN_NUM_VARS. Do not assume every numeric 24 is a field count. The generated
compile-time enum assertions protect the expected field order, not every
arbitrary local source edit.

The archive contains all new EKG implementation files and the adapter. It does
not contain a vendored public BSSN tree. To own the expanded native source as
ordinary version-controlled files, follow the EKG_NATIVE_SOURCE_DIR workflow in
README.md.
