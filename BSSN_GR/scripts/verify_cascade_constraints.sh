#!/usr/bin/env bash
# Gate for the IR polynomial-cascade constraint/psi4 kernel: builds the scalar
# reference and every vector variant the target arch can run, then checks that
# each agrees with the reference and that no block silently fell back to scalar.
set -uo pipefail

REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
PARFILE=${PARFILE:?set PARFILE to a BSSN parameter file}
WORK=${WORK:-$(mktemp -d -t cascade-cons-XXXXXX)}
ARCH=${ARCH:-native}
NPROC=${NPROC:-4}
JOBS=${JOBS:-8}
DENDROLIB=${DENDROLIB:-$REPO/../dendrolib_dfvk_copy}
MPIRUN=${MPIRUN:-mpirun}

PARFILE=$(readlink -f "$PARFILE")
mkdir -p "$WORK"
echo "verify_cascade_constraints: arch=$ARCH nproc=$NPROC work=$WORK"

CMAKE_COMMON=(
    -G Ninja
    -DDENDRO_dendrolib_DIR="$DENDROLIB"
    -DCMAKE_BUILD_TYPE=Release
    -DCPU_ARCH="$ARCH"
    -DBSSN_COMPUTE_CONSTRAINTS=ON
    -DBSSN_EXTRACT_GRAVITATIONAL_WAVES=ON
    -DBSSN_PHYSCON_BLOCK_HISTO=ON
    -DDVEC_ZERO_ALLOC=ON
    -DOCT2BLK_COARSEST_LEV=0
    -DBSSN_ENABLE_SSL_HD=ON
    -DBSSN_USE_6TH_ORDER_DERIVS=ON
)

variants=(scalar avx2)
# The intrinsics do not compile without the ISA, so ask the compiler, not the host.
if echo | ${CXX:-g++} -march="$ARCH" -dM -E - 2>/dev/null | grep -q __AVX512F__; then
    variants+=(avx512)
else
    echo "  note: $ARCH has no AVX-512; skipping the 8-wide variant"
fi

flags_for() {
    case $1 in
    scalar) : ;;
    avx2) echo -DBSSN_USE_CASCADE_CONSTRAINTS_AVX=ON ;;
    avx512) echo -DBSSN_USE_CASCADE_CONSTRAINTS_AVX512=ON ;;
    esac
}

fail=0
for v in "${variants[@]}"; do
    b=$WORK/build_$v
    read -r -a extra <<<"$(flags_for "$v")"
    cmake -S "$REPO" -B "$b" "${CMAKE_COMMON[@]}" "${extra[@]}" \
        >"$WORK/$v.cfg.log" 2>&1 || {
        echo "  FAIL configure $v"; tail -5 "$WORK/$v.cfg.log"; exit 1
    }
    cmake --build "$b" -j"$JOBS" --target bssnSolver >"$WORK/$v.build.log" 2>&1
    # cmake can exit 0 with no binary if the link was interrupted.
    [ -f "$b/BSSN_GR/bssnSolver" ] || {
        echo "  FAIL build $v"; grep -m5 "error:" "$WORK/$v.build.log"; exit 1
    }
    grep -m1 "BSSN constraint kernel:" "$WORK/$v.cfg.log" | sed 's/^-- /  /'
done

for v in "${variants[@]}"; do
    for w in "${variants[@]}"; do
        [ "$v" \< "$w" ] || continue
        a=$(sha256sum "$WORK/build_$v/BSSN_GR/bssnSolver" | cut -d' ' -f1)
        c=$(sha256sum "$WORK/build_$w/BSSN_GR/bssnSolver" | cut -d' ' -f1)
        [ "$a" = "$c" ] && {
            echo "  FAIL $v and $w built the same binary; the flag did nothing"
            fail=1
        }
    done
done

for v in "${variants[@]}"; do
    d=$WORK/run_$v
    rm -rf "$d"; mkdir -p "$d/aeh"
    cp "$PARFILE" "$d/par.toml"
    (cd "$d" && OMP_NUM_THREADS=1 "$MPIRUN" -np "$NPROC" \
        "$WORK/build_$v/BSSN_GR/bssnSolver" par.toml >run.log 2>err.log) || {
        echo "  FAIL run $v"; tail -20 "$d/err.log"; exit 1
    }
done

for v in "${variants[@]}"; do
    [ "$v" = scalar ] && continue
    tally=$(grep -h '^\[PBH\]' "$WORK/run_$v/err.log" | awk '{
        split($2,a,"="); split($3,b,"="); split($4,c,"=");
        B[a[2]]+=b[2]; S[a[2]]+=c[2] }
        END { for (k in B) printf "nx=%s blocks=%d scalar=%d\n", k, B[k], S[k] }' \
        | sort -k1.4 -n)
    [ -z "$tally" ] && { echo "  FAIL $v: no dispatch tally"; fail=1; continue; }
    echo "  $v dispatch:"; echo "$tally" | sed 's/^/    /'
    if echo "$tally" | grep -qv 'scalar=0$'; then
        echo "    FAIL $v fell back to the scalar loop on some blocks"
        fail=1
    fi
done

for v in "${variants[@]}"; do
    [ "$v" = scalar ] && continue
    echo "  $v vs scalar:"
    python3 "$REPO/BSSN_GR/scripts/compare_constraint_output.py" \
        "$WORK/run_scalar" "$WORK/run_$v" || fail=1
done

if [ "$fail" -eq 0 ]; then
    echo "PASS  ($WORK)"
else
    echo "FAIL  ($WORK)"
fi
exit "$fail"
