#!/usr/bin/env bash
# End-to-end BSSN solver profiling.
# Builds each RHS variant FRESH with identical config (only the kernel CMake flag
# differs), runs the SAME short evolution under MPI (real ghost exchange), and
# tabulates the per-step RK-loop phase breakdown: ETS / rkStep / rhs / deriv /
# unzip_sync(=ghost exchange) / zip, with ETS speed-up vs production.
#
#   ./run_profile.sh
#   NP=8 VARIANTS="production cascade_ir_avx512_fused" ./run_profile.sh
#
# Knobs (env): NP (ranks), VARIANTS, CPU_ARCH, JOBS (build -j), PARFILE,
#              ENV_SETUP (sourced for MKL; e.g. oneAPI setvars, or `module load` wrapper).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"                       # dendrogr_upstream

NP=${NP:-4}
PARFILE=${PARFILE:-"$HERE/q1.profile.par.toml"}
CPU_ARCH=${CPU_ARCH:-native}
JOBS=${JOBS:-$(nproc)}
ENV_SETUP=${ENV_SETUP:-/opt/intel/oneapi/setvars.sh} # MKL (+icpx). On HPC: a script that `module load`s.
# Local dendrolib with dendro_bh::BHHistory (bh_history.h); empty => CMake FetchContents master.
DENDROLIB_DIR=${DENDROLIB_DIR:-$HOME/research/dendrolib_dfvk_copy}
BUILD_ROOT=${BUILD_ROOT:-"$HERE"}   # set to a shared scratch path to keep build dirs off home
OUT=${OUT:-"$HERE/profile_results.csv"}
VARIANTS=${VARIANTS:-"production cascade_ir_avx512 cascade_ir_avx512_fused cascade_ir_avx cascade_ir_avx_fused"}

# variant label -> CMake flag (production = none). The full kernel catalog:
declare -A FLAG=(
  [production]=""
  [naive]="-DBSSN_USE_NAIVE=ON"
  [cascade]="-DBSSN_USE_CASCADE=ON"
  [cascade_avx]="-DBSSN_USE_CASCADE_AVX=ON"
  [cascade_avx_fused]="-DBSSN_USE_CASCADE_AVX_FUSED=ON"
  [cascade_avx512]="-DBSSN_USE_CASCADE_AVX512=ON"
  [cascade_avx512_fused]="-DBSSN_USE_CASCADE_AVX512_FUSED=ON"
  [cascade_ir_avx]="-DBSSN_USE_CASCADE_IR_AVX=ON"
  [cascade_ir_avx_fused]="-DBSSN_USE_CASCADE_IR_AVX_FUSED=ON"
  [cascade_ir_avx512]="-DBSSN_USE_CASCADE_IR_AVX512=ON"
  [cascade_ir_avx512_fused]="-DBSSN_USE_CASCADE_IR_AVX512_FUSED=ON"
)

[[ -f "$ENV_SETUP" ]] && { echo "## env: sourcing $ENV_SETUP"; set +u; source "$ENV_SETUP" >/dev/null 2>&1 || true; set -u; }
command -v mpirun >/dev/null || { echo "ERROR: mpirun not found"; exit 1; }
command -v cmake  >/dev/null || { echo "ERROR: cmake not found";  exit 1; }
[[ -f "$PARFILE" ]] || { echo "ERROR: parfile $PARFILE missing"; exit 1; }

# pull a "   <label>   : <value>" timer line from a run log
tmr() { grep -E "^[[:space:]]+$2[[:space:]]+:" "$1" | tail -1 | awk -F: '{gsub(/[[:space:]]/,"",$2); print $2}'; }

echo "variant,ets_s,rkstep_s,rhs_s,deriv_s,unzip_ghost_s,zip_s,ets_speedup" > "$OUT"
base_ets=""
for v in $VARIANTS; do
  flag="${FLAG[$v]-__UNSET__}"
  [[ "$flag" == "__UNSET__" ]] && { echo "## unknown variant '$v' — skipping"; continue; }
  mkdir -p "$BUILD_ROOT"; bdir="$BUILD_ROOT/build_$v"; bin="$bdir/BSSN_GR/bssnSolver"
  echo ""; echo "===== $v  (${flag:-production CSE}) ====="
  if [[ ! -x "$bin" ]]; then
    echo "## building (fresh, identical config + $flag) ..."
    cmake -S "$REPO" -B "$bdir" -DCMAKE_BUILD_TYPE=Release -DCPU_ARCH="$CPU_ARCH" \
          ${DENDROLIB_DIR:+-DDENDRO_dendrolib_DIR="$DENDROLIB_DIR"} $flag >"$bdir.cfg.log" 2>&1 \
      && cmake --build "$bdir" -j"$JOBS" >"$bdir.build.log" 2>&1 \
      || { echo "## BUILD FAILED ($v) — see $bdir.build.log"; tail -4 "$bdir.build.log"; continue; }
  fi
  log="$HERE/run_$v.log"
  echo "## running: mpirun -np $NP bssnSolver $(basename "$PARFILE")"
  ( cd "$HERE" && timeout 1800 mpirun --oversubscribe -np "$NP" "$bin" "$PARFILE" >"$log" 2>&1 )
  grep -q 'phase breakdown' "$log" || { echo "## RUN FAILED ($v) — tail:"; tail -6 "$log"; continue; }
  ets=$(grep 'ETS time' "$log" | tail -1 | awk -F: '{gsub(/[[:space:]]/,"",$2);print $2}')
  rk=$(tmr "$log" rkStep); rhs=$(tmr "$log" rhs); dv=$(tmr "$log" deriv)
  uz=$(tmr "$log" unzip_sync); zp=$(tmr "$log" zip)
  [[ -z "$base_ets" ]] && base_ets="$ets"
  sp=$(awk -v b="$base_ets" -v e="$ets" 'BEGIN{printf (e+0>0)?"%.2f":"-", b/e}')
  printf "%s,%s,%s,%s,%s,%s,%s,%s\n" "$v" "$ets" "$rk" "$rhs" "$dv" "$uz" "$zp" "$sp" >> "$OUT"
  printf "  ETS=%ss  rkStep=%ss  rhs=%ss  deriv=%ss  ghost(unzip)=%ss  -> %sx vs production\n" "$ets" "$rk" "$rhs" "$dv" "$uz" "$sp"
done
echo ""; echo "## results -> $OUT"; column -t -s, "$OUT" 2>/dev/null || cat "$OUT"
