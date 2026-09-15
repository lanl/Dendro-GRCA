#!/usr/bin/env bash
#SBATCH --job-name=bssn_cascade_e2e
#SBATCH --account=sundar-np
#SBATCH --partition=soc-np
#SBATCH --constraint=csl              # Intel Cascadelake (AVX-512). Use 'rom' for AMD Zen2 (AVX2-only).
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40          # full csl node (Xeon Gold 6230 = 2x20 cores). Drop if the octree over-decomposes.
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=slurm_e2e_%j.out
#
# Multi-node end-to-end DendroGR BSSN profiling: production CSE vs the deployed
# cascade kernels, real MPI ghost exchange across nodes. Builds each variant
# FRESH with identical config (only the kernel CMake flag differs), runs the
# SAME short BBH evolution under srun, and tabulates the per-step RK-loop phase
# breakdown (ETS / rkStep / rhs / deriv / unzip_sync=ghost / zip) + ETS speedup.
#
# This is the cluster companion to ./run_profile.sh (single-node, validated on
# the laptop). Same engine, SLURM-native launch + module setup.
#
#   sbatch run_profile.slurm
#   # or interactively on an allocated node:
#   salloc --account=sundar-np --partition=soc-np --constraint=csl --nodes=4 --exclusive --time=2:00:00
#   ./run_profile.slurm
#
# Knobs (env, all optional):
#   VARIANTS   space list of kernel labels (default: production + the two deployed fused kernels)
#   PARFILE    parfile to evolve (default q1.profile.4node.par.toml — MAXDEPTH=11)
#   CPU_ARCH   -march target (default native; on csl = cascadelake)
#   NTASKS     total MPI ranks (default: $SLURM_NTASKS from the allocation)
#   MPI_LAUNCH 'srun' (default, SLURM-native) or 'mpirun'
#   JOBS       build parallelism (default: cores on the build node)
set -uo pipefail
# Under sbatch, BASH_SOURCE points at the spooled script copy; SLURM_SUBMIT_DIR is
# the dir you launched from (profile_e2e). Fall back to BASH_SOURCE when interactive.
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO="$(cd "$HERE/.." && pwd)"                       # dendrogr_upstream

# ---- modules ------------------------------------------------------------
# gcc is the COMPILER; MKL + MPI are libraries. gcc/15.1.0 (NOT 16 — it
# regresses the no-CSE baseline). intel-oneapi-mkl is a *library* module and
# does NOT trigger the compiler-family clash that intel-oneapi-COMPILERS would.
# Adjust these names to your site (`module spider <name>` to find versions).
MODULES=${MODULES:-"gcc/15.1.0 intel-oneapi-mkl/2025.3.1 openmpi/5.0.3"}
if command -v module >/dev/null 2>&1; then
  echo "## module load $MODULES"
  module purge 2>/dev/null || true
  # shellcheck disable=SC2086
  module load $MODULES || { echo "ERROR: module load failed — fix MODULES at top of script"; exit 1; }
fi

NTASKS=${NTASKS:-${SLURM_NTASKS:-4}}
PARFILE=${PARFILE:-"$HERE/q1.profile.4node.par.toml"}
[[ "$PARFILE" = /* ]] || PARFILE="$HERE/$PARFILE"   # resolve a bare filename relative to the submit dir
CPU_ARCH=${CPU_ARCH:-native}
JOBS=${JOBS:-$(nproc)}
MPI_LAUNCH=${MPI_LAUNCH:-mpirun}  # OpenMPI: mpirun is SLURM-aware + needs no srun MPI plugin.
                                  # Use srun only if your SLURM has the matching --mpi plugin (pmix/pmi2).
# Local dendrolib with dendro_bh::BHHistory (bh_history.h); REQUIRED on dfvk —
# vanilla FetchContent pulls an old dendrolib master that lacks it. Point at the
# CHPC checkout of dendrolib_dfvk_copy.
DENDROLIB_DIR=${DENDROLIB_DIR:-$HOME/research/dendrolib_dfvk_copy}
# Where the (large, ~1-2GB each) build dirs go. Set to a SHARED scratch path to
# keep them off home quota, e.g. /scratch/general/vast/$USER/dendrogr_builds.
# Must be shared (not node-local /scratch/local) so multi-node MPI sees the binary.
BUILD_ROOT=${BUILD_ROOT:-"$HERE"}
OUT=${OUT:-"$HERE/profile_results_${SLURM_JOB_ID:-local}.csv"}
# Default decomposes the speedup: production -> non-fused (kernel restructuring)
# -> fused (+ stencil fusion), both widths. For strong scaling, override with just
# the deployment kernel: VARIANTS="production cascade_ir_avx512_fused".
VARIANTS=${VARIANTS:-"production cascade_ir_avx512 cascade_ir_avx512_fused cascade_ir_avx cascade_ir_avx_fused"}

# variant label -> CMake flag (production = none). Full kernel catalog:
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

echo "## $(g++ --version | head -1)"
echo "## ranks=$NTASKS  launch=$MPI_LAUNCH  arch=$CPU_ARCH  parfile=$(basename "$PARFILE")"
command -v cmake >/dev/null || { echo "ERROR: cmake not found (module load cmake?)"; exit 1; }
[[ -f "$PARFILE" ]] || { echo "ERROR: parfile $PARFILE missing"; exit 1; }

# launch wrapper: srun is SLURM-native (no host/oversubscribe flags needed);
# mpirun fallback reads the allocation from the SLURM env.
run_mpi() {  # $1 = binary, $2 = parfile, runs in $HERE
  if [[ "$MPI_LAUNCH" == srun ]]; then
    srun --mpi=pmix -n "$NTASKS" "$1" "$2"
  else
    mpirun -np "$NTASKS" "$1" "$2"
  fi
}

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
      || { echo "## BUILD FAILED ($v) — see $bdir.build.log"; tail -6 "$bdir.build.log"; continue; }
  fi
  log="$HERE/run_${v}_${SLURM_JOB_ID:-local}.log"
  echo "## running: $MPI_LAUNCH -n $NTASKS bssnSolver $(basename "$PARFILE")"
  ( cd "$HERE" && run_mpi "$bin" "$PARFILE" >"$log" 2>&1 )
  grep -q 'phase breakdown' "$log" || { echo "## RUN FAILED ($v) — tail:"; tail -8 "$log"; continue; }
  ets=$(grep 'ETS time' "$log" | tail -1 | awk -F: '{gsub(/[[:space:]]/,"",$2);print $2}')
  rk=$(tmr "$log" rkStep); rhs=$(tmr "$log" rhs); dv=$(tmr "$log" deriv)
  uz=$(tmr "$log" unzip_sync); zp=$(tmr "$log" zip)
  [[ -z "$base_ets" ]] && base_ets="$ets"
  sp=$(awk -v b="$base_ets" -v e="$ets" 'BEGIN{printf (e+0>0)?"%.2f":"-", b/e}')
  printf "%s,%s,%s,%s,%s,%s,%s,%s\n" "$v" "$ets" "$rk" "$rhs" "$dv" "$uz" "$zp" "$sp" >> "$OUT"
  printf "  ETS=%ss  rkStep=%ss  rhs=%ss  deriv=%ss  ghost(unzip)=%ss  -> %sx vs production\n" "$ets" "$rk" "$rhs" "$dv" "$uz" "$sp"
done
echo ""; echo "## results -> $OUT"; column -t -s, "$OUT" 2>/dev/null || cat "$OUT"
