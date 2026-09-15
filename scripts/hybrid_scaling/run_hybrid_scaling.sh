#!/usr/bin/env bash
# Pure-MPI vs hybrid MPI+OpenMP scaling for DendroGR BSSN. See README.md.
# Builds bssnScalingBench once, then sweeps OpenMP threads/rank at fixed total
# cores (T=1 is pure-MPI); parses the per-step JSONL into a CSV. Run multi-node.

#SBATCH --job-name=bssn_hybrid_scaling
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=hybrid_scaling_%j.out
# account/partition/constraint differ per allocation -> pass on the sbatch line.

set -uo pipefail

# --- config (all overridable via env) ----------------------------------------
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO="$HERE"; while [[ "$REPO" != / && ! -f "$REPO/CMakeLists.txt" ]]; do REPO="$(dirname "$REPO")"; done

# Site preset: modules + arch + cores/node + dendrolib path + MPI fabric. Default
# stampede3 (TACC); 'chpc' for Utah csl nodes; 'none' if you load modules yourself.
# Any individual value below is still overridable by its own env var. See PORTING.md.
SITE="${SITE:-stampede3}"
case "$SITE" in
  stampede3)                                       # Sapphire Rapids, Omni-Path fabric
    # sapphirerapids needs GCC >= 11 / Clang >= 12 / ICX; CMake now errors with the
    # remedy if the loaded compiler is older. Do NOT "fix" that by switching to
    # native: Stampede3 has SPR/ICX/SKX partitions and builds usually happen on a
    # login node, so native bakes in the build host's ISA under the benchmark.
    : "${CPU_ARCH:=sapphirerapids}"; SITE_CORES=112
    : "${DENDROLIB_DIR:=$HOME/Projects/dendro/DVK_experimental/Dendro-5.01}"
    export OMPI_MCA_mtl="${OMPI_MCA_mtl:-ofi}" FI_PROVIDER="${FI_PROVIDER:-opx}" ;;  # select Omni-Path, not TCP
  chpc)                                            # Utah CHPC, csl (Cascade Lake)
    : "${CPU_ARCH:=cascadelake}"; SITE_CORES=40
    : "${DENDROLIB_DIR:=$HOME/research/dendrolib_dfvk_copy}" ;;
  none|skip) SITE_CORES="" ;;                      # no presets; load modules yourself
  *) echo "ERROR: unknown SITE=$SITE (stampede3|chpc|none)"; exit 1 ;;
esac
CPU_ARCH="${CPU_ARCH:-native}"
DENDROLIB_DIR="${DENDROLIB_DIR:-}"                 # from site preset or env; must have the hybrid path (checked below)
CORES_PER_NODE="${CORES_PER_NODE:-${SLURM_CPUS_ON_NODE:-${SITE_CORES:-$(nproc)}}}"

GRID="${GRID:-bbh}"; LEV="${LEV:-7}"     # uniform blocks ~= 8^(LEV-3); use it to saturate many nodes
# Block fusion (dendrolib -DOCT2BLK_COARSEST_LEV): 0 = fuse same-level octants under a common
# ancestor into big blocks (production RHS default); 31 = no fusion, block-per-element. Uniform
# needs many small blocks to fill ranks, so default it to 31; bbh keeps the fused default.
OCT2BLK_LEV="${OCT2BLK_LEV:-$([[ "$GRID" == uniform ]] && echo 31 || echo 0)}"

BUILD_ROOT="${BUILD_ROOT:-$HERE}"                  # put on shared scratch for multi-node
BUILD_DIR="${BUILD_DIR:-$BUILD_ROOT/build_hybrid_ob${OCT2BLK_LEV}}"
CASCADE_FLAG="${CASCADE_FLAG:-}"                   # optional, e.g. -DBSSN_USE_CASCADE_AVX512_FUSED=ON
JOBS="${JOBS:-$(nproc)}"
NNODES="${NNODES:-${SLURM_NNODES:-1}}"
TOTAL_CORES=$(( CORES_PER_NODE * NNODES ))
MPI_LAUNCH="${MPI_LAUNCH:-mpirun}"       # or srun (set SRUN_MPI plugin)
SRUN_MPI="${SRUN_MPI:-pmix}"

THREADS_LIST="${THREADS_LIST:-1 2 4}"    # threads/rank; 1 = hybrid-T1 (NOT pure MPI). Per-socket divisors.

# PROFILE=roofline: instead of the sweep, measure one single-node config under likwid
# (MEM_DP counters = achieved memory bandwidth + DP GFLOP/s + operational intensity) and a
# likwid-bench STREAM ceiling -> where the workload sits vs the memory-bandwidth wall.
PROFILE="${PROFILE:-}"
ROOFLINE_T="${ROOFLINE_T:-1}"                       # threads/rank for the roofline run (sweep to probe splits)
ROOFLINE_STREAM_KERNEL="${ROOFLINE_STREAM_KERNEL:-stream}"  # likwid-bench kernel; stream_avx512 for widest SIMD

# Default = the self-contained BBH grid here. NOT BSSN_GR/pars/q1.par.toml: it needs
# an external TwoPunctures (tp_*) file and segfaults without one.
PARFILE="${PARFILE:-$HERE/q1.scaling.par.toml}"
if [[ "$PARFILE" != /* ]]; then
  [[ -f "$HERE/$PARFILE" ]] && PARFILE="$HERE/$PARFILE" || PARFILE="$REPO/BSSN_GR/pars/$PARFILE"
fi
STEPS="${STEPS:-10}"; WARMUP="${WARMUP:-2}"

STAMP="${SLURM_JOB_ID:-local}"
OUTDIR="${OUTDIR:-$HERE/results_${STAMP}}"
OUT_CSV="${OUT_CSV:-$OUTDIR/hybrid_scaling_${STAMP}.csv}"

# --- preflight ---------------------------------------------------------------
mkdir -p "$OUTDIR"
if [[ "$SITE" != none && "$SITE" != skip ]] && command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  case "$SITE" in
    stampede3)   # openmpi is off the default module tree -> needs 'module use' first
      module load autotools/1.4 cmake/4.1.1 python/3.12.11 gcc/15.1.0 mkl/25.1 gsl/2.8 \
        && module use /scratch/projects/compilers/modulefiles \
        && module load openmpi/5.0.9 ;;
    chpc)
      module load gcc/15.1.0 intel-oneapi-mkl/2025.3.1 openmpi/5.0.3 ;;
  esac || { echo "ERROR: module load failed for SITE=$SITE (try SITE=none if you load modules yourself)"; exit 1; }
  # Use module-provided paths; TACC_GSL_DIR is set by the gsl module to the correct
  # compiler-matched build. GSL_ROOT_DIR may already be in the env pointing at a
  # different build, so override it explicitly here.
  export GSL_ROOT_DIR="${TACC_GSL_DIR:-/home1/apps/gcc15/gsl/2.8}"  # override any pre-set intel path
fi
command -v cmake   >/dev/null || { echo "ERROR: cmake not found";   exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
[[ -f "$PARFILE" ]]     || { echo "ERROR: parfile not found: $PARFILE"; exit 1; }
[[ -d "$DENDROLIB_DIR" ]] || { echo "ERROR: DENDROLIB_DIR not set/found: '$DENDROLIB_DIR' (set DENDROLIB_DIR, or SITE with a preset)"; exit 1; }

echo "site=$SITE  nodes=$NNODES cores/node=$CORES_PER_NODE total=$TOTAL_CORES  launch=$MPI_LAUNCH arch=$CPU_ARCH"
echo "sweep T=[$THREADS_LIST]  parfile=$(basename "$PARFILE") grid=$GRID oct2blk=$OCT2BLK_LEV steps=$STEPS  -> $OUTDIR"

# --- memory preflight --------------------------------------------------------
# Uniform + no-fusion keeps 8^LEV blocks, each with its own ghost zone, so the
# unzipped arrays OOM long before the compute is interesting. Estimate per-node
# RAM up front and fail with a node count instead of a cryptic bad_alloc/SIGSEGV
# (a bad_alloc thrown in an OMP region crashes as a corrupt-pointer segfault).
if [[ "$GRID" == uniform && "${MEM_CHECK:-on}" != off ]]; then
  ord=$(grep -oiE 'BSSN_ELE_ORDER[[:space:]]*=[[:space:]]*[0-9]+' "$PARFILE" | grep -oE '[0-9]+$'); ord="${ord:-6}"
  pad=$(( (2*ord + 1) ** 3 ))                              # padded pts per unzipped block
  fuse=$([[ "$OCT2BLK_LEV" -ge 20 ]] && echo 1 || echo 512)  # ob31 ~ block/elem; ob0 fuses ~8^3
  nblk=$(( 8 ** LEV / fuse )); (( nblk < 1 )) && nblk=1
  gb=$(( nblk * pad * 85 * 8 / 1073741824 ))              # ~85 arrays/block (24 fields x RK stages + unzip/deriv peak)
  per_node=$(( gb / NNODES ))
  ram=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0) / 1048576 ))
  (( ram < 1 )) && ram="${NODE_RAM_GB:-192}"
  budget=$(( ram * 85 / 100 ))
  echo "mem: ~${nblk} blocks -> ~${per_node} GB/node (usable ~${budget} of ${ram} GB)"
  if (( per_node > budget )); then
    echo "ERROR: est ~${per_node} GB/node > ~${budget} GB usable. Need >= $(( gb / budget + 1 )) nodes, a smaller LEV, or OCT2BLK_LEV=0 (fused)."
    echo "       Override with MEM_CHECK=off if you know the run fits."
    exit 1
  fi
fi

# --- builds -------------------------------------------------------------------
# TWO binaries, not one, because "run the hybrid binary at OMP_NUM_THREADS=1" is NOT
# pure MPI. -DDENDRO_HYBRID_OMP=ON forces DENDRO_UNZIP_OMP_EFFECTIVE=ON
# (dendrolib/CMakeLists.txt:295), and that path allocates + zero-fills a per-rank
# all_dg buffer (mesh.tcc:11295) the serial path never touches. Measured cost at T=1
# with zero threading benefit: ~13% on unzip, ~6% on rk_step (8 cores, depth 9,
# 2026-07-16). So T=1 is "hybrid-T1" and the CSV labels it that way.
#
# The baseline MUST also set -DDENDRO_UNZIP_SPEEDUP=ON. dendrolib/CMakeLists.txt:302
# says:
#     if(DENDRO_UNZIP_SPEEDUP OR DENDRO_UNZIP_OMP_EFFECTIVE)
#       set(DENDRO_UNZIP_SCATTER_FAST_EFFECTIVE ON) ... TENSOR_SIMD ... SPEEDUP ...
# i.e. DENDRO_HYBRID_OMP=ON silently ALSO enables the integer-index scatter fast path
# and the SIMD tensor kernels, which are MPI-safe and have nothing to do with
# threading. A bare -DDENDRO_HYBRID_OMP=OFF baseline therefore gives up those too
# (measured: unzip 193ms vs 117ms, ~39% slower) and would FLATTER hybrid by crediting
# it with unrelated optimizations. UNZIP_SPEEDUP=ON isolates the threading variable.
# BASELINE_OFF=off skips the baseline if you only want the T-trend.
BASELINE_OFF="${BASELINE_OFF:-on}"
BUILD_DIR_OFF="${BUILD_DIR}_puremip"

# Pinned on EVERY variant so the only difference between binaries is the flag under
# test. These are not defaults you can rely on: a build tree configured in an earlier
# session keeps its own cached values, and DENDRO_USE_NEW_DERIVS in particular swaps
# the whole RHS translation unit (BSSN_GR/CMakeLists.txt:92 picks rhs_experimental.cpp
# over rhs.cpp), while DVEC_ZERO_ALLOC switches calloc/malloc for the unzip padding.
# On 2026-07-16 an unpinned pair differed in both, which produced a bogus "not
# bit-exact" verdict and a bogus all_dg timing -- the binaries were running different
# physics. Pin them; do not remove.
PIN_FLAGS=(-DDENDRO_USE_NEW_DERIVS=OFF -DDVEC_ZERO_ALLOC=ON)

build_variant() {  # $1=builddir  $2=tag  $3..=extra cmake args
  local bdir="$1" tag="$2"; shift 2
  local bin="$bdir/BSSN_GR/bssnScalingBench"
  local want="$*"
  if [[ -x "$bin" ]]; then
    # NEVER reuse blind. A stale tree can hold different flags than we asked for, and
    # the binary looks identical from the outside. Record what each tree was built
    # with and refuse to reuse a mismatch.
    if [[ -f "$bdir/.sweep_flags" ]] && [[ "$(cat "$bdir/.sweep_flags")" == "$want ${PIN_FLAGS[*]}" ]]; then
      echo "## reusing build $bin"; return 0
    fi
    echo "## STALE BUILD: $bdir was configured with"
    echo "##     $(cat "$bdir/.sweep_flags" 2>/dev/null || echo '<unknown - no .sweep_flags>')"
    echo "##   but this sweep needs"
    echo "##     $want ${PIN_FLAGS[*]}"
    echo "##   Reusing it would compare binaries that differ in more than the flag under test."
    echo "##   Remove it and re-run:  rm -rf $bdir"
    exit 1
  fi
  echo "## building $tag ($want) ..."
  cmake -S "$REPO" -B "$bdir" -DCMAKE_BUILD_TYPE=Release -DCPU_ARCH="$CPU_ARCH" \
        -DDENDRO_dendrolib_DIR="$DENDROLIB_DIR" \
        -DOCT2BLK_COARSEST_LEV="$OCT2BLK_LEV" \
        -DENABLE_DENDRO_PROFILE_COUNTERS=ON $CASCADE_FLAG "${PIN_FLAGS[@]}" "$@" \
        >"$OUTDIR/configure_${tag}.log" 2>&1 \
    || { echo "## CONFIGURE FAILED ($tag) -- $OUTDIR/configure_${tag}.log"; tail -20 "$OUTDIR/configure_${tag}.log"; exit 1; }
  grep -E "unzip/zip speedups|hybrid OpenMP" "$OUTDIR/configure_${tag}.log" | sed 's/^/   /'
  cmake --build "$bdir" --target bssnScalingBench -j"$JOBS" >"$OUTDIR/build_${tag}.log" 2>&1 \
    || { echo "## BUILD FAILED ($tag) -- $OUTDIR/build_${tag}.log"; tail -20 "$OUTDIR/build_${tag}.log"; exit 1; }
  echo "$want ${PIN_FLAGS[*]}" > "$bdir/.sweep_flags"
}

# UNZIP_BATCH: routes Ctx::unzip through Mesh::unzip_scatter_batch, which reuses a
# dof=1-sized DG scratch instead of materializing and zero-filling
# (numTotalElements * dof * nPe) doubles per call. Verified bit-exact against
# BATCH=OFF (R=2, T=4; gate confirmed non-vacuous against a known-different build),
# and measured 2.0x faster unzip / -18% rk_step at R=8 T=1 on 2026-07-16. ON by
# default here; UNZIP_BATCH=off to measure the old path.
UNZIP_BATCH="${UNZIP_BATCH:-on}"
BATCH_FLAG="-DDENDRO_UNZIP_BATCH=$([[ "$UNZIP_BATCH" == on ]] && echo ON || echo OFF)"

BIN="$BUILD_DIR/BSSN_GR/bssnScalingBench"
BIN_OFF="$BUILD_DIR_OFF/BSSN_GR/bssnScalingBench"
build_variant "$BUILD_DIR" hybrid -DDENDRO_HYBRID_OMP=ON "$BATCH_FLAG"
[[ "$BASELINE_OFF" == on ]] && \
  build_variant "$BUILD_DIR_OFF" puremip -DDENDRO_HYBRID_OMP=OFF -DDENDRO_UNZIP_SPEEDUP=ON

# --- launch one split: R ranks x T threads. Pinning is what makes hybrid work. -
run_split() {
  local R="$1" T="$2" RPN="$3" PFX="$4" EXE="${5:-$BIN}"
  export OMP_NUM_THREADS="$T" OMP_PROC_BIND=close OMP_PLACES=cores
  local args=("$PARFILE" --grid "$GRID" --lev "$LEV" --steps "$STEPS" --warmup "$WARMUP" --prefix "$PFX")
  if [[ "$MPI_LAUNCH" == srun ]]; then
    srun --mpi="$SRUN_MPI" --nodes="$NNODES" --ntasks="$R" --ntasks-per-node="$RPN" \
         --cpus-per-task="$T" --cpu-bind=cores --distribution=block:block "$EXE" "${args[@]}"
  else
    mpirun --np "$R" --map-by "ppr:${RPN}:node:pe=${T}" --bind-to core "$EXE" "${args[@]}"
  fi
}

# --- roofline: one single-node config under likwid MEM_DP vs a STREAM ceiling -
run_roofline() {
  command -v likwid-perfctr >/dev/null 2>&1 && command -v likwid-bench >/dev/null 2>&1 \
    || { echo "ERROR: PROFILE=roofline needs likwid (likwid-perfctr + likwid-bench). Try 'module load likwid'."; return 1; }
  local T="$ROOFLINE_T"; (( CORES_PER_NODE % T == 0 )) || { echo "ERROR: ROOFLINE_T=$T must divide $CORES_PER_NODE"; return 1; }
  local R=$(( CORES_PER_NODE / T )) hi=$(( CORES_PER_NODE - 1 ))
  local PFX="$OUTDIR/roofline_t${T}"
  export OMP_NUM_THREADS="$T" OMP_PROC_BIND=close OMP_PLACES=cores
  echo "roofline: 1 node, $R ranks x $T threads, MEM_DP over cores 0-$hi"
  # IMC counters see all node memory traffic; keep all ranks on one node (ppr:R:node)
  likwid-perfctr -g MEM_DP -C "0-$hi" -f \
    mpirun --np "$R" --map-by "ppr:${R}:node:pe=${T}" --bind-to core \
      "$BIN" "$PARFILE" --grid "$GRID" --lev "$LEV" --steps "$STEPS" --warmup "$WARMUP" --prefix "$PFX" \
    >"${PFX}.likwid" 2>&1 \
    || { echo "## likwid-perfctr failed (counter access? need perf_event_paranoid<=0 or the likwid accessdaemon)"; tail -15 "${PFX}.likwid"; return 1; }
  local nsock; nsock=$(likwid-topology 2>/dev/null | awk -F: '/Sockets:/{gsub(/ /,"",$2);print $2}'); nsock="${nsock:-1}"
  likwid-bench -t "$ROOFLINE_STREAM_KERNEL" -W "S0:2GB:$(( CORES_PER_NODE / nsock ))" \
    >"$OUTDIR/roofline_stream.likwid" 2>&1 || echo "## likwid-bench $ROOFLINE_STREAM_KERNEL failed; ceiling unknown"
  python3 - "${PFX}.likwid" "$OUTDIR/roofline_stream.likwid" "$nsock" "$T" "$R" "$OUTDIR/roofline_${STAMP}.csv" <<'PYEOF'
import re, sys
run, stream, nsock, T, R, out_csv = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]
def maxfloat(path, needle):                     # STAT rows put the node-wide Sum as the largest value on the line
    best = None
    try: lines = open(path).read().splitlines()
    except OSError: return None
    for ln in lines:
        if needle in ln:
            for m in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", ln.split(needle,1)[1]):
                v = float(m); best = v if best is None or v > best else best
    return best
bw   = maxfloat(run, "Memory bandwidth [MBytes/s]")     # achieved, node-wide
flop = maxfloat(run, "DP [MFLOP/s]")
oi   = maxfloat(run, "Operational intensity")
# STREAM ceiling: likwid-bench reports one socket's MByte/s; scale to the node
sbw = None
for ln in open(stream).read().splitlines() if stream else []:
    if "MByte/s" in ln:
        m = re.findall(r"[\d.]+", ln.split("MByte/s",1)[1]);  sbw = float(m[0]) if m else sbw
peak = sbw * nsock if sbw else None
bw_g   = bw/1000 if bw else None
peak_g = peak/1000 if peak else None
pct    = 100*bw/peak if (bw and peak) else None
flop_g = flop/1000 if flop else None
f = lambda x: f"{x:.2f}" if isinstance(x,float) else ("" if x is None else str(x))
print("\n== roofline ==")
print(f"  achieved memory BW : {f(bw_g)} GB/s (node)")
print(f"  STREAM ceiling     : {f(peak_g)} GB/s (node, {nsock}x socket)")
print(f"  % of BW wall       : {f(pct)} %")
print(f"  DP compute         : {f(flop_g)} GFLOP/s")
print(f"  operational intens.: {f(oi)} FLOP/byte")
with open(out_csv, "w") as fh:
    fh.write("threads,ranks,achieved_bw_gbs,stream_bw_gbs,pct_of_wall,dp_gflops,op_intensity\n")
    fh.write(",".join(f(x) for x in [T, R, bw_g, peak_g, pct, flop_g, oi]) + "\n")
print(f"  wrote {out_csv}")
PYEOF
}

if [[ "$PROFILE" == roofline ]]; then run_roofline; exit $?; fi

# --- sweep -------------------------------------------------------------------
# The real pure-MPI reference: DENDRO_HYBRID_OMP=OFF + UNZIP_SPEEDUP=ON, one rank per
# core. This is the only row that is honestly "pure MPI"; every T in THREADS_LIST
# (including T=1) runs the hybrid binary and its OMP unzip path.
if [[ "$BASELINE_OFF" == on ]]; then
  RPN=$CORES_PER_NODE; R=$(( RPN * NNODES ))
  PFX="$OUTDIR/run_N${NNODES}_r${R}_t1_off"; rm -f "${PFX}_steps.jsonl"
  echo "===== baseline: pure-MPI (HYBRID_OMP=OFF, UNZIP_SPEEDUP=ON) ranks=$R ranks/node=$RPN ====="
  run_split "$R" 1 "$RPN" "$PFX" "$BIN_OFF" >"${PFX}.log" 2>&1
  [[ -s "${PFX}_steps.jsonl" ]] && echo "## ok" || { echo "## FAILED (no JSONL):"; tail -15 "${PFX}.log"; }
fi

for T in $THREADS_LIST; do
  (( CORES_PER_NODE % T == 0 )) || { echo "## skip T=$T (does not divide $CORES_PER_NODE cores/node)"; continue; }
  RPN=$(( CORES_PER_NODE / T )); R=$(( RPN * NNODES ))
  PFX="$OUTDIR/run_N${NNODES}_r${R}_t${T}"; rm -f "${PFX}_steps.jsonl"
  echo "===== T=$T (hybrid$([[ $T == 1 ]] && echo '-T1, NOT pure MPI')) ranks=$R ranks/node=$RPN ====="
  run_split "$R" "$T" "$RPN" "$PFX" >"${PFX}.log" 2>&1
  [[ -s "${PFX}_steps.jsonl" ]] && echo "## ok" || { echo "## FAILED (no JSONL):"; tail -15 "${PFX}.log"; }
done

# --- parse -> CSV ------------------------------------------------------------
echo ""; echo "## parsing -> $OUT_CSV"
python3 - "$OUTDIR" "$OUT_CSV" <<'PYEOF'
import glob, json, os, re, sys
outdir, out_csv = sys.argv[1], sys.argv[2]
PHASES = ["rk_step", "rhs_wall", "unzip", "unzip_wcomm", "zip",
          "ghost_pack", "ghost_wait", "ghost_unpack"]
mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")

rows = []
for path in sorted(glob.glob(os.path.join(outdir, "run_N*_steps.jsonl"))):
    m = re.search(r"run_N(\d+)_r(\d+)_t(\d+)(_off)?_steps\.jsonl$", os.path.basename(path))
    if not m: continue
    nodes, ranks, threads = map(int, m.groups()[:3])
    is_off = m.group(4) is not None   # DENDRO_HYBRID_OMP=OFF build => real pure MPI
    # per phase, mean over timed steps of the per-step max-across-ranks (the critical path)
    acc = {p: [] for p in PHASES}; nblocks = None
    # Load imbalance = rhs_wall max/mean across ranks (1.0 = perfect). Measured on rhs_wall,
    # NOT rk_step: the ghost exchange synchronizes ranks, so a rank that finishes early just
    # blocks longer in comm and rk_step max/mean is ~1.000 regardless of how skewed the work is.
    # rhs_wall is the comm-free compute phase, so it exposes the skew the dynamic schedule
    # exists to fix.
    # rhs_wall, NOT the old "rhs" key: that was a thread-0 sample (~blocks_per_rank/T), so its
    # max/mean mixed rank skew with thread-0's luck under schedule(dynamic,1).
    # Active ranks only (the emitter skips inactive ones) => understated if the mesh spans fewer
    # ranks than exist.
    imb = []
    for line in open(path):
        line = line.strip()
        if not line: continue
        try: rec = json.loads(line)
        except json.JSONDecodeError: continue
        ph = rec.get("phase", {})
        for p in PHASES:
            if p in ph and "max" in ph[p]: acc[p].append(float(ph[p]["max"]))
        r = ph.get("rhs_wall")
        if r and r.get("mean", 0) > 0: imb.append(float(r["max"]) / float(r["mean"]))
        if nblocks is None: nblocks = rec.get("num_blocks")
    v = {p: mean(acc[p]) for p in PHASES}
    # mode: only the flag-OFF build is "pure-MPI". The hybrid binary at T=1 is "hybrid-T1":
    # -DDENDRO_HYBRID_OMP=ON forces DENDRO_UNZIP_OMP_EFFECTIVE=ON (dendrolib/CMakeLists.txt:296),
    # so even at T=1 it runs the OMP unzip path and pays its all_dg cost (0.84x vs flag-OFF at
    # t=1, findings/unzip_openmp.txt).
    mode = "pure-MPI" if is_off else ("hybrid-T1" if threads == 1 else "hybrid")
    rows.append(dict(mode=mode, nodes=nodes, ranks=ranks,
                     threads=threads, ranks_per_node=ranks // nodes if nodes else 0,
                     ghost_pack_wait=v["unzip_wcomm"] - v["unzip"], num_blocks=nblocks,
                     blocks_per_rank=(nblocks / ranks) if (nblocks and ranks) else float("nan"),
                     imbalance=mean(imb), **v))

# Two baselines, deliberately. speedup_vs_hybrid_T1 is the T-trend (does threading help?).
# speedup_vs_pure_mpi is the one that answers "should we ship this?" -- and the gap between
# the two IS the all_dg cost the hybrid binary pays even at T=1. If the pure-MPI row is
# missing (BASELINE_OFF=off), that column is blank rather than silently falling back to
# hybrid-T1 and calling it pure MPI.
base_hyb1 = {r["nodes"]: r["rk_step"] for r in rows if r["threads"] == 1 and r["mode"] != "pure-MPI"}
base_pure = {r["nodes"]: r["rk_step"] for r in rows if r["mode"] == "pure-MPI"}
for r in rows:
    b, p = base_hyb1.get(r["nodes"]), base_pure.get(r["nodes"])
    r["speedup"]      = (b / r["rk_step"]) if (b and r["rk_step"] > 0) else float("nan")
    r["speedup_pure"] = (p / r["rk_step"]) if (p and r["rk_step"] > 0) else float("nan")
rows.sort(key=lambda r: (r["nodes"], r["mode"] != "pure-MPI", r["threads"]))

cols = ["mode", "nodes", "ranks", "threads", "ranks_per_node", "rk_step", "rhs_wall",
        "unzip", "unzip_wcomm", "ghost_pack_wait", "ghost_pack", "ghost_wait", "ghost_unpack",
        "zip", "num_blocks", "blocks_per_rank", "imbalance", "speedup", "speedup_pure"]
# ghost_pack_wait = unzip_wcomm - unzip, the whole exchange lumped together. The three
# ghost_* columns break it down and sum back to it within ~1% (on means; do NOT check the
# identity on the max columns -- different ranks are the max for different sub-timers, so
# summing maxima overshoots by ~40%).
#
# Measured 2026-07-16 (8 cores, depth 9, B/rank 36, properly bound): ghost_wait is ~62% of
# the exchange at T=1/T=2. The exchange is dominated by Waitall, and Waitall is mostly ranks
# waiting on slower NEIGHBOURS -- load imbalance billed to comm, not wire time. So a rising
# ghost_pack_wait_s is usually an imbalance signal, not a network signal, and arguments about
# ghost VOLUME are arguing about a minority of the cost.
hdr = {"rk_step": "rk_step_s", "rhs_wall": "rhs_wall_s", "unzip": "unzip_s",
       "unzip_wcomm": "unzip_wcomm_s", "ghost_pack_wait": "ghost_pack_wait_s",
       "ghost_pack": "ghost_pack_s", "ghost_wait": "ghost_wait_s",
       "ghost_unpack": "ghost_unpack_s", "zip": "zip_s",
       "imbalance": "rhs_wall_imbalance_max_over_mean", "speedup": "speedup_vs_hybrid_T1",
       "speedup_pure": "speedup_vs_pure_mpi"}
fmt = lambda x: (f"{x:.4f}" if x == x else "") if isinstance(x, float) else ("" if x is None else str(x))
with open(out_csv, "w") as fh:
    fh.write(",".join(hdr.get(c, c) for c in cols) + "\n")
    for r in rows: fh.write(",".join(fmt(r[c]) for c in cols) + "\n")
print(f"wrote {len(rows)} rows -> {out_csv}")

# saturation check: too few blocks/rank => imbalance dominates, results are noise
bprs = [r["blocks_per_rank"] for r in rows if r["blocks_per_rank"] == r["blocks_per_rank"]]
if bprs and min(bprs) < 4:
    print(f"\n  !! UNDER-SATURATED: as few as {min(bprs):.1f} blocks/rank -- results unreliable."
          f"\n     Use a bigger mesh: GRID=uniform LEV=7 (~4096 blocks) or LEV=8, or fewer nodes.")
elif bprs and min(bprs) < 16:
    print(f"\n  ! marginal: as few as {min(bprs):.1f} blocks/rank (prefer >=16). Try a larger GRID=uniform LEV.")
PYEOF

echo ""; echo "== $OUT_CSV =="
column -t -s, "$OUT_CSV" 2>/dev/null || cat "$OUT_CSV"
echo "(headline rk_step_s. speedup_vs_pure_mpi is the ship/no-ship number, measured against the"
echo " DENDRO_HYBRID_OMP=OFF build; speedup_vs_hybrid_T1 is only the T-trend within the hybrid"
echo " binary. The gap between the two at T=1 is what the hybrid unzip path costs before any"
echo " threading happens (all_dg, mesh.tcc:11295)."
echo " ghost_pack_wait_s = pack+Isend/Irecv+Waitall+unpack, NOT wire time -- fewer ranks each pack"
echo " MORE bytes (per-rank ghost ~ (V/R)^(2/3) RISES with T even as job-aggregate volume falls),"
echo " so this rising is expected and is not evidence against hybrid."
echo " phase times carry profiler-barrier overhead -> read them relative, not absolute.)"
