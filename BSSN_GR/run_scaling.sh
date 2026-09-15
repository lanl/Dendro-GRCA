#!/usr/bin/env bash
#
# Strong / weak scaling sweep driver for bssnScalingBench.
#
# Runs the benchmark across a list of (MPI ranks x OpenMP threads) configs and
# collects one <outdir>/<tag>_steps.jsonl per config, ready for the closure /
# scaling-analysis tooling. The benchmark itself does a fixed number of RK steps
# with NO remesh and NO IO, so the JSONL is pure per-kernel timing.
#
# Usage:
#   ./run_scaling.sh <paramFile> [strong|weak] [bbh|uniform] [steps] [warmup]
#
# Override the binary path with BIN=..., and the sweep with CONFIGS="1 1,1 4,...".
# Each CONFIGS entry is "RANKS THREADS".
set -u

BIN=${BIN:-./bssnScalingBench}
PARAM=${1:?usage: run_scaling.sh <paramFile> [strong|weak] [bbh|uniform] [steps] [warmup]}
MODE=${2:-strong}
GRID=${3:-bbh}
STEPS=${4:-20}
WARMUP=${5:-3}

# (ranks threads) configs -- edit for your node, or pass CONFIGS="1 1,2 4,...".
IFS=',' read -r -a CONFIGS <<< "${CONFIGS:-1 1,1 2,1 4,2 2,4 1}"

OUTDIR="scaling_${MODE}_${GRID}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTDIR"
echo "[sweep] param=$PARAM mode=$MODE grid=$GRID steps=$STEPS warmup=$WARMUP -> $OUTDIR"

for cfg in "${CONFIGS[@]}"; do
  read -r NP NT <<< "$cfg"
  tag="${MODE}_${GRID}_np${NP}_t${NT}"
  echo ">>> $tag  (ranks=$NP threads=$NT)"
  OMP_NUM_THREADS="$NT" OMP_PROC_BIND=close OMP_PLACES=cores \
    mpirun -np "$NP" "$BIN" "$PARAM" --mode "$MODE" --grid "$GRID" \
      --steps "$STEPS" --warmup "$WARMUP" --prefix "$OUTDIR/$tag" \
    > "$OUTDIR/$tag.log" 2>&1 \
    && echo "    ok -> $OUTDIR/${tag}_steps.jsonl" \
    || echo "    FAILED (see $OUTDIR/$tag.log)"
done

echo "[sweep] done. records in $OUTDIR/*_steps.jsonl"
