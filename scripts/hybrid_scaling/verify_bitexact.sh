#!/usr/bin/env bash
# Bit-exactness gate for mesh-construction changes (Phase B: threading
# buildE2EMap / buildE2NMap / performBlocksSetup / OCT2BLK).
#
# The invariant: the mesh is a deterministic function of the octree, so at a
# FIXED rank count its structure must not depend on the thread count. Any
# divergence is a race or a reordering -- and a reordered E2N map silently
# changes every numeric result forever, which is the failure this exists to stop.
#
# WHY RANKS ARE HELD FIXED: the mesh is legitimately rank-count dependent. The
# SFC partition gives each rank different elements, node numbering is rank-local,
# and getAllElements() includes ghosts. So comparing across a fixed-CORE sweep
# (R = cores/T) would show differences that are entirely correct. We vary ONLY
# OMP_NUM_THREADS, oversubscribing if need be -- this is a correctness gate, so
# timing does not matter.
#
# Digests are per construction stage, so a failure localizes: e2n differing while
# e2e matches points straight at buildE2NMap.
#
#   ./verify_bitexact.sh                  # default: R=4, T=1,2,4, depth 9
#   RANKS=8 THREADS_LIST="1 4" ./verify_bitexact.sh
#   BIN=<other build>/bssnScalingBench ./verify_bitexact.sh   # e.g. flag OFF vs ON
#   REQUIRE_MESH_THREADS=1 ./verify_bitexact.sh   # ALSO fail unless the Mesh ctor
#                                                 # really threaded (needs a
#                                                 # DENDRO_MESH_OMP=ON build)
#
# A PASS from this gate is necessary, not sufficient. It runs a handful of ranks
# on one node; "thread-invariant here" does not imply "thread-invariant at 768
# ranks", where the ghost and hanging-node structure differs. Treat a local PASS
# as a smoke test and re-run at scale before trusting a mesh change in production.
set -uo pipefail

REPO="${REPO:-/home/denv/research/dendrogr_dfvk}"
BIN="${BIN:-$REPO/build_hybrid/BSSN_GR/bssnScalingBench}"
BASE_PAR="${BASE_PAR:-$REPO/scripts/hybrid_scaling/q1.scaling.par.toml}"
OUTDIR="${OUTDIR:-$PWD/bitexact}"
RANKS="${RANKS:-4}"                       # FIXED across the sweep -- see above
THREADS_LIST="${THREADS_LIST:-1 2 4}"
DEPTH="${DEPTH:-9}"                       # small: correctness, not performance
STEPS="${STEPS:-2}"
WARMUP="${WARMUP:-1}"

mkdir -p "$OUTDIR"
[[ -x "$BIN" ]] || { echo "no binary at $BIN (build target bssnScalingBench)"; exit 2; }

PAR="$OUTDIR/verify.d${DEPTH}.par.toml"
sed "s/^BSSN_MAXDEPTH *=.*/BSSN_MAXDEPTH = $DEPTH/" "$BASE_PAR" > "$PAR"
# An unchecked sed silently leaves the BASE depth in place if the pattern ever
# moves -- the gate then runs at an unintended grid size and still says PASS.
# (And BSSN_MAXDEPTH<=6 on q1 bbh aborts mid-run, so a silent depth change is not
# harmless.) Assert the substitution actually took.
if ! grep -qE "^BSSN_MAXDEPTH *= *${DEPTH}\b" "$PAR"; then
  echo "sed did not apply: $PAR has no 'BSSN_MAXDEPTH = $DEPTH'"
  grep -nE "^BSSN_MAXDEPTH" "$PAR" | sed 's/^/    got: /'
  exit 2
fi

echo "bit-exactness gate: R=$RANKS (fixed)  T=[$THREADS_LIST]  depth=$DEPTH  bin=$BIN"

ref=""; rc=0
for T in $THREADS_LIST; do
  f="$OUTDIR/fp_t${T}.txt"
  OMP_NUM_THREADS="$T" OMP_PROC_BIND=close OMP_PLACES=cores \
    mpirun --np "$RANKS" --oversubscribe "$BIN" "$PAR" \
      --grid bbh --steps "$STEPS" --warmup "$WARMUP" --fingerprint \
      --prefix "$OUTDIR/be_t${T}" > "$OUTDIR/run_t${T}.log" 2>&1
  # cols: tag field hash count  (count = items hashed, summed over ranks)
  grep "^\[fingerprint\]" "$OUTDIR/run_t${T}.log" | awk '{print $2,$3,$4,$5}' > "$f"
  if [[ ! -s "$f" ]]; then
    echo "  T=$T: FAILED to produce digests (see $OUTDIR/run_t${T}.log)"; rc=1; continue
  fi

  # NON-VACUITY. The old check here grepped the digests for the un-hashed FNV
  # offset basis and could never fire, twice over: the constant it looked for
  # wasn't the one grUtils.cpp used (that one was missing a digit), and even
  # fixed it tests a pre-combine_ranks value against post-combine_ranks output,
  # so the bare basis never appears at all. Digest-sniffing cannot answer this.
  # The counts can: a stage that hashed zero items is vacuous no matter how
  # respectable its hash looks. Assert on the stages this gate exists to protect.
  for stage in elements e2e e2n blocks; do
    n=$(awk -v s="$stage" '$2==s {print $4}' "$f")
    if [[ -z "$n" ]]; then
      echo "  T=$T: SUSPECT -- no '$stage' digest emitted at all"; rc=1
    elif [[ "$n" == "0" ]]; then
      echo "  T=$T: SUSPECT -- '$stage' hashed 0 items (vacuous digest)"; rc=1
    fi
  done

  # THREADING CANARY. The digests cannot tell "threaded and correct" from
  # "threading never happened" -- the latter matches trivially at every T and
  # reports PASS. That is the 9bb8f45 failure, which was bit-exact and silent for
  # two days. omp_ctor_threads is read at Mesh-ctor entry inside dendrolib, so it
  # answers specifically whether the CTOR threaded (bssn::BSSN_HYBRID_NTHREADS
  # cannot: it is the RHS thread count and is set after the Mesh already exists).
  ctor_thr=$(grep -m1 "^\[meshcanary\]" "$OUTDIR/run_t${T}.log" | sed 's/.*omp_ctor_threads=//')
  if [[ -z "$ctor_thr" ]]; then
    echo "  T=$T: no [meshcanary] line -- cannot tell whether the ctor threaded"
    if [[ "${REQUIRE_MESH_THREADS:-0}" == "1" ]]; then rc=1; fi
  elif [[ "${REQUIRE_MESH_THREADS:-0}" == "1" && "$ctor_thr" != "$T" ]]; then
    echo "  T=$T: *** VACUOUS *** Mesh ctor saw omp_ctor_threads=$ctor_thr, expected $T"
    echo "        -> the ctor did NOT thread; matching digests here prove nothing."
    rc=1
  else
    echo "  T=$T: mesh ctor threads=$ctor_thr"
  fi
  if [[ -z "$ref" ]]; then
    ref="$f"; echo "  T=$T: reference"; sed 's/^/      /' "$f"
  elif diff -q "$ref" "$f" >/dev/null; then
    echo "  T=$T: MATCH"
  else
    echo "  T=$T: *** MISMATCH ***  (the differing stage names the culprit)"
    diff "$ref" "$f" | sed 's/^/      /'
    rc=1
  fi
done

echo
if (( rc == 0 )); then
  echo "PASS - mesh + evolved state are bit-identical across threads at R=$RANKS"
else
  echo "FAIL - see above. A stage listed here is not thread-invariant."
fi
exit $rc
