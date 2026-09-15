# Hybrid MPI+OpenMP scaling — handoff

Status as of 2026-07-13. Branch `experimental` on both repos, all pushed.

## Goal
Characterize **pure-MPI vs hybrid MPI+OpenMP** for DendroGR BSSN across HPC
architectures, multi-node, with real inter-node communication. Same binary both
modes — only `OMP_NUM_THREADS` differs (T=1 pure-MPI, T>1 hybrid) at fixed total
cores. Deliverable is shaping up as a **characterization** (speedup vs arch, vs
node count, and *why*), not a new optimization.

## Headline result (production `ob0` BBH, q1)
Hybrid wins, and it's architecture-dependent. All at N=2, T=4 vs each arch's own
pure-MPI baseline; reproduced 4× on SPR:

| Arch | cores/node (sockets) | T=4 speedup | notes |
|------|----------------------|-------------|-------|
| SKX  | 48 (2×24)  | 1.53× | saturates at T=2 already (24-core socket) |
| ICX  | 80 (2×40)  | 2.18× | sits between |
| SPR  | 112 (2×56) | 2.77× | best; big socket keeps T=4 NUMA-local |

**T=4 is the ceiling** on SPR — T=8 regresses (2.47×) because ghost_comm climbs
back up. Node sweep: N=2 T=4 is the clean point; N≥4 is under-saturated on q1
(1–2 blocks/rank) — that's what q8 fixes.

## The mechanism (important — it's not what we first thought)
The hybrid win is **OpenMP dynamic-schedule load-balancing of heterogeneous large
blocks**, NOT memory bandwidth. Proof: the `ob31` control (uniform, no fusion, all
blocks tiny+equal, well-saturated) is **flat, 1.03×**. If the win were bandwidth
contention, ob31 (more memory-bound) would win *more*; it wins nothing. With few
big size-varying blocks (`ob0` at high rank count → ~5 blocks/rank), pure-MPI
load-balances terribly across ranks; the threaded `schedule(dynamic)` block loop
rebalances within a rank. Uniform tiny blocks are already balanced → nothing to fix.

Consequence: the hybrid path is already mature (RHS, unzip via `unzip_scatter`,
constraints, ghost pack/unpack all threaded; DVectors NUMA first-touched). There's
no "unthreaded hot loop" to grab. Further speedup would come from better *load
balance* (weighted MPI partitioning, block-loop chunking), not more threading.

## Files (`scripts/hybrid_scaling/`)
- `run_hybrid_scaling.sh` — build-once, sweep threads, parse JSONL→CSV. SITE presets
  (stampede3 default / chpc / none), memory preflight, `PROFILE=roofline`.
- `q1.scaling.par.toml` — self-contained q=1 BBH (~1100 blocks). Good for N=1–2.
- `q8.scaling.par.toml` — q=8 (MAXDEPTH=17), several× more blocks for N=4/N=8.
- `README.md` — how to run + reading the CSV + sizing + roofline.

## Commits
**Dendro-GR** `experimental`: `de0842c` OCT2BLK knob · `3a07a56` memory preflight ·
`8421134` roofline mode · `0b8b877` q8 parfile.
**Dendro-5.01** `experimental`: `87429da` `getDegOfFreedomUnZip` 64-bit fix (latent).

## How to run
```bash
# Stampede3 (default site: modules, sapphirerapids, 112 cores, Omni-Path all preset)
sbatch -A ACCT -p spr --nodes=2 run_hybrid_scaling.sh                 # q1, N=2
PARFILE=q8.scaling.par.toml sbatch -A ACCT -p spr --nodes=4 run_hybrid_scaling.sh
SITE=chpc sbatch -A ACCT -p soc-np -C csl --nodes=2 run_hybrid_scaling.sh
```
Output → `results_<jobid>/hybrid_scaling_<jobid>.csv`. Headline col `rk_step_s`;
`speedup_vs_purempi` >1 = hybrid wins; watch `blocks_per_rank` (want ≥16).

## Open items / next steps
1. **Confirm q8 block count on the cluster.** Verified locally only that it parses/
   builds/runs (at reduced depth); the count at MAXDEPTH=17 is *unverified* (didn't
   run full-depth locally — memory). First run prints `blocks_per_rank`; if short of
   ~16 at target N, bump `BSSN_MAXDEPTH` to 18 in the parfile. Then rerun N=4/N=8.
2. Optional `q4` parfile as a lighter fallback if q8@17 is heavier than wanted.
3. **Fabric + pinning confirmation** was recommended but never explicitly verified
   (`FI_PROVIDER=opx` active? `mpirun --report-bindings` sane?). SPR runs used the
   stampede3 preset so fabric is *probably* right; worth one confirming run.
4. If the characterization is the deliverable, it's nearly done — SPR/ICX/SKX at N=2
   plus clean N=4/N=8 on q8, the ob0-vs-ob31 control, and the T-sweep ceiling.

## Gotchas / lessons
- **uniform grid is a dead end for the hybrid comparison**: `ob0` uniform OOMs (every
  block full-size), `ob31` uniform is flat (no win). Use higher mass ratio (q8) for
  more blocks instead. `MEM_CHECK` guards uniform against OOM.
- **Memory**: no-fusion deep uniform explodes (each tiny block carries a full ghost
  zone); a `bad_alloc` inside an OMP region surfaces as a corrupt-pointer SIGSEGV, not
  a clean OOM. Do the memory math; run big meshes on the cluster, not a laptop.
- **`OCT2BLK_COARSEST_LEV`**: `0` = fused big blocks (production/CPU default, shows the
  win), `31` = no fusion (CUDA setting; many tiny blocks). Baked per-build into
  `build_hybrid_ob<N>/`. `GRID=uniform` auto-picks 31, `bbh` picks 0.
- **Don't** point `PARFILE` at `BSSN_GR/pars/q1.par.toml` — ID_TYPE=0 needs an external
  `tp_*` TwoPunctures file and segfaults without it. The shipped parfiles use ID_TYPE=1.
- **likwid is not on Stampede3** (`module spider likwid` empty), so `PROFILE=roofline`
  can't run there. It's also premature — the bottleneck is load balance, not bandwidth.
  CHPC may have likwid. A counter-free STREAM+derived fallback was scoped but not built.
- Collaborator's first builds grabbed **system gcc 11.5**, not the module gcc 15 — same
  binary both modes so it doesn't bias the comparison, but rebuild with 15 for
  representative absolute numbers.
