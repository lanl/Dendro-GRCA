# Porting testing scripts to HPC sites

When a runner under `scripts/` fails on a new machine it's almost always modules,
a hardcoded path, or the MPI fabric:

- **Modules** — write the sequence explicitly per site (a `module use <path>` step
  can't fit in `module load $VAR`, and names differ: `mkl` vs `intel-oneapi-mkl`).
  Give an escape hatch (e.g. `SITE=none`) for users who load modules elsewhere.
- **Paths** — no personal defaults (`$HOME/research/...`); derive from the repo or
  leave empty and let the preflight fail loudly. A wrong-but-present path is worst.
- **MPI fabric** — the wrong provider silently falls back to TCP and inflates
  inter-node comm times. Set it per site, make it overridable (`FI_PROVIDER=...`).

`hybrid_scaling/run_hybrid_scaling.sh` bakes these into a `SITE` switch
(`stampede3` default, `chpc`, `none`) — the reference for new scripts. What each
preset does:

## Stampede3 (TACC) — Omni-Path, `SITE=stampede3`

```bash
module load autotools/1.4 cmake/4.1.1 python/3.12.11 gcc/15.1.0 mkl/25.1 gsl/2.8
module use /scratch/projects/compilers/modulefiles   # openmpi is off the default tree
module load openmpi/5.0.9
export OMPI_MCA_mtl=ofi FI_PROVIDER=opx               # REQUIRED: Omni-Path, not TCP
# CPU_ARCH=sapphirerapids, CORES_PER_NODE=112, DENDROLIB_DIR=$HOME/Projects/dendro/DVK_experimental/Dendro-5.01
```

## CHPC (Utah) — csl, `SITE=chpc`

```bash
module load gcc/15.1.0 intel-oneapi-mkl/2025.3.1 openmpi/5.0.3
# CPU_ARCH=cascadelake, CORES_PER_NODE=40, DENDROLIB_DIR=$HOME/research/dendrolib_dfvk_copy
```
