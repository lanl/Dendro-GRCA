#include "sommerfeld_bcs.h"

#include <cmath>

#include "dendro.h"  // OCT_DIR_* face flags

namespace dgr {

// Index into a zip-block laid out as i + nx*(j + ny*k).
#define IDX3(i, j, k, nx, ny) ((i) + (nx) * ((j) + (ny) * (k)))

static inline void apply_sommerfeld_cell(const SommerfeldVar *__restrict vars,
                                         unsigned int nVars, unsigned int pp,
                                         double x, double y, double z,
                                         double inv_r) {
    for (unsigned int v = 0; v < nVars; v++) {
        const SommerfeldVar &V = vars[v];
        const double advect =
            x * V.dxf[pp] + y * V.dyf[pp] + z * V.dzf[pp];
        V.f_rhs[pp] = -inv_r * (advect + V.falloff * (V.f[pp] - V.asymptotic));
    }
}

void sommerfeld_outer_bcs(const SommerfeldVar *vars, unsigned int nVars,
                          const double *pmin, const double *pmax,
                          const unsigned int *sz, unsigned int pad,
                          unsigned int bflag) {
    if (bflag == 0u || nVars == 0u) return;

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const double hx       = (pmax[0] - pmin[0]) / (double)(nx - 1);
    const double hy       = (pmax[1] - pmin[1]) / (double)(ny - 1);
    const double hz       = (pmax[2] - pmin[2]) / (double)(nz - 1);

    const unsigned int ib = pad;
    const unsigned int jb = pad;
    const unsigned int kb = pad;
    const unsigned int ie = nx - pad;
    const unsigned int je = ny - pad;
    const unsigned int ke = nz - pad;

    // LEFT/RIGHT: x fixed, inner loop along j is stride-nx (no SIMD).
    if (bflag & (1u << OCT_DIR_LEFT)) {
        const double x  = pmin[0] + ib * hx;
        const double x2 = x * x;
        for (unsigned int k = kb; k < ke; k++) {
            const double z  = pmin[2] + k * hz;
            const double z2 = z * z;
            for (unsigned int j = jb; j < je; j++) {
                const double y     = pmin[1] + j * hy;
                const double inv_r = 1.0 / std::sqrt(x2 + y * y + z2);
                const unsigned int pp = IDX3(ib, j, k, nx, ny);
                apply_sommerfeld_cell(vars, nVars, pp, x, y, z, inv_r);
            }
        }
    }

    if (bflag & (1u << OCT_DIR_RIGHT)) {
        const double x  = pmin[0] + (ie - 1) * hx;
        const double x2 = x * x;
        for (unsigned int k = kb; k < ke; k++) {
            const double z  = pmin[2] + k * hz;
            const double z2 = z * z;
            for (unsigned int j = jb; j < je; j++) {
                const double y     = pmin[1] + j * hy;
                const double inv_r = 1.0 / std::sqrt(x2 + y * y + z2);
                const unsigned int pp = IDX3(ie - 1, j, k, nx, ny);
                apply_sommerfeld_cell(vars, nVars, pp, x, y, z, inv_r);
            }
        }
    }

    // DOWN/UP/BACK/FRONT: inner i-loop is unit stride -> SIMD.
    if (bflag & (1u << OCT_DIR_DOWN)) {
        const double y  = pmin[1] + jb * hy;
        const double y2 = y * y;
        for (unsigned int k = kb; k < ke; k++) {
            const double z  = pmin[2] + k * hz;
            const double z2 = z * z;
#pragma omp simd
            for (unsigned int i = ib; i < ie; i++) {
                const double x     = pmin[0] + i * hx;
                const double inv_r = 1.0 / std::sqrt(x * x + y2 + z2);
                const unsigned int pp = IDX3(i, jb, k, nx, ny);
                apply_sommerfeld_cell(vars, nVars, pp, x, y, z, inv_r);
            }
        }
    }

    if (bflag & (1u << OCT_DIR_UP)) {
        const double y  = pmin[1] + (je - 1) * hy;
        const double y2 = y * y;
        for (unsigned int k = kb; k < ke; k++) {
            const double z  = pmin[2] + k * hz;
            const double z2 = z * z;
#pragma omp simd
            for (unsigned int i = ib; i < ie; i++) {
                const double x     = pmin[0] + i * hx;
                const double inv_r = 1.0 / std::sqrt(x * x + y2 + z2);
                const unsigned int pp = IDX3(i, je - 1, k, nx, ny);
                apply_sommerfeld_cell(vars, nVars, pp, x, y, z, inv_r);
            }
        }
    }

    if (bflag & (1u << OCT_DIR_BACK)) {
        const double z  = pmin[2] + kb * hz;
        const double z2 = z * z;
        for (unsigned int j = jb; j < je; j++) {
            const double y  = pmin[1] + j * hy;
            const double y2 = y * y;
#pragma omp simd
            for (unsigned int i = ib; i < ie; i++) {
                const double x     = pmin[0] + i * hx;
                const double inv_r = 1.0 / std::sqrt(x * x + y2 + z2);
                const unsigned int pp = IDX3(i, j, kb, nx, ny);
                apply_sommerfeld_cell(vars, nVars, pp, x, y, z, inv_r);
            }
        }
    }

    if (bflag & (1u << OCT_DIR_FRONT)) {
        const double z  = pmin[2] + (ke - 1) * hz;
        const double z2 = z * z;
        for (unsigned int j = jb; j < je; j++) {
            const double y  = pmin[1] + j * hy;
            const double y2 = y * y;
#pragma omp simd
            for (unsigned int i = ib; i < ie; i++) {
                const double x     = pmin[0] + i * hx;
                const double inv_r = 1.0 / std::sqrt(x * x + y2 + z2);
                const unsigned int pp = IDX3(i, j, ke - 1, nx, ny);
                apply_sommerfeld_cell(vars, nVars, pp, x, y, z, inv_r);
            }
        }
    }
}

}  // namespace dgr
