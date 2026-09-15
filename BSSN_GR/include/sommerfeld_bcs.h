// Solver-agnostic Sommerfeld outer-boundary BCs.
//
// Hoists inv_r = 1/sqrt(x^2+y^2+z^2) out of the per-variable apply: one
// sqrt per cell drives all N variables (vs the legacy bssn_bcs() that
// redid the sqrt per variable).
//
// Reuse rule: add a variable -> append to the SommerfeldVar array. Add a
// new BC formula -> write a sibling function. Don't generalize this into
// a god-function with per-cell callbacks (kills SIMD inlining).
#pragma once

namespace dgr {

struct SommerfeldVar {
    double *f_rhs;
    const double *f;
    const double *dxf;
    const double *dyf;
    const double *dzf;
    double falloff;     // outgoing-wave decay rate
    double asymptotic;  // far-field limit
};

// Applies Sommerfeld BC to faces flagged in bflag for all nVars.
// pad = padding width; faces start at index `pad`.
// bflag bits: OCT_DIR_LEFT/RIGHT/DOWN/UP/BACK/FRONT.
void sommerfeld_outer_bcs(const SommerfeldVar *vars, unsigned int nVars,
                          const double *pmin, const double *pmax,
                          const unsigned int *sz, unsigned int pad,
                          unsigned int bflag);

}  // namespace dgr
