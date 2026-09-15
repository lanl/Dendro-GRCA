// bssn_cascade_ir_avx512_interior.inc.cpp -- interior RHS triple loop using
// the IR-generated polynomial-cascade body (bssneqs_cascade_ir_avx.cpp),
// 8-wide AVX-512. Unlike the hand-ported cascade, the body is emitted by
// codegen/cascade_emit.py from the symbolic spec and carries the SSL+CAHD
// gauge terms via runtime coefficients (ssl_fac/cahd_coef): h_ssl = 0 or
// BSSN_CAHD_C = 0 disables the term with no rebuild.
//
// Expects in scope (all present in rhs_experimental.cpp at the include
// point): state/output/derivative pointers with DendroGR names, nx/ny/nz,
// PW, pmin, hx/hy/hz, lambda[4], lambda_f[2], t, h_ssl, sig_ssl, dt, dx_i,
// BSSN_CAHD_C. Narrow blocks (interior x-extent < 8) are handled by the
// caller via the 4-wide twin (bssn_cascade_ir_avx2_interior.inc.cpp).
{
#define VEC __m512d
#define VLOAD(p)      _mm512_loadu_pd(p)
#define VSTORE(p, v)  _mm512_storeu_pd((p), (v))
#define VSET(x)       _mm512_set1_pd(x)
#define VADD(a, b)    _mm512_add_pd((a), (b))
#define VSUB(a, b)    _mm512_sub_pd((a), (b))
#define VMUL(a, b)    _mm512_mul_pd((a), (b))
#define VDIV(a, b)    _mm512_div_pd((a), (b))
#define VFMA(a, b, c) _mm512_fmadd_pd((a), (b), (c))
#define VFNMADD(a, b, c) _mm512_fnmadd_pd((a), (b), (c))
#define VSQRT(a)      _mm512_sqrt_pd(a)

    // Per-call gauge coefficients (constant across the block).
    const VEC ssl_fac =
        VSET(-h_ssl * std::exp(-0.5 * t * t / (sig_ssl * sig_ssl)));
    const VEC cahd_coef =
        VSET(BSSN_CAHD_C * dx_i * dx_i / (1.0 + 10.0 * dx_i * dx_i) / dt);

    const unsigned int i_lo = PW;
    const unsigned int i_hi = nx - PW;

    if (i_hi > i_lo + 7) {  // at least one 8-batch fits
        for (unsigned int k = PW; k < nz - PW; k++) {
            const double zc = pmin[2] + k * hz;
            for (unsigned int j = PW; j < ny - PW; j++) {
                const double yc = pmin[1] + j * hy;

                auto run_batch = [&](unsigned int i_start) {
                    const unsigned int pp = i_start + nx * (j + ny * k);
                    // eta: per-lane RIT formula, packed to VEC
                    alignas(64) double eta_arr[8];
                    for (int l = 0; l < 8; l++) {
                        const double xl = pmin[0] + (i_start + l) * hx;
                        const double r  = std::sqrt(xl*xl + yc*yc + zc*zc);
                        const double w  = r / bssn::RIT_ETA_WIDTH;
                        const double a  = -w*w*w*w;
                        eta_arr[l] = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER)
                                     * std::exp(a) + bssn::RIT_ETA_OUTER;
                    }
                    const VEC eta = _mm512_load_pd(eta_arr);

                    // IR-generated body (VLOADs its own leaves at offset pp)
                    {
#include "bssneqs_cascade_ir_avx.cpp"
                    }
                };

                unsigned int i = i_lo;
                for (; i + 8 <= i_hi; i += 8) run_batch(i);
                if (i < i_hi) run_batch(i_hi - 8);  // shift-back tail
            }
        }
    }

#undef VEC
#undef VLOAD
#undef VSTORE
#undef VSET
#undef VADD
#undef VSUB
#undef VMUL
#undef VDIV
#undef VFMA
#undef VFNMADD
#undef VSQRT
}
