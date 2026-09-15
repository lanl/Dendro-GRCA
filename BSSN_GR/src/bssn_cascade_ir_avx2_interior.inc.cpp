// bssn_cascade_ir_avx2_interior.inc.cpp -- 4-wide twin of
// bssn_cascade_ir_avx512_interior.inc.cpp (same IR-generated body; AVX2
// macros and 4-lane batches). Used for blocks whose interior x-extent is
// too narrow for an 8-batch.
{
#define VEC __m256d
#define VLOAD(p)      _mm256_loadu_pd(p)
#define VSTORE(p, v)  _mm256_storeu_pd((p), (v))
#define VSET(x)       _mm256_set1_pd(x)
#define VADD(a, b)    _mm256_add_pd((a), (b))
#define VSUB(a, b)    _mm256_sub_pd((a), (b))
#define VMUL(a, b)    _mm256_mul_pd((a), (b))
#define VDIV(a, b)    _mm256_div_pd((a), (b))
#define VFMA(a, b, c) _mm256_fmadd_pd((a), (b), (c))
#define VFNMADD(a, b, c) _mm256_fnmadd_pd((a), (b), (c))
#define VSQRT(a)      _mm256_sqrt_pd(a)

    const VEC ssl_fac =
        VSET(-h_ssl * std::exp(-0.5 * t * t / (sig_ssl * sig_ssl)));
    const VEC cahd_coef =
        VSET(BSSN_CAHD_C * dx_i * dx_i / (1.0 + 10.0 * dx_i * dx_i) / dt);

    const unsigned int i_lo = PW;
    const unsigned int i_hi = nx - PW;

    if (i_hi > i_lo + 3) {  // at least one 4-batch fits
        for (unsigned int k = PW; k < nz - PW; k++) {
            const double zc = pmin[2] + k * hz;
            for (unsigned int j = PW; j < ny - PW; j++) {
                const double yc = pmin[1] + j * hy;

                auto run_batch = [&](unsigned int i_start) {
                    const unsigned int pp = i_start + nx * (j + ny * k);
                    alignas(32) double eta_arr[4];
                    for (int l = 0; l < 4; l++) {
                        const double xl = pmin[0] + (i_start + l) * hx;
                        const double r  = std::sqrt(xl*xl + yc*yc + zc*zc);
                        const double w  = r / bssn::RIT_ETA_WIDTH;
                        const double a  = -w*w*w*w;
                        eta_arr[l] = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER)
                                     * std::exp(a) + bssn::RIT_ETA_OUTER;
                    }
                    const VEC eta = _mm256_load_pd(eta_arr);

                    {
#include "bssneqs_cascade_ir_avx.cpp"
                    }
                };

                unsigned int i = i_lo;
                for (; i + 4 <= i_hi; i += 4) run_batch(i);
                if (i < i_hi) run_batch(i_hi - 4);  // shift-back tail
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
