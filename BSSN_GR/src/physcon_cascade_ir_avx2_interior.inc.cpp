// 4-wide AVX2 driver for physconeqs_cascade_ir_avx.cpp; the AVX-512 twin runs
// the same body 8 lanes wide. psi4's tetrad needs the point coordinates: y and
// z are fixed within an i-batch, x steps with the lane.
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
#define PHYSCON_LANES 4

    const unsigned int i_lo = PW;
    const unsigned int i_hi = nx - PW;

    if (i_hi >= i_lo + PHYSCON_LANES) {
        for (unsigned int k = PW; k < nz - PW; k++) {
            const double zc = pmin[2] + k * hz;
            for (unsigned int j = PW; j < ny - PW; j++) {
                const double yc = pmin[1] + j * hy;

                auto run_batch = [&](unsigned int i_start) {
                    const unsigned int pp = i_start + nx * (j + ny * k);
                    alignas(32) double x_arr[PHYSCON_LANES];
                    for (int l = 0; l < PHYSCON_LANES; l++)
                        x_arr[l] = pmin[0] + (i_start + l) * hx;

                    const VEC x = _mm256_load_pd(x_arr);
                    const VEC y = VSET(yc);
                    const VEC z = VSET(zc);

                    {
#include "physconeqs_cascade_ir_avx.cpp"
                    }

                    // The tetrad is degenerate on the z axis; psi4 is 0 there.
                    for (int l = 0; l < PHYSCON_LANES; l++) {
                        if (fabs(x_arr[l]) <= 1e-7 && fabs(yc) <= 1e-7) {
                            psi4_real[pp + l] = 0.0;
                            psi4_img[pp + l]  = 0.0;
                        }
                    }
                };

                unsigned int i = i_lo;
                for (; i + PHYSCON_LANES <= i_hi; i += PHYSCON_LANES)
                    run_batch(i);
                // Shift-back tail: re-writes lanes; every output is pure.
                if (i < i_hi) run_batch(i_hi - PHYSCON_LANES);
            }
        }
        cascade_done = true;
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
#undef PHYSCON_LANES
}
