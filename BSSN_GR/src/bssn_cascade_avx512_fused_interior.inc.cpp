// bssn_cascade_avx512_fused_interior.inc.cpp -- 8-wide AVX-512 fused cascade.
// Only used when block is wide enough (interior width >= 8). Narrow blocks
// are handled by the AVX2 non-fused path (full precompute) in the else branch
// of rhs.cpp.

// <immintrin.h> and <cmath> are included at top of rhs.cpp under the
// BSSN_USE_CASCADE_AVX512* flag.
{
    // Pack deriv pointers. Same layout as AVX2 fused: only mixed 2nd derivs
    // need arrays; 1st and pure 2nd derivs are computed inline.
    const double *const d2_al_p[6] = {
        grad2_0_0_alpha, grad2_0_1_alpha, grad2_0_2_alpha,
        grad2_1_1_alpha, grad2_1_2_alpha, grad2_2_2_alpha};
    const double *const d2_ch_p[6] = {
        grad2_0_0_chi, grad2_0_1_chi, grad2_0_2_chi,
        grad2_1_1_chi, grad2_1_2_chi, grad2_2_2_chi};
    const double *const d2_be_p[3][6] = {
        {grad2_0_0_beta0, grad2_0_1_beta0, grad2_0_2_beta0,
         grad2_1_1_beta0, grad2_1_2_beta0, grad2_2_2_beta0},
        {grad2_0_0_beta1, grad2_0_1_beta1, grad2_0_2_beta1,
         grad2_1_1_beta1, grad2_1_2_beta1, grad2_2_2_beta1},
        {grad2_0_0_beta2, grad2_0_1_beta2, grad2_0_2_beta2,
         grad2_1_1_beta2, grad2_1_2_beta2, grad2_2_2_beta2}};
    const double *const d2_gt_p[6][6] = {
        {grad2_0_0_gt0, grad2_0_1_gt0, grad2_0_2_gt0,
         grad2_1_1_gt0, grad2_1_2_gt0, grad2_2_2_gt0},
        {grad2_0_0_gt1, grad2_0_1_gt1, grad2_0_2_gt1,
         grad2_1_1_gt1, grad2_1_2_gt1, grad2_2_2_gt1},
        {grad2_0_0_gt2, grad2_0_1_gt2, grad2_0_2_gt2,
         grad2_1_1_gt2, grad2_1_2_gt2, grad2_2_2_gt2},
        {grad2_0_0_gt3, grad2_0_1_gt3, grad2_0_2_gt3,
         grad2_1_1_gt3, grad2_1_2_gt3, grad2_2_2_gt3},
        {grad2_0_0_gt4, grad2_0_1_gt4, grad2_0_2_gt4,
         grad2_1_1_gt4, grad2_1_2_gt4, grad2_2_2_gt4},
        {grad2_0_0_gt5, grad2_0_1_gt5, grad2_0_2_gt5,
         grad2_1_1_gt5, grad2_1_2_gt5, grad2_2_2_gt5}};

    auto cascade_sym_idx = [](int i, int j) {
        if (i > j) { int t = i; i = j; j = t; }
        static const int tbl[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        return tbl[i][j];
    };

    const unsigned int i_lo = PW;
    const unsigned int i_hi = nx - PW;

    for (unsigned int k = PW; k < nz - PW; k++) {
        const double zc = pmin[2] + k * hz;
        for (unsigned int j = PW; j < ny - PW; j++) {
            const double yc = pmin[1] + j * hy;

            auto run_batch = [&](unsigned int i_start) {
                const unsigned int pp = i_start + nx * (j + ny * k);
                alignas(64) double eta_arr[8];
                for (int l = 0; l < 8; l++) {
                    const double xl = pmin[0] + (i_start + l) * hx;
                    const double r = std::sqrt(xl*xl + yc*yc + zc*zc);
                    const double w = r / bssn::RIT_ETA_WIDTH;
                    const double a = -w*w*w*w;
                    eta_arr[l] = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER) * std::exp(a) + bssn::RIT_ETA_OUTER;
                }
                const __m512d eta = _mm512_load_pd(eta_arr);

                #include "bssneqs_cascade_avx512_fused.cpp"
            };

            unsigned int i = i_lo;
            for (; i + 8 <= i_hi; i += 8) run_batch(i);
            if (i < i_hi) run_batch(i_hi - 8);  // shift-back tail
        }
    }
}
