// bssn_cascade_avx_interior.inc.cpp -- drop-in replacement for the interior
// RHS triple loop using the AVX2-batched cascade (4 lanes per SIMD op).
//
// Included in BSSN_GR/src/rhs.cpp under #ifdef BSSN_USE_CASCADE_AVX.
// Expects the following to be in scope (all present in rhs.cpp at this point):
//   - State pointers: alpha, chi, K, gt0..gt5, beta0..beta2, At0..At5,
//     Gt0..Gt2, B0..B2  (const double * const)
//   - Output pointers: a_rhs, chi_rhs, K_rhs, gt_rhs00..gt_rhs22,
//     b_rhs0..b_rhs2, At_rhs00..At_rhs22, Gt_rhs0..Gt_rhs2, B_rhs0..B_rhs2
//   - Derivative pointers: grad_0_alpha .. grad2_2_2_beta2 (allocated in
//     bssnrhs_derivs.h)
//   - Block metadata: nx, ny, nz, hx, hy, hz, pmin[3], PW, bflag
//   - BSSN params: lambda[4], lambda_f[2]
//
// Writes interior (PW..n-PW) RHS values. Boundary layer is handled by
// bssn_bcs() as usual after this loop. Assumes BSSN_USE_CASCADE_AVX is
// compatible with the standard gauge (no Rochester/SSL/CAHD).

// <immintrin.h> and <cmath> are included at top of rhs.cpp under
// #ifdef BSSN_USE_CASCADE_AVX.
{
    // Pack deriv pointers for the generated cascade body
    const double *const d_al_p[3] = {grad_0_alpha, grad_1_alpha, grad_2_alpha};
    const double *const d_ch_p[3] = {grad_0_chi,   grad_1_chi,   grad_2_chi};
    const double *const d_K_p[3]  = {grad_0_K,     grad_1_K,     grad_2_K};
    const double *const d_be_p[3][3] = {
        {grad_0_beta0, grad_1_beta0, grad_2_beta0},
        {grad_0_beta1, grad_1_beta1, grad_2_beta1},
        {grad_0_beta2, grad_1_beta2, grad_2_beta2}};
    const double *const d_Gt_p[3][3] = {
        {grad_0_Gt0, grad_1_Gt0, grad_2_Gt0},
        {grad_0_Gt1, grad_1_Gt1, grad_2_Gt1},
        {grad_0_Gt2, grad_1_Gt2, grad_2_Gt2}};
    const double *const d_B_p[3][3] = {
        {grad_0_B0, grad_1_B0, grad_2_B0},
        {grad_0_B1, grad_1_B1, grad_2_B1},
        {grad_0_B2, grad_1_B2, grad_2_B2}};
    const double *const d_gt_p[6][3] = {
        {grad_0_gt0, grad_1_gt0, grad_2_gt0},
        {grad_0_gt1, grad_1_gt1, grad_2_gt1},
        {grad_0_gt2, grad_1_gt2, grad_2_gt2},
        {grad_0_gt3, grad_1_gt3, grad_2_gt3},
        {grad_0_gt4, grad_1_gt4, grad_2_gt4},
        {grad_0_gt5, grad_1_gt5, grad_2_gt5}};
    const double *const d_At_p[6][3] = {
        {grad_0_At0, grad_1_At0, grad_2_At0},
        {grad_0_At1, grad_1_At1, grad_2_At1},
        {grad_0_At2, grad_1_At2, grad_2_At2},
        {grad_0_At3, grad_1_At3, grad_2_At3},
        {grad_0_At4, grad_1_At4, grad_2_At4},
        {grad_0_At5, grad_1_At5, grad_2_At5}};
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

    // sym(i,j) -> 0..5 for the 3x3 symmetric index (inline lambda)
    auto cascade_sym_idx = [](int i, int j) {
        if (i > j) { int t = i; i = j; j = t; }
        static const int tbl[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        return tbl[i][j];
    };

    const unsigned int i_lo = PW;
    const unsigned int i_hi = nx - PW;

    if (i_hi > i_lo + 3) {  // at least one 4-batch fits
        for (unsigned int k = PW; k < nz - PW; k++) {
            const double zc = pmin[2] + k * hz;
            for (unsigned int j = PW; j < ny - PW; j++) {
                const double yc = pmin[1] + j * hy;

                auto run_batch = [&](unsigned int i_start) {
                    const unsigned int pp = i_start + nx * (j + ny * k);
                    // eta: compute scalar per-lane (RIT formula), pack to VEC
                    alignas(32) double eta_arr[4];
                    for (int l = 0; l < 4; l++) {
                        const double xl = pmin[0] + (i_start + l) * hx;
                        const double r = std::sqrt(xl*xl + yc*yc + zc*zc);
                        const double w = r / bssn::RIT_ETA_WIDTH;
                        const double a = -w*w*w*w;
                        eta_arr[l] = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER) * std::exp(a) + bssn::RIT_ETA_OUTER;
                    }
                    const __m256d eta = _mm256_load_pd(eta_arr);

                    // Generated cascade body: opens its own { ... } scope
                    #include "bssneqs_cascade_avx.cpp"
                };

                unsigned int i = i_lo;
                for (; i + 4 <= i_hi; i += 4) run_batch(i);
                if (i < i_hi) run_batch(i_hi - 4);  // shift-back tail
            }
        }
    }
}
