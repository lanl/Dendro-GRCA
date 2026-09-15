// Compute ONLY the mixed 2nd derivatives needed by the fused AVX cascade.
// The 1st derivs and pure 2nd derivs (xx/yy/zz) are computed inline by the
// fused cascade kernel; the mixed 2nds (xy/xz/yz) require 2D stencils and
// stay precomputed.
//
// For each variable v in {alpha, chi, beta0, beta1, beta2, gt0..gt5}, we need:
//   - grad_0_v, grad_1_v  (as intermediate steps for the mixed 2nds)
//   - grad2_0_1_v = d/dy (grad_0_v),  grad2_0_2_v = d/dz (grad_0_v),
//     grad2_1_2_v = d/dz (grad_1_v)
//
// Total: 11 vars × (2 + 3) = 55 deriv passes, vs 138 in the full workspace.
// grad_2_v is NOT computed here (only needed for 3rd-axis 1st-deriv and pure
// 2nd xx/yy/zz; both inlined by the fused cascade).

// alpha
deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
deriv_y(grad_1_alpha, alpha, hy, sz, bflag);
deriv_y(grad2_0_1_alpha, grad_0_alpha, hy, sz, bflag);
deriv_z(grad2_0_2_alpha, grad_0_alpha, hz, sz, bflag);
deriv_z(grad2_1_2_alpha, grad_1_alpha, hz, sz, bflag);

// chi
deriv_x(grad_0_chi, chi, hx, sz, bflag);
deriv_y(grad_1_chi, chi, hy, sz, bflag);
deriv_y(grad2_0_1_chi, grad_0_chi, hy, sz, bflag);
deriv_z(grad2_0_2_chi, grad_0_chi, hz, sz, bflag);
deriv_z(grad2_1_2_chi, grad_1_chi, hz, sz, bflag);

// beta0
deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
deriv_y(grad_1_beta0, beta0, hy, sz, bflag);
deriv_y(grad2_0_1_beta0, grad_0_beta0, hy, sz, bflag);
deriv_z(grad2_0_2_beta0, grad_0_beta0, hz, sz, bflag);
deriv_z(grad2_1_2_beta0, grad_1_beta0, hz, sz, bflag);

// beta1
deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
deriv_y(grad_1_beta1, beta1, hy, sz, bflag);
deriv_y(grad2_0_1_beta1, grad_0_beta1, hy, sz, bflag);
deriv_z(grad2_0_2_beta1, grad_0_beta1, hz, sz, bflag);
deriv_z(grad2_1_2_beta1, grad_1_beta1, hz, sz, bflag);

// beta2
deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
deriv_y(grad_1_beta2, beta2, hy, sz, bflag);
deriv_y(grad2_0_1_beta2, grad_0_beta2, hy, sz, bflag);
deriv_z(grad2_0_2_beta2, grad_0_beta2, hz, sz, bflag);
deriv_z(grad2_1_2_beta2, grad_1_beta2, hz, sz, bflag);

// gt0..gt5 (6 symmetric components)
deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
deriv_y(grad_1_gt0, gt0, hy, sz, bflag);
deriv_y(grad2_0_1_gt0, grad_0_gt0, hy, sz, bflag);
deriv_z(grad2_0_2_gt0, grad_0_gt0, hz, sz, bflag);
deriv_z(grad2_1_2_gt0, grad_1_gt0, hz, sz, bflag);

deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
deriv_y(grad_1_gt1, gt1, hy, sz, bflag);
deriv_y(grad2_0_1_gt1, grad_0_gt1, hy, sz, bflag);
deriv_z(grad2_0_2_gt1, grad_0_gt1, hz, sz, bflag);
deriv_z(grad2_1_2_gt1, grad_1_gt1, hz, sz, bflag);

deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
deriv_y(grad_1_gt2, gt2, hy, sz, bflag);
deriv_y(grad2_0_1_gt2, grad_0_gt2, hy, sz, bflag);
deriv_z(grad2_0_2_gt2, grad_0_gt2, hz, sz, bflag);
deriv_z(grad2_1_2_gt2, grad_1_gt2, hz, sz, bflag);

deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
deriv_y(grad_1_gt3, gt3, hy, sz, bflag);
deriv_y(grad2_0_1_gt3, grad_0_gt3, hy, sz, bflag);
deriv_z(grad2_0_2_gt3, grad_0_gt3, hz, sz, bflag);
deriv_z(grad2_1_2_gt3, grad_1_gt3, hz, sz, bflag);

deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
deriv_y(grad_1_gt4, gt4, hy, sz, bflag);
deriv_y(grad2_0_1_gt4, grad_0_gt4, hy, sz, bflag);
deriv_z(grad2_0_2_gt4, grad_0_gt4, hz, sz, bflag);
deriv_z(grad2_1_2_gt4, grad_1_gt4, hz, sz, bflag);

deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
deriv_y(grad_1_gt5, gt5, hy, sz, bflag);
deriv_y(grad2_0_1_gt5, grad_0_gt5, hy, sz, bflag);
deriv_z(grad2_0_2_gt5, grad_0_gt5, hz, sz, bflag);
deriv_z(grad2_1_2_gt5, grad_1_gt5, hz, sz, bflag);
