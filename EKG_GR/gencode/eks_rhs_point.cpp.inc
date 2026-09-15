// Separate CSE scope: no DENDRO_* collisions with the vacuum kernel.
{
    [[maybe_unused]] double K_scalar_src = 0.0;
    [[maybe_unused]] double At_scalar_src00 = 0.0;
    [[maybe_unused]] double At_scalar_src01 = 0.0;
    [[maybe_unused]] double At_scalar_src02 = 0.0;
    [[maybe_unused]] double At_scalar_src11 = 0.0;
    [[maybe_unused]] double At_scalar_src12 = 0.0;
    [[maybe_unused]] double At_scalar_src22 = 0.0;
    [[maybe_unused]] double Gt_scalar_src0 = 0.0;
    [[maybe_unused]] double Gt_scalar_src1 = 0.0;
    [[maybe_unused]] double Gt_scalar_src2 = 0.0;
    [[maybe_unused]] double scalar_rho = 0.0;
    [[maybe_unused]] double scalar_S = 0.0;
    [[maybe_unused]] double scalar_V = 0.0;
    [[maybe_unused]] double scalar_grad_sq = 0.0;
    #include "eks_matter_rhs.cpp.inc"
    if (eks::config.mode == "coupled") {
        K_rhs[pp] += K_scalar_src;
        At_rhs00[pp] += At_scalar_src00;
        At_rhs01[pp] += At_scalar_src01;
        At_rhs02[pp] += At_scalar_src02;
        At_rhs11[pp] += At_scalar_src11;
        At_rhs12[pp] += At_scalar_src12;
        At_rhs22[pp] += At_scalar_src22;
        Gt_rhs0[pp] += Gt_scalar_src0;
        B_rhs0[pp] += Gt_scalar_src0;
        Gt_rhs1[pp] += Gt_scalar_src1;
        B_rhs1[pp] += Gt_scalar_src1;
        Gt_rhs2[pp] += Gt_scalar_src2;
        B_rhs2[pp] += Gt_scalar_src2;
        // B correction is for the CHECKED-IN standard Gamma-driver kernel.
        // Re-derive this if you change to Rochester/first-order shift.
    }
}
