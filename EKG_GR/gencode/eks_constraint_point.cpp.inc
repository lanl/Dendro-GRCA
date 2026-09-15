if(eks::config.mode=="coupled") {
    // CodeGen/constraints.py stores the covariant momentum residual M_i.
    // No extra chi is attached to the -8*pi*S_i term in this convention.
    constexpr double pi_const=3.141592653589793238462643383279502884;
    const double mu=eks::config.mu, fa=eks::config.fa, phi=eks_phi[pp];
    const double V=std::string(EKS_GENERATED_POTENTIAL)=="axion"
        ? 2.0*mu*mu*fa*fa*std::pow(std::sin(phi/(2.0*fa)),2)
        : 0.5*mu*mu*phi*phi;
    const double rho=eks::rho_scalar(
        {gt0[pp],gt1[pp],gt2[pp],gt3[pp],gt4[pp],gt5[pp]},chi[pp],
        {eks_dx[pp],eks_dy[pp],eks_dz[pp]},eks_pi[pp],V);
    ham[pp]-=16.0*pi_const*rho;
    mom0[pp]-=8.0*pi_const*eks_pi[pp]*eks_dx[pp];
    mom1[pp]-=8.0*pi_const*eks_pi[pp]*eks_dy[pp];
    mom2[pp]-=8.0*pi_const*eks_pi[pp]*eks_dz[pp];
    // The inherited Psi4 expression has NOT been matter-corrected here.
    // GW extraction must remain disabled until that separate audit is done.
}
