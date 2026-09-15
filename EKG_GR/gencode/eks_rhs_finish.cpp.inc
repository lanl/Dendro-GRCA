// Included at the END of bssnrhs: after the original metric BC/KO code.
// Scalar gradients in eks_derivs have NOT been overwritten by metric KO.
if (eks::config.boundary == "sommerfeld" && bflag) {
    bssn_bcs(scalarPhi_rhs,scalarPhi,grad_0_scalarPhi,grad_1_scalarPhi,
              grad_2_scalarPhi,pmin,pmax,1.0,0.0,sz,bflag);
    bssn_bcs(scalarPi_rhs,scalarPi,grad_0_scalarPi,grad_1_scalarPi,
              grad_2_scalarPi,pmin,pmax,1.0,0.0,sz,bflag);
}
if (eks::config.ko_sigma > 0.0) {
    // Reuse scalar derivative scratch ONLY AFTER all source terms and BCs.
    for (unsigned int v=24; v<26; ++v) {
        const double* q=uZipVars[v]+offset;
        double* out=unzipVarsRHS[v]+offset;
        ko_deriv_x(grad_0_scalarPhi,q,hx,sz,bflag);
        ko_deriv_y(grad_1_scalarPhi,q,hy,sz,bflag);
        ko_deriv_z(grad_2_scalarPhi,q,hz,sz,bflag);
        for(unsigned int k=PW;k<nz-PW;++k)
            for(unsigned int j=PW;j<ny-PW;++j)
                for(unsigned int i=PW;i<nx-PW;++i) {
                    const auto p=i+nx*(j+ny*k);
                    out[p]+=eks::config.ko_sigma*(grad_0_scalarPhi[p]+
                         grad_1_scalarPhi[p]+grad_2_scalarPhi[p]);
                }
    }
}
if (eks::config.boundary == "analytic" && bflag) {
    // Exact time-dependent boundary RHS; no periodic-grid assumption.
    // Use RK STAGE time t, never the beginning-of-step global time.
    for(unsigned int k=PW;k<nz-PW;++k)
        for(unsigned int j=PW;j<ny-PW;++j)
            for(unsigned int i=PW;i<nx-PW;++i) {
                const bool boundary =
                  ((bflag&(1u<<OCT_DIR_LEFT)) && i==PW) ||
                  ((bflag&(1u<<OCT_DIR_RIGHT)) && i==nx-PW-1) ||
                  ((bflag&(1u<<OCT_DIR_DOWN)) && j==PW) ||
                  ((bflag&(1u<<OCT_DIR_UP)) && j==ny-PW-1) ||
                  ((bflag&(1u<<OCT_DIR_BACK)) && k==PW) ||
                  ((bflag&(1u<<OCT_DIR_FRONT)) && k==nz-PW-1);
                if(!boundary) continue;
                const auto p=i+nx*(j+ny*k);
                const std::array<double,3> x={
                  pmin[0]+i*hx-eks::config.center[0],
                  pmin[1]+j*hy-eks::config.center[1],
                  pmin[2]+k*hz-eks::config.center[2]};
                eks::ExactState q;
                if(bssn::BSSN_ID_TYPE==100) {
                    const double r=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
                    q=eks::gaussian_exact(r,t,eks::config.amplitude,eks::config.width);
                } else {
                    q=eks::plane_exact(x,t,eks::config.amplitude,eks::config.mu,
                                       eks::config.wavevector);
                }
                scalarPhi_rhs[p]=-q.pi;
                scalarPi_rhs[p]=q.pi_t;
            }
}
if (eks::config.mode == "frozen") {
    // Freeze AFTER every metric boundary/dissipation contribution.
    for(unsigned int v=0;v<24;++v)
        std::fill_n(unzipVarsRHS[v]+offset,n,0.0);
}
