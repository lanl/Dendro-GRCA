const double* eks_phi = uZipVars[bssn::U_SCALARPHI] + offset;
const double* eks_pi = uZipVars[bssn::U_SCALARPI] + offset;
std::vector<double> eks_cgrad(3 * static_cast<std::size_t>(n),0.0);
double* eks_dx=eks_cgrad.data();
double* eks_dy=eks_dx+n;
double* eks_dz=eks_dy+n;
if(eks::config.mode=="coupled") {
    deriv_x(eks_dx,eks_phi,hx,sz,bflag);
    deriv_y(eks_dy,eks_phi,hy,sz,bflag);
    deriv_z(eks_dz,eks_phi,hz,sz,bflag);
}
