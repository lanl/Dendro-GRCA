// one planned grad_set per variable (hand-generated from bssnrhs_derivs.h);
// interior blocks only, puncture blocks keep bssnrhs_derivs.h
using DD = dendroderivs::DendroDerivatives;
{
    DD::DerivSet o;
    o.x = grad_0_alpha;
    o.y = grad_1_alpha;
    o.z = grad_2_alpha;
    o.xx = grad2_0_0_alpha;
    o.yy = grad2_1_1_alpha;
    o.zz = grad2_2_2_alpha;
    o.xy = grad2_0_1_alpha;
    o.xz = grad2_0_2_alpha;
    o.yz = grad2_1_2_alpha;
    bssn::active_derivs()->grad_set(o, alpha, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_beta0;
    o.y = grad_1_beta0;
    o.z = grad_2_beta0;
    o.xx = grad2_0_0_beta0;
    o.yy = grad2_1_1_beta0;
    o.zz = grad2_2_2_beta0;
    o.xy = grad2_0_1_beta0;
    o.xz = grad2_0_2_beta0;
    o.yz = grad2_1_2_beta0;
    bssn::active_derivs()->grad_set(o, beta0, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_beta1;
    o.y = grad_1_beta1;
    o.z = grad_2_beta1;
    o.xx = grad2_0_0_beta1;
    o.yy = grad2_1_1_beta1;
    o.zz = grad2_2_2_beta1;
    o.xy = grad2_0_1_beta1;
    o.xz = grad2_0_2_beta1;
    o.yz = grad2_1_2_beta1;
    bssn::active_derivs()->grad_set(o, beta1, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_beta2;
    o.y = grad_1_beta2;
    o.z = grad_2_beta2;
    o.xx = grad2_0_0_beta2;
    o.yy = grad2_1_1_beta2;
    o.zz = grad2_2_2_beta2;
    o.xy = grad2_0_1_beta2;
    o.xz = grad2_0_2_beta2;
    o.yz = grad2_1_2_beta2;
    bssn::active_derivs()->grad_set(o, beta2, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_B0;
    o.y = grad_1_B0;
    o.z = grad_2_B0;
    bssn::active_derivs()->grad_set(o, B0, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_B1;
    o.y = grad_1_B1;
    o.z = grad_2_B1;
    bssn::active_derivs()->grad_set(o, B1, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_B2;
    o.y = grad_1_B2;
    o.z = grad_2_B2;
    bssn::active_derivs()->grad_set(o, B2, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_chi;
    o.y = grad_1_chi;
    o.z = grad_2_chi;
    o.xx = grad2_0_0_chi;
    o.yy = grad2_1_1_chi;
    o.zz = grad2_2_2_chi;
    o.xy = grad2_0_1_chi;
    o.xz = grad2_0_2_chi;
    o.yz = grad2_1_2_chi;
    bssn::active_derivs()->grad_set(o, chi, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_Gt0;
    o.y = grad_1_Gt0;
    o.z = grad_2_Gt0;
    bssn::active_derivs()->grad_set(o, Gt0, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_Gt1;
    o.y = grad_1_Gt1;
    o.z = grad_2_Gt1;
    bssn::active_derivs()->grad_set(o, Gt1, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_Gt2;
    o.y = grad_1_Gt2;
    o.z = grad_2_Gt2;
    bssn::active_derivs()->grad_set(o, Gt2, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_K;
    o.y = grad_1_K;
    o.z = grad_2_K;
    bssn::active_derivs()->grad_set(o, K, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_gt0;
    o.y = grad_1_gt0;
    o.z = grad_2_gt0;
    o.xx = grad2_0_0_gt0;
    o.yy = grad2_1_1_gt0;
    o.zz = grad2_2_2_gt0;
    o.xy = grad2_0_1_gt0;
    o.xz = grad2_0_2_gt0;
    o.yz = grad2_1_2_gt0;
    bssn::active_derivs()->grad_set(o, gt0, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_gt1;
    o.y = grad_1_gt1;
    o.z = grad_2_gt1;
    o.xx = grad2_0_0_gt1;
    o.yy = grad2_1_1_gt1;
    o.zz = grad2_2_2_gt1;
    o.xy = grad2_0_1_gt1;
    o.xz = grad2_0_2_gt1;
    o.yz = grad2_1_2_gt1;
    bssn::active_derivs()->grad_set(o, gt1, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_gt2;
    o.y = grad_1_gt2;
    o.z = grad_2_gt2;
    o.xx = grad2_0_0_gt2;
    o.yy = grad2_1_1_gt2;
    o.zz = grad2_2_2_gt2;
    o.xy = grad2_0_1_gt2;
    o.xz = grad2_0_2_gt2;
    o.yz = grad2_1_2_gt2;
    bssn::active_derivs()->grad_set(o, gt2, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_gt3;
    o.y = grad_1_gt3;
    o.z = grad_2_gt3;
    o.xx = grad2_0_0_gt3;
    o.yy = grad2_1_1_gt3;
    o.zz = grad2_2_2_gt3;
    o.xy = grad2_0_1_gt3;
    o.xz = grad2_0_2_gt3;
    o.yz = grad2_1_2_gt3;
    bssn::active_derivs()->grad_set(o, gt3, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_gt4;
    o.y = grad_1_gt4;
    o.z = grad_2_gt4;
    o.xx = grad2_0_0_gt4;
    o.yy = grad2_1_1_gt4;
    o.zz = grad2_2_2_gt4;
    o.xy = grad2_0_1_gt4;
    o.xz = grad2_0_2_gt4;
    o.yz = grad2_1_2_gt4;
    bssn::active_derivs()->grad_set(o, gt4, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_gt5;
    o.y = grad_1_gt5;
    o.z = grad_2_gt5;
    o.xx = grad2_0_0_gt5;
    o.yy = grad2_1_1_gt5;
    o.zz = grad2_2_2_gt5;
    o.xy = grad2_0_1_gt5;
    o.xz = grad2_0_2_gt5;
    o.yz = grad2_1_2_gt5;
    bssn::active_derivs()->grad_set(o, gt5, DD::DM_X | DD::DM_Y | DD::DM_Z | DD::DM_XX | DD::DM_YY | DD::DM_ZZ | DD::DM_XY | DD::DM_XZ | DD::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_At0;
    o.y = grad_1_At0;
    o.z = grad_2_At0;
    bssn::active_derivs()->grad_set(o, At0, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_At1;
    o.y = grad_1_At1;
    o.z = grad_2_At1;
    bssn::active_derivs()->grad_set(o, At1, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_At2;
    o.y = grad_1_At2;
    o.z = grad_2_At2;
    bssn::active_derivs()->grad_set(o, At2, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_At3;
    o.y = grad_1_At3;
    o.z = grad_2_At3;
    bssn::active_derivs()->grad_set(o, At3, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_At4;
    o.y = grad_1_At4;
    o.z = grad_2_At4;
    bssn::active_derivs()->grad_set(o, At4, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
{
    DD::DerivSet o;
    o.x = grad_0_At5;
    o.y = grad_1_At5;
    o.z = grad_2_At5;
    bssn::active_derivs()->grad_set(o, At5, DD::DM_X | DD::DM_Y | DD::DM_Z, hx, hy, hz, sz, bflag);
}
