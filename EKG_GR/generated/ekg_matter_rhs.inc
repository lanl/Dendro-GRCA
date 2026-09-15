// GENERATED: edit CodeGen/ekg_equations.py, not this file.
// potential=massive; gauge=puncture; backreaction=on; advection=centered
// Dendro: C++ Equation Code Generation {{{{ 
// Dendro: TEMPORARY VARIABLES
const double DENDRO_0000 = alpha[pp]*scalarPi[pp];
const double DENDRO_0001 = scalar_mu*scalar_mu;
const double DENDRO_0002 = gt4[pp]*gt4[pp];
const double DENDRO_0003 = gt3[pp]*gt5[pp];
const double DENDRO_0004 = -DENDRO_0002 + DENDRO_0003;
const double DENDRO_0005 = DENDRO_0004*grad_0_scalarPhi[pp];
const double DENDRO_0006 = gt1[pp]*gt4[pp];
const double DENDRO_0007 = DENDRO_0006 - gt2[pp]*gt3[pp];
const double DENDRO_0008 = DENDRO_0007*grad_2_scalarPhi[pp];
const double DENDRO_0009 = DENDRO_0007*grad_0_scalarPhi[pp];
const double DENDRO_0010 = gt2[pp]*gt2[pp];
const double DENDRO_0011 = -DENDRO_0010 + gt0[pp]*gt5[pp];
const double DENDRO_0012 = DENDRO_0011*grad_1_scalarPhi[pp];
const double DENDRO_0013 = gt1[pp]*gt1[pp];
const double DENDRO_0014 = -DENDRO_0013 + gt0[pp]*gt3[pp];
const double DENDRO_0015 = DENDRO_0014*grad_2_scalarPhi[pp];
const double DENDRO_0016 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
const double DENDRO_0017 = DENDRO_0016*grad_1_scalarPhi[pp];
const double DENDRO_0018 = DENDRO_0016*grad_0_scalarPhi[pp];
const double DENDRO_0019 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
const double DENDRO_0020 = DENDRO_0019*grad_2_scalarPhi[pp];
const double DENDRO_0021 = DENDRO_0019*grad_1_scalarPhi[pp];
const double DENDRO_0022 = 1.0/((DENDRO_0002*gt0[pp] - DENDRO_0003*gt0[pp] - 2*DENDRO_0006*gt2[pp] + DENDRO_0010*gt3[pp] + DENDRO_0013*gt5[pp]));
const double DENDRO_0023 = DENDRO_0022*chi[pp];
const double DENDRO_0024 = 1.0/(chi[pp]);
const double DENDRO_0025 = -DENDRO_0007*grad_0_chi[pp] - DENDRO_0014*grad_2_chi[pp] + DENDRO_0019*grad_1_chi[pp];
const double DENDRO_0026 = DENDRO_0022*DENDRO_0025;
const double DENDRO_0027 = DENDRO_0014*DENDRO_0022;
const double DENDRO_0028 = grad_1_gt5[pp] - 2*grad_2_gt4[pp];
const double DENDRO_0029 = DENDRO_0019*DENDRO_0022;
const double DENDRO_0030 = grad_0_gt5[pp] - 2*grad_2_gt2[pp];
const double DENDRO_0031 = -DENDRO_0004*grad_0_chi[pp] - DENDRO_0007*grad_2_chi[pp] + DENDRO_0016*grad_1_chi[pp];
const double DENDRO_0032 = DENDRO_0022*grad_0_scalarPhi[pp];
const double DENDRO_0033 = -DENDRO_0011*grad_1_chi[pp] + DENDRO_0016*grad_0_chi[pp] + DENDRO_0019*grad_2_chi[pp];
const double DENDRO_0034 = DENDRO_0022*grad_1_scalarPhi[pp];
const double DENDRO_0035 = DENDRO_0022*DENDRO_0033;
const double DENDRO_0036 = DENDRO_0011*DENDRO_0022;
const double DENDRO_0037 = grad_0_gt3[pp] - 2*grad_1_gt1[pp];
const double DENDRO_0038 = DENDRO_0016*DENDRO_0022;
const double DENDRO_0039 = 2*grad_1_gt4[pp] - grad_2_gt3[pp];
const double DENDRO_0040 = DENDRO_0024*gt3[pp];
const double DENDRO_0041 = DENDRO_0022*grad_2_scalarPhi[pp];
const double DENDRO_0042 = DENDRO_0022*DENDRO_0031;
const double DENDRO_0043 = DENDRO_0004*DENDRO_0022;
const double DENDRO_0044 = 2*grad_0_gt2[pp] - grad_2_gt0[pp];
const double DENDRO_0045 = DENDRO_0007*DENDRO_0022;
const double DENDRO_0046 = 2*grad_0_gt1[pp] - grad_1_gt0[pp];
const double DENDRO_0047 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_0048 = 2*DENDRO_0019;
const double DENDRO_0049 = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_0050 = 2*DENDRO_0016;
const double DENDRO_0051 = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_0052 = 2*DENDRO_0007;
const double DENDRO_0053 = scalarPi[pp]*scalarPi[pp];
const double DENDRO_0054 = DENDRO_0001*scalarPhi[pp]*scalarPhi[pp];
const double DENDRO_0055 = alpha[pp]*ekg_pi;
const double DENDRO_0056 = grad_0_scalarPhi[pp]*grad_0_scalarPhi[pp];
const double DENDRO_0057 = grad_1_scalarPhi[pp]*grad_1_scalarPhi[pp];
const double DENDRO_0058 = grad_2_scalarPhi[pp]*grad_2_scalarPhi[pp];
const double DENDRO_0059 = grad_0_scalarPhi[pp]*grad_1_scalarPhi[pp];
const double DENDRO_0060 = grad_1_scalarPhi[pp]*grad_2_scalarPhi[pp];
const double DENDRO_0061 = grad_0_scalarPhi[pp]*grad_2_scalarPhi[pp];
const double DENDRO_0062 = DENDRO_0004*DENDRO_0056 + DENDRO_0011*DENDRO_0057 + DENDRO_0014*DENDRO_0058 - DENDRO_0048*DENDRO_0060 - DENDRO_0050*DENDRO_0059 + DENDRO_0052*DENDRO_0061;
const double DENDRO_0063 = (1.0/3.0)*DENDRO_0022*DENDRO_0062;
const double DENDRO_0064 = 8*DENDRO_0055*chi[pp];
const double DENDRO_0065 = 16*DENDRO_0000*DENDRO_0022*ekg_pi;
const double DENDRO_0066 = DENDRO_0065*(DENDRO_0005 + DENDRO_0008 - DENDRO_0017);
const double DENDRO_0067 = -DENDRO_0065*(-DENDRO_0012 + DENDRO_0018 + DENDRO_0020);
const double DENDRO_0068 = DENDRO_0065*(DENDRO_0009 + DENDRO_0015 - DENDRO_0021);
const double DENDRO_0069 = -DENDRO_0023*DENDRO_0062;
// Dendro: END TEMPORARY VARIABLES

// Dendro: MAIN VARIABLES
//--
out_scalarPhi[pp] = -DENDRO_0000 + beta0[pp]*grad_0_scalarPhi[pp] + beta1[pp]*grad_1_scalarPhi[pp] + beta2[pp]*grad_2_scalarPhi[pp];

//--
out_scalarPi[pp] = DENDRO_0000*K[pp] + DENDRO_0001*alpha[pp]*scalarPhi[pp] - 1.0/2.0*DENDRO_0023*alpha[pp]*(DENDRO_0004*(DENDRO_0034*(-DENDRO_0011*DENDRO_0046 + DENDRO_0016*grad_0_gt0[pp] + DENDRO_0019*DENDRO_0044 + DENDRO_0024*DENDRO_0033*gt0[pp]) + DENDRO_0041*(-DENDRO_0007*grad_0_gt0[pp] - DENDRO_0014*DENDRO_0044 + DENDRO_0019*DENDRO_0046 + DENDRO_0024*DENDRO_0025*gt0[pp]) - 2*grad2_0_0_scalarPhi[pp] + grad_0_scalarPhi[pp]*(DENDRO_0016*DENDRO_0022*DENDRO_0046 - DENDRO_0024*(-DENDRO_0042*gt0[pp] + 2*grad_0_chi[pp]) - DENDRO_0043*grad_0_gt0[pp] - DENDRO_0044*DENDRO_0045)) + DENDRO_0011*(DENDRO_0032*(DENDRO_0004*DENDRO_0037 - DENDRO_0007*DENDRO_0039 + DENDRO_0016*grad_1_gt3[pp] + DENDRO_0031*DENDRO_0040) + DENDRO_0041*(DENDRO_0007*DENDRO_0037 - DENDRO_0014*DENDRO_0039 + DENDRO_0019*grad_1_gt3[pp] + DENDRO_0025*DENDRO_0040) - 2*grad2_1_1_scalarPhi[pp] + grad_1_scalarPhi[pp]*(DENDRO_0019*DENDRO_0022*DENDRO_0039 - DENDRO_0024*(-DENDRO_0035*gt3[pp] + 2*grad_1_chi[pp]) - DENDRO_0036*grad_1_gt3[pp] - DENDRO_0037*DENDRO_0038)) + DENDRO_0014*(DENDRO_0032*(DENDRO_0004*DENDRO_0030 - DENDRO_0007*grad_2_gt5[pp] - DENDRO_0016*DENDRO_0028 + DENDRO_0024*DENDRO_0031*gt5[pp]) + DENDRO_0034*(DENDRO_0011*DENDRO_0028 - DENDRO_0016*DENDRO_0030 + DENDRO_0019*grad_2_gt5[pp] + DENDRO_0024*DENDRO_0033*gt5[pp]) - 2*grad2_2_2_scalarPhi[pp] + grad_2_scalarPhi[pp]*(DENDRO_0007*DENDRO_0022*DENDRO_0030 - DENDRO_0024*(-DENDRO_0026*gt5[pp] + 2*grad_2_chi[pp]) - DENDRO_0027*grad_2_gt5[pp] - DENDRO_0028*DENDRO_0029)) - DENDRO_0048*(DENDRO_0032*(-DENDRO_0004*DENDRO_0047 - DENDRO_0007*grad_1_gt5[pp] + DENDRO_0016*grad_2_gt3[pp] + DENDRO_0024*DENDRO_0031*gt4[pp]) - 2*grad2_1_2_scalarPhi[pp] + grad_1_scalarPhi[pp]*(-DENDRO_0024*(-DENDRO_0035*gt4[pp] + grad_2_chi[pp]) + DENDRO_0029*grad_1_gt5[pp] - DENDRO_0036*grad_2_gt3[pp] + DENDRO_0038*DENDRO_0047) + grad_2_scalarPhi[pp]*(DENDRO_0019*DENDRO_0022*grad_2_gt3[pp] - DENDRO_0024*(-DENDRO_0026*gt4[pp] + grad_1_chi[pp]) - DENDRO_0027*grad_1_gt5[pp] - DENDRO_0045*DENDRO_0047)) - DENDRO_0050*(DENDRO_0041*(-DENDRO_0007*grad_1_gt0[pp] - DENDRO_0014*DENDRO_0049 + DENDRO_0019*grad_0_gt3[pp] + DENDRO_0024*DENDRO_0025*gt1[pp]) - 2*grad2_0_1_scalarPhi[pp] + grad_0_scalarPhi[pp]*(DENDRO_0016*DENDRO_0022*grad_0_gt3[pp] - DENDRO_0024*(-DENDRO_0042*gt1[pp] + grad_1_chi[pp]) - DENDRO_0043*grad_1_gt0[pp] - DENDRO_0045*DENDRO_0049) + grad_1_scalarPhi[pp]*(DENDRO_0016*DENDRO_0022*grad_1_gt0[pp] + DENDRO_0019*DENDRO_0022*DENDRO_0049 - DENDRO_0024*(-DENDRO_0035*gt1[pp] + grad_0_chi[pp]) - DENDRO_0036*grad_0_gt3[pp])) + DENDRO_0052*(DENDRO_0034*(-DENDRO_0011*DENDRO_0051 + DENDRO_0016*grad_2_gt0[pp] + DENDRO_0019*grad_0_gt5[pp] + DENDRO_0024*DENDRO_0033*gt2[pp]) - 2*grad2_0_2_scalarPhi[pp] + grad_0_scalarPhi[pp]*(DENDRO_0016*DENDRO_0022*DENDRO_0051 - DENDRO_0024*(-DENDRO_0042*gt2[pp] + grad_2_chi[pp]) - DENDRO_0043*grad_2_gt0[pp] - DENDRO_0045*grad_0_gt5[pp]) + grad_2_scalarPhi[pp]*(DENDRO_0019*DENDRO_0022*DENDRO_0051 - DENDRO_0024*(-DENDRO_0026*gt2[pp] + grad_0_chi[pp]) - DENDRO_0027*grad_0_gt5[pp] - DENDRO_0045*grad_2_gt0[pp]))) + DENDRO_0023*(DENDRO_0005*grad_0_alpha[pp] + DENDRO_0008*grad_0_alpha[pp] + DENDRO_0009*grad_2_alpha[pp] + DENDRO_0012*grad_1_alpha[pp] + DENDRO_0015*grad_2_alpha[pp] - DENDRO_0017*grad_0_alpha[pp] - DENDRO_0018*grad_1_alpha[pp] - DENDRO_0020*grad_1_alpha[pp] - DENDRO_0021*grad_2_alpha[pp]) + beta0[pp]*grad_0_scalarPi[pp] + beta1[pp]*grad_1_scalarPi[pp] + beta2[pp]*grad_2_scalarPi[pp];

//--
const double ekg_srcK = -4*DENDRO_0055*(-2*DENDRO_0053 + DENDRO_0054);

//--
const double ekg_srcAt0 = -DENDRO_0064*(DENDRO_0056 + DENDRO_0063*gt0[pp]);

//--
const double ekg_srcAt1 = -DENDRO_0064*(DENDRO_0059 + DENDRO_0063*gt1[pp]);

//--
const double ekg_srcAt2 = -DENDRO_0064*(DENDRO_0061 + DENDRO_0063*gt2[pp]);

//--
const double ekg_srcAt3 = -DENDRO_0064*(DENDRO_0057 + DENDRO_0063*gt3[pp]);

//--
const double ekg_srcAt4 = -DENDRO_0064*(DENDRO_0060 + DENDRO_0063*gt4[pp]);

//--
const double ekg_srcAt5 = -DENDRO_0064*(DENDRO_0058 + DENDRO_0063*gt5[pp]);

//--
const double ekg_srcGt0 = DENDRO_0066;

//--
const double ekg_srcGt1 = DENDRO_0067;

//--
const double ekg_srcGt2 = DENDRO_0068;

//--
const double ekg_srcB0 = DENDRO_0066;

//--
const double ekg_srcB1 = DENDRO_0067;

//--
const double ekg_srcB2 = DENDRO_0068;

//--
const double ekg_rho = (1.0/2.0)*DENDRO_0053 + (1.0/2.0)*DENDRO_0054 + (1.0/2.0)*DENDRO_0069;

//--
const double ekg_S = (1.0/2.0)*DENDRO_0023*DENDRO_0062 + (3.0/2.0)*DENDRO_0053 - 3.0/2.0*DENDRO_0054;

//--
const double ekg_V = (1.0/2.0)*DENDRO_0054;

//--
const double ekg_gradphi2 = DENDRO_0069;

//--
const double ekg_Scov0 = grad_0_scalarPhi[pp]*scalarPi[pp];

//--
const double ekg_Scov1 = grad_1_scalarPhi[pp]*scalarPi[pp];

//--
const double ekg_Scov2 = grad_2_scalarPhi[pp]*scalarPi[pp];
// Dendro: END MAIN VARIABLES

// Dendro: INFORMATION
// Dendro: number of original operations: 4684 
// Dendro: number of reduced operations: 432 
// Dendro: preprocessing reduced the number of operations by 4252
// Dendro: a 90.77711%% reduction
// Dendro: }}}} End Code Generation 
