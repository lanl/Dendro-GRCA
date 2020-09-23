//
// Created by milinda on 10/8/18.
//

#ifndef DENDRO_5_0_GRDEF_H
#define DENDRO_5_0_GRDEF_H

#define Rx (massgrav::MASSGRAV_COMPD_MAX[0]-massgrav::MASSGRAV_COMPD_MIN[0])
#define Ry (massgrav::MASSGRAV_COMPD_MAX[1]-massgrav::MASSGRAV_COMPD_MIN[1])
#define Rz (massgrav::MASSGRAV_COMPD_MAX[2]-massgrav::MASSGRAV_COMPD_MIN[2])

#define RgX (massgrav::MASSGRAV_OCTREE_MAX[0]-massgrav::MASSGRAV_OCTREE_MIN[0])
#define RgY (massgrav::MASSGRAV_OCTREE_MAX[1]-massgrav::MASSGRAV_OCTREE_MIN[1])
#define RgZ (massgrav::MASSGRAV_OCTREE_MAX[2]-massgrav::MASSGRAV_OCTREE_MIN[2])

#define GRIDX_TO_X(xg) (((Rx/RgX)*(xg-massgrav::MASSGRAV_OCTREE_MIN[0]))+massgrav::MASSGRAV_COMPD_MIN[0])
#define GRIDY_TO_Y(yg) (((Ry/RgY)*(yg-massgrav::MASSGRAV_OCTREE_MIN[1]))+massgrav::MASSGRAV_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg) (((Rz/RgZ)*(zg-massgrav::MASSGRAV_OCTREE_MIN[2]))+massgrav::MASSGRAV_COMPD_MIN[2])

#define X_TO_GRIDX(xc) (((RgX/Rx)*(xc-massgrav::MASSGRAV_COMPD_MIN[0]))+massgrav::MASSGRAV_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc) (((RgY/Ry)*(yc-massgrav::MASSGRAV_COMPD_MIN[1]))+massgrav::MASSGRAV_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc) (((RgZ/Rz)*(zc-massgrav::MASSGRAV_COMPD_MIN[2]))+massgrav::MASSGRAV_OCTREE_MIN[2])

// type of the rk method.
enum RKType{RK3,RK4,RK45};


namespace massgrav
{

    enum VAR {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,
              U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,
              U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};

    //enum VAR_PSI4 {C_PSI4_REAL, C_PSI4_IMG};
    enum VAR_CONSTRAINT {C_HAM=0, C_MOM0, C_MOM1, C_MOM2, C_PSI4_REAL, C_PSI4_IMG};

    static const char * MASSGRAV_VAR_NAMES[]={"U_ALPHA","U_CHI","U_K","U_GT0","U_GT1","U_GT2",
                                          "U_BETA0","U_BETA1","U_BETA2","U_B0","U_B1","U_B2",
                                          "U_SYMGT0","U_SYMGT1","U_SYMGT2","U_SYMGT3","U_SYMGT4","U_SYMGT5",
                                          "U_SYMAT0","U_SYMAT1","U_SYMAT2","U_SYMAT3","U_SYMAT4","U_SYMAT5"};


    static const char * MASSGRAV_CONSTRAINT_VAR_NAMES[]={"C_HAM","C_MOM0","C_MOM1","C_MOM2","C_PSI4_REAL","C_PSI4_IMG"};

    /**
     * @brief Refinement mode types. 
     * WAMR : Wavelet based refinement. 
     * EH : black hole event horizon based refinement. 
     * EH_WAMR: both even horizon as well as WAMR based refinement. 
     */
    enum RefinementMode{WAMR=0, EH, EH_WAMR};



}// end of namespace massgrav


#endif //DENDRO_5_0_GRDEF_H