// Created by milinda on 9/29/16.
/**
 * @author Milinda Fernando.
 * School of Computing University of Utah.
 * @breif Contains the data structure to store blackhole properties.
 * */
#ifndef SFCSORTBENCH_BH_H
#define SFCSORTBENCH_BH_H

#include "TreeNode.h"
#include <math.h>
#include "point.h"


/**
 *

@brief : Data to initialize the black holes.

mass1 :=  0.4824
bh1x  :=  6.0
bh1y  :=  0.0
bh1z  :=  0.00123
vx1   :=  -0.02088
vy1   :=  -0.512349
vz1   :=   0.0
spin1 :=   0.0
spin1_th := 0.0
spin1_ph := 0.0

mass2 := 0.4824
bh2x  := -6.0
bh2y  :=  0.0
bh2z  := -0.00123
vx2   :=  0.02088
vy2   :=  0.512349
vz2   :=  0.0
spin2 :=  0.0
spin2_th := 0.0
spin2_ph := 0.0
 * */


/**
 *
 *
 mass1 = 0.0086895
 bh1x =  4.9526
 bh1y =  0
 bh1z =  sqrt(3)*1.0e-6
 vp1[0] = -0.000010265 / mass1
 vp1[1] =  0.0067226 / mass1
 vp1[2] = 0.0
 spin1 = 0.0;
 spin1_th = 0.0;
 spin1_phi = 0.0;

 mass2 =  0.98962
 bh2x =  -0.047437
 bh2y =  0
 bh2z =  sqrt(3)*1.0e-6
 vp2[0] = -0.000010265 / mass2
 vp2[1] = -0.0067226 / mass2
 vp2[2] = 0.0
 spin2 = 0.0;
 spin2_th = 0.0;
 spin2_phi = 0.0;

 * */




/*
 *
 * Here are the preprocessor derectives for the 24 function values calculated in the ini_puncture function.
 *
 * */

#define X_MIN -400.0
#define X_MAX  400.0
#define RANGE1 (X_MAX-X_MIN)
#define RANGE2 1024.0

#define NUM_FUNC_VALUES 24


#define U_ALPHA 0
#define U_CHI 1
#define U_TRK 2
#define U_SHIFTX 3
#define U_SHIFTY 4
#define U_SHIFTZ 5

#define U_GAMTX 6
#define U_GAMTY 7
#define U_GAMTZ 8

#define U_GBX 9
#define U_GBY 10
#define U_GBZ 11

#define U_GTXX 12
#define U_GTXY 13
#define U_GTXZ 14
#define U_GTYY 15
#define U_GTYZ 16
#define U_GTZZ 17

#define U_ATXX 18
#define U_ATXY 19
#define U_ATXZ 20
#define U_ATYY 21
#define U_ATYZ 22
#define U_ATZZ 23


struct BH
{

    double mass;

    double bhx;
    double bhy;
    double bhz;

    double vx;
    double vy;
    double vz;

    double spin;
    double spin_th;
    double spin_phi;

    BH(double pmass, double pbhx, double pbhy, double pbhz, double pvx,double pvy, double pvz, double pspin, double pspin_th, double pspin_ph) // Note that parameters should be passed in the order that we have defined the properties of the black holes.
    {
        mass=pmass;

        bhx=pbhx;
        bhy=pbhy;
        bhz=pbhz;

        vx=pvx;
        vy=pvy;
        vz=pvz;

        spin=pspin;
        spin_th=pspin_th;
        spin_phi=pspin_ph;

    }


};

extern BH bh1;
extern BH bh2;


void init_data_puncture(Point pts, double *u);



#endif //SFCSORTBENCH_BH_H
