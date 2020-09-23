//
// Created by milinda on 10/5/16.
//

/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 *
 * Contains the basis functions to perform inteplotation of the child bh_init value form the parent octants.
 *
 * */

#ifndef SFCSORTBENCH_INTERPOLATE_H
#define SFCSORTBENCH_INTERPOLATE_H

#include "TreeNode.h"
#include "bh.h"

/**
 * @author: Milinda Fernando.
 * @breif: This is a simple linear interpolation fucntion to decide whether we should split a corresponding octant or not according to the bh data.
 * */
template<typename T>
bool isSplit_linear(T octant, double tolerance);

template<typename T>
bool isSplit_linear(T octant, double tolerance)
{


    unsigned int myX,myY,myZ,myLev,mySz;
    myLev=octant.getLevel();
    bool state=false;
    if(myLev==m_uiMaxDepth)
        return false; // since we have reached the max depth allowed we cannot split this anymore.


    myX=octant.getX();
    myY=octant.getY();
    myZ=octant.getZ();
    mySz=1u<<(m_uiMaxDepth-myLev);

    unsigned int num_children=(1u<<m_uiDim);

    double* u=new double[num_children*NUM_FUNC_VALUES];
    Point p[num_children];
    if(m_uiDim==2) {
        p[0] = Point(myX, myY, myZ);
        p[1] = Point(myX+mySz, myY, myZ);
        p[2] = Point(myX, myY+mySz, myZ);
        p[3] = Point(myX+mySz, myY+mySz, myZ);
    }else if(m_uiDim==3)
    {
        p[0] = Point(myX, myY, myZ);
        p[1] = Point(myX+mySz, myY, myZ);
        p[2] = Point(myX, myY+mySz, myZ);
        p[3] = Point(myX+mySz, myY+mySz, myZ);

        p[4] = Point(myX, myY, myZ+mySz);
        p[5] = Point(myX+mySz, myY, myZ+mySz);
        p[6] = Point(myX, myY+mySz, myZ+mySz);
        p[7] = Point(myX+mySz, myY+mySz, myZ+mySz);

    }

    for(unsigned int i=0;i<num_children;i++)
        init_data_puncture(p[i],(u+i*NUM_FUNC_VALUES));

    /*for(unsigned int j=0;j<num_children;j++)
        for(unsigned int i=0;i<NUM_FUNC_VALUES;i++)
          std::cout<<i<<" func value  of child : "<<j<<" : "<<u[j*NUM_FUNC_VALUES+i]<<std::endl;*/


    double* u_interp=new double[NUM_FUNC_VALUES];
    double* u_computed=new double[NUM_FUNC_VALUES];
    double sum=0;
    for(unsigned int i=0;i<NUM_FUNC_VALUES;i++)
    {
        for(unsigned int j=0;j<num_children;j++)
            sum+=u[j*NUM_FUNC_VALUES+i];
        u_interp[i]=sum/(double)(num_children);
        sum=0;
    }

    // compute the function values at the mid point.
    //T tmp_oct = T(myX+(mySz>>1),myY+(mySz>>1),myZ+(mySz>>1),(myLev+1),m_uiDim,m_uiMaxDepth);
    Point tmpPnt(myX+(mySz>>1),myY+(mySz>>1),myZ+(mySz>>1));
    //std::cout<<" computing refined: "<<tmpPnt.x()<<" "<<tmpPnt.y()<<" "<<tmpPnt.z()<<std::endl;
    init_data_puncture(tmpPnt,u_computed);
    double waveletCoeff=0;
    for(unsigned int i=0;i<NUM_FUNC_VALUES;i++)
    {

        waveletCoeff=fabs(u_interp[i]-u_computed[i]);
        //std::cout<<"computed ["<<i<<"]: "<<u_computed[i]<<" , interpolated: "<<u_interp[i]<<" tolerance: "<<(tolerance/num_children)<<" diff: "<<waveletCoeff<<std::endl;
        if(waveletCoeff>(tolerance/num_children)) {
            //std::cout<<"Splitting Octant: "<<octant<<std::endl;
            state=true;
            break;
        }
        //std::cout<<"Error : "<< abs(u_interp[i]-u_computed[i])<<std::endl;
    }

    delete [] u;
    delete [] u_interp;
    delete [] u_computed;

    return state;



}



#endif //SFCSORTBENCH_INTERPOLATE_H
