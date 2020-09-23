#include "enutsOp.h"
namespace ts
{
    
    ENUTSOp::ENUTSOp(ETSType type)
    {
        m_uiType = type;

        if(type == ETSType::RK3)
        {
            m_uiNumStages = 3;
            const unsigned int d1 = m_uiNumStages;
            const unsigned int d2 = (m_uiNumStages*m_uiNumStages);


            const DendroScalar ETS_U[] =        { 0.0     ,0.0      , 0.0, 
                                                  1.0     ,0.0      , 0.0,
                                                  1.0/4.0 ,1.0/4.0  , 0.0 };

            // coefficient matrix, result in by Taylor expansion of the stages. 
            const DendroScalar ETS_C[] =        { 1.0     , 0.0      , 0.0, 
                                                  1.0     , 1.0      , 0.0,
                                                  1.0     , 1.0/2.0  , 1.0/4.0 };

            const DendroScalar ETS_InvC[] =     { 1.0     , 0.0      , 0.0, 
                                                  1.0     , 1.0      , 0.0,
                                                  1.0     , 2.0      , 4.0 };


            
            m_uiAij.resize(d2);
            m_uiBij.resize(d2);
            m_uiCij.resize(d2);
            m_uiInvCij.resize(d2);
            m_uiPi.resize(d1);
            m_uiInvPi.resize(d1);
            m_uiMij.resize(d2);
            m_uiVin.resize(d1);
            m_uiVout.resize(d1);

            std::memcpy(m_uiAij.data(), ETS_U, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiCij.data(), ETS_C, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiInvCij.data(), ETS_InvC, sizeof(DendroScalar)*d2);
            

            
        
        }else if( type == ETSType::RK4)
        {
            m_uiNumStages = 4; 
            const unsigned int d1 = m_uiNumStages;
            const unsigned int d2 = (m_uiNumStages*m_uiNumStages);
            
            const DendroScalar ETS_U[] =        { 0.0   , 0.0       , 0.0   , 0.0,
                                                1.0/2.0 , 0.0       , 0.0   , 0.0,
                                                0.0     , 1.0/2.0   , 0.0   , 0.0,
                                                0.0     , 0.0       , 1.0   , 0.0 };

            // coefficient matrix, result in by Taylor expansion of the stages. 
            const DendroScalar ETS_C[] =        { 1.0     , 0.0          , 0.0,      0.0,
                                                  1.0     , 1.0/2.0      , 0.0,      0.0,
                                                  1.0     , 1.0/2.0      , 1.0/4.0,  0.0,
                                                  1.0     , 1.0          , 1.0/2.0,  1.0/4.0};


            const DendroScalar ETS_InvC[] =     { 1.0     , 0.0          , 0.0,      0.0,
                                                  -2.0    , 2.0          , 0.0,      0.0,
                                                  0.0     , -4.0         , 4.0,      0.0,
                                                  4.0     , 0.0          , -8.0,     4.0};


            
            m_uiAij.resize(d2);
            m_uiBij.resize(d2);
            m_uiCij.resize(d2);
            m_uiInvCij.resize(d2);
            m_uiPi.resize(d1);
            m_uiInvPi.resize(d1);
            m_uiMij.resize(d2);
            m_uiVin.resize(d1);
            m_uiVout.resize(d1);

            std::memcpy(m_uiAij.data(), ETS_U, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiCij.data(), ETS_C, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiInvCij.data(), ETS_InvC, sizeof(DendroScalar)*d2);

        }else if ( type == ETSType::RK5)
        {
            // these coefficients looks wrong. - Milinda. (need to fix those and enable the rk5 method. )
            return;
            /*m_uiNumStages = 5;
            const unsigned int d1 = m_uiNumStages;
            const unsigned int d2 = (m_uiNumStages*m_uiNumStages);
            
            const DendroScalar ETS_U[] =        {   0.0       ,  0.0        , 0.0       , 0.0        , 0.0,
                                                    1.0/8.0   ,  1.0/8.0    , 0.0       , 0.0        , 0.0,
                                                    0.0       , -1.0/2.0    , 1.0       , 0.0        , 0.0,
                                                    3.0/16.0  ,  0.0        , 0.0       , 9.0/16.0   , 0.0,
                                                    -3.0/7.0  ,  2.0/7.0    , 12.0/7.0  , -12.0/7.0  , 8.0/7.0};
            
            m_uiAij.resize(d2);
            m_uiBij.resize(d2);
            m_uiCij.resize(d2);
            m_uiPi.resize(d1);

            std::memcpy(m_uiAij.data(), ETS_U, sizeof(DendroScalar)*(d2));*/

        }

        std::cout<<" m_uinum stages : "<<m_uiNumStages<<std::endl;
        return;

    }

    void ENUTSOp::Pdt(DendroScalar dt, unsigned int rk_s)
    {
        
        m_uiPi[0]    = 1.0;
        m_uiInvPi[0] = 1.0;
        for(unsigned int s=2; s <= rk_s; s++)
        {
            m_uiPi[s-1] = m_uiPi[s-2]*dt;
            m_uiInvPi[s-1] = 1.0/m_uiPi[s-1];
        }
            

        
    }



    void ENUTSOp::Bdt(DendroScalar dt, unsigned int rk_s)
    {
        
        if(m_uiType == ETSType::RK3)
        {

            /*const DendroScalar ETS_B[] =        { 1.0     , 1.0      , 1.0/2.0, 
                                                  0.0     , 1.0      , 1.0,
                                                  0.0     , 0.0      , 1.0  };*/

            m_uiBij[0 * m_uiNumStages + 0 ] = 1.0;
            m_uiBij[0 * m_uiNumStages + 1 ] = 1.0 * dt;
            m_uiBij[0 * m_uiNumStages + 2 ] = (1.0/2.0) *dt*dt;

            m_uiBij[1 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[1 * m_uiNumStages + 1 ] = 1.0;
            m_uiBij[1 * m_uiNumStages + 2 ] = 1.0 * dt;


            m_uiBij[2 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 1 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 2 ] = 1.0;




        }else if(m_uiType == ETSType::RK4)
        {

            /*const DendroScalar ETS_B[] =       { 1.0     , 1.0      , 1.0/2.0 , 1.0/6.0,
                                                  0.0     , 1.0      , 1.0     , 1.0/2.0,    
                                                  0.0     , 0.0      , 1.0     , 1.0    ,
                                                  0.0     , 0.0      , 0.0     , 1.0    };*/
            
            m_uiBij[0 * m_uiNumStages + 0 ] = 1.0;
            m_uiBij[0 * m_uiNumStages + 1 ] = 1.0 * dt;
            m_uiBij[0 * m_uiNumStages + 2 ] = (1.0/2.0) * dt*dt;
            m_uiBij[0 * m_uiNumStages + 3 ] = (1.0/6.0) * dt*dt*dt;

            m_uiBij[1 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[1 * m_uiNumStages + 1 ] = 1.0;
            m_uiBij[1 * m_uiNumStages + 2 ] = 1.0 *dt;
            m_uiBij[1 * m_uiNumStages + 3 ] = (1.0/2.0)*dt*dt;


            m_uiBij[2 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 1 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 2 ] = 1.0;
            m_uiBij[2 * m_uiNumStages + 3 ] = 1.0 * dt;

            m_uiBij[3 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[3 * m_uiNumStages + 1 ] = 0.0;
            m_uiBij[3 * m_uiNumStages + 2 ] = 0.0;
            m_uiBij[3 * m_uiNumStages + 3 ] = 1.0;



        }
        
    }


    void ENUTSOp::Cfc(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, unsigned int rk_s, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof)
    {

        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];

        const unsigned int NN = nx*ny*nz;
        const unsigned int d = m_uiNumStages;


        Bdt(dt,rk_s);

        // B x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                for(unsigned int k=0; k < d; k++)
                    m_uiMij[i*d + j] += m_uiBij[i*d+k] * m_uiCij[k*d + j]; 
            }
                
        std::swap(m_uiMij,m_uiBij);

        // C x B x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                for(unsigned int k=0; k < d; k++)
                    m_uiMij[i*d + j] += m_uiCij[i*d+k] * m_uiBij[k*d + j]; 
            }
        
        // now M = C x B x C^{-1}
        
        for(unsigned int v=0; v < dof; v++)
        {
            for(unsigned int n=0; n < NN; n++)
            {
                
                Pdt(dt_f,rk_s);

                for(unsigned int s=0; s <rk_s; s++)
                    m_uiVin[s] = in[s][ v*NN +  n] * m_uiInvPi[s];
                

                for(unsigned int i=0; i < d; i++)
                {
                    m_uiVout[i]=0;
                    for(unsigned int k=0; k < d; k++)
                        m_uiVout[i]+= m_uiMij[i*d + k ] * m_uiVin[k];
                }

                Pdt(dt_c,rk_s);
                for(unsigned int s=0; s <rk_s; s++)
                    m_uiVout[s] *= m_uiPi[s];

                for(unsigned int s=0; s <rk_s; s++)
                    out[s][v*NN +  n] = m_uiVout[s];

            }

        }


    }


    void ENUTSOp::Ccf(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, unsigned int rk_s, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof)
    {

        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];

        const unsigned int NN = nx*ny*nz;
        const unsigned int d = m_uiNumStages;


        Bdt(dt,rk_s);

        // B x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                for(unsigned int k=0; k < d; k++)
                    m_uiMij[i*d + j] += m_uiBij[i*d+k] * m_uiCij[k*d + j]; 
            }
                
        std::swap(m_uiMij,m_uiBij);

        // C x B x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                for(unsigned int k=0; k < d; k++)
                    m_uiMij[i*d + j] += m_uiCij[i*d+k] * m_uiBij[k*d + j]; 
            }
        
        // now M = C x B x C^{-1}
        
        for(unsigned int v=0; v < dof; v++)
        {
            for(unsigned int n=0; n < NN; n++)
            {
                
                Pdt(dt_c,rk_s);

                for(unsigned int s=0; s <rk_s; s++)
                    m_uiVin[s] = in[s][ v*NN +  n] * m_uiInvPi[s];
                

                for(unsigned int i=0; i < d; i++)
                {
                    m_uiVout[i]=0;
                    for(unsigned int k=0; k < d; k++)
                        m_uiVout[i]+= m_uiMij[i*d + k ] * m_uiVin[k];
                }

                Pdt(dt_f,rk_s);
                for(unsigned int s=0; s <rk_s; s++)
                    m_uiVout[s] *= m_uiPi[s];

                for(unsigned int s=0; s <rk_s; s++)
                    out[s][v*NN +  n] = m_uiVout[s];

            }

        }


    }



}