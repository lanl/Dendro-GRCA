#include "physcon.h"
#include "gr.h"


//Define massive term for dRGT gravity
//TODO : fix it as constexpr and as parameter
#define M_dRGT 0

using namespace massgrav;


/*----------------------------------------------------------------------
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void physical_constraints(double **uZipConVars, const double **uZipVars,
                       const unsigned int& offset,
                       const double *pmin, const double *pmax,
                       const unsigned int *sz, const unsigned int& bflag)
{

  const unsigned int nx = sz[0];
  const unsigned int ny = sz[1];
  const unsigned int nz = sz[2];
  const unsigned int n = nx * ny * nz;

  double hx = (pmax[0] - pmin[0]) / (nx - 1);
  double hy = (pmax[1] - pmin[1]) / (ny - 1);
  double hz = (pmax[2] - pmin[2]) / (nz - 1);

  double *ham = &uZipConVars[VAR_CONSTRAINT::C_HAM][offset];
  double *mom0 = &uZipConVars[VAR_CONSTRAINT::C_MOM0][offset];
  double *mom1 = &uZipConVars[VAR_CONSTRAINT::C_MOM1][offset];
  double *mom2 = &uZipConVars[VAR_CONSTRAINT::C_MOM2][offset];
  double *psi4_real = &uZipConVars[VAR_CONSTRAINT::C_PSI4_REAL][offset];
  double *psi4_img = &uZipConVars[VAR_CONSTRAINT::C_PSI4_IMG][offset];

  const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
  const double *chi = &uZipVars[VAR::U_CHI][offset];
  const double *K = &uZipVars[VAR::U_K][offset];
  const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
  const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
  const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
  const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
  const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
  const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
  const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
  const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
  const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
  const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
  const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
  const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
  const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
  const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
  const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
  const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
  const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
  const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
  const double *B0 = &uZipVars[VAR::U_B0][offset];
  const double *B1 = &uZipVars[VAR::U_B1][offset];
  const double *B2 = &uZipVars[VAR::U_B2][offset];


#include "constraint_memalloc.h"
#include "constraint_derivs.h"

// enforce hamiltonian and momentum constraints
  for (unsigned int k = 3; k < nz - 3; k++) {
    double z = pmin[2] + k*hz;
    for (unsigned int j = 3; j < ny - 3; j++) {
      double y = pmin[1] + j*hy;
      for (unsigned int i = 3; i < nx - 3; i++) {
        double x = pmin[0] + i*hx;
        unsigned int pp = i + nx * (j + ny * k);
#include "physconeqs.cpp"

      }
    }
  }
#if 0 
// We don't need it but keep until we have confident to delete it
// computes the psi4 function.
for (unsigned int k = 3; k < nz - 3; k++) {
  double z = pmin[2] + k*hz;
  for (unsigned int j = 3; j < ny - 3; j++) {
    double y = pmin[1] + j*hy;
    for (unsigned int i = 3; i < nx - 3; i++) {
      double x = pmin[0] + i*hx;
      unsigned int pp = i + nx * (j + ny * k);
#include "psi4eqs.cpp"

    }
  }
}
#endif


#include "constraint_dealloc.h"

}
