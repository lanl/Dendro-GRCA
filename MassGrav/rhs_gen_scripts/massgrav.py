####################################################################
# Jun.10.2020
# Evolution equation generator for massive gravity
#####################################################################

import dendro
import math
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# Additional parameters for damping term
R0 = symbols('MASSGRAV_ETA_R0')
ep1, ep2 = symbols('MASSGRAV_ETA_POWER[0] MASSGRAV_ETA_POWER[1]')

PI = 3.14159265358979323846

# Mass from dRGT
M_dRGT = symbols('M_dRGT')

# declare variables
a   = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")

Gt  = dendro.vec3("Gt", "[pp]")
b   = dendro.vec3("beta", "[pp]")
#B   = dendro.vec3("B", "[pp]") # Just comment out

gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

# declare reference metric related vars
# TODO : this is not really evolution variables... but somewhat need to be defined
# I comment out them.
#a_ref   = dendro.scalar("alpha_ref", "[pp]")
#b_ref   = dendro.vec3("beta_ref", "[pp]")
#f_ref  = dendro.sym_3x3("f_ref", "[pp]")

#Alternative ways to have this. Define here since they are not evolution variables
#Define flat spacetime for reference metric
a_ref = 1
b_ref = Matrix([[0,0,0]])
f_ref = eye(3)

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')

d2 = dendro.d2

#f = Function('f')

# generate metric related quantities
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

dendro.set_ref_metric(f_ref)
if_ref = dendro.get_inverse_ref_metric()


C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)

# Compute some energy momentum tensor related values
# Precomputation for sqrt(g^-1 f) matrix
# We use power serise to compute sqrt of matrix. 
# Let A = g^-1 f then we have:
#  sqrt(A) = I - |Binomial(1/2,1)| (I-A) - |Binomial(1/2,2)| (I-A)^2 - ...
#     where I is identity matrix

# Define mass square
M_dRGT_sq = M_dRGT*M_dRGT

# Define A = g^-1 f
AMat_00 = a_ref*a_ref + sum([b[l]*b_ref[l] for l in dendro.e_i])
AMat_01 = Matrix([-b_ref[j] + sum([b[l]*f_ref[l,j] for l in dendro.e_i]) for j in dendro.e_i])
AMat_10 = Matrix([a*a*sum([igt[i,l]*b_ref[l] for l in dendro.e_i]) - b[i]*(a_ref*a_ref+sum([b[l]*b_ref[l] for l in dendro.e_i])) for i in dendro.e_i])
AMat_11 = Matrix([a*a*sum([igt[i,l]*f_ref[l,j] for l in dendro.e_i]) - b[i] * (sum([b[l]*f_ref[l,j] for l in dendro.e_i]) - b_ref[j]) for i,j in dendro.e_ij])
AMat_11 = AMat_11.reshape(3,3)

# Define A^2
AMatsq_00 = (a_ref*a_ref + sum([b[l]*b_ref[l] for l in dendro.e_i])) * \
            (a_ref*a_ref + sum([b[l]*b_ref[l] for l in dendro.e_i])) + \
            sum([(-b_ref[k] + sum([b[l]*f_ref[l,k] for l in dendro.e_i])) * \
                (a*a*sum([igt[k,l]*b_ref[l] for l in dendro.e_i]) - \
                b[k]*(a_ref*a_ref+sum([b[l]*b_ref[l] for l in dendro.e_i]))) for k in dendro.e_i])
AMatsq_01 = (a_ref*a_ref + sum([b[l]*b_ref[l] for l in dendro.e_i])) * \
            Matrix([-b_ref[j] + sum([b[l]*f_ref[l,j] for l in dendro.e_i]) for j in dendro.e_i]) +\
            Matrix([sum([(-b_ref[k] + sum([b[l]*f_ref[l,k] for l in dendro.e_i])) * \
            (a*a*sum([igt[k,l]*f_ref[l,j] for l in dendro.e_i]) - \
            b[k] * (sum([b[l]*f_ref[l, j] for l in dendro.e_i]) - b_ref[k])) \
            for k in dendro.e_i]) for j in dendro.e_i])
AMatsq_10 = (a_ref*a_ref + sum([b[l]*b_ref[l] for l in dendro.e_i])) * \
            Matrix([-b_ref[i] + sum([b[l]*f_ref[l,i] for l in dendro.e_i]) for i in dendro.e_i]) +\
            Matrix([sum([(a*a*sum([igt[i,l]*f_ref[l,k] for l in dendro.e_i]) - \
                         b[i] * (sum([b[l]*f_ref[l, k] for l in dendro.e_i]) - b_ref[k]))* \
                         (a*a*sum([igt[k,l]*b_ref[l] for l in dendro.e_i]) - \
                         b[k]*(a_ref*a_ref+sum([b[l]* b_ref[l] for l in dendro.e_i]))) for k in dendro.e_i]) for i in dendro.e_i])
AMatsq_11 = Matrix([(a*a*sum([igt[i,l]*b_ref[l] for l in dendro.e_i]) - \
                    b[i]*(a_ref*a_ref+sum([b[l]* b_ref[l] for l in dendro.e_i]))) * \
                    (-b_ref[j] + sum([b[l]*f_ref[l,j] for l in dendro.e_i])) \
                    for i,j in dendro.e_ij]) + \
            Matrix([sum([(a*a*sum([igt[i,l]*f_ref[l,k] for l in dendro.e_i]) - \
                    b[i] * (sum([b[l]*f_ref[l, k] for l in dendro.e_i]) - b_ref[k])) * (a*a*sum([igt[k,l]*f_ref[l,j] for l in dendro.e_i]) - \
                    b[k] * (sum([b[l]*f_ref[l, j] for l in dendro.e_i]) - b_ref[j])) \
                    for k in dendro.e_i]) \
                    for i,j in dendro.e_ij])
AMatsq_11 = AMatsq_11.reshape(3,3)
# Define sqrt(A) from power series (upto 2nd order)
sqrtA_00 = Rational(3,4)*AMat_00 - Rational(1,8)*AMatsq_00 + Rational(3,8)
sqrtA_01 = Rational(3,4)*AMat_01 - Rational(1,8)*AMatsq_01
sqrtA_10 = Rational(3,4)*AMat_10 - Rational(1,8)*AMatsq_10
sqrtA_11 = Rational(3,4)*AMat_11 - Rational(1,8)*AMatsq_11 + Rational(3,8)*eye(3)

# Trace of sqrt(A)
Tr_sqrtA = sqrtA_00+sqrtA_11[0,0]+sqrtA_11[1,1]+sqrtA_11[2,2]

#Define other var
x_var = a_ref*a_ref + sum([sum([(b_ref[i]*b_ref[j]-b[i]*b[j])*if_ref[i,j] for i in dendro.e_i]) for j in dendro.e_i])


# General potential computations
V_alpha = 2*M_dRGT_sq*(sum([sum([b[i]*b[j]*f_ref[i,j] for i in dendro.e_i]) for j in dendro.e_i]) - \
          a_ref*a_ref - 2*sum([b[i]*b_ref[i] for i in dendro.e_i]))/(a*a) + \
          2*M_dRGT_sq*(Tr_sqrtA-3)
V_beta_i = Matrix([2*M_dRGT_sq*sum([b[j]*f_ref[j,i] for j in dendro.e_i])/sqrt(x_var) for i in dendro.e_i])
calgt = sum([sum([gt[i,j]*igt[i,j] for i in dendro.e_i]) for j in dendro.e_i])
gtjk = Matrix([sum([igt[j,l]*gt[l,k] for l in dendro.e_i]) for j,k in dendro.e_ij])
gtjk = gtjk.reshape(3,3)
V_pi_ij = Matrix([ Tr_sqrtA*igt[i,j] for i,j in dendro.e_ij]) + \
          Matrix([ sqrtA_01[i]*b_ref[j] for i,j in dendro.e_ij]) + \
          Matrix([ sum([2*sqrtA_11[i,k]*gtjk[j,k] for k in dendro.e_i]) for i,j in dendro.e_ij])

V_pi_ij = 2*M_dRGT_sq*a*sqrt(calgt)/2 * V_pi_ij

V_pi_ij = V_pi_ij.reshape(3,3)

#define energy momentum quantities
rhot = -V_alpha/2
Sit = -V_beta_i
Sijt = - V_pi_ij/a
St = sum([sum([Sijt[i,j]*igt[i,j] for i in dendro.e_i]) for j in dendro.e_i])

###################################################################
# evolution equations
###################################################################

a_rhs = l1*dendro.lie(b, a) - 2*a*K + 0*dendro.kodiss(a)


gt_rhs = dendro.lie(b, gt, weight) - 2*a*At + 0*dendro.kodiss(gt)

chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) + 0*dendro.kodiss(chi)

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*(R-8*PI*Sijt) - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) + 0*dendro.kodiss(At)

K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) + 0*dendro.kodiss(K) + 4*a*PI*(rhot+St)

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
         Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
         Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
         Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i]) - \
         Matrix([sum([dendro.inv_metric[i,j]*Sit[j] for j in dendro.e_i]) for i in dendro.e_i])
         # + kod(i,Gt[i])

Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

eta_func = R0*sqrt(sum([igt[i,j]*d(i,chi)*d(j,chi) for i,j in dendro.e_ij]))/((1-chi**ep1)**ep2)

#B_rhs = [Gt_rhs[i] - eta_func * B[i] +
#         l3 * dendro.vec_j_ad_j(b, B[i]) -
#         l4 * dendro.vec_j_ad_j(b, Gt[i]) + 0*kod(i,B[i])
#         for i in dendro.e_i]

#Compute Rti

At_rhs_aux = At_rhs.reshape(3,3)

Rti = Matrix([sum([igt[j,k]*(  d(k,At[i,j]) - \
              sum(dendro.C2[m,k,i]*At[j,m] for m in dendro.e_i)) \
                  for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
      Matrix([sum([Gt[j]*At[i,j] for j in dendro.e_i]) for i in dendro.e_i]) -\
      Rational(3,2)*Matrix([ \
            sum([igt[j,k]*At[k,i]*d(j,chi)/chi for j,k in dendro.e_ij])  \
            for i in dendro.e_i]) -\
      Rational(2,3)*Matrix([d(i,K) for i in dendro.e_i]) 
Rti= [item for sublist in Rti.tolist() for item in sublist]

#TODO : Check this term. First order system?
Rti_dt = 6*Matrix([sum([d(j,chi)*At_rhs[i,j] for j in dendro.e_i ]) for i in dendro.e_i]) \
        + 6*Matrix([d(i,chi_rhs) for i in dendro.e_i]) + \
         Rational(2,3)*Matrix([d(i,K_rhs) for i in dendro.e_i]) 
#         Matrix([ d(j,At_rhs_aux[i,j]) for i,j in dendro.e_ij])
Rti_dt = [item for sublist in Rti_dt.tolist() for item in sublist]

# Some prefactor

det_gamma = sum([sum([gt[i,j]*igt[i,j] for i in dendro.e_i]) for j in dendro.e_i])

shift_fac = 1/sqrt(4*M_dRGT_sq*M_dRGT_sq*det_gamma+sum([sum([Rti[k]*f_ref[k,l] for k in dendro.e_i]) for l in dendro.e_i]))
shift_fac_sq = shift_fac*shift_fac

# Theory dependent shift condition
b_rhs = a_ref*shift_fac*(Matrix([sum([f_ref[i,j]*Rti_dt[j] for j in dendro.e_i]) for i in dendro.e_i]) - Matrix([sum([f_ref[i,j]*Rti[j] for j in dendro.e_i])*sum([sum([Rti[m]*f_ref[m,n]*Rti_dt[n] for m in dendro.e_i]) for n in dendro.e_i])/shift_fac_sq for i in dendro.e_i]))
b_rhs = [item for sublist in b_rhs.tolist() for item in sublist]

#_I = gt*igt
#print(simplify(_I))

#_I = gt*dendro.inv_metric
#print(simplify(_I))


###
# Substitute ...
#for expr in [a_rhs, b_rhs[0], b_rhs[1], b_rhs[2], B_rhs[0], B_rhs[1], B_rhs[2], K_rhs, chi_rhs, Gt_rhs[0], Gt_rhs[1], Gt_rhs[2], gt_rhs[0], gt_rhs[0,0], gt_rhs[1,1], gt_rhs[2,2], gt_rhs[0,1], gt_rhs[0,2], gt_rhs[1,2], At_rhs[0,0], At_rhs[0,1], At_rhs[0,2], At_rhs[1,1], At_rhs[1,2], At_rhs[2,2]]:
#    for var in [a, b[0], b[1], b[2], B[0], B[1], B[2], chi, K, gt[0,0], gt[0,1], gt[0,2], gt[1,1], gt[1,2], gt[2,2], Gt[0], Gt[1], Gt[2], At[0,0], At[0,1], At[0,2], At[1,1], At[1,2], At[2,2]]:
#        expr.subs(d2(1,0,var), d2(0,1,var))
#        expr.subs(d2(2,1,var), d2(1,2,var))
#        expr.subs(d2(2,0,var), d2(0,2,var))
#
#print (a_rhs)
#print (G_rhs)


###################################################################
# generate code
###################################################################

#outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
#vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs']
#dendro.generate_debug(outs, vnames)
dendro.generate(outs, vnames, '[pp]')
#numVars=len(outs)
#for i in range(0,numVars):
#    dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')
