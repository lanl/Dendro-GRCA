####################################################################
#
# Date : Dec.12.2017
# Updated : Jun.10.2020
# Python script that generates Psi4 for gravitational waves and
# the momentum and Hamiltonian constraint equations for massive
# gravity. 
# 
####################################################################

#!/usr/bin/env/ python3

import dendro
from sympy import *

###################################################################
# initialize
###################################################################

# Declare variables.
# These include the BSSN variables that we need for the Psi4 
# calculation. 
chi = dendro.scalar("chi","[pp]")
K   = dendro.scalar("K","[pp]")
Gt  = dendro.vec3("Gt","[pp]")
gt  = dendro.sym_3x3("gt","[pp]")
At  = dendro.sym_3x3("At","[pp]")
b = dendro.vec3("beta","[pp]")
a = dendro.scalar("alpha","[pp]")

# declare reference metric related vars
# TODO : this is not really evolution variables... but somewhat need to be defined
#a_ref   = dendro.scalar("alpha_ref", "[pp]")
#b_ref   = dendro.vec3("beta_ref", "[pp]")
#f_ref  = dendro.sym_3x3("f_ref", "[pp]")

# Alternative way (same as in evolution eqs part)
a_ref = 1
b_ref = Matrix([[0,0,0]])
f_ref = eye(3)

# Specify the operators needed for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2 = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction

# Metric related quantities, i.e. the metric and its inverse  
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

dendro.set_ref_metric(f_ref)
if_ref = dendro.get_inverse_ref_metric()

# Some constants
PI = 3.14159265358979323846

# Mass from dRGT
M_dRGT = symbols('M_dRGT')

# Christoffels, Ricci, et al  
C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails 
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)

# Compute some energy momentum tensor related values
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
AMat_10 = Matrix([a*a*sum([igt[i,l]*b_ref[l] for l in dendro.e_i]) - b[i]*(a_ref*a_ref+sum([b[l]* b_ref[l] for l in dendro.e_i])) for i in dendro.e_i])
AMat_11 = Matrix([a*a*sum([igt[i,l]*f_ref[l,j] for l in dendro.e_i]) - b[i] * (sum([b[l]*f_ref[l, j] for l in dendro.e_i]) - b_ref[j]) for i,j in dendro.e_ij])
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
                         b[k]*(a_ref*a_ref+sum([b[l]* b_ref[l] for l in dendro.e_i]))) for k in   dendro.e_i]) for i in dendro.e_i])
AMatsq_11 = Matrix([(a*a*sum([igt[i,l]*b_ref[l] for l in dendro.e_i]) - \
                    b[i]*(a_ref*a_ref+sum([b[l]* b_ref[l] for l in dendro.e_i]))) * \
                    (-b_ref[j] + sum([b[l]*f_ref[l,j] for l in dendro.e_i])) \
                    for i,j in dendro.e_ij]) + \
            Matrix([sum([(a*a*sum([igt[i,l]*f_ref[l,k] for l in dendro.e_i]) - \
                    b[i] * (sum([b[l]*f_ref[l, k] for l in dendro.e_i]) - b_ref[k])) * (a*a*      sum([igt[k,l]*f_ref[l,j] for l in dendro.e_i]) - \
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
x_var = a_ref*a_ref + sum([sum([(b_ref[i]*b_ref[j]-b[i]*b[j])*if_ref[i,j] for i in dendro.e_i])   for j in dendro.e_i])

# General potential computations
V_alpha = 2*M_dRGT_sq*(sum([sum([b[i]*b[j]*f_ref[i,j] for i in dendro.e_i]) for j in dendro.      e_i]) - \
          a_ref*a_ref - 2*sum([b[i]*b_ref[i] for i in dendro.e_i]))/(a*a) + \
          2*M_dRGT_sq*(Tr_sqrtA-3)
V_beta_i = Matrix([2*M_dRGT_sq*sum([b[j]*f_ref[j,i] for j in dendro.e_i])/sqrt(x_var) for i in    dendro.e_i])
#define energy momentum quantities
rhot = -V_alpha/2
Sit = -V_beta_i

###################################################################
# Calculate the tetrad used in the Psi4 calculation
####################################################################

# Define coordinates
x, y, z = symbols('x, y, z')

# Some other values
invsqrt2 = 0.7071067811865475244
inv_chi = 1/chi

# Define the original spatial vectors in our tetrad 
r_vec = Matrix([[x,y,z]])
theta = Matrix([[x*z,y*z,-(x*x+y*y)]])
phi = Matrix([[-y,x,0.0]])

# We use Gram-Schmidt to make the basis orthonormal. 
# Note that we use the original (not conformally rescaled) metric to define
# the tetrad and correspondingly Psi4.  
gd = gt*inv_chi

# For r_vec
inner_product = 0.0
inner_product = sum([sum([gd[i,j] * r_vec[i] * r_vec[j] for i in dendro.e_i]) for j in dendro.e_i])

r_vec /= sqrt(inner_product)

# For theta
inner_product_1 = 0.0
inner_product_2 = 0.0

inner_product_1 = sum([sum([gd[i,j] * theta[i] * theta[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * theta[j] for i in dendro.e_i]) for j in dendro.e_i])

theta -= inner_product_2 * r_vec
theta /= sqrt(inner_product_1 - inner_product_2 * inner_product_2)

# For phi
inner_product_1 = 0.0
inner_product_2 = 0.0
inner_product_3 = 0.0

inner_product_1 = sum([sum([gd[i,j] * phi[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_3 = sum([sum([gd[i,j] * theta[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])

phi -= inner_product_2 * r_vec + inner_product_3 * theta
phi /= sqrt(inner_product_1 - inner_product_2 * inner_product_2 - inner_product_3 * inner_product_3)

# This completes the tetrad construction. 

###################################################################
# Calculate the Weyl scalar, Psi4, for graviational wave extraction
###################################################################

# Rename the tetrad quantities for calculating Psi4 
r_np = Matrix([[r_vec[0],r_vec[1],r_vec[2]]])
m_np_real = Matrix([[theta[0],theta[1],theta[2]]])*invsqrt2
m_np_img = Matrix([[phi[0],phi[1],phi[2]]])*invsqrt2

# Some auxilary variables
# MM and NN are symmetric 2nd rank objects and
# MR and NR are anti-symmetric 2nd rank objects 

MM = Matrix([m_np_real[i]*m_np_real[j] - m_np_img[i]*m_np_img[j] for i,j in dendro.e_ij]) 
MM = MM.reshape(3,3)
NN = Matrix([m_np_real[i]*m_np_img[j] + m_np_real[j]*m_np_img[i] for i,j in dendro.e_ij])
NN = NN.reshape(3,3)  
MR = Matrix([m_np_real[i]*r_np[j] - m_np_real[j]*r_np[i] for i,j in dendro.e_ij])
MR = MR.reshape(3,3) 
NR = Matrix([m_np_img[i]*r_np[j] - m_np_img[j]*r_np[i] for i,j in dendro.e_ij])
NR = NR.reshape(3,3)  

# Additional intermediate variables
#A_vec = Matrix([[sum([At[j,0]*r_np[j] for j in dendro.e_i]), sum([At[j,1]*r_np[j] for j in dendro.e_i]),sum([At[j,2]*r_np[j] for j in dendro.e_i])]])
A_vec = [ sum([At[i,j]*r_np[j] for j in dendro.e_i]) for i in dendro.e_i ] 

Uu = Matrix([sum([m_np_real[k] * (d(j, At[k,i]) + sum([C2[m,k,i] * At[m,j] for m in dendro.e_i])) for k in dendro.e_i]) for i,j in dendro.e_ij])
Uu = Uu.reshape(3,3)
Vv = Matrix([sum([m_np_img[k] * (d(j, At[k,i]) + sum([C2[m,k,i] * At[m,j] for m in dendro.e_i])) for k in dendro.e_i]) for i,j in dendro.e_ij])
Vv = Vv.reshape(3,3)

r_d_chi = sum([r_np[i] * d(i, chi) for i in dendro.e_i]) 

A_temp = inv_chi * inv_chi * ( sum([A_vec[i] * r_np[i] for i in dendro.e_i]) + K * chi/3 + 0.5 * r_d_chi ) 

m_real_d_chi = sum([m_np_real[i] * d(i, chi) for i in dendro.e_i])  
m_img_d_chi  = sum([m_np_img [i] * d(i, chi) for i in dendro.e_i]) 

m_real_A_vec = sum([m_np_real[i] * A_vec[i] for i in dendro.e_i]) 
m_img_A_vec  = sum([m_np_img [i] * A_vec[i] for i in dendro.e_i])  


# Calculate Psi4

psi4_1_real = sum([R[i,i] * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([R[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 
psi4_1_img  = sum([R[i,i] * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([R[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 

psi4_2_real = A_temp * (sum([At[i,i] * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]))
psi4_2_img = A_temp * (sum([At[i,i] * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]))  

psi4_3_real = inv_chi * sum([sum([MR[i,j]* Uu[i,j] - NR[i,j]*Vv[i,j] for i in dendro.e_i]) for j in dendro.e_i])  
psi4_3_img = inv_chi * sum([sum([NR[i,j]* Uu[i,j] + MR[i,j]*Vv[i,j] for i in dendro.e_i]) for j in dendro.e_i]) 

psi4_4_real = inv_chi * inv_chi * (m_real_A_vec * (m_real_A_vec + 0.5 * m_real_d_chi) - m_img_A_vec * (m_img_A_vec + 0.5 * m_img_d_chi))  
psi4_4_img = inv_chi * inv_chi * (m_real_A_vec * (m_img_A_vec - 0.5 * m_img_d_chi ) + m_img_A_vec * (m_real_A_vec - 0.5 * m_real_d_chi))  

# Adding previous auxilary Psi4 calculations

psi4_real =     psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real
psi4_img  = - ( psi4_1_img  + psi4_2_img  - psi4_3_img  - psi4_4_img  )

###################################################################
# Constraint Equations
###################################################################

# The Hamiltonian constraint
ham = sum(chi*igt[j,k]*R[j,k] for j,k in dendro.e_ij) - dendro.sqr(At) + Rational(2,3)*K**2 - 2*rhot

# The momentum  constraints 
mom = Matrix([sum([igt[j,k]*(  d(k,At[i,j]) - \
              sum(dendro.C2[m,k,i]*At[j,m] for m in dendro.e_i)) \
                  for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
      Matrix([sum([Gt[j]*At[i,j] for j in dendro.e_i]) for i in dendro.e_i]) -\
      Rational(3,2)*Matrix([ \
            sum([igt[j,k]*At[k,i]*d(j,chi)/chi for j,k in dendro.e_ij])  \
            for i in dendro.e_i]) -\
      Rational(2,3)*Matrix([d(i,K) for i in dendro.e_i]) - \
    Matrix([Sit[i] for i in dendro.e_i])
mom = [item for sublist in mom.tolist() for item in sublist]

# Output for this should be included psi4_real and psi4_img as double precision  
###################################################################
# generate code
###################################################################

outs = [psi4_real, psi4_img, ham, mom]
vnames = ['psi4_real', 'psi4_img', 'ham', 'mom']
dendro.generate(outs, vnames, '[pp]')
