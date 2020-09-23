####################################################################
# May.8.2018
# Adding gamma driver into shift equation to optimize BBH behavior
# with large mass ratio
#####################################################################


import dendro
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

#QG related constants
a_const = symbols('a_const')
b_const = symbols('b_const')
qg_ho_coup = symbols('qg_ho_coup')

PI = 3.14159265358979323846
kappa = 1/(16*PI)

# Additional parameters for damping term
R0 = symbols('QUADGRAV_ETA_R0')
ep1, ep2 = symbols('QUADGRAV_ETA_POWER[0] QUADGRAV_ETA_POWER[1]')

# declare variables
a   = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")

Gt  = dendro.vec3("Gt", "[pp]")
b   = dendro.vec3("beta", "[pp]")
B   = dendro.vec3("B", "[pp]")

gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

Rsc = dendro.scalar("Rsc", "[pp]")
Rsch = dendro.scalar("Rsch", "[pp]")
Rtt = dendro.sym_3x3("Rtt", "[pp]")
Vat = dendro.sym_3x3("Vat", "[pp]")

Qsc = dendro.scalar("sc", "[pp]")
Yabnp = dendro.sym_3x3("Yabp", "[pp]")
Yabp = dendro.sym_3x3("Yabnp", "[pp]")

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

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)
###################################################################
# evolution equations
###################################################################

a_rhs = l1*dendro.lie(b, a) - 2*a*K + 0*dendro.kodiss(a)

b_rhs = [ S(3)/4 * (lf0 + lf1*a) * B[i] +
        l2 * dendro.vec_j_ad_j(b, b[i])
         for i in dendro.e_i ] + 0*dendro.kodiss(b)

gt_rhs = dendro.lie(b, gt, weight) - 2*a*At + 0*dendro.kodiss(gt)

chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) + 0*dendro.kodiss(chi)

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) + 0*dendro.kodiss(At)

K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) + 0*dendro.kodiss(K)

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
         Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
         Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
         Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])
         # + kod(i,Gt[i])

Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

eta_func = R0*sqrt(sum([igt[i,j]*d(i,chi)*d(j,chi) for i,j in dendro.e_ij]))/((1-chi**ep1)**ep2)

B_rhs = [Gt_rhs[i] - eta_func * B[i] +
         l3 * dendro.vec_j_ad_j(b, B[i]) -
         l4 * dendro.vec_j_ad_j(b, Gt[i]) + 0*kod(i,B[i])
         for i in dendro.e_i]

# Additional Equations

Rsc_rhs = dendro.lie(b, Rsc) - a*Rsch

Rsch_rhs = dendro.lie(b, Rsch) - a*chi*sum(igt[i,j]*d2(i,j,Rsc) for i,j in dendro.e_ij) - \
           chi*sum(igt[i,j]*d(i,a)*d(j,Rsc) for i,j in dendro.e_ij) + \
           a*chi*sum(Gt[i]*d(i,Rsc) for i in dendro.e_i) + \
           1/2*a*sum(igt[i,j]*d(i,Rsc)*d(j,chi) for i,j in dendro.e_ij) + \
           a*K*Rsch + a*Rsc/(32*PI*(3*b_const-2*a_const))

#TODO: not sure if I understand but in the implementation Rtt is a 3D object while in the notes R_ab is still 4D? Same for Vat.
Rtt_rhs = dendro.lie(b, Rtt, weight) - a * Vat 

#TODO : Check eqns.. somewhat not working

#TODO: gt refers to 3D metric according to above assignment. in the notes these terms are still 4D metrics.
#TODO: implement Q terms correctly (in which there appear 1st order time derivatives again) ... presume this works along the lines of usual BSSN
Yabnp = - kappa*Rtt - \
	kappa*(2*b_const - a_const)/(8*(3*b_const-a_const))*gt*Rsc
#TODO: implement problematic terms correctly ... probably requires rewriting of other RHS bits too
Yabp = Rtt

#TODO: changed to current form in the notes ... terms in line 1 of (4.30) are ommitted for now because there may be a direct implementation of the spatial covariant derivative ... otherwise, I can still generate the expressions via xAct but they are quite huge
#TODO; confirm that dendro.laplacian(*) works for non-scalar *
Vat_rhs = Matrix([sum(b[k]*d(k,Vat[i,j]) for k in dendro.e_i) for i,j in dendro.e_ij]) + \
          Matrix([sum(Vat[i,k]*d(j,b[k]) for k in dendro.e_i) for i,j in dendro.e_ij]) + \
          Matrix([sum(Vat[k,j]*d(i,b[k]) for k in dendro.e_i) for i,j in dendro.e_ij]) - \
          a*Matrix([dendro.laplacian(Rtt[i,j],chi) for i,j in dendro.e_ij]) + \
	  a*Matrix([K*Vat[i,j] for i,j in dendro.e_ij])
# TODO: commented out Yabnp and Yabp terms since these require rewriting before implementation
#+ \
#	  a*Matrix([Yabnp[i,j] for i,j in dendro.e_ij]) + \
#	  a*qg_ho_coup*Matrix([Yabp[i,j] for i,j in dendro.e_ij])

          
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

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs, Rsc_rhs, Rsch_rhs, Rtt_rhs, Vat_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs', 'Rsc_rhs','Rsch_rhs','Rtt_rhs','Vat_rhs']
#dendro.generate_debug(outs, vnames)
dendro.generate(outs, vnames, '[pp]')
#numVars=len(outs)
#for i in range(0,numVars):
#    dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')
