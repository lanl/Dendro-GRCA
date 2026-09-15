"""Symbolic, second-order-in-space BSSN + real Klein-Gordon equations.

All evolution algebra lives here, not in a hand-written C++ point kernel.
DendroSym's low-level CSE/CPU emitter accepts these ordinary SymPy expressions.
Derivative inputs are explicit array symbols; their declarations AND finite-
difference calls are generated from the same dependency list.

Signature (-,+,+,+), Pi=-n^a partial_a phi, gamma_ij=gt_ij/chi,
K_ij=-1/2 L_n gamma_ij, G=c=1. Standard puncture gauge follows the
bssn_puncture_gauge structure in the supplied bssn_eqns_config.py.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import sympy as s

I = range(3)
PAIRS = [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
FIELDS = (['alpha','chi','K'] + [f'Gt{i}' for i in I]
          + [f'beta{i}' for i in I] + [f'B{i}' for i in I]
          + [f'gt{i}' for i in range(6)] + [f'At{i}' for i in range(6)]
          + ['scalarPhi','scalarPi'])
ENUMS = (['U_ALPHA','U_CHI','U_K'] + [f'U_GT{i}' for i in I]
         + [f'U_BETA{i}' for i in I] + [f'U_B{i}' for i in I]
         + [f'U_SYMGT{i}' for i in range(6)] + [f'U_SYMAT{i}' for i in range(6)]
         + ['U_SCALARPHI','U_SCALARPI'])
R = s.Rational

@dataclass
class System:
    vacuum: dict
    matter: dict
    full: dict
    constraints: dict
    sources: dict
    derivatives: dict
    fields: dict
    tensors: dict
    parameters: list


def build_system(potential='massive', gauge='puncture', advection='centered',
                 backreaction=True):
    if potential not in ('massive','axion'):
        raise ValueError('potential must be massive or axion')
    if gauge not in ('puncture','synchronous'):
        raise ValueError('gauge must be puncture or synchronous')
    if advection not in ('centered','upwind'):
        raise ValueError('advection must be centered or upwind')
    q = {n:s.Symbol(n+'[pp]', real=True) for n in FIELDS}
    derivative_map = {}
    name_of = {v:k for k,v in q.items()}
    def d(i, v, kind='grad'):
        name = f'{kind}_{i}_{name_of[v]}'
        out = s.Symbol(name+'[pp]', real=True)
        derivative_map[out] = (kind, i, -1, name_of[v])
        return out
    def dd(i,j,v):
        i,j = sorted((i,j))
        name = f'grad2_{i}_{j}_{name_of[v]}'
        out = s.Symbol(name+'[pp]', real=True)
        derivative_map[out] = ('grad2',i,j,name_of[v])
        return out
    def mat(prefix):
        out = s.zeros(3)
        for n,(i,j) in enumerate(PAIRS):
            out[i,j] = out[j,i] = q[f'{prefix}{n}']
        return out
    a,chi,K = (q[n] for n in ('alpha','chi','K'))
    phi,Pi = q['scalarPhi'],q['scalarPi']
    beta = [q[f'beta{i}'] for i in I]
    B = [q[f'B{i}'] for i in I]
    Gt = [q[f'Gt{i}'] for i in I]
    g,A = mat('gt'),mat('At')
    # Actual matrix inverse: do not assume det(gt)=1 inside algebraic contractions.
    ig = g.adjugate()/g.det()
    dg = [[[d(k,g[i,j]) for j in I] for i in I] for k in I]
    @lru_cache(None)
    def dig(k,i,j):
        return -sum(ig[i,l]*dg[k][l][m]*ig[m,j] for l,m in product(I,I))
    C = [[[sum(ig[i,l]*(dg[j][l][k]+dg[k][l][j]-dg[l][j][k])/2
                 for l in I) for k in I] for j in I] for i in I]
    C1 = [[[sum(g[i,l]*C[l][j][k] for l in I) for k in I] for j in I] for i in I]
    calG = [sum(ig[j,k]*C[i][j][k] for j,k in product(I,I)) for i in I]
    @lru_cache(None)
    def dC(e,i,j,k):
        return sum(dig(e,i,l)*(dg[j][l][k]+dg[k][l][j]-dg[l][j][k])/2
             + ig[i,l]*(dd(e,j,g[l,k])+dd(e,k,g[l,j])-dd(e,l,g[j,k]))/2 for l in I)
    @lru_cache(None)
    def dcalG(e,i):
        return sum(dig(e,j,k)*C[i][j][k]+ig[j,k]*dC(e,i,j,k)
                   for j,k in product(I,I))
    def conformal_ricci(use_connection_constraint):
        out=s.zeros(3)
        for i,j in PAIRS:
            def gradG(e,k):
                return dcalG(e,k) if use_connection_constraint else d(e,Gt[k])
            v=-sum(ig[l,m]*dd(l,m,g[i,j])/2 for l,m in product(I,I))
            v+=sum((g[k,i]*gradG(j,k)+g[k,j]*gradG(i,k))/2 for k in I)
            v+=sum(calG[k]*(C1[i][j][k]+C1[j][i][k])/2 for k in I)
            v+=sum(ig[l,m]*(C[k][l][i]*C1[j][k][m]
                           +C[k][l][j]*C1[i][k][m]
                           +C[k][i][m]*C1[k][l][j])
                   for l,m,k in product(I,I,I))
            out[i,j]=out[j,i]=v
        return out
    def direct_ricci():
        return s.Matrix(3,3,lambda i,j:
            sum(dC(k,k,i,j)-dC(j,k,i,k)
                +sum(C[k][i][j]*C[l][k][l]-C[l][i][k]*C[k][j][l] for l in I)
                for k in I))
    Hchi=s.Matrix(3,3,lambda i,j:dd(i,j,chi)-sum(C[k][i][j]*d(k,chi) for k in I))
    lapchi=sum(ig[i,j]*Hchi[i,j] for i,j in product(I,I))
    gradchi2=sum(ig[i,j]*d(i,chi)*d(j,chi) for i,j in product(I,I))
    Rchi=s.Matrix(3,3,lambda i,j:
        (Hchi[i,j]+g[i,j]*lapchi)/(2*chi)
        -(d(i,chi)*d(j,chi)+3*g[i,j]*gradchi2)/(4*chi**2))
    Rt=conformal_ricci(False)
    Rij=Rt+Rchi
    # Constraints use geometric Ricci, not derivatives of independently evolved Gt.
    Rgeom=conformal_ricci(True)+Rchi
    Cp=[[[C[k][i][j]-(s.KroneckerDelta(k,i)*d(j,chi)
          +s.KroneckerDelta(k,j)*d(i,chi)
          -g[i,j]*sum(ig[k,l]*d(l,chi) for l in I))/(2*chi)
          for j in I] for i in I] for k in I]
    def hess(v):
        return s.Matrix(3,3,lambda i,j:dd(i,j,v)-sum(Cp[k][i][j]*d(k,v) for k in I))
    def trace(v):
        return sum(ig[i,j]*v[i,j] for i,j in product(I,I))
    def tf(v): return v-g*trace(v)/3
    def adv(v):
        return sum(beta[i]*d(i,v,'agrad' if advection=='upwind' else 'grad') for i in I)
    divb=sum(d(i,beta[i]) for i in I)
    def lie_tensor(v):
        return s.Matrix(3,3,lambda i,j:adv(v[i,j])
            +sum(v[i,k]*d(j,beta[k])+v[k,j]*d(i,beta[k]) for k in I)
            -R(2,3)*v[i,j]*divb)
    mu,fa,eta,pi=s.symbols('scalar_mu scalar_fa eta ekg_pi', real=True)
    lam=s.symbols('lambda0:4', real=True)
    lf=s.symbols('lambda_f0:2', real=True)
    if potential=='massive':
        V=mu**2*phi**2/2
        Vp=mu**2*phi
    else:
        V=2*mu**2*fa**2*s.sin(phi/(2*fa))**2
        Vp=mu**2*fa*s.sin(phi/fa)
    gradphi2=chi*sum(ig[i,j]*d(i,phi)*d(j,phi) for i,j in product(I,I))
    rho=(Pi**2+gradphi2)/2+V
    Sc=s.Matrix([Pi*d(i,phi) for i in I])
    Sdd=s.Matrix(3,3,lambda i,j:d(i,phi)*d(j,phi)
          +g[i,j]/chi*((Pi**2-gradphi2)/2-V))
    S=R(3,2)*Pi**2-gradphi2/2-3*V
    scalar_rhs={'scalarPhi':adv(phi)-a*Pi,
        'scalarPi':adv(Pi)+a*K*Pi-a*chi*trace(hess(phi))
        -chi*sum(ig[i,j]*d(i,a)*d(j,phi) for i,j in product(I,I))+a*Vp}
    srcK=4*pi*a*(rho+S)
    # Use the exactly equivalent gradient-only form of the TF source.
    srcA=-8*pi*a*chi*tf(s.Matrix(3,3,lambda i,j:d(i,phi)*d(j,phi)))
    srcGt=s.Matrix([-16*pi*a*sum(ig[i,j]*Sc[j] for j in I) for i in I])
    AU=ig*A*ig
    Asq=sum(AU[i,j]*A[i,j] for i,j in product(I,I))
    gt_rhs=lie_tensor(g)-2*a*A
    chi_rhs=adv(chi)+R(2,3)*chi*(a*K-divb)
    At_vac=lie_tensor(A)+chi*tf(a*Rij-hess(a))+a*(K*A-2*A*ig*A)
    K_vac=adv(K)-chi*trace(hess(a))+a*(K*K/3+Asq)
    G_vac=s.Matrix([adv(Gt[i])-sum(calG[j]*d(j,beta[i]) for j in I)
        +R(2,3)*calG[i]*divb
        +sum(ig[j,k]*dd(j,k,beta[i])+ig[i,j]*dd(j,k,beta[k])/3
             for j,k in product(I,I))
        -2*sum(AU[i,j]*d(j,a) for j in I)
        +2*a*sum(C[i][j][k]*AU[j,k] for j,k in product(I,I))
        -a*sum(3*AU[i,j]*d(j,chi)/chi+R(4,3)*ig[i,j]*d(j,K) for j in I)
        for i in I])
    def assemble(coupled):
        GG=G_vac+srcGt if coupled else G_vac
        AA=At_vac+srcA if coupled else At_vac
        KK=K_vac+srcK if coupled else K_vac
        ans={'alpha':lam[0]*adv(a)-2*a*K if gauge=='puncture' else s.S.Zero,
             'chi':chi_rhs,'K':KK}
        for i in I:
            ans[f'Gt{i}']=GG[i]
            ans[f'beta{i}']=(R(3,4)*(lf[0]+lf[1]*a)*B[i]+lam[1]*adv(beta[i])
                              if gauge=='puncture' else s.S.Zero)
            ans[f'B{i}']=(GG[i]-eta*B[i]+lam[2]*adv(B[i])-lam[3]*adv(Gt[i])
                           if gauge=='puncture' else s.S.Zero)
        for n,(i,j) in enumerate(PAIRS):
            ans[f'gt{n}']=gt_rhs[i,j]
            ans[f'At{n}']=AA[i,j]
        return {name:ans[name] for name in FIELDS[:24]}
    vacuum=assemble(False)
    full=assemble(backreaction)
    full.update(scalar_rhs)
    sources={'srcK':srcK}
    sources.update({f'srcAt{n}':srcA[i,j] for n,(i,j) in enumerate(PAIRS)})
    sources.update({f'srcGt{i}':srcGt[i] for i in I})
    sources.update({f'srcB{i}':srcGt[i] if gauge=='puncture' else s.S.Zero for i in I})
    matter=dict(scalar_rhs)
    matter.update(sources)
    matter.update({'rho':rho,'S':S,'V':V,'gradphi2':gradphi2})
    matter.update({f'Scov{i}':Sc[i] for i in I})
    # Momentum in covariant-index form: D_j A^j_i - 2/3 partial_i K.
    Amix=ig*A
    mom=[]
    for i in I:
        mv=sum(sum(dig(j,j,k)*A[k,i]+ig[j,k]*d(j,A[k,i]) for k in I)
               +sum(Cp[j][j][l]*Amix[l,i]-Cp[l][j][i]*Amix[j,l] for l in I)
               for j in I)-R(2,3)*d(i,K)
        if backreaction: mv-=8*pi*Sc[i]
        mom.append(mv)
    ham=chi*trace(Rgeom)+R(2,3)*K**2-Asq
    if backreaction: ham-=16*pi*rho
    constraints={'ham':ham,**{f'mom{i}':mom[i] for i in I}}
    return System(vacuum,matter,full,constraints,sources,derivative_map,q,
        {'rho':rho,'S':S,'Sdd':Sdd,'g':g,'ig':ig,'Asq':Asq,'Rtilde':Rt,
         'Rtilde_geom':conformal_ricci(True),'Rtilde_direct':direct_ricci(),
         'calG':calG,'C':C,'V':V,'dV':Vp,'srcA':srcA,'srcGt':srcGt},
        [mu,fa,eta,pi,*lam,*lf])
