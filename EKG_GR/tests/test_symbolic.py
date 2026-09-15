#!/usr/bin/env python3
"""Algebra tests; these do not stand in for Dendro integration/AMR tests."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'CodeGen'))
import numpy as np
import sympy as s
from ekg_equations import build_system, PAIRS

x=build_system()
rng=np.random.default_rng(6127)
# Compare the BSSN conformal Ricci identity with direct Christoffel Ricci.
diff=list(x.tensors['Rtilde_geom']-x.tensors['Rtilde_direct'])
syms=sorted(set().union(*(v.free_symbols for v in diff)),key=str)
f=s.lambdify(syms,diff,'numpy',cse=True,docstring_limit=0)
for sample in range(12):
    z=rng.normal(size=(3,3));g=z@z.T+2*np.eye(3)
    values={str(q):rng.normal(scale=.04) for q in syms}
    for n,(i,j) in enumerate(PAIRS):values[f'gt{n}[pp]']=g[i,j]
    err=np.max(np.abs(f(*[values[str(q)] for q in syms])))
    assert err<2e-12,(sample,err)
print('PASS: conformal Ricci identity on 12 general SPD metric jets')
# Matter identities and the coupled Gamma-driver dependency.
q=x.fields
pi=s.Symbol('ekg_pi',real=True)
assert s.simplify(x.sources['srcK']-8*pi*q['alpha']*(q['scalarPi']**2-x.tensors['V']))==0
for i in range(3):
    assert s.simplify(x.full[f'B{i}']-x.vacuum[f'B{i}']-x.sources[f'srcGt{i}'])==0
print('PASS: K source and all three Gamma-driver source identities')
# Flat-space KG on a general nonzero first/second derivative jet.
flat={q[n]:0 for n in q}
for n in ('alpha','chi','gt0','gt3','gt5'):flat[q[n]]=1
flat[q['scalarPhi']]=s.Symbol('P');flat[q['scalarPi']]=s.Symbol('Q')
for v,spec in x.derivatives.items():
    if spec[3] not in ('scalarPhi','scalarPi'):flat[v]=0
phi_rhs=x.matter['scalarPhi'].subs(flat,simultaneous=True)
pi_rhs=x.matter['scalarPi'].subs(flat,simultaneous=True)
expected=-sum(s.Symbol(f'grad2_{i}_{i}_scalarPhi[pp]',real=True) for i in range(3))+s.Symbol('scalar_mu',real=True)**2*s.Symbol('P')
assert s.simplify(phi_rhs+s.Symbol('Q'))==0
assert s.simplify(pi_rhs-expected)==0
print('PASS: flat massive KG signs/dispersion operator')
