#!/usr/bin/env python3
"""Emit matter-only, vacuum, fused BSSN+matter and constraint CPU kernels.

Normal use:
  python generate_ekg.py --out-dir ../generated --backreaction on

Default emitter: DendroSym. No source installation/patching happens here.
--reference-emitter is ONLY for algebra/syntax testing without DendroSym;
CMake refuses that output for the production ekgSolver target.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import sys
import sympy as s
from ekg_equations import build_system, FIELDS, ENUMS


def emit(expressions, names, reference=False):
    count=sum(s.count_ops(x) for x in expressions)
    if reference:
        pairs,values=s.cse(expressions, symbols=s.numbered_symbols('EKG_TMP_'),
                           optimizations=None, order='canonical')
        code='\n'.join(f'const double {v} = {s.ccode(e)};' for v,e in pairs)
        code+='\n'+'\n'.join(f'{n}[pp] = {s.ccode(e)};' for n,e in zip(names,values))+'\n'
    else:
        try:
            import dendrosym
            from dendrosym import codegen
        except ImportError as exc:
            raise RuntimeError('This Python cannot import dendrosym. Use the interpreter '
                               'where you installed it; set Python3_EXECUTABLE in CMake.') from exc
        for name in ('construct_cse_from_list','generate_cpu_preextracted'):
            if not callable(getattr(codegen,name,None)):
                raise RuntimeError(f'DendroSym codegen.{name} is unavailable in this checkout')
        cse=codegen.construct_cse_from_list(expressions)
        code=codegen.generate_cpu_preextracted(cse,names,'[pp]',count,
                                               dtype='double',use_const=True)
        if not isinstance(code,str):
            raise RuntimeError('generate_cpu_preextracted did not return a C++ string; '
                               'check your installed API rather than silently using a fallback')
    return code,int(count)


def make_local(code,name):
    pattern=rf'(?m)^\s*{re.escape(name)}\[pp\]\s*='
    out,n=re.subn(pattern,f'const double {name} =',code)
    if n!=1:
        raise RuntimeError(f'Expected exactly one emitted assignment for {name}; found {n}')
    return out


def generate(args):
    dest=args.out_dir.resolve();dest.mkdir(parents=True,exist_ok=True)
    # Make stale fragments visibly unusable until generation completes.
    manifest_path=dest/'manifest.json'
    manifest_path.unlink(missing_ok=True)
    coupled=args.backreaction=='on'
    system=build_system(args.potential,args.gauge,args.advection,coupled)
    groups=[('ekg_full_rhs.inc',system.full,'out_',False),
            ('ekg_vacuum_rhs.inc',system.vacuum,'out_',False),
            ('ekg_matter_rhs.inc',system.matter,'ekg_',True),
            ('ekg_constraints.inc',system.constraints,'constraint_',False)]
    metrics={}
    used_symbols=set()
    for filename,group,prefix,is_matter in groups:
        names=[]
        for name in group:
            names.append('out_'+name if is_matter and name in ('scalarPhi','scalarPi')
                         else prefix+name)
        print(f'Emitting {filename} ({len(group)} outputs)',flush=True)
        code,ops=emit(list(group.values()),names,args.reference_emitter)
        if is_matter:
            for name in group:
                if name not in ('scalarPhi','scalarPi'):code=make_local(code,'ekg_'+name)
        header=('// GENERATED: edit CodeGen/ekg_equations.py, not this file.\n'
                f'// potential={args.potential}; gauge={args.gauge}; '
                f'backreaction={args.backreaction}; advection={args.advection}\n')
        (dest/filename).write_text(header+code)
        metrics[filename]={'outputs':len(group),'original_ops':ops}
        for expr in group.values(): used_symbols |= expr.free_symbols
    # Include centered gradients for boundaries and mixed Hessian assembly.
    derivs=dict(system.derivatives)
    for field in FIELDS:
        for i in range(3):
            symbol=s.Symbol(f'grad_{i}_{field}[pp]',real=True)
            derivs[symbol]=('grad',i,-1,field)
            used_symbols.add(symbol)
    selected=[(v,spec) for v,spec in derivs.items() if v in used_symbols]
    kind_order={'grad':0,'grad2':1,'agrad':2}
    selected.sort(key=lambda item:(kind_order[item[1][0]],FIELDS.index(item[1][3]),
                                  item[1][1],item[1][2]))
    names={v.name[:-4]:i for i,(v,_) in enumerate(selected)}
    specs=[]
    for symbol,(kind,i,j,field) in selected:
        parent=names[f'grad_{i}_{field}'] if kind=='grad2' and i!=j else -1
        specs.append(({'grad':0,'grad2':1,'agrad':2}[kind],FIELDS.index(field),i,j,parent))
    h=['#pragma once','#include <array>','#include <cstddef>',
       'namespace ekg::generated {',
       'struct DerivativeSpec { int kind,field,i,j,parent; };',
       f'inline constexpr std::array<DerivativeSpec,{len(specs)}> derivatives = {{{{']
    h += ['    {'+','.join(map(str,v))+'},' for v in specs]
    h+=['}};','inline constexpr unsigned fields = 26;',
        'inline constexpr unsigned geometry_fields = 24;','}']
    (dest/'ekg_derivative_specs.h').write_text('\n'.join(h)+'\n')
    def write(name,lines): (dest/name).write_text('\n'.join(lines)+'\n')
    write('ekg_input_aliases.inc',[f'const double* {f} = in[{i}]+offset;' for i,f in enumerate(FIELDS)])
    write('ekg_output_aliases.inc',[f'double* out_{f} = out[{i}]+offset;' for i,f in enumerate(FIELDS)])
    write('ekg_derivative_aliases.inc',[f'const double* {v.name[:-4]} = work[{i}].data();'
                                      for i,(v,_) in enumerate(selected)])
    write('ekg_apply_sources.inc',[
        '#if EKG_BACKREACTION',
        'out_K[pp] += ekg_srcK;',
        *[f'out_At{i}[pp] += ekg_srcAt{i};' for i in range(6)],
        *[f'out_Gt{i}[pp] += ekg_srcGt{i};' for i in range(3)],
        *[f'out_B{i}[pp] += ekg_srcB{i};' for i in range(3)],
        '#endif'])
    write('ekg_field_asserts.inc',[
        'static_assert(bssn::BSSN_NUM_VARS == 26, "Native context must use 26 fields");',
        *[f'static_assert(bssn::{e} == {i}, "EKG field layout mismatch: {e}");'
          for i,e in enumerate(ENUMS)]])
    write('ekg_generated_config.h',[
        '#pragma once',f'#define EKG_GENERATED_BACKREACTION {int(coupled)}',
        f'#define EKG_GENERATED_REFERENCE_EMITTER {int(args.reference_emitter)}',
        f'#define EKG_GENERATED_GAUGE_SYNCHRONOUS {int(args.gauge=="synchronous")}',
        f'#define EKG_GENERATED_POTENTIAL_AXION {int(args.potential=="axion")}',
        f'#define EKG_GENERATED_ADVECTION_UPWIND {int(args.advection=="upwind")}',
        '#define EKG_GENERATED_SCHEMA 2'])
    # Symbolic expressions used by the elliptic initial-data solver as well.
    x,mu,fa=s.symbols('phi mu fa')
    V=mu**2*x*x/2 if args.potential=='massive' else 2*mu**2*fa**2*s.sin(x/(2*fa))**2
    write('ekg_potential.h',[
        '#pragma once','#include <cmath>','namespace ekg {',
        'inline double potential(double phi,double mu,double fa) {',
        '    return '+s.ccode(V)+';','}',
        'inline double potential_prime(double phi,double mu,double fa) {',
        '    return '+s.ccode(s.diff(V,x))+';','}','}'])
    digest=hashlib.sha256()
    for f in ('ekg_equations.py','generate_ekg.py'):
        digest.update((Path(__file__).parent/f).read_bytes())
    try: ds_version=importlib.metadata.version('dendrosym')
    except importlib.metadata.PackageNotFoundError: ds_version='not installed'
    manifest={'schema':2,'backreaction':coupled,'potential':args.potential,
              'gauge':args.gauge,'advection':args.advection,'field_order':FIELDS,
              'derivative_arrays':len(specs),'outputs':metrics,
              'generator_sha256':digest.hexdigest(),'sympy_version':s.__version__,
              'dendrosym_version':ds_version,
              'emitter':'sympy-reference' if args.reference_emitter else 'dendrosym',
              'full_solver_run_validated':False}
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Wrote {dest}; derivatives={len(specs)}; emitter={manifest["emitter"]}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out-dir',type=Path,required=True)
    p.add_argument('--backreaction',choices=['on','off'],default='on')
    p.add_argument('--potential',choices=['massive','axion'],default='massive')
    p.add_argument('--gauge',choices=['puncture','synchronous'],default='puncture')
    p.add_argument('--advection',choices=['centered','upwind'],default='centered')
    p.add_argument('--reference-emitter',action='store_true',help='TESTS ONLY; not accepted by ekgSolver')
    args=p.parse_args()
    try:generate(args)
    except (RuntimeError,ValueError,OSError) as exc:
        print('Code generation failed:',exc,file=sys.stderr);raise SystemExit(2)
if __name__=='__main__':main()
