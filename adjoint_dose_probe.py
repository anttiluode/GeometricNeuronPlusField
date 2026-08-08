"""Conductance-step radius of validity for the soma adjoint.

Corrects the tiny exterior-boundary omission in the first adjoint surrogate, then
holds the base adjoint fixed while each structural event is applied gradually.
See ADJOINT_DOSE_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
from transfer_decomposition_probe import safe_corr
from structural_interference_probe import body_probe as nonlinear_event_probe

ALPHAS = (1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1.0)


def exact_weighted_lap(u, wh, wv):
    """Match FunctionalArbor _lap including weak coupling to zero outside grid."""
    u=np.asarray(u)
    out=np.zeros_like(u,dtype=np.result_type(u,np.complex128))
    dh=u[:,1:]-u[:,:-1]
    out[:,:-1]+=wh*dh; out[:,1:]-=wh*dh
    dv=u[1:,:]-u[:-1,:]
    out[:-1,:]+=wv*dv; out[1:,:]-=wv*dv
    kb=float(min(np.min(wh),np.min(wv)))
    out[:,0]-=kb*u[:,0]; out[:,-1]-=kb*u[:,-1]
    out[0,:]-=kb*u[0,:]; out[-1,:]-=kb*u[-1,:]
    return out

# All imported forward/adjoint functions look up this module-global function at runtime.
ae.weighted_lap = exact_weighted_lap


def apply_event_fraction(base_wh,base_wv,body,event,dk,alpha):
    wh=base_wh.copy();wv=base_wv.copy()
    for p,q,sgn in ae.event_edges(body,event):
        w0=ae.edge_lookup(wh,wv,p,q)
        ae.set_edge(wh,wv,p,q,w0+float(sgn)*float(alpha)*float(dk))
    return wh,wv


def percentile_sign(pred,actual,eps=1e-12):
    pairs=[(p,a) for p,a in zip(pred,actual) if np.isfinite(p) and np.isfinite(a) and abs(a)>eps]
    if not pairs:return float('nan')
    return float(np.mean([np.sign(p)==np.sign(a) for p,a in pairs]))


def body_probe(m,lag,steps,max_each):
    observed=nonlinear_event_probe(m,lag,steps,max_each)
    seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    base=ae.contrast_adjoint(m,m.body,seqT,seqD)
    fd=ae.finite_diff_check(m,base,seqT,seqD)

    rng=np.random.default_rng(int(m.cfg.seed)+992)
    u=rng.normal(size=m.body.shape)+1j*rng.normal(size=m.body.shape)
    lap0=m._lap(u,True);lap1=exact_weighted_lap(u,base['wh'],base['wv'])
    laprel=float(np.linalg.norm(lap1-lap0)/(np.linalg.norm(lap0)+1e-30))
    dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)

    events=[]
    for e in observed['events']:
        pred_full=ae.event_score(base['gh'],base['gv'],m.body,e,dk)
        row=dict(kind=e['kind'],cell=e['cell'],dist_soma=e['dist_soma'],
                 base_derivative_full_step=float(pred_full),
                 dC_int=float(e['dCint']),dC_peak=float(e['dCpeak']),dose={})
        for alpha in ALPHAS:
            wh,wv=apply_event_fraction(base['wh'],base['wv'],m.body,e,dk,alpha)
            ca,_,_=ae.linear_contrast(m,wh,wv,seqT,seqD)
            d=float(ca-base['C'])
            row['dose'][str(alpha)]=dict(
                exact_dC_lin=d,
                adjoint_prediction=float(alpha*pred_full),
            )
        events.append(row)

    dose_summary={}
    useful=[]
    for alpha in ALPHAS:
        key=str(alpha)
        pred=np.asarray([e['dose'][key]['adjoint_prediction'] for e in events],float)
        act=np.asarray([e['dose'][key]['exact_dC_lin'] for e in events],float)
        r=safe_corr(pred,act)
        dose_summary[key]=dict(r=float(r),sign_agreement=percentile_sign(pred,act,1e-12),
                               mean_abs_exact=float(np.mean(np.abs(act))))
        if np.isfinite(r) and r>=.70:useful.append(alpha)
    useful_max=float(max(useful)) if useful else 0.0

    full=np.asarray([e['dose']['1.0']['exact_dC_lin'] for e in events],float)
    dint=np.asarray([e['dC_int'] for e in events],float)
    dpeak=np.asarray([e['dC_peak'] for e in events],float)
    return dict(seed=int(m.cfg.seed),cells=int(m.body.sum()),n_events=len(events),
                laplacian_relative_error=laprel,fd_check=fd,
                base_C_lin=float(base['C']),base_C_int=float(observed['base']['Cint']),
                dose_summary=dose_summary,useful_alpha_max=useful_max,
                r_full_lin_int=float(safe_corr(full,dint)),
                r_full_lin_peak=float(safe_corr(full,dpeak)),events=events)


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=192)
    ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-each',type=int,default=6)
    ap.add_argument('--out',default='runs/adjoint_dose_discovery/adjoint_dose.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    wh=np.ones((3,2));wv=np.ones((2,3));u=np.ones((3,3),complex)
    L=exact_weighted_lap(u,wh,wv)
    # Interior constant field has zero lap; boundaries leak to zero exterior.
    assert L[1,1]==0
    assert L[0,0].real==-2
    print('selftest ok')


def main():
    a=parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor

    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,a.lag,a.steps,a.max_each);rows.append(r)
        print(f"seed {seed}: r1e-4={r['dose_summary']['0.0001']['r']:+.3f} "
              f"r1e-2={r['dose_summary']['0.01']['r']:+.3f} r1={r['dose_summary']['1.0']['r']:+.3f} "
              f"useful={r['useful_alpha_max']:.4g} lin/int={r['r_full_lin_int']:+.3f} "
              f"lap={r['laplacian_relative_error']:.2e}",flush=True)

    if not rows:raise SystemExit('No valid bodies')
    summary=dict(bodies=len(rows),alphas=list(ALPHAS),
                 mean_laplacian_relative_error=float(np.mean([r['laplacian_relative_error'] for r in rows])),
                 max_laplacian_relative_error=float(np.max([r['laplacian_relative_error'] for r in rows])),
                 mean_fd_relative_error=float(np.mean([r['fd_check']['relative_error'] for r in rows])),
                 max_fd_relative_error=float(np.max([r['fd_check']['relative_error'] for r in rows])))
    curves={}
    for alpha in ALPHAS:
        k=str(alpha);rs=np.asarray([r['dose_summary'][k]['r'] for r in rows],float)
        ss=np.asarray([r['dose_summary'][k]['sign_agreement'] for r in rows],float)
        curves[k]=dict(mean_r=float(np.nanmean(rs)),median_r=float(np.nanmedian(rs)),positive_bodies=int(np.sum(rs>0)),
                       mean_sign_agreement=float(np.nanmean(ss)))
    summary['dose_curve']=curves
    useful=np.asarray([r['useful_alpha_max'] for r in rows],float)
    summary['median_useful_alpha_max']=float(np.median(useful))
    summary['mean_useful_alpha_max']=float(np.mean(useful))
    fullint=np.asarray([r['r_full_lin_int'] for r in rows],float)
    summary['mean_r_full_lin_int']=float(np.nanmean(fullint));summary['positive_full_lin_int_bodies']=int(np.sum(fullint>0))
    summary['mean_r_full_lin_peak']=float(np.nanmean([r['r_full_lin_peak'] for r in rows]))
    summary['D0_pass']=bool(summary['mean_laplacian_relative_error']<1e-10 and summary['mean_fd_relative_error']<1e-3)
    summary['D1_pass']=bool(curves['0.0001']['mean_r']>.95 and curves['0.0001']['positive_bodies']==len(rows))
    summary['D2_pass']=bool(curves['0.0001']['mean_r']-curves['1.0']['mean_r']>.60)
    summary['D3_pass']=bool(summary['median_useful_alpha_max']>=1e-3)
    summary['D4_pass']=bool(summary['mean_r_full_lin_int']>.90 and summary['positive_full_lin_int_bodies']==len(rows))

    payload=dict(experiment='adjoint_dose_discovery_v01',prereg='ADJOINT_DOSE_DISCOVERY_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,max_each=a.max_each,
                 summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nADJOINT DOSE DISCOVERY RECEIPT')
    for k,v in summary.items():print(f' {k}: {v}')

if __name__=='__main__':main()
