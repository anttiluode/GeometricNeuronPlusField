"""Predict return-gate gradient damage from boundary spectral sideband contamination.

Held-out controls established that equal 50% duty is not equal information:
fast regular interleaving preserves the gradient map, random masks are worse, and
contiguous/slow masks are much worse.

This probe tests a boundary-accessible mechanism.

For an exact soma return g(t) and a 50% mask m(t), write

    m(t) = 1/2 + r(t)

so

    FFT[m g] = 1/2 G + FFT[r g].

The first term is a scaled copy of the desired return.  The second term is spectral
sideband contamination created by the non-DC part of the mask.

Previous held-out work showed that a small common set of frequency bins selected from
source/return PORT spectra predicts most of the internal gradient direction.  Therefore
we ask whether contamination landing specifically in those K=8/K=16 boundary-selected
bins predicts gradient-map damage across periodic, random, and contiguous masks.

Important: bin selection and contamination score use only boundary signals and the mask.
The gated internal gradient is evaluated afterward as the dependent variable.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # mature-boundary patch
from biological_return_code_probe import build_exact, gradient_from_codes, map_metrics, l2_match
from reciprocal_adjoint_probe import retro_source_sequence
from spectral_correlation_compression_probe import port_spectrum_score
from transfer_decomposition_probe import safe_corr

PERIODS=(2,6,10,14,30,42,70)
KS=(8,16)
N_RANDOM=12
N_BLOCK=12


def periodic_mask(T,period,offset):
    p=int(period); on=p//2; t=np.arange(T)
    return ((t+int(offset))%p)<on


def block_mask(T,offset):
    n=T//2; idx=(np.arange(n)+int(offset))%T
    m=np.zeros(T,bool); m[idx]=True; return m


def random_mask(T,rng):
    n=T//2; idx=rng.choice(T,size=n,replace=False)
    m=np.zeros(T,bool); m[idx]=True; return m


def source_scalar_sequence(m,seq):
    """Collapse each external source frame to total complex injected amplitude."""
    return np.asarray([np.sum(np.asarray(x,np.complex128)) for x in seq],np.complex128)


def port_importance(m,seqT,seqD,gT,gD):
    rT=retro_source_sequence(m,gT,reverse=True)
    rD=retro_source_sequence(m,gD,reverse=True)
    # Reuse the already-validated boundary ranking from spectral compression.
    return np.asarray(port_spectrum_score(seqT,rT),float)+np.asarray(port_spectrum_score(seqD,rD),float)


def reverse_fft(g):
    return np.fft.fft(np.asarray(g,np.complex128)[::-1])


def contamination(g,mask):
    """Non-DC sideband term FFT[(m-1/2)g] in physical reverse-source coordinates."""
    r=np.asarray(mask,float)-0.5
    return reverse_fft(r*np.asarray(g,np.complex128))


def selected_contam_score(gT,gD,mask,kk,weights=None):
    GT=reverse_fft(gT); GD=reverse_fft(gD)
    CT=contamination(gT,mask); CD=contamination(gD,mask)
    kk=np.asarray(kk,int)
    if weights is None:
        w=np.ones(len(kk),float)
    else:
        w=np.asarray(weights,float)[kk]
        w=w/(np.mean(w)+1e-30)
    desired=np.concatenate([0.5*GT[kk],0.5*GD[kk]])
    contam=np.concatenate([CT[kk],CD[kk]])
    ww=np.concatenate([w,w])
    ratio=float(np.sqrt(np.sum(ww*np.abs(contam)**2)/(np.sum(ww*np.abs(desired)**2)+1e-30)))
    # Coherent residual after allowing one complex scalar: useful descriptive control.
    gated=desired+contam
    num=np.sum(ww*np.conj(desired)*gated)
    den=np.sum(ww*np.abs(desired)**2)+1e-30
    c=num/den
    resid=float(np.sqrt(np.sum(ww*np.abs(gated-c*desired)**2)/(np.sum(ww*np.abs(c*desired)**2)+1e-30)))
    cos=float(np.abs(np.sum(ww*np.conj(desired)*gated))/(np.sqrt(np.sum(ww*np.abs(desired)**2)*np.sum(ww*np.abs(gated)**2))+1e-30))
    return dict(contamination_ratio=ratio,regressed_residual=resid,complex_cosine=cos)


def evaluate_mask(m,z,mask,orders,importance):
    gt=l2_match(z['gT']*mask,z['gT']); gd=l2_match(z['gD']*mask,z['gD'])
    h,v=gradient_from_codes(m,z['wh'],z['wv'],z['pT'],z['pD'],gt,gd)
    mm=map_metrics(z['exact_h'],z['exact_v'],h,v)
    mm['spectral']={}
    for K,kk in orders.items():
        mm['spectral'][str(K)]=selected_contam_score(z['gT'],z['gD'],mask,kk,importance)
    return mm


def one(m,lag,steps):
    z=build_exact(m,lag,steps)
    seqT=ae.source_sequence(m,True,lag,steps)
    seqD=ae.source_sequence(m,False,lag,steps)
    imp=port_importance(m,seqT,seqD,z['gT'],z['gD'])
    order=np.argsort(imp)[::-1]
    orders={K:np.asarray(order[:K],int) for K in KS}
    T=len(z['gT'])
    rng=np.random.default_rng(int(m.cfg.seed)+990011)
    points=[]

    for P in PERIODS:
        for off in range(P):
            mask=periodic_mask(T,P,off)
            q=evaluate_mask(m,z,mask,orders,imp)
            q.update(kind='periodic',period=int(P),offset=int(off))
            points.append(q)

    for j in range(N_RANDOM):
        mask=random_mask(T,rng)
        q=evaluate_mask(m,z,mask,orders,imp)
        q.update(kind='random',draw=int(j))
        points.append(q)

    offs=np.linspace(0,T-1,N_BLOCK,dtype=int)
    for j,off in enumerate(offs):
        mask=block_mask(T,int(off))
        q=evaluate_mask(m,z,mask,orders,imp)
        q.update(kind='block',draw=int(j),offset=int(off))
        points.append(q)

    return dict(seed=int(m.cfg.seed),C=float(z['C']),bins={str(K):[int(x) for x in orders[K]] for K in KS},points=points)


def summarize(rows):
    allp=[(r['seed'],p) for r in rows for p in r['points']]
    out=dict(bodies=len(rows),K={},periodic={})
    for K in KS:
        k=str(K)
        contam=np.asarray([p['spectral'][k]['contamination_ratio'] for _,p in allp],float)
        resid=np.asarray([p['spectral'][k]['regressed_residual'] for _,p in allp],float)
        cos=np.asarray([p['spectral'][k]['complex_cosine'] for _,p in allp],float)
        map_err=np.asarray([1.0-p['corr'] for _,p in allp],float)
        rel=np.asarray([p['relative_l2'] for _,p in allp],float)
        out['K'][k]=dict(
            contamination_vs_one_minus_corr=float(safe_corr(contam,map_err)),
            contamination_vs_relative_l2=float(safe_corr(contam,rel)),
            regressed_residual_vs_one_minus_corr=float(safe_corr(resid,map_err)),
            spectral_cosine_vs_map_corr=float(safe_corr(cos,np.asarray([p['corr'] for _,p in allp],float))),
        )

    for P in PERIODS:
        q=[p for _,p in allp if p['kind']=='periodic' and p['period']==P]
        out['periodic'][str(P)]=dict(
            mean_map_corr=float(np.mean([x['corr'] for x in q])),
            K8_mean_contamination=float(np.mean([x['spectral']['8']['contamination_ratio'] for x in q])),
            K16_mean_contamination=float(np.mean([x['spectral']['16']['contamination_ratio'] for x in q])),
        )
    for kind in ('random','block'):
        q=[p for _,p in allp if p['kind']==kind]
        out[kind]=dict(
            mean_map_corr=float(np.mean([x['corr'] for x in q])),
            K8_mean_contamination=float(np.mean([x['spectral']['8']['contamination_ratio'] for x in q])),
            K16_mean_contamination=float(np.mean([x['spectral']['16']['contamination_ratio'] for x in q])),
        )
    return out


def selftest():
    rng=np.random.default_rng(1); T=210
    g=rng.normal(size=T)+1j*rng.normal(size=T)
    m=periodic_mask(T,14,3)
    lhs=reverse_fft(m*g)
    rhs=0.5*reverse_fft(g)+contamination(g,m)
    assert np.max(np.abs(lhs-rhs))<1e-10
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=534)
    ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/return_gate_spectral/dev.json')
    ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps);rows.append(r)
        # compact per-body diagnostic
        pp={P:[x for x in r['points'] if x['kind']=='periodic' and x['period']==P] for P in PERIODS}
        print('seed',seed,
              'P6',round(np.mean([x['corr'] for x in pp[6]]),4),round(np.mean([x['spectral']['8']['contamination_ratio'] for x in pp[6]]),4),
              'P42',round(np.mean([x['corr'] for x in pp[42]]),4),round(np.mean([x['spectral']['8']['contamination_ratio'] for x in pp[42]]),4),flush=True)
    if not rows:raise SystemExit('No valid bodies')
    s=summarize(rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(dict(experiment='return_gate_spectral_mechanism_dev_v01',summary=s,rows=rows),indent=2))
    print('\nRETURN GATE SPECTRAL MECHANISM DEV')
    print(json.dumps(s,indent=2))

if __name__=='__main__':main()
