"""Is the soma a good TASK projection even if it is a poor full modal observer?

`local_observability_probe.py` asks a strong question: can a local soma ball
independently reconstruct the confirmed three-mode task band?  It cannot do that
robustly if the restricted mode matrix is ill-conditioned.  But a neuron does not
need to reconstruct all modal coordinates.  It may only need one consequential
combination.

This probe therefore returns to the actual frozen wave.  For A->B versus B->A at
lag 20 it maps the ordinary point-power contrast over EVERY occupied cell and asks
where the soma lies in that empirical task-selectivity distribution.  It also tests
uniform coherent graph balls around the soma against equal-size random apertures.

Same body, same wave, no growth, credit or learning.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np


def n4(y,x,h,w):
    for dy,dx in ((-1,0),(1,0),(0,-1),(0,1)):
        yy,xx=y+dy,x+dx
        if 0<=yy<h and 0<=xx<w:
            yield yy,xx


def dists(body,start):
    b=body.astype(bool); d={tuple(start):0}; q=deque([tuple(start)])
    while q:
        p=q.popleft()
        for r in n4(*p,*b.shape):
            if b[r] and r not in d:
                d[r]=d[p]+1; q.append(r)
    return d


def reset(m):
    try:m.reset_fast(clear_traces=True)
    except TypeError:m.reset_fast(True)


def addsrc(a,b):
    if isinstance(a,(float,int,np.floating)): return b
    if isinstance(b,(float,int,np.floating)): return a
    return a+b


def trace_fields(m,lag,target,steps,coords):
    reset(m)
    first,second=(0,1) if target else (1,0)
    Z=np.zeros((steps,len(coords)),np.complex128)
    for t in range(steps):
        a=m.pulse_source(first,t,False); b=m.pulse_source(second,t-lag,False)
        m.advance(addsrc(a,b),False,True,'none')
        Z[t]=[m.psi[p] for p in coords]
    return Z


def contrast(a,b):
    return (a-b)/(a+b+1e-12)


def pct(x,arr):
    arr=np.asarray(arr,float)
    return float((np.sum(arr<x)+.5*np.sum(arr==x))/len(arr))


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='../FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=0)
    ap.add_argument('--seeds',type=int,default=24)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=150)
    ap.add_argument('--radii',default='1,2,3,4,5,6')
    ap.add_argument('--random-controls',type=int,default=250)
    ap.add_argument('--out',default='runs/task_bottleneck/task_bottleneck.json')
    return ap.parse_args()


def main():
    a=parse_args(); fa=Path(a.functional_arbors).resolve()
    if not fa.exists(): raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    radii=[int(x) for x in a.radii.split(',') if x.strip()]
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True
        coords=[tuple(map(int,p)) for p in np.argwhere(m.body>0)]
        idx={p:i for i,p in enumerate(coords)}
        tg=trace_fields(m,a.lag,True,a.steps,coords)
        ds=trace_fields(m,a.lag,False,a.steps,coords)
        pT=np.max(np.abs(tg)**2,axis=0); pD=np.max(np.abs(ds)**2,axis=0)
        C=contrast(pT,pD); absC=np.abs(C)
        si=idx[tuple(m.soma)]
        somaC=float(C[si]); soma_abs=float(absC[si]); soma_pct=pct(soma_abs,absC)
        gd=dists(m.body,m.soma); rng=np.random.default_rng(seed+5150)
        balls=[]
        for R in radii:
            ids=np.asarray([idx[p] for p in coords if gd[p]<=R],int); k=len(ids)
            # Uniform coherent aperture: square magnitude of the mean complex field.
            tball=np.max(np.abs(tg[:,ids].mean(axis=1))**2)
            dball=np.max(np.abs(ds[:,ids].mean(axis=1))**2)
            cball=float(contrast(tball,dball)); ac=abs(cball)
            rc=[]
            for _ in range(a.random_controls):
                rid=rng.choice(len(coords),k,replace=False)
                rt=np.max(np.abs(tg[:,rid].mean(axis=1))**2)
                rd=np.max(np.abs(ds[:,rid].mean(axis=1))**2)
                rc.append(abs(float(contrast(rt,rd))))
            balls.append(dict(radius=R,n=k,contrast=cball,abs_contrast=ac,
                              random_absC_mean=float(np.mean(rc)),
                              percentile=pct(ac,rc)))
        rows.append(dict(seed=seed,cells=len(coords),soma_contrast=somaC,
                         soma_absC=soma_abs,soma_cell_percentile=soma_pct,
                         cell_absC_mean=float(absC.mean()),cell_absC_max=float(absC.max()),balls=balls))
        print(f'seed {seed:2d}: soma |C|={soma_abs:.3f} cell-percentile={soma_pct:.2f}')
    if not rows: raise SystemExit('No bodies')
    summary=dict(
        bodies=len(rows),
        soma_absC_mean=float(np.mean([r['soma_absC'] for r in rows])),
        soma_cell_percentile_median=float(np.median([r['soma_cell_percentile'] for r in rows])),
        soma_top_quartile_count=int(sum(r['soma_cell_percentile']>=.75 for r in rows)),
        soma_top_decile_count=int(sum(r['soma_cell_percentile']>=.90 for r in rows)),
        cell_absC_mean=float(np.mean([r['cell_absC_mean'] for r in rows])),
        balls=[])
    for R in radii:
        rr=[next(q for q in r['balls'] if q['radius']==R) for r in rows]
        summary['balls'].append(dict(radius=R,n_mean=float(np.mean([q['n'] for q in rr])),
            absC_mean=float(np.mean([q['abs_contrast'] for q in rr])),
            random_absC_mean=float(np.mean([q['random_absC_mean'] for q in rr])),
            percentile_median=float(np.median([q['percentile'] for q in rr])),
            above_random_count=int(sum(q['percentile']>.5 for q in rr))))
    payload=dict(experiment='task_bottleneck_v01',lag=a.lag,summary=summary,rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(payload,indent=2))
    print('\nTASK BOTTLENECK RECEIPT')
    print(f'  soma mean |C|             {summary["soma_absC_mean"]:.4f}')
    print(f'  all-cell mean |C|         {summary["cell_absC_mean"]:.4f}')
    print(f'  soma median cell pct      {summary["soma_cell_percentile_median"]:.3f}')
    print(f'  soma top quartile         {summary["soma_top_quartile_count"]}/{len(rows)}')
    for q in summary['balls']:
        print(f'  R={q["radius"]}: |C|={q["absC_mean"]:.3f} random={q["random_absC_mean"]:.3f} pct={q["percentile_median"]:.3f}')
    print(f'  wrote {out}')

if __name__=='__main__': main()
