"""Is the soma a good TASK projection even if it is a poor full modal observer?

Maps A->B versus B->A temporal-order contrast over every occupied cell and over
coherent graph-radius apertures. Soma-centered balls are compared both with
scattered equal-count controls and, crucially, with same-radius contiguous balls
centered on every other body cell.

Same frozen body, same wave, no growth, credit or learning.
"""
from __future__ import annotations
import argparse, json, sys
from collections import deque
from pathlib import Path
import numpy as np

def n4(y,x,h,w):
    for dy,dx in ((-1,0),(1,0),(0,-1),(0,1)):
        yy,xx=y+dy,x+dx
        if 0<=yy<h and 0<=xx<w: yield yy,xx

def dists(body,start):
    b=body.astype(bool); d={tuple(start):0}; q=deque([tuple(start)])
    while q:
        p=q.popleft()
        for r in n4(*p,*b.shape):
            if b[r] and r not in d: d[r]=d[p]+1; q.append(r)
    return d

def reset(m):
    try:m.reset_fast(clear_traces=True)
    except TypeError:m.reset_fast(True)

def addsrc(a,b):
    if isinstance(a,(float,int,np.floating)): return b
    if isinstance(b,(float,int,np.floating)): return a
    return a+b

def trace_fields(m,lag,target,steps,coords):
    reset(m); first,second=(0,1) if target else (1,0)
    Z=np.zeros((steps,len(coords)),np.complex128)
    for t in range(steps):
        a=m.pulse_source(first,t,False); b=m.pulse_source(second,t-lag,False)
        m.advance(addsrc(a,b),False,True,'none'); Z[t]=[m.psi[p] for p in coords]
    return Z

def contrast(a,b): return (a-b)/(a+b+1e-12)
def pct(x,arr):
    arr=np.asarray(arr,float)
    return float((np.sum(arr<x)+.5*np.sum(arr==x))/len(arr)) if len(arr) else float('nan')
def aperture_absC(tg,ds,ids):
    ids=np.asarray(ids,int)
    if not len(ids): return 0.0
    a=np.max(np.abs(tg[:,ids].mean(axis=1))**2)
    b=np.max(np.abs(ds[:,ids].mean(axis=1))**2)
    return abs(float(contrast(a,b)))

def parse_args():
    ap=argparse.ArgumentParser(); ap.add_argument('--functional-arbors',default='../FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=0); ap.add_argument('--seeds',type=int,default=24)
    ap.add_argument('--lag',type=int,default=20); ap.add_argument('--steps',type=int,default=150)
    ap.add_argument('--radii',default='1,2,3,4,5,6'); ap.add_argument('--random-controls',type=int,default=250)
    ap.add_argument('--out',default='runs/task_bottleneck/task_bottleneck.json'); return ap.parse_args()

def main():
    a=parse_args(); fa=Path(a.functional_arbors).resolve()
    if not fa.exists(): raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa)); from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    radii=[int(x) for x in a.radii.split(',') if x.strip()]; rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True; coords=[tuple(map(int,p)) for p in np.argwhere(m.body>0)]; idx={p:i for i,p in enumerate(coords)}
        tg=trace_fields(m,a.lag,True,a.steps,coords); ds=trace_fields(m,a.lag,False,a.steps,coords)
        pT=np.max(np.abs(tg)**2,axis=0); pD=np.max(np.abs(ds)**2,axis=0); C=contrast(pT,pD); absC=np.abs(C)
        si=idx[tuple(m.soma)]; soma_abs=float(absC[si]); soma_pct=pct(soma_abs,absC)
        allgd={p:dists(m.body,p) for p in coords}; gd=allgd[tuple(m.soma)]; rng=np.random.default_rng(seed+5150); balls=[]
        for R in radii:
            ids=np.asarray([idx[p] for p in coords if gd[p]<=R],int); k=len(ids); ac=aperture_absC(tg,ds,ids)
            rc=[]
            for _ in range(a.random_controls): rc.append(aperture_absC(tg,ds,rng.choice(len(coords),k,replace=False)))
            local=[]; local_n=[]
            for center in coords:
                cid=[idx[p] for p in coords if allgd[center].get(p,999)<=R]
                local.append(aperture_absC(tg,ds,cid)); local_n.append(len(cid))
            balls.append(dict(radius=R,n=k,abs_contrast=ac,random_absC_mean=float(np.mean(rc)),random_percentile=pct(ac,rc),
                              localball_absC_mean=float(np.mean(local)),localball_percentile=pct(ac,local),
                              localball_n_mean=float(np.mean(local_n))))
        rows.append(dict(seed=seed,cells=len(coords),soma_absC=soma_abs,soma_cell_percentile=soma_pct,
                         cell_absC_mean=float(absC.mean()),cell_absC_max=float(absC.max()),balls=balls))
        print(f'seed {seed:2d}: soma |C|={soma_abs:.3f} cell-pct={soma_pct:.2f}')
    if not rows: raise SystemExit('No bodies')
    summary=dict(bodies=len(rows),soma_absC_mean=float(np.mean([r['soma_absC'] for r in rows])),
        soma_cell_percentile_median=float(np.median([r['soma_cell_percentile'] for r in rows])),
        soma_top_quartile_count=int(sum(r['soma_cell_percentile']>=.75 for r in rows)),
        soma_top_decile_count=int(sum(r['soma_cell_percentile']>=.90 for r in rows)),
        cell_absC_mean=float(np.mean([r['cell_absC_mean'] for r in rows])),balls=[])
    for R in radii:
        rr=[next(q for q in r['balls'] if q['radius']==R) for r in rows]
        summary['balls'].append(dict(radius=R,n_mean=float(np.mean([q['n'] for q in rr])),absC_mean=float(np.mean([q['abs_contrast'] for q in rr])),
            random_absC_mean=float(np.mean([q['random_absC_mean'] for q in rr])),random_percentile_median=float(np.median([q['random_percentile'] for q in rr])),
            localball_absC_mean=float(np.mean([q['localball_absC_mean'] for q in rr])),localball_percentile_median=float(np.median([q['localball_percentile'] for q in rr])),
            localball_above_median_count=int(sum(q['localball_percentile']>.5 for q in rr))))
    payload=dict(experiment='task_bottleneck_v02',lag=a.lag,summary=summary,rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(payload,indent=2))
    print('\nTASK BOTTLENECK RECEIPT v0.2'); print(f'  soma mean |C|        {summary["soma_absC_mean"]:.4f}')
    print(f'  all-cell mean |C|    {summary["cell_absC_mean"]:.4f}'); print(f'  soma median cell pct {summary["soma_cell_percentile_median"]:.3f}')
    for q in summary['balls']: print(f'  R={q["radius"]}: |C|={q["absC_mean"]:.3f} local={q["localball_absC_mean"]:.3f} local-pct={q["localball_percentile_median"]:.3f}')
    print(f'  wrote {out}')
if __name__=='__main__': main()
