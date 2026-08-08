"""Corrected radial-vs-full material control.

This supersedes the development optimizer in `material_radial_vs_full.py` before
any held-out branch-residual bodies are examined.

The original shell learner computed the exact derivative with respect to one
shared shell density,

    dJ / d x_s = sum_{i in shell s} dJ / d d_i,

but projected its step using the expanded per-cell Euclidean metric.  Those two
choices are inconsistent: the exact shell-parameter gradient belongs with the
ordinary Euclidean metric in shell-parameter space.

For the constraint

    sum_s n_s x_s = material_budget,

ordinary Euclidean projection has KKT form

    x_s = clip(v_s - lambda * n_s, 0, cap),

not a common scalar shift.  This file installs that corrected projection and
reuses the otherwise frozen development comparison.

No held-out seeds 628+ are touched by this file/workflow.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import material_radial_vs_full as rvf


def shell_euclidean_project(v, weights, total, cap):
    """Euclidean projection in shell-parameter space.

    Minimize 0.5*sum_s (x_s-v_s)^2 subject to
        sum_s weights_s*x_s = total
        0 <= x_s <= cap.

    KKT gives x_s = clip(v_s - lambda*weights_s, 0, cap).
    """
    v=np.asarray(v,float)
    w=np.asarray(weights,float)
    total=float(total);cap=float(cap)
    if total < -1e-12 or total > cap*float(w.sum())+1e-12:
        raise ValueError('infeasible weighted bounded simplex')

    def residual(lam):
        x=np.clip(v-float(lam)*w,0.0,cap)
        return float(np.dot(w,x)-total)

    lo=-1.0
    while residual(lo)<0:
        lo*=2.0
    hi=1.0
    while residual(hi)>0:
        hi*=2.0
    for _ in range(140):
        mid=(lo+hi)/2.0
        if residual(mid)>0:
            lo=mid
        else:
            hi=mid
    x=np.clip(v-((lo+hi)/2.0)*w,0.0,cap)
    return x


# `rvf.train_shells` looks up this function in its own module namespace.
rvf.weighted_bounded_project=shell_euclidean_project


def selftest():
    v=np.array([-.1,.2,.8,.3])
    w=np.array([1.,4.,2.,7.])
    total=0.37
    cap=.08
    x=shell_euclidean_project(v,w,total,cap)
    assert np.all(x>=-1e-12) and np.all(x<=cap+1e-12)
    assert abs(float(np.dot(w,x))-total)<1e-10
    print('selftest ok',x,float(np.dot(w,x)))


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=622)
    ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2.0)
    ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005)
    ap.add_argument('--ratio',type=float,default=10.0)
    ap.add_argument('--full-steps',type=int,default=50)
    ap.add_argument('--radial-steps',type=int,default=160)
    ap.add_argument('--step-fraction',type=float,default=.10)
    ap.add_argument('--out',default='runs/material_radial_vs_full_v02/dev_622_627.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def main():
    a=parse_args()
    if a.selftest:
        selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    oms=[float(x) for x in a.omegas.split(',')]
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=rvf.body_probe(m,oms,a.tau,a.mu,a.g0,a.ratio,
                         a.full_steps,a.radial_steps,a.step_fraction)
        rows.append(r)
        e=r['evaluation']
        print(f"seed {seed}: uniform={e['uniform']['coherence2']:.4f} "
              f"radial={e['radial']['coherence2']:.4f} full={e['full']['coherence2']:.4f} "
              f"delta={e['full']['coherence2']-e['radial']['coherence2']:+.4f} shells={r['shells']}")
    out=dict(config=vars(a),projection='ordinary Euclidean shell-parameter projection',
             rows=rows,summary=rvf.summarize(rows))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True)
    q.write_text(json.dumps(out,indent=2))
    print('\nRADIAL VS FULL MATERIAL V02')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
