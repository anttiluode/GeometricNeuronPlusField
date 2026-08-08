"""Compare full local material learning with an optimized distance-shell learner.

Exploratory mechanism control motivated by dendritic-democracy prior art.

The confirmed full learner has one quasi-active material density per occupied
cell.  Its dominant emergent coordinate is graph distance from the readout, but
smaller branch-specific corrections remain.

Here the control is given every advantage short of branch identity:

* same frozen arbor and quasi-active material physics;
* same uniform start;
* same total material budget and per-cell cap;
* same .03/.04 phase-coherence objective;
* exact analytic gradient;
* but all cells at the same integer graph distance must share one density.

Thus this is not a hand-drawn linear gradient.  It is the best distance-shell
material profile the same gradient machinery can discover.  The comparison
asks how much objective value the extra branch/cell degrees of freedom buy.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import hcn_material_learning as hml
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites


def weighted_bounded_project(v,weights,total,cap):
    """Project shell densities under sum_s weights[s]*x[s]=total, 0<=x<=cap.

    The metric is the expanded per-cell Euclidean metric
        sum_s weights[s] (x_s-v_s)^2,
    for which the KKT shift is common across shells.
    """
    v=np.asarray(v,float);w=np.asarray(weights,float)
    total=float(total);cap=float(cap)
    if total < -1e-12 or total > cap*w.sum()+1e-12:
        raise ValueError('infeasible weighted bounded simplex')
    lo=float(np.min(v-cap))-1.0;hi=float(np.max(v))+1.0
    for _ in range(120):
        mid=(lo+hi)/2
        x=np.clip(v-mid,0.0,cap)
        if float(np.dot(w,x))>total:lo=mid
        else:hi=mid
    return np.clip(v-(lo+hi)/2,0.0,cap)


def shell_structure(m,var_cells,var_ids):
    dmap=m.graph_distance_from_soma().astype(int)
    dist=np.asarray([int(dmap[p]) for p in var_cells],int)
    levels=np.unique(dist)
    groups=[np.where(dist==q)[0] for q in levels]
    counts=np.asarray([len(g) for g in groups],float)
    return dmap,dist,levels,groups,counts


def shell_to_grid(shape,var_ids,groups,x):
    vals=np.zeros(len(var_ids),float)
    for g,a in zip(groups,x):vals[g]=float(a)
    z=np.zeros(int(np.prod(shape)),float);z[var_ids]=vals
    return z.reshape(shape)


def train_shells(m,L,initial,omegas,tau,mu,sites,var_ids,groups,counts,total,cap,steps,step_fraction):
    # uniform initial state, represented by one value per distance shell
    initvals=np.asarray(initial,float).ravel()[var_ids]
    x=np.asarray([np.mean(initvals[g]) for g in groups],float)
    best_x=x.copy();best=-np.inf;hist=[]
    for it in range(int(steps)+1):
        d=shell_to_grid(m.body.shape,var_ids,groups,x)
        obj,g,det=hml.coherence_and_gradient(m,L,d,omegas,tau,mu,sites,var_ids,True)
        # Exact derivative wrt one shared shell density is sum of its cell derivatives.
        gs=np.asarray([np.sum(g[gidx]) for gidx in groups],float)
        hist.append(dict(iteration=it,objective=float(obj),grad_max=float(np.max(np.abs(gs))),
                         shell_min=float(x.min()),shell_max=float(x.max())))
        if obj>best:best=float(obj);best_x=x.copy()
        if it==steps:break
        gm=float(np.max(np.abs(gs)))
        if gm<1e-14:break
        base=float(step_fraction)*float(cap)/gm
        accepted=False
        for bt in range(14):
            trial=weighted_bounded_project(x+base*(.5**bt)*gs,counts,total,cap)
            td=shell_to_grid(m.body.shape,var_ids,groups,trial)
            tobj,_=hml.coherence_and_gradient(m,L,td,omegas,tau,mu,sites,var_ids,False)
            if tobj>=obj-1e-12:
                x=trial;accepted=True;break
        if not accepted:break
    return shell_to_grid(m.body.shape,var_ids,groups,best_x),best,hist,best_x


def body_probe(m,omegas,tau,mu,g0,ratio,full_steps,radial_steps,step_fraction):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    body=m.body.astype(bool);h,w=body.shape
    cells=[tuple(map(int,p)) for p in np.argwhere(body)]
    ids=np.asarray([p[0]*w+p[1] for p in cells],int)
    dmap,dist,levels,groups,counts=shell_structure(m,cells,ids)
    dmax=max(float(dist.max()),1.0)
    hand=np.zeros_like(dmap,float);hand[body]=float(g0)*(1+(float(ratio)-1)*dmap[body]/dmax)
    total=float(hand[body].sum());cap=float(g0)*float(ratio)
    uniform=np.zeros_like(hand,float);uniform[body]=total/len(ids)
    sites=injection_sites(m)

    full,full_best,full_hist=hml.train_density(m,L,uniform,omegas,tau,mu,sites,ids,total,cap,full_steps,step_fraction)
    radial,rad_best,rad_hist,shell_vals=train_shells(m,L,uniform,omegas,tau,mu,sites,ids,groups,counts,total,cap,radial_steps,step_fraction)

    ev=dict(
        uniform=hml.eval_profile(m,L,uniform,omegas,tau,mu,sites),
        hand=hml.eval_profile(m,L,hand,omegas,tau,mu,sites),
        radial=hml.eval_profile(m,L,radial,omegas,tau,mu,sites),
        full=hml.eval_profile(m,L,full,omegas,tau,mu,sites),
    )
    return dict(seed=int(m.cfg.seed),cells=int(body.sum()),shells=len(levels),
                budget=total,cap=cap,evaluation=ev,
                full_best=float(full_best),radial_best=float(rad_best),
                full_history=full_hist,radial_history=rad_hist,
                radial_profile=[dict(distance=int(q),cells=int(n),density=float(v)) for q,n,v in zip(levels,counts,shell_vals)])


def summarize(rows):
    def mean(name,key='coherence2'):
        return float(np.mean([r['evaluation'][name][key] for r in rows]))
    dif=[r['evaluation']['full']['coherence2']-r['evaluation']['radial']['coherence2'] for r in rows]
    return dict(
        bodies=len(rows),
        uniform_R2=mean('uniform'),hand_R2=mean('hand'),radial_R2=mean('radial'),full_R2=mean('full'),
        full_minus_radial=float(np.mean(dif)),full_beats_radial=int(sum(x>0 for x in dif)),
        radial_minus_uniform=float(np.mean([r['evaluation']['radial']['coherence2']-r['evaluation']['uniform']['coherence2'] for r in rows])),
        full_minus_uniform=float(np.mean([r['evaluation']['full']['coherence2']-r['evaluation']['uniform']['coherence2'] for r in rows])),
        full_phase_rms=mean('full','mean_phase_rms'),radial_phase_rms=mean('radial','mean_phase_rms'),
        full_amp=mean('full','mean_median_amp'),radial_amp=mean('radial','mean_median_amp'),
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=622)
    ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2.0);ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10.0)
    ap.add_argument('--full-steps',type=int,default=50);ap.add_argument('--radial-steps',type=int,default=120)
    ap.add_argument('--step-fraction',type=float,default=.10)
    ap.add_argument('--out',default='runs/material_radial_vs_full/dev_622_627.json')
    return ap.parse_args()


def main():
    a=parse_args();fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    oms=[float(x) for x in a.omegas.split(',')];rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,oms,a.tau,a.mu,a.g0,a.ratio,a.full_steps,a.radial_steps,a.step_fraction);rows.append(r)
        e=r['evaluation']
        print(f"seed {seed}: uniform={e['uniform']['coherence2']:.4f} radial={e['radial']['coherence2']:.4f} full={e['full']['coherence2']:.4f} delta={e['full']['coherence2']-e['radial']['coherence2']:+.4f} shells={r['shells']}")
    out=dict(config=vars(a),rows=rows,summary=summarize(rows))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nRADIAL VS FULL MATERIAL')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
