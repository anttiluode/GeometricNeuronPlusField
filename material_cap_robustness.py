"""Constraint robustness for emergent material-distance organization.

Reuse the exact phase-coherence objective and material adjoint from
hcn_material_learning.py, but vary the local density cap while keeping:

- the frozen arbor,
- total material budget,
- frequencies,
- quasi-active kinetics,
- uniform initialization,
- optimizer and step rule

fixed.

The total budget remains the budget of the confirmed hand profile.  Graph
distance never enters training.  We ask after training whether positive density
vs readout-distance correlation and objective gain survive when the bang-bang
boundary is made tighter or looser.
"""
from __future__ import annotations

import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import weighted_laplacian_sparse,injection_sites
from hcn_material_learning import train_density,coherence_and_gradient,eval_profile


def probe_body(m,caps,omegas,tau,mu,g0,ratio,steps,step_fraction):
    body=m.body.astype(bool);h,w=body.shape
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    cells=[tuple(map(int,p)) for p in np.argwhere(body)]
    ids=np.asarray([p[0]*w+p[1] for p in cells],int)
    dist=m.graph_distance_from_soma().astype(float);dmax=max(float(dist[body].max()),1)
    hand=np.zeros_like(dist,float);hand[body]=float(g0)*(1+(float(ratio)-1)*dist[body]/dmax)
    total=float(hand[body].sum());sites=injection_sites(m)
    rows=[]
    for cap in caps:
        if float(cap)*len(ids)<total-1e-12:
            rows.append(dict(cap=float(cap),feasible=False));continue
        uniform=np.zeros_like(dist,float);uniform[body]=total/len(ids)
        start,_=coherence_and_gradient(m,L,uniform,omegas,tau,mu,sites,ids,False)
        learned,best,hist=train_density(m,L,uniform,omegas,tau,mu,sites,ids,total,float(cap),steps,step_fraction)
        den=learned.ravel()[ids];xd=np.asarray([dist[p] for p in cells],float)
        rho=float(spearmanr(xd,den).statistic)
        pear=float(np.corrcoef(xd,den)[0,1])
        q25,q75=np.quantile(xd,[.25,.75]);prox=xd<=q25;far=xd>=q75
        rows.append(dict(cap=float(cap),feasible=True,uniform_coherence2=float(start),learned_coherence2=float(best),
                         gain=float(best-start),spearman_distance=rho,pearson_distance=pear,
                         proximal_mean=float(np.mean(den[prox])),distal_mean=float(np.mean(den[far])),
                         zero_fraction=float(np.mean(den<1e-8)),cap_fraction=float(np.mean(den>float(cap)-1e-8)),
                         steps_completed=len(hist)-1))
    return dict(seed=int(m.cfg.seed),budget=total,uniform_density=total/len(ids),rows=rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=580);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--caps',default='0.035,0.05,0.075,0.10');ap.add_argument('--omegas',default='0.03,0.04');ap.add_argument('--tau',type=float,default=2);ap.add_argument('--mu',type=float,default=.5);ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10)
    ap.add_argument('--steps',type=int,default=50);ap.add_argument('--step-fraction',type=float,default=.10);ap.add_argument('--out',default='runs/material_cap_robustness/dev.json');a=ap.parse_args()
    caps=[float(x) for x in a.caps.split(',')];omegas=[float(x) for x in a.omegas.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=probe_body(m,caps,omegas,a.tau,a.mu,a.g0,a.ratio,a.steps,a.step_fraction);rows.append(r)
        print('seed',seed,[(x['cap'],round(x.get('gain',0),4),round(x.get('spearman_distance',0),3)) for x in r['rows'] if x.get('feasible')],flush=True)
    bycap={}
    for cap in caps:
        z=[x for r in rows for x in r['rows'] if x.get('feasible') and abs(x['cap']-cap)<1e-12]
        if not z:continue
        bycap[str(cap)]=dict(n=len(z),mean_gain=float(np.mean([x['gain'] for x in z])),mean_spearman=float(np.mean([x['spearman_distance'] for x in z])),
                             min_spearman=float(np.min([x['spearman_distance'] for x in z])),positive_spearman=int(sum(x['spearman_distance']>0 for x in z)),
                             mean_zero_fraction=float(np.mean([x['zero_fraction'] for x in z])),mean_cap_fraction=float(np.mean([x['cap_fraction'] for x in z])))
    summary=dict(bodies=len(rows),bycap=bycap,all_feasible_positive_distance=bool(all(v['positive_spearman']==v['n'] for v in bycap.values())))
    payload=dict(experiment='material_cap_robustness_dev_v01',development_only=True,summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
