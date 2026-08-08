"""Descriptive ablation on the already-opened matched-tuner benchmark bodies.

Which oracle coordinates make F8 beat local bonds: pole frequencies, residues, or both?
No fresh confirmatory claim is made here.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
import analog_frontier_learning as afl
from matched_tuner_audit import modal_setup,aggregate_modal,train_graph_any,eval_graph


def train_allowed(cfg,lam,bA,bB,c,lags,steps,P,eta,iters,allowed):
    N=len(lam);names='FABC';T=[np.zeros(N) for _ in range(4)]
    z=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,*T,True);C0=z[0];G=list(z[1:5])
    cand=[]
    for j,t in enumerate(names):
        if t not in allowed:continue
        for i,g in enumerate(G[j]):cand.append((abs(float(g)),t,i,float(g)))
    cand.sort(reverse=True,key=lambda q:q[0]);coords=[(t,i,g) for _,t,i,g in cand[:P]]
    traj=[C0]
    for _ in range(iters):
        z=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,*T,True);G=list(z[1:5])
        gg=np.asarray([G[names.index(t)][i] for t,i,_ in coords],float);ss=afl.normalized_step(gg,eta)
        for (t,i,_),s in zip(coords,ss):
            j=names.index(t);lo,hi=(-.5,.5) if t=='F' else (-.9,1.0);T[j][i]=np.clip(T[j][i]+s,lo,hi)
        C2,_=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,*T,False);traj.append(C2)
    return T,coords,traj


def one(m,train,test,steps,P,eta,iters):
    G=train_graph_any(m,train,steps,P,eta,iters);lam,bA,bB,c=modal_setup(m,G['base_wh'],G['base_wv']);N=len(lam)
    base,_=eval_graph(m,G['base_wh'],G['base_wv'],test,steps);gtest,_=eval_graph(m,G['wh'],G['wv'],test,steps)
    arms={}
    for label,allowed in [('F','F'),('ABC','ABC'),('FC','FC'),('FABC','FABC')]:
        T,coords,traj=train_allowed(m.cfg,lam,bA,bB,c,train,steps,P,eta,iters,set(allowed))
        te,_=aggregate_modal(m.cfg,lam,bA,bB,c,test,steps,*T,False)
        tr,_=aggregate_modal(m.cfg,lam,bA,bB,c,train,steps,*T,False)
        arms[label]=dict(train=float(tr),test=float(te),delta_test=float(te-base),coords=[dict(type=t,index=int(i),initial_gradient=float(g)) for t,i,g in coords],traj=[float(x) for x in traj])
    return dict(seed=int(m.cfg.seed),base_test=float(base),graph_test=float(gtest),graph_delta_test=float(gtest-base),arms=arms)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=288);ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--test-lags',default='14,18,22,26');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=40);ap.add_argument('--out',default='runs/free_coordinate_ablation/ablation.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        assert set('FC').issubset(set('FABC'));print('selftest ok');return
    train=[int(x) for x in a.train_lags.split(',')];test=[int(x) for x in a.test_lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,train,test,a.steps,a.P,a.eta,a.iterations);rows.append(r);print('seed',seed,'G',round(r['graph_test'],3),{k:round(v['test'],3) for k,v in r['arms'].items()},flush=True)
    summary=dict(bodies=len(rows),mean_base_test=float(np.mean([r['base_test'] for r in rows])),mean_graph_test=float(np.mean([r['graph_test'] for r in rows])))
    for label in ['F','ABC','FC','FABC']:
        summary[label]=dict(mean_test=float(np.mean([r['arms'][label]['test'] for r in rows])),mean_delta_test=float(np.mean([r['arms'][label]['delta_test'] for r in rows])),
                            beats_graph=int(np.sum([r['arms'][label]['test']>r['graph_test'] for r in rows])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='free_coordinate_ablation_descriptive_v01',summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
