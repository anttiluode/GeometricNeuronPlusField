"""Development probe: can smaller local trust steps plus more relinearization passes
recover the performance of a larger-step spectral optimizer at equal cumulative
nominal task-space trust budget?

All graph schedules use the SAME 8 coordinates selected for delta=.01 at base.
Free reference is F8 at delta=.01 for 40 iterations. Reused bodies only.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from matched_tuner_audit import modal_setup,weights_from_edge_rho,eval_graph
from matched_tuner_trust import (choose_graph,graph_feature_jacobian,scale_invariant_step,
                                 train_free)

SCHEDULES=((.01,40),(.005,80),(.0025,160))


def train_graph_fixed(m,train_lags,test_lags,steps,edges,rho0,delta,iters):
    bwh,bwv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb);rho=np.asarray(rho0,float).copy();traj=[]
    for it in range(iters+1):
        wh,wv=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
        C,J=graph_feature_jacobian(m,wh,wv,train_lags,steps,edges);traj.append(float(np.mean(C)))
        if it==iters:break
        d,_=scale_invariant_step(J,delta);rho=np.clip(rho+d,0.,1.)
    wh,wv=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
    tr,_=eval_graph(m,wh,wv,train_lags,steps);te,_=eval_graph(m,wh,wv,test_lags,steps)
    return dict(delta=float(delta),iterations=int(iters),nominal_budget=float(delta*iters),train=float(tr),test=float(te),rho=[float(x) for x in rho],traj=traj)


def one(m,train_lags,test_lags,steps,P):
    bwh,bwv=ae.bond_weights(m,m.body);base_train,_=eval_graph(m,bwh,bwv,train_lags,steps);base_test,_=eval_graph(m,bwh,bwv,test_lags,steps)
    edges,rho0,meta,pool=choose_graph(m,bwh,bwv,train_lags,steps,P,.01)
    if len(edges)<P:raise RuntimeError('insufficient graph coordinates')
    graph=[train_graph_fixed(m,train_lags,test_lags,steps,edges,rho0,d,it) for d,it in SCHEDULES]
    lam,bA,bB,c=modal_setup(m,bwh,bwv);free=train_free(m,lam,bA,bB,c,train_lags,test_lags,steps,P,.01,40)
    return dict(seed=int(m.cfg.seed),base_train=float(base_train),base_test=float(base_test),graph_coords=[list(e) for e in edges],graph_meta=meta,graph=graph,
                free_reference=dict(delta=.01,iterations=40,nominal_budget=.4,train=float(free['train']),test=float(free['test'])))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--test-lags',default='14,18,22,26');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--out',default='runs/pass_budget_dev/pass_budget.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        assert all(abs(d*i-.4)<1e-12 for d,i in SCHEDULES);print('selftest ok');return
    train=[int(x) for x in a.train_lags.split(',')];test=[int(x) for x in a.test_lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,train,test,a.steps,a.P);rows.append(r)
        print('seed',seed,'base',round(r['base_test'],3),'F40',round(r['free_reference']['test'],3),'G',[(q['iterations'],round(q['test'],3)) for q in r['graph']],flush=True)
    summary=dict(bodies=len(rows),P=a.P,train_lags=train,test_lags=test,free_reference=dict(delta=.01,iterations=40,nominal_budget=.4,
                 mean_test=float(np.mean([r['free_reference']['test'] for r in rows])),mean_delta_test=float(np.mean([r['free_reference']['test']-r['base_test'] for r in rows]))))
    summary['graph_schedules']=[]
    for j,(d,it) in enumerate(SCHEDULES):
        te=np.asarray([r['graph'][j]['test'] for r in rows]);base=np.asarray([r['base_test'] for r in rows]);fr=np.asarray([r['free_reference']['test'] for r in rows])
        summary['graph_schedules'].append(dict(delta=d,iterations=it,nominal_budget=d*it,mean_test=float(te.mean()),mean_delta_test=float(np.mean(te-base)),mean_graph_minus_free=float(np.mean(te-fr)),graph_beats_free=int(np.sum(te>fr))))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='pass_budget_dev_v01',summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
