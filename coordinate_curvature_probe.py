"""Finite-step task-space curvature of local bond vs direct spectral coordinates.

Select the same scale-invariant G8/F8 coordinates at the base state, then move one
coordinate at a time by a parameter amount whose *training-lag* first-order tangent
predicts a fixed RMS trust step. Compare the actual seven-lag output change with
that initial tangent prediction.

Development probe; see MATCHED_TUNER_TRUST_CONFIRM_V01.md for motivation.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from matched_tuner_audit import modal_setup, weights_from_edge_rho, eval_graph
from matched_tuner_trust import (rms, choose_graph, choose_free, graph_feature_jacobian,
                                 free_feature_jacobian, aggregate_modal, NAMES, free_bounds)

DELTAS=(.001,.0025,.005,.01)


def metrics(pred,actual):
    pred=np.asarray(pred,float);actual=np.asarray(actual,float)
    npred=float(np.linalg.norm(pred));nact=float(np.linalg.norm(actual))
    rel=float(np.linalg.norm(actual-pred)/(npred+1e-30))
    cos=float(np.dot(actual,pred)/(nact*npred+1e-30))
    return dict(relative_linearization_error=rel,cosine=cos,
                actual_to_predicted_norm=float(nact/(npred+1e-30)),
                predicted_rms=rms(pred),actual_rms=rms(actual))


def graph_curvature(m,train_lags,eval_lags,steps,P):
    wh,wv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb)
    edges,rho0,meta,pool=choose_graph(m,wh,wv,train_lags,steps,P,max(DELTAS))
    if len(edges)<P:raise RuntimeError('insufficient graph coordinates')
    base_mean,base_eval=eval_graph(m,wh,wv,eval_lags,steps);base_eval=np.asarray(base_eval,float)
    _,Jeval=graph_feature_jacobian(m,wh,wv,eval_lags,steps,edges)
    rows=[]
    for j,(e,r0,md) in enumerate(zip(edges,rho0,meta)):
        sign=float(md['initial_sign']);strain=float(md['rms_sensitivity'])
        rr=[]
        for delta in DELTAS:
            dr=sign*float(delta)/(strain+1e-30)
            r1=float(np.clip(r0+dr,0.,1.));actual_dr=r1-r0
            wh1,wv1=weights_from_edge_rho(wh,wv,[e],np.asarray([r1]),kb,dk)
            _,c1=eval_graph(m,wh1,wv1,eval_lags,steps);actual=np.asarray(c1,float)-base_eval
            pred=Jeval[:,j]*actual_dr
            rr.append(dict(delta=float(delta),parameter_step=float(actual_dr),clipped=bool(abs(actual_dr-dr)>1e-12),**metrics(pred,actual)))
        rows.append(dict(edge=list(e),rho0=float(r0),selection=md,steps=rr))
    return dict(pool=pool,selected=rows)


def free_curvature(m,train_lags,eval_lags,steps,P):
    wh,wv=ae.bond_weights(m,m.body);lam,bA,bB,c=modal_setup(m,wh,wv);N=len(lam)
    coords,T0,meta,pool=choose_free(m,lam,bA,bB,c,train_lags,steps,P,max(DELTAS))
    if len(coords)<P:raise RuntimeError('insufficient free coordinates')
    _,base_eval=aggregate_modal(m.cfg,lam,bA,bB,c,eval_lags,steps,*T0,False);base_eval=np.asarray(base_eval,float)
    _,Jeval=free_feature_jacobian(m,lam,bA,bB,c,eval_lags,steps,T0,coords)
    rows=[]
    for j,((typ,i),md) in enumerate(zip(coords,meta)):
        sign=float(md['initial_sign']);strain=float(md['rms_sensitivity']);lo,hi=free_bounds(typ)
        rr=[]
        for delta in DELTAS:
            dp=sign*float(delta)/(strain+1e-30);T=[x.copy() for x in T0]
            old=float(T[NAMES.index(typ)][i]);new=float(np.clip(old+dp,lo,hi));actual_dp=new-old;T[NAMES.index(typ)][i]=new
            _,c1=aggregate_modal(m.cfg,lam,bA,bB,c,eval_lags,steps,*T,False);actual=np.asarray(c1,float)-base_eval
            pred=Jeval[:,j]*actual_dp
            rr.append(dict(delta=float(delta),parameter_step=float(actual_dp),clipped=bool(abs(actual_dp-dp)>1e-12),**metrics(pred,actual)))
        rows.append(dict(type=typ,index=int(i),selection=md,steps=rr))
    return dict(pool=pool,selected=rows)


def summarize(rows,arm):
    out={}
    for delta in DELTAS:
        q=[]
        for r in rows:
            for c in r[arm]['selected']:
                q.append(next(x for x in c['steps'] if abs(x['delta']-delta)<1e-12))
        err=np.asarray([x['relative_linearization_error'] for x in q],float);cos=np.asarray([x['cosine'] for x in q],float);rat=np.asarray([x['actual_to_predicted_norm'] for x in q],float)
        out[str(delta)]=dict(n=len(q),mean_rel_error=float(err.mean()),median_rel_error=float(np.median(err)),
                             mean_cosine=float(cos.mean()),median_cosine=float(np.median(cos)),
                             mean_norm_ratio=float(rat.mean()),median_norm_ratio=float(np.median(rat)),
                             clipped_n=int(np.sum([x['clipped'] for x in q])))
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--eval-lags',default='14,16,18,20,22,24,26');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8)
    ap.add_argument('--out',default='runs/coordinate_curvature_dev/curvature.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        m=metrics(np.array([1.,2.]),np.array([1.,2.]));assert m['relative_linearization_error']<1e-12 and abs(m['cosine']-1)<1e-12;print('selftest ok');return
    train=[int(x) for x in a.train_lags.split(',')];ev=[int(x) for x in a.eval_lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=dict(seed=int(seed),graph=graph_curvature(m,train,ev,a.steps,a.P),free=free_curvature(m,train,ev,a.steps,a.P));rows.append(r)
        print('seed',seed,'done',flush=True)
    summary=dict(bodies=len(rows),P=a.P,train_lags=train,eval_lags=ev,deltas=list(DELTAS),graph=summarize(rows,'graph'),free=summarize(rows,'free'))
    for d in DELTAS:
        g=summary['graph'][str(d)]['median_rel_error'];f=summary['free'][str(d)]['median_rel_error'];print(f'delta {d:.4g}: median err G={g:.4f} F={f:.4f} ratio={(g/(f+1e-30)):.2f}',flush=True)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='coordinate_curvature_dev_v01',summary=summary,rows=rows),indent=2))
if __name__=='__main__':main()
