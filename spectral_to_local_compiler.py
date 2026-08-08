"""Task-space compiler audit: how many feasible local bond directions approximate one
winning free spectral coordinate across a family of temporal conditions?

This is descriptive on already-opened benchmark bodies. It improves on the raw
operator-support audit by comparing *functional tangent vectors* rather than
matrix entries.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import nnls

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from matched_tuner_audit import modal_setup,aggregate_modal,modal_C_grad,all_edges,edge_value,edge_grad

KS=(1,2,4,8,16,24)


def local_feature_jacobian(m,wh,wv,lags,steps):
    """Rows=lags, cols=all 1860 local bonds, oriented into each bond's feasible
    one-sided direction at the binary base state (bath:+rho, arbor:-rho)."""
    edges=all_edges(m.body.shape);kb=float(m.cfg.k_mature_bath);ka=float(m.cfg.k_arbor);dk=ka-kb
    J=np.zeros((len(lags),len(edges)),float)
    for r,lag in enumerate(lags):
        st=ae.contrast_adjoint_weights(m,wh,wv,ae.source_sequence(m,True,lag,steps),ae.source_sequence(m,False,lag,steps))
        for j,e in enumerate(edges):
            raw=dk*edge_grad(st['gh'],st['gv'],e)
            rho=(edge_value((wh,wv),e)-kb)/dk
            J[r,j]=raw if rho<.5 else -raw
    return edges,J


def free_feature_jacobian(m,lam,bA,bB,c,lags,steps):
    N=len(lam);Z=[np.zeros(N) for _ in range(4)];names='FABC'
    J=np.zeros((len(lags),4*N),float)
    for r,lag in enumerate(lags):
        z=modal_C_grad(m.cfg,lam,bA,bB,c,lag,steps,*Z,True)
        for j in range(4):J[r,j*N:(j+1)*N]=z[j+1]
    return names,J


def select_free_on_train(m,lam,bA,bB,c,train_lags,steps,P):
    N=len(lam);Z=[np.zeros(N) for _ in range(4)];z=aggregate_modal(m.cfg,lam,bA,bB,c,train_lags,steps,*Z,True)
    G=list(z[1:5]);names='FABC';flat=np.concatenate([np.abs(g) for g in G]);pick=np.argsort(flat)[::-1][:P]
    out=[]
    for q in pick:
        j=int(q//N);i=int(q%N);g=float(G[j][i]);out.append(dict(type=names[j],index=i,flat=int(q),train_gradient=g,improve_sign=1.0 if g>=0 else -1.0))
    return out


def greedy_nnls(A,target,ks=KS):
    """Greedy nonnegative sparse approximation. A columns are already oriented
    feasible local directions. Column normalization is used only for candidate
    selection; NNLS fits the original physical derivatives."""
    y=np.asarray(target,float);yn=float(np.linalg.norm(y))
    if yn<1e-15:return {str(k):dict(rel_residual=float('nan'),selected=0,l1_coeff=float('nan'),max_coeff=float('nan')) for k in ks}
    chosen=[];coef=np.zeros(0);res=y.copy();out={};maxk=max(ks)
    norms=np.linalg.norm(A,axis=0)+1e-30
    for t in range(1,maxk+1):
        corr=(A.T@res)/norms
        corr[np.asarray([j in chosen for j in range(A.shape[1])])]=-np.inf
        j=int(np.argmax(corr))
        if not np.isfinite(corr[j]) or corr[j]<=0:break
        chosen.append(j)
        coef,_=nnls(A[:,chosen],y)
        res=y-A[:,chosen]@coef
        if t in ks:
            out[str(t)]=dict(rel_residual=float(np.linalg.norm(res)/yn),selected=len(chosen),
                             l1_coeff=float(np.sum(coef)),max_coeff=float(np.max(coef) if len(coef) else 0.0),indices=[int(x) for x in chosen])
    last=dict(rel_residual=float(np.linalg.norm(res)/yn),selected=len(chosen),l1_coeff=float(np.sum(coef)) if len(coef) else 0.0,
              max_coeff=float(np.max(coef)) if len(coef) else 0.0,indices=[int(x) for x in chosen])
    for k in ks:
        out.setdefault(str(k),last.copy())
    return out


def one(m,train_lags,eval_lags,steps,P):
    wh,wv=ae.bond_weights(m,m.body);lam,bA,bB,c=modal_setup(m,wh,wv)
    selected=select_free_on_train(m,lam,bA,bB,c,train_lags,steps,P)
    edges,JL=local_feature_jacobian(m,wh,wv,eval_lags,steps);names,JF=free_feature_jacobian(m,lam,bA,bB,c,eval_lags,steps)
    rows=[]
    for q in selected:
        y=q['improve_sign']*JF[:,q['flat']]
        comp=greedy_nnls(JL,y)
        rows.append(dict(type=q['type'],index=q['index'],train_gradient=q['train_gradient'],target_norm=float(np.linalg.norm(y)),compiler=comp))
    return dict(seed=int(m.cfg.seed),local_candidates=len(edges),eval_features=len(eval_lags),selected=rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=288);ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--eval-lags',default='8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32')
    ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--out',default='runs/spectral_to_local_compiler/compiler.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        A=np.eye(3);r=greedy_nnls(A,np.array([1.,2.,0.]),(1,2));assert r['2']['rel_residual']<1e-12;print('selftest ok');return
    train=[int(x) for x in a.train_lags.split(',')];ev=[int(x) for x in a.eval_lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,train,ev,a.steps,a.P);rows.append(r)
        med={k:round(float(np.median([x['compiler'][str(k)]['rel_residual'] for x in r['selected']])),3) for k in KS};print('seed',seed,med,flush=True)
    sel=[x for r in rows for x in r['selected']]
    summary=dict(bodies=len(rows),selected_coordinates=len(sel),eval_features=len(ev),local_candidates=rows[0]['local_candidates'] if rows else 0)
    for k in KS:
        rr=np.asarray([x['compiler'][str(k)]['rel_residual'] for x in sel],float);cc=np.asarray([x['compiler'][str(k)]['l1_coeff'] for x in sel],float)
        summary[str(k)]=dict(mean_rel_residual=float(np.nanmean(rr)),median_rel_residual=float(np.nanmedian(rr)),fraction_residual_lt_0_1=float(np.nanmean(rr<.1)),
                             fraction_residual_lt_0_2=float(np.nanmean(rr<.2)),median_l1_coeff=float(np.nanmedian(cc)))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='spectral_to_local_compiler_dev_v01',train_lags=train,eval_lags=ev,summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
