"""Coordinate-scale-invariant development benchmark.

Both local-bond and free-modal arms are selected by normalized task-space alignment
and updated with the same predicted RMS change of the multi-lag training-output
vector per iteration. This removes the arbitrary 'eta=0.01 means the same thing in
different coordinates' flaw exposed by SPECTRAL_TO_LOCAL_COMPILER_V01.md.
"""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
import analog_frontier_learning as afl
from full_pole_containment import stiffness_matrix
from matched_tuner_audit import (modal_setup, modal_C_grad, aggregate_modal,
    all_edges, edge_value, edge_grad, set_edge_value, weights_from_edge_rho, eval_graph)

NAMES='FABC'


def rms(x):
    x=np.asarray(x,float);return float(np.sqrt(np.mean(x*x)))


def normalized_alignment(v):
    s=rms(v)
    return float(np.mean(v)/s) if s>1e-14 else -np.inf


def graph_feature_jacobian(m,wh,wv,lags,steps,edges=None):
    if edges is None: edges=all_edges(m.body.shape)
    dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    J=np.zeros((len(lags),len(edges)),float);C=[]
    for r,lag in enumerate(lags):
        st=afl.contrast_adjoint_weights(m,wh,wv,ae.source_sequence(m,True,lag,steps),ae.source_sequence(m,False,lag,steps))
        C.append(float(st['C']))
        for j,e in enumerate(edges):J[r,j]=dk*edge_grad(st['gh'],st['gv'],e)
    return np.asarray(C,float),J


def free_feature_jacobian(m,lam,bA,bB,c,lags,steps,T,coords=None):
    C=[]
    if coords is None:
        J=np.zeros((len(lags),4*len(lam)),float)
    else:
        J=np.zeros((len(lags),len(coords)),float)
    for r,lag in enumerate(lags):
        z=modal_C_grad(m.cfg,lam,bA,bB,c,lag,steps,*T,True);C.append(float(z[0]));G=list(z[1:5])
        if coords is None:
            N=len(lam)
            for j in range(4):J[r,j*N:(j+1)*N]=G[j]
        else:
            for j,(typ,i) in enumerate(coords):J[r,j]=G[NAMES.index(typ)][i]
    return np.asarray(C,float),J


def choose_graph(m,wh,wv,lags,steps,P,delta):
    edges=all_edges(m.body.shape);C,J=graph_feature_jacobian(m,wh,wv,lags,steps,edges)
    kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb)
    scored=[]
    for j,e in enumerate(edges):
        rho=(edge_value((wh,wv),e)-kb)/dk
        # feasible direction at the binary base state
        sign=+1.0 if rho<.5 else -1.0
        v=sign*J[:,j];s=rms(v);avail=(1-rho) if sign>0 else rho
        score=normalized_alignment(v)
        if np.isfinite(score) and score>0 and s*avail>=delta:
            scored.append((score,j,e,float(rho),sign,s))
    scored.sort(reverse=True,key=lambda q:q[0])
    out=scored[:P]
    return [q[2] for q in out],np.asarray([q[3] for q in out],float),[dict(score=q[0],edge=list(q[2]),rho0=q[3],initial_sign=q[4],rms_sensitivity=q[5]) for q in out],len(scored)


def free_bounds(typ):
    return (-.5,.5) if typ=='F' else (-.9,1.0)


def choose_free(m,lam,bA,bB,c,lags,steps,P,delta):
    N=len(lam);T=[np.zeros(N) for _ in range(4)];C,J=free_feature_jacobian(m,lam,bA,bB,c,lags,steps,T,None)
    scored=[]
    for q in range(J.shape[1]):
        typ=NAMES[int(q//N)];i=int(q%N);v=J[:,q];s=rms(v)
        if s<1e-14:continue
        sg=1.0 if np.mean(v)>=0 else -1.0;lo,hi=free_bounds(typ);avail=hi if sg>0 else -lo
        score=normalized_alignment(sg*v)
        if np.isfinite(score) and score>0 and s*avail>=delta:
            scored.append((score,typ,i,sg,s))
    scored.sort(reverse=True,key=lambda q:q[0]);out=scored[:P]
    coords=[(q[1],q[2]) for q in out]
    meta=[dict(score=q[0],type=q[1],index=q[2],initial_sign=q[3],rms_sensitivity=q[4]) for q in out]
    return coords,T,meta,len(scored)


def scale_invariant_step(J,delta):
    """Return dtheta with predicted RMS(J dtheta)=delta before box clipping.
    Invariant to independent positive rescaling of coordinate units."""
    J=np.asarray(J,float);s=np.sqrt(np.mean(J*J,axis=0));ok=s>1e-14
    U=np.zeros_like(J);U[:,ok]=J[:,ok]/s[ok]
    gz=np.mean(U,axis=0);gz[~ok]=0.0
    ng=float(np.linalg.norm(gz))
    if ng<1e-14:return np.zeros(J.shape[1]),0.0
    zdir=gz/ng;pred=U@zdir;rp=rms(pred)
    if rp<1e-14:return np.zeros(J.shape[1]),0.0
    alpha=float(delta/rp);z=alpha*zdir
    d=np.zeros_like(z);d[ok]=z[ok]/s[ok]
    return d,float(rms(J@d))


def train_graph(m,train_lags,test_lags,steps,P,delta,iters):
    bwh,bwv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb)
    edges,rho,meta,pool=choose_graph(m,bwh,bwv,train_lags,steps,P,delta)
    if len(edges)<P:raise RuntimeError(f'only {len(edges)} feasible graph coords')
    traj=[];pred=[]
    for it in range(iters+1):
        wh,wv=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
        C,J=graph_feature_jacobian(m,wh,wv,train_lags,steps,edges);traj.append(float(np.mean(C)))
        if it==iters:break
        d,pr=scale_invariant_step(J,delta);pred.append(pr);rho=np.clip(rho+d,0.,1.)
    wh,wv=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
    train,_=eval_graph(m,wh,wv,train_lags,steps);test,_=eval_graph(m,wh,wv,test_lags,steps)
    return dict(train=train,test=test,edges=[list(e) for e in edges],rho=[float(x) for x in rho],meta=meta,pool=pool,traj=traj,predicted_rms=pred,wh=wh,wv=wv)


def train_free(m,lam,bA,bB,c,train_lags,test_lags,steps,P,delta,iters):
    coords,T,meta,pool=choose_free(m,lam,bA,bB,c,train_lags,steps,P,delta)
    if len(coords)<P:raise RuntimeError(f'only {len(coords)} feasible free coords')
    traj=[];pred=[]
    for it in range(iters+1):
        C,J=free_feature_jacobian(m,lam,bA,bB,c,train_lags,steps,T,coords);traj.append(float(np.mean(C)))
        if it==iters:break
        d,pr=scale_invariant_step(J,delta);pred.append(pr)
        for j,(typ,i) in enumerate(coords):
            ti=NAMES.index(typ);lo,hi=free_bounds(typ);T[ti][i]=np.clip(T[ti][i]+d[j],lo,hi)
    train,_=aggregate_modal(m.cfg,lam,bA,bB,c,train_lags,steps,*T,False);test,_=aggregate_modal(m.cfg,lam,bA,bB,c,test_lags,steps,*T,False)
    return dict(train=float(train),test=float(test),coords=[dict(type=t,index=i) for t,i in coords],T=[[float(x) for x in q] for q in T],meta=meta,pool=pool,traj=traj,predicted_rms=pred)


def one(m,train_lags,test_lags,steps,P,delta,iters):
    bwh,bwv=ae.bond_weights(m,m.body);base_train,_=eval_graph(m,bwh,bwv,train_lags,steps);base_test,_=eval_graph(m,bwh,bwv,test_lags,steps)
    lam,bA,bB,c=modal_setup(m,bwh,bwv);N=len(lam);Z=[np.zeros(N) for _ in range(4)];mbase,_=aggregate_modal(m.cfg,lam,bA,bB,c,train_lags,steps,*Z,False)
    G=train_graph(m,train_lags,test_lags,steps,P,delta,iters);F=train_free(m,lam,bA,bB,c,train_lags,test_lags,steps,P,delta,iters)
    # strip heavy arrays not needed in JSON
    G.pop('wh',None);G.pop('wv',None)
    return dict(seed=int(m.cfg.seed),base_train=float(base_train),base_test=float(base_test),base_identity_error=float(abs(base_train-mbase)),graph=G,free=F,
                graph_delta_test=float(G['test']-base_test),free_delta_test=float(F['test']-base_test),graph_minus_free_test=float(G['test']-F['test']))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--test-lags',default='14,18,22,26');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8)
    ap.add_argument('--delta',type=float,default=.005);ap.add_argument('--iterations',type=int,default=40);ap.add_argument('--out',default='runs/matched_tuner_trust/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        J=np.array([[1.,2.],[2.,4.],[3.,6.]])
        d1,p1=scale_invariant_step(J,.01);d2,p2=scale_invariant_step(J@np.diag([.1,10.]),.01)
        # Mapping second coordinates back to first units reproduces the same task step.
        assert abs(p1-.01)<1e-10 and abs(p2-.01)<1e-10
        assert np.linalg.norm(J@d1-(J@np.diag([.1,10.]))@d2)<1e-10
        print('selftest ok');return
    train=[int(x) for x in a.train_lags.split(',')];test=[int(x) for x in a.test_lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,train,test,a.steps,a.P,a.delta,a.iterations);rows.append(r)
        print(f"seed {seed}: base={r['base_test']:+.3f} G={r['graph']['test']:+.3f} F={r['free']['test']:+.3f} G-F={r['graph_minus_free_test']:+.3f}",flush=True)
    summary=dict(bodies=len(rows),P=a.P,delta=a.delta,iterations=a.iterations,train_lags=train,test_lags=test,max_identity_error=float(max(r['base_identity_error'] for r in rows)),
                 mean_base_test=float(np.mean([r['base_test'] for r in rows])),mean_graph_test=float(np.mean([r['graph']['test'] for r in rows])),mean_free_test=float(np.mean([r['free']['test'] for r in rows])),
                 mean_graph_delta_test=float(np.mean([r['graph_delta_test'] for r in rows])),mean_free_delta_test=float(np.mean([r['free_delta_test'] for r in rows])),mean_graph_minus_free_test=float(np.mean([r['graph_minus_free_test'] for r in rows])),graph_beats_free=int(np.sum([r['graph_minus_free_test']>0 for r in rows])),
                 mean_graph_train=float(np.mean([r['graph']['train'] for r in rows])),mean_free_train=float(np.mean([r['free']['train'] for r in rows])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='matched_tuner_trust_dev_v01',summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
