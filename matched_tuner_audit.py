"""Development symmetry audit for the matched-tuner benchmark.

Strengthens BOTH arms relative to v0.1:
- graph: choose 8 feasible tuners from every local bond in the exact 31x31 medium;
- free modal: choose 8 from pole, source-A, source-B, and soma/output residues.
Reused seeds only. No held-out claim.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
import analog_frontier_learning as afl
from full_pole_containment import stiffness_matrix


def pulse(cfg,q):
    q=int(q)
    if not (0 <= q < int(cfg.pulse_frames)): return 0j
    env=math.sin(math.pi*(q+1)/(cfg.pulse_frames+1))**2
    return complex(cfg.source_amp*env*np.exp(1j*cfg.carrier_omega*q))


def aggregate_graph(m,wh,wv,lags,steps):
    Cs=[];gh=None;gv=None
    for lag in lags:
        st=afl.contrast_adjoint_weights(m,wh,wv,ae.source_sequence(m,True,lag,steps),ae.source_sequence(m,False,lag,steps))
        Cs.append(float(st['C']))
        gh=st['gh'].copy() if gh is None else gh+st['gh']
        gv=st['gv'].copy() if gv is None else gv+st['gv']
    n=len(lags)
    return float(np.mean(Cs)),gh/n,gv/n,Cs


def all_edges(shape):
    h,w=shape;out=[]
    for y in range(h):
        for x in range(w-1):out.append(('h',y,x))
    for y in range(h-1):
        for x in range(w):out.append(('v',y,x))
    return out


def edge_value(A,e):
    typ,y,x=e;return float(A[0][y,x] if typ=='h' else A[1][y,x])


def set_edge_value(wh,wv,e,v):
    typ,y,x=e
    if typ=='h':wh[y,x]=v
    else:wv[y,x]=v


def edge_grad(gh,gv,e):
    typ,y,x=e;return float(gh[y,x] if typ=='h' else gv[y,x])


def weights_from_edge_rho(base_wh,base_wv,edges,rho,kb,dk):
    wh=base_wh.copy();wv=base_wv.copy()
    for e,r in zip(edges,rho):set_edge_value(wh,wv,e,kb+float(r)*dk)
    return wh,wv


def train_graph_any(m,lags,steps,P,eta,iters):
    bwh,bwv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);ka=float(m.cfg.k_arbor);dk=ka-kb
    C0,gh,gv,_=aggregate_graph(m,bwh,bwv,lags,steps)
    pool=all_edges(m.body.shape)
    rho0=np.asarray([(edge_value((bwh,bwv),e)-kb)/dk for e in pool],float)
    g=np.asarray([dk*edge_grad(gh,gv,e) for e in pool],float)
    # projected first-step improvement: weak edges can only increase, strong only decrease
    score=np.where(rho0<.5,np.maximum(g,0.0),np.maximum(-g,0.0))
    pick=np.argsort(score)[::-1][:P]
    edges=[pool[int(i)] for i in pick];rho=rho0[pick].copy();initial_g=g[pick].copy();traj=[C0]
    for _ in range(iters):
        wh,wv=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
        C,gh,gv,_=aggregate_graph(m,wh,wv,lags,steps)
        gg=np.asarray([dk*edge_grad(gh,gv,e) for e in edges],float)
        rho=np.clip(rho+afl.normalized_step(gg,eta),0.,1.)
        wh2,wv2=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
        C2,_,_,_=aggregate_graph(m,wh2,wv2,lags,steps);traj.append(C2)
    wh,wv=weights_from_edge_rho(bwh,bwv,edges,rho,kb,dk)
    return dict(base_wh=bwh,base_wv=bwv,wh=wh,wv=wv,edges=edges,rho0=rho0[pick],rho=rho,
                initial_g=initial_g,traj=np.asarray(traj),pool_n=len(pool))


def modal_setup(m,wh,wv):
    K=stiffness_matrix(wh,wv);lam,Phi=np.linalg.eigh(K);w=m.body.shape[1]
    ia=int(m.source_terminal(0)[0])*w+int(m.source_terminal(0)[1]);ib=int(m.source_terminal(1)[0])*w+int(m.source_terminal(1)[1]);s=int(m.soma[0])*w+int(m.soma[1])
    return lam,Phi[ia,:].copy(),Phi[ib,:].copy(),Phi[s,:].copy()


def modal_order(cfg,lam,bA,bB,c,lag,target,steps,tF,tA,tB,tC,grad=True):
    N=len(lam);dt=float(cfg.dt);damp=float(cfg.damping);kap=(float(cfg.restoring)+float(cfg.stiffness)*lam)*np.exp(tF)
    mA=1+tA;mB=1+tB;ce=c*(1+tC)
    q=np.zeros(N,complex);v=np.zeros(N,complex);E=0.0
    if grad:
        fq=np.zeros(N,complex);fv=np.zeros(N,complex);aq=np.zeros(N,complex);av=np.zeros(N,complex);bq=np.zeros(N,complex);bv=np.zeros(N,complex)
        gF=np.zeros(N);gA=np.zeros(N);gB=np.zeros(N);gC=np.zeros(N)
    for t in range(steps):
        if target:sA=pulse(cfg,t);sB=pulse(cfg,t-lag)
        else:sB=pulse(cfg,t);sA=pulse(cfg,t-lag)
        src=bA*mA*sA+bB*mB*sB
        vn=v+dt*(-kap*q-damp*v+src);qn=q+dt*vn
        if grad:
            fvn=fv+dt*(-kap*fq-kap*q-damp*fv);fqn=fq+dt*fvn
            avn=av+dt*(-kap*aq-damp*av+bA*sA);aqn=aq+dt*avn
            bvn=bv+dt*(-kap*bq-damp*bv+bB*sB);bqn=bq+dt*bvn
        z=np.dot(ce,qn);E+=float(abs(z)**2)
        if grad:
            cz=np.conj(z);gF+=2*np.real(cz*(ce*fqn));gA+=2*np.real(cz*(ce*aqn));gB+=2*np.real(cz*(ce*bqn));gC+=2*np.real(cz*(c*qn))
            fq,fv,aq,av,bq,bv=fqn,fvn,aqn,avn,bqn,bvn
        q,v=qn,vn
    return (E,gF,gA,gB,gC) if grad else E


def modal_C_grad(cfg,lam,bA,bB,c,lag,steps,tF,tA,tB,tC,grad=True):
    if grad:
        T=modal_order(cfg,lam,bA,bB,c,lag,True,steps,tF,tA,tB,tC,True);D=modal_order(cfg,lam,bA,bB,c,lag,False,steps,tF,tA,tB,tC,True)
        ET,ED=T[0],D[0];S=ET+ED+1e-30;wT=2*ED/S**2;wD=-2*ET/S**2
        return ((ET-ED)/S,)+tuple(wT*T[i]+wD*D[i] for i in range(1,5))
    ET=modal_order(cfg,lam,bA,bB,c,lag,True,steps,tF,tA,tB,tC,False);ED=modal_order(cfg,lam,bA,bB,c,lag,False,steps,tF,tA,tB,tC,False)
    return float((ET-ED)/(ET+ED+1e-30))


def aggregate_modal(cfg,lam,bA,bB,c,lags,steps,tF,tA,tB,tC,grad=True):
    Cs=[]
    if grad:
        G=[np.zeros(len(lam)) for _ in range(4)]
        for lag in lags:
            z=modal_C_grad(cfg,lam,bA,bB,c,lag,steps,tF,tA,tB,tC,True);Cs.append(float(z[0]))
            for j in range(4):G[j]+=z[j+1]
        return (float(np.mean(Cs)),)+tuple(g/len(lags) for g in G)+(Cs,)
    for lag in lags:Cs.append(modal_C_grad(cfg,lam,bA,bB,c,lag,steps,tF,tA,tB,tC,False))
    return float(np.mean(Cs)),Cs


def train_free4(cfg,lam,bA,bB,c,lags,steps,P,eta,iters):
    N=len(lam);T=[np.zeros(N) for _ in range(4)]
    z=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,*T,True);C0=z[0];G=list(z[1:5]);names='FABC'
    vals=np.concatenate([np.abs(g) for g in G]);pick=np.argsort(vals)[::-1][:P]
    coords=[]
    for z0 in pick:
        j=int(z0//N);i=int(z0%N);coords.append((names[j],i,float(G[j][i])))
    traj=[C0]
    for _ in range(iters):
        z=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,*T,True);G=list(z[1:5]);gg=np.asarray([G[names.index(t)][i] for t,i,_ in coords])
        ss=afl.normalized_step(gg,eta)
        for (t,i,_),s in zip(coords,ss):
            j=names.index(t);lo,hi=(-.5,.5) if t=='F' else (-.9,1.0);T[j][i]=np.clip(T[j][i]+s,lo,hi)
        C2,_=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,*T,False);traj.append(C2)
    return dict(T=T,coords=coords,traj=np.asarray(traj))


def eval_graph(m,wh,wv,lags,steps):
    Cs=[]
    for lag in lags:
        C,_,_=ae.linear_contrast(m,wh,wv,ae.source_sequence(m,True,lag,steps),ae.source_sequence(m,False,lag,steps));Cs.append(float(C))
    return float(np.mean(Cs)),Cs


def body(m,train,test,steps,P,eta,iters):
    G=train_graph_any(m,train,steps,P,eta,iters);lam,bA,bB,c=modal_setup(m,G['base_wh'],G['base_wv']);Z=[np.zeros(len(lam)) for _ in range(4)]
    bgt,_=eval_graph(m,G['base_wh'],G['base_wv'],train,steps);bmt,_=aggregate_modal(m.cfg,lam,bA,bB,c,train,steps,*Z,False)
    F=train_free4(m.cfg,lam,bA,bB,c,train,steps,P,eta,iters)
    gtr,_=eval_graph(m,G['wh'],G['wv'],train,steps);gte,_=eval_graph(m,G['wh'],G['wv'],test,steps);bte,_=eval_graph(m,G['base_wh'],G['base_wv'],test,steps)
    ftr,_=aggregate_modal(m.cfg,lam,bA,bB,c,train,steps,*F['T'],False);fte,_=aggregate_modal(m.cfg,lam,bA,bB,c,test,steps,*F['T'],False)
    return dict(seed=int(m.cfg.seed),base_identity_error=float(abs(bgt-bmt)),base_test=bte,graph_train=gtr,graph_test=gte,free_train=ftr,free_test=fte,
                graph_delta_test=float(gte-bte),free_delta_test=float(fte-bte),graph_minus_free_test=float(gte-fte),graph_pool=G['pool_n'],
                graph_edges=[list(e) for e in G['edges']],graph_rho0=[float(x) for x in G['rho0']],graph_rho=[float(x) for x in G['rho']],
                graph_initial_g=[float(x) for x in G['initial_g']],free_coords=[dict(type=t,index=i,initial_gradient=g) for t,i,g in F['coords']],
                graph_traj=[float(x) for x in G['traj']],free_traj=[float(x) for x in F['traj']])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=40)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--test-lags',default='14,18,22,26');ap.add_argument('--out',default='runs/matched_tuner_audit/audit.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        assert len(all_edges((31,31)))==1860;print('selftest ok');return
    train=[int(x) for x in a.train_lags.split(',')];test=[int(x) for x in a.test_lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=body(m,train,test,a.steps,a.P,a.eta,a.iterations);rows.append(r);print(f"seed {seed}: G={r['graph_test']:+.3f} F={r['free_test']:+.3f} G-F={r['graph_minus_free_test']:+.3f}",flush=True)
    summary=dict(bodies=len(rows),P=a.P,train_lags=train,test_lags=test,max_base_identity_error=float(max(r['base_identity_error'] for r in rows)),
                 mean_base_test=float(np.mean([r['base_test'] for r in rows])),mean_graph_test=float(np.mean([r['graph_test'] for r in rows])),mean_free_test=float(np.mean([r['free_test'] for r in rows])),
                 mean_graph_delta_test=float(np.mean([r['graph_delta_test'] for r in rows])),mean_free_delta_test=float(np.mean([r['free_delta_test'] for r in rows])),
                 mean_graph_minus_free_test=float(np.mean([r['graph_minus_free_test'] for r in rows])),graph_beats_free=int(np.sum([r['graph_minus_free_test']>0 for r in rows])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='matched_tuner_audit_dev_v01',summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
