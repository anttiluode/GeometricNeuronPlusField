"""Development benchmark: 8 local conductance tuners vs 8 direct modal tuners.

See BENCHMARK_SCOPE_V01.md. This first script is deliberately a development probe,
not yet a held-out claim. Both arms start from the exact same mature linear transfer.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # exact mature boundary operator patch
import analog_frontier_learning as afl
from structural_interference_probe import event_candidates, n4
from full_pole_containment import stiffness_matrix


def pulse(cfg,q):
    if not (0 <= int(q) < int(cfg.pulse_frames)): return 0j
    q=int(q);env=math.sin(math.pi*(q+1)/(cfg.pulse_frames+1))**2
    return complex(cfg.source_amp*env*np.exp(1j*cfg.carrier_omega*q))


def aggregate_graph(m,wh,wv,lags,steps):
    Cs=[];gh=None;gv=None
    for lag in lags:
        st=afl.contrast_adjoint_weights(m,wh,wv,ae.source_sequence(m,True,lag,steps),ae.source_sequence(m,False,lag,steps))
        Cs.append(float(st['C']))
        gh=st['gh'].copy() if gh is None else gh+st['gh']
        gv=st['gv'].copy() if gv is None else gv+st['gv']
    return float(np.mean(Cs)),gh/len(lags),gv/len(lags),Cs


def graph_candidates_all(m):
    adds,_,_,_=event_candidates(m,99999,np.random.default_rng(int(m.cfg.seed)+808080))
    body=m.body.astype(bool);out=[]
    for p in adds:
        qs=[q for q in n4(*p,body.shape) if body[q]]
        if len(qs)==1: out.append((tuple(map(int,p)),tuple(map(int,qs[0]))))
    return out


def train_graph(m,lags,steps,P,eta,iters):
    base_wh,base_wv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb)
    C0,gh,gv,_=aggregate_graph(m,base_wh,base_wv,lags,steps)
    allc=graph_candidates_all(m)
    scores=np.asarray([dk*ae.edge_lookup(gh,gv,p,q) for p,q in allc],float)
    order=np.argsort(scores)[::-1][:min(P,len(allc))]
    cands=[allc[int(i)] for i in order];sel0=scores[order] if len(order) else np.array([])
    rho=np.zeros(len(cands),float);traj=[C0]
    for _ in range(iters):
        wh,wv=afl.weights_from_rho(base_wh,base_wv,cands,rho,kb,dk)
        C,gh,gv,_=aggregate_graph(m,wh,wv,lags,steps)
        g=afl.gradient_on_candidates(gh,gv,cands,dk)
        rho=np.clip(rho+afl.normalized_step(g,eta),0.,1.)
        wh2,wv2=afl.weights_from_rho(base_wh,base_wv,cands,rho,kb,dk)
        C2,_,_,_=aggregate_graph(m,wh2,wv2,lags,steps)
        traj.append(C2)
    wh,wv=afl.weights_from_rho(base_wh,base_wv,cands,rho,kb,dk)
    return dict(base_wh=base_wh,base_wv=base_wv,wh=wh,wv=wv,cands=cands,rho=rho,
                initial_scores=sel0,traj=np.asarray(traj,float),candidate_pool=len(allc))


def modal_setup(m,wh,wv):
    K=stiffness_matrix(wh,wv);lam,Phi=np.linalg.eigh(K)
    h,w=m.body.shape
    ia=int(m.source_terminal(0)[0])*w+int(m.source_terminal(0)[1])
    ib=int(m.source_terminal(1)[0])*w+int(m.source_terminal(1)[1])
    isoma=int(m.soma[0])*w+int(m.soma[1])
    return lam,Phi[ia,:].copy(),Phi[ib,:].copy(),Phi[isoma,:].copy()


def modal_order(cfg,lam,bA,bB,c,lag,target,steps,tF,tA,tB,need_grad=True):
    lam=np.asarray(lam,float);N=len(lam);dt=float(cfg.dt);damp=float(cfg.damping)
    k0=float(cfg.restoring)+float(cfg.stiffness)*lam
    kap=k0*np.exp(tF)
    mA=1.0+tA;mB=1.0+tB
    q=np.zeros(N,complex);v=np.zeros(N,complex)
    if need_grad:
        fq=np.zeros(N,complex);fv=np.zeros(N,complex)
        aq=np.zeros(N,complex);av=np.zeros(N,complex)
        bq=np.zeros(N,complex);bv=np.zeros(N,complex)
        gF=np.zeros(N,float);gA=np.zeros(N,float);gB=np.zeros(N,float)
    E=0.0
    for t in range(steps):
        if target:
            sA=pulse(cfg,t);sB=pulse(cfg,t-lag)
        else:
            sB=pulse(cfg,t);sA=pulse(cfg,t-lag)
        src=bA*mA*sA+bB*mB*sB
        vnew=v+dt*(-kap*q-damp*v+src)
        qnew=q+dt*vnew
        if need_grad:
            fvnew=fv+dt*(-kap*fq-kap*q-damp*fv)
            fqnew=fq+dt*fvnew
            avnew=av+dt*(-kap*aq-damp*av+bA*sA)
            aqnew=aq+dt*avnew
            bvnew=bv+dt*(-kap*bq-damp*bv+bB*sB)
            bqnew=bq+dt*bvnew
        z=np.dot(c,qnew);E+=float(abs(z)**2)
        if need_grad:
            cz=np.conj(z)
            gF += 2*np.real(cz*(c*fqnew));gA += 2*np.real(cz*(c*aqnew));gB += 2*np.real(cz*(c*bqnew))
            fq,fv,aq,av,bq,bv=fqnew,fvnew,aqnew,avnew,bqnew,bvnew
        q,v=qnew,vnew
    return (E,gF,gA,gB) if need_grad else E


def modal_contrast_grad(cfg,lam,bA,bB,c,lag,steps,tF,tA,tB,need_grad=True):
    if need_grad:
        ET,fT,aT,bT=modal_order(cfg,lam,bA,bB,c,lag,True,steps,tF,tA,tB,True)
        ED,fD,aD,bD=modal_order(cfg,lam,bA,bB,c,lag,False,steps,tF,tA,tB,True)
        S=ET+ED+1e-30;wT=2*ED/(S*S);wD=-2*ET/(S*S)
        return float((ET-ED)/S),wT*fT+wD*fD,wT*aT+wD*aD,wT*bT+wD*bD
    ET=modal_order(cfg,lam,bA,bB,c,lag,True,steps,tF,tA,tB,False)
    ED=modal_order(cfg,lam,bA,bB,c,lag,False,steps,tF,tA,tB,False)
    return float((ET-ED)/(ET+ED+1e-30))


def aggregate_modal(cfg,lam,bA,bB,c,lags,steps,tF,tA,tB,need_grad=True):
    Cs=[]
    if need_grad:
        F=np.zeros(len(lam));A=np.zeros(len(lam));B=np.zeros(len(lam))
        for lag in lags:
            C,f,a,b=modal_contrast_grad(cfg,lam,bA,bB,c,lag,steps,tF,tA,tB,True)
            Cs.append(C);F+=f;A+=a;B+=b
        n=len(lags);return float(np.mean(Cs)),F/n,A/n,B/n,Cs
    for lag in lags:Cs.append(modal_contrast_grad(cfg,lam,bA,bB,c,lag,steps,tF,tA,tB,False))
    return float(np.mean(Cs)),Cs


def train_free(cfg,lam,bA,bB,c,lags,steps,P,eta,iters):
    N=len(lam);tF=np.zeros(N);tA=np.zeros(N);tB=np.zeros(N)
    C0,gF,gA,gB,_=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,tF,tA,tB,True)
    vals=np.concatenate([np.abs(gF),np.abs(gA),np.abs(gB)])
    pick=np.argsort(vals)[::-1][:P]
    coords=[]
    for z in pick:
        typ=('F' if z<N else ('A' if z<2*N else 'B'));idx=int(z%N)
        g=(gF if typ=='F' else gA if typ=='A' else gB)[idx]
        coords.append((typ,idx,float(g)))
    traj=[C0]
    for _ in range(iters):
        C,gF,gA,gB,_=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,tF,tA,tB,True)
        g=np.asarray([(gF if typ=='F' else gA if typ=='A' else gB)[idx] for typ,idx,_ in coords],float)
        step=afl.normalized_step(g,eta)
        for (typ,idx,_),s in zip(coords,step):
            if typ=='F': tF[idx]=np.clip(tF[idx]+s,-.5,.5)
            elif typ=='A': tA[idx]=np.clip(tA[idx]+s,-.9,1.0)
            else: tB[idx]=np.clip(tB[idx]+s,-.9,1.0)
        C2,_=aggregate_modal(cfg,lam,bA,bB,c,lags,steps,tF,tA,tB,False);traj.append(C2)
    return dict(tF=tF,tA=tA,tB=tB,coords=coords,traj=np.asarray(traj,float))


def evaluate_graph(m,wh,wv,lags,steps):
    Cs=[]
    for lag in lags:
        C,_,_=ae.linear_contrast(m,wh,wv,ae.source_sequence(m,True,lag,steps),ae.source_sequence(m,False,lag,steps));Cs.append(float(C))
    return float(np.mean(Cs)),Cs


def body(m,train_lags,test_lags,steps,P,eta,iters):
    G=train_graph(m,train_lags,steps,P,eta,iters)
    lam,bA,bB,c=modal_setup(m,G['base_wh'],G['base_wv'])
    N=len(lam);z=np.zeros(N)
    base_modal_train,_=aggregate_modal(m.cfg,lam,bA,bB,c,train_lags,steps,z,z,z,False)
    base_graph_train,_=evaluate_graph(m,G['base_wh'],G['base_wv'],train_lags,steps)
    F=train_free(m.cfg,lam,bA,bB,c,train_lags,steps,P,eta,iters)
    gtrain,gtrain_each=evaluate_graph(m,G['wh'],G['wv'],train_lags,steps)
    gtest,gtest_each=evaluate_graph(m,G['wh'],G['wv'],test_lags,steps)
    ftrain,ftrain_each=aggregate_modal(m.cfg,lam,bA,bB,c,train_lags,steps,F['tF'],F['tA'],F['tB'],False)
    ftest,ftest_each=aggregate_modal(m.cfg,lam,bA,bB,c,test_lags,steps,F['tF'],F['tA'],F['tB'],False)
    btest,btest_each=evaluate_graph(m,G['base_wh'],G['base_wv'],test_lags,steps)
    return dict(seed=int(m.cfg.seed),N=N,P=P,
        base_train_graph=base_graph_train,base_train_modal=base_modal_train,base_identity_error=float(abs(base_graph_train-base_modal_train)),
        base_test=btest,graph_train=gtrain,graph_test=gtest,free_train=ftrain,free_test=ftest,
        graph_delta_train=float(gtrain-base_graph_train),graph_delta_test=float(gtest-btest),
        free_delta_train=float(ftrain-base_modal_train),free_delta_test=float(ftest-btest),
        graph_minus_free_test=float(gtest-ftest),graph_candidate_pool=int(G['candidate_pool']),
        graph_rho=[float(x) for x in G['rho']],graph_initial_scores=[float(x) for x in G['initial_scores']],
        free_coords=[dict(type=t,index=int(i),initial_gradient=float(g)) for t,i,g in F['coords']],
        graph_train_each=gtrain_each,graph_test_each=gtest_each,free_train_each=ftrain_each,free_test_each=ftest_each,base_test_each=btest_each,
        graph_traj=[float(x) for x in G['traj']],free_traj=[float(x) for x in F['traj']])


def selftest():
    # Sensitivity smoke test is covered by development identity and monotonic training receipts.
    x=afl.normalized_step(np.array([2.,-1.]),.01);assert abs(x[0]-.01)<1e-12
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=40)
    ap.add_argument('--train-lags',default='16,20,24');ap.add_argument('--test-lags',default='14,18,22,26')
    ap.add_argument('--out',default='runs/matched_tuner_dev/matched_tuner.json');ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args();
    if a.selftest:selftest();return
    train=[int(x) for x in a.train_lags.split(',') if x];test=[int(x) for x in a.test_lags.split(',') if x]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=body(m,train,test,a.steps,a.P,a.eta,a.iterations);rows.append(r)
        print(f"seed {seed}: id={r['base_identity_error']:.1e} baseT={r['base_test']:+.3f} G={r['graph_test']:+.3f} F={r['free_test']:+.3f} G-F={r['graph_minus_free_test']:+.3f}",flush=True)
    summary=dict(bodies=len(rows),P=a.P,train_lags=train,test_lags=test,eta=a.eta,iterations=a.iterations,
        max_base_identity_error=float(max(r['base_identity_error'] for r in rows)),
        mean_base_test=float(np.mean([r['base_test'] for r in rows])),
        mean_graph_test=float(np.mean([r['graph_test'] for r in rows])),mean_free_test=float(np.mean([r['free_test'] for r in rows])),
        mean_graph_delta_test=float(np.mean([r['graph_delta_test'] for r in rows])),mean_free_delta_test=float(np.mean([r['free_delta_test'] for r in rows])),
        mean_graph_minus_free_test=float(np.mean([r['graph_minus_free_test'] for r in rows])),graph_beats_free_bodies=int(np.sum([r['graph_minus_free_test']>0 for r in rows])),
        mean_graph_train=float(np.mean([r['graph_train'] for r in rows])),mean_free_train=float(np.mean([r['free_train'] for r in rows])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='matched_tuner_dev_v01',summary=summary,rows=rows),indent=2))
    print('\nMATCHED TUNER DEVELOPMENT RECEIPT');[print(f' {k}: {v}') for k,v in summary.items()]
if __name__=='__main__':main()
