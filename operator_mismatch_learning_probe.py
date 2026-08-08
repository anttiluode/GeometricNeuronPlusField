"""Closed-loop broadband learning when the retro pass uses a mismatched operator.

Development probe after operator_mismatch_probe.py located a map-degradation wall.
Each arm has one fixed fractional edge-mismatch pattern for the whole training run.
The forward computation uses the current nominal conductances.  The retro physical
adjoint replay uses current_conductance * (1 + fixed_edge_error), clipped positive.
Thus each pass remains reciprocal internally, but the forward and backward operators
are not identical.

Arms:
  exact       exact discrete adjoint
  K8_ideal    8-bin boundary-selected compressed physical gradient, no mismatch
  K8_m20      K=8, sigma=0.20 fixed fractional retro-pass edge mismatch
  K8_m30      K=8, sigma=0.30
  K16_m30     K=16, sigma=0.30
  K8_m50      K=8, sigma=0.50 stress arm
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from analog_frontier_learning import candidate_bonds,weights_from_rho,gradient_on_candidates,normalized_step,eval_C
from reciprocal_adjoint_probe import retro_source_sequence
from device_error_probe import complex_spectral_maps
from spectral_correlation_compression_probe import port_spectrum_score
from transfer_decomposition_probe import safe_corr

ARMS={
 'K8_ideal':dict(K=8,sigma=0.0),
 'K8_m20':dict(K=8,sigma=.20),
 'K8_m30':dict(K=8,sigma=.30),
 'K16_m30':dict(K=16,sigma=.30),
 'K8_m50':dict(K=8,sigma=.50),
}


def exact_state(m,wh,wv,sT,sD):
    ET=ae.linear_forward(m,wh,wv,sT,store=False);ED=ae.linear_forward(m,wh,wv,sD,store=False);S=ET+ED+1e-30
    aT=2*ED/S**2;aD=-2*ET/S**2
    pT,vT,_=ae.linear_forward(m,wh,wv,sT,store=True);pD,vD,_=ae.linear_forward(m,wh,wv,sD,store=True)
    ghT,gvT=ae.adjoint_grad(m,wh,wv,pT,vT,aT);ghD,gvD=ae.adjoint_grad(m,wh,wv,pD,vD,aD)
    return dict(C=float((ET-ED)/S),gh=ghT+ghD,gv=gvT+gvD,map_corr=1.,map_rel=0.)


def one_order(m,wh,wv,whr,wvr,seq,weight):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    g=weight*np.asarray(p[1:,m.soma[0],m.soma[1]],np.complex128)
    rseq=retro_source_sequence(m,g,reverse=True)
    rp,rv,_=ae.linear_forward(m,whr,wvr,rseq,store=True)
    ZH,ZV=complex_spectral_maps(m,p[:-1],rp[1:])
    eh,ev=ae.adjoint_grad(m,wh,wv,p,v,weight)
    return dict(ZH=ZH,ZV=ZV,exact=(eh,ev),score=port_spectrum_score(seq,rseq))


def physical_state(m,wh,wv,sT,sD,pars,dh,dv):
    whr=np.maximum(1e-12,wh*(1.+dh));wvr=np.maximum(1e-12,wv*(1.+dv))
    ET=ae.linear_forward(m,wh,wv,sT,store=False);ED=ae.linear_forward(m,wh,wv,sD,store=False);S=ET+ED+1e-30
    aT=2*ED/S**2;aD=-2*ET/S**2
    T=one_order(m,wh,wv,whr,wvr,sT,aT);D=one_order(m,wh,wv,whr,wvr,sD,aD)
    ZH=T['ZH']+D['ZH'];ZV=T['ZV']+D['ZV'];eh=T['exact'][0]+D['exact'][0];ev=T['exact'][1]+D['exact'][1]
    score=T['score']+D['score'];order=np.argsort(score)[::-1];kk=np.asarray(order[:int(pars['K'])],int)
    gh=np.real(np.sum(ZH[kk],axis=0));gv=np.real(np.sum(ZV[kk],axis=0))
    ex=np.concatenate([eh.ravel(),ev.ravel()]);ap=np.concatenate([gh.ravel(),gv.ravel()])
    return dict(C=float((ET-ED)/S),gh=gh,gv=gv,map_corr=float(safe_corr(ex,ap)),map_rel=float(np.linalg.norm(ap-ex)/(np.linalg.norm(ex)+1e-30)),bins=[int(x) for x in kk])


def mismatch_pattern(seed,name,shapeh,shapev,sigma):
    if sigma==0:return np.zeros(shapeh),np.zeros(shapev)
    rng=np.random.default_rng(int(seed)*1000003 + sum(map(ord,name))*1009)
    return rng.normal(0.,sigma,size=shapeh),rng.normal(0.,sigma,size=shapev)


def run_arm(m,bwh,bwv,cands,sT,sD,eta,iters,kbath,dk,name,pars=None):
    rho=np.zeros(len(cands),float);traj=[];corr=[];rel=[]
    if pars is not None:dh,dv=mismatch_pattern(m.cfg.seed,name,bwh.shape,bwv.shape,float(pars['sigma']))
    for it in range(int(iters)):
        wh,wv=weights_from_rho(bwh,bwv,cands,rho,kbath,dk)
        st=exact_state(m,wh,wv,sT,sD) if pars is None else physical_state(m,wh,wv,sT,sD,pars,dh,dv)
        if it==0:traj.append(float(st['C']))
        g=gradient_on_candidates(st['gh'],st['gv'],cands,dk);rho=np.clip(rho+normalized_step(g,eta),0.,1.)
        wh2,wv2=weights_from_rho(bwh,bwv,cands,rho,kbath,dk);traj.append(eval_C(m,wh2,wv2,sT,sD));corr.append(st['map_corr']);rel.append(st['map_rel'])
    tr=np.asarray(traj,float)
    return dict(start_C=float(tr[0]),final_C=float(tr[-1]),delta_C=float(tr[-1]-tr[0]),improved=bool(tr[-1]>tr[0]+1e-10),
                monotone_fraction=float(np.mean(np.diff(tr)>=-1e-8)),sum_rho=float(np.sum(rho)),rho=[float(x) for x in rho],trajectory=[float(x) for x in tr],
                mean_map_corr=float(np.mean(corr)),mean_map_relative_l2=float(np.mean(rel)),min_map_corr=float(np.min(corr)),max_map_relative_l2=float(np.max(rel)))


def one(m,lag,steps,max_candidates,eta,iters):
    rng=np.random.default_rng(int(m.cfg.seed)+240240);cands=candidate_bonds(m,max_candidates,rng)
    if len(cands)<2:return None
    bwh,bwv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb)
    sT=ae.source_sequence(m,True,lag,steps);sD=ae.source_sequence(m,False,lag,steps)
    exact=run_arm(m,bwh,bwv,cands,sT,sD,eta,iters,kb,dk,'exact',None);arms={}
    for name,pars in ARMS.items():
        q=run_arm(m,bwh,bwv,cands,sT,sD,eta,iters,kb,dk,name,pars);q['minus_exact_delta']=float(q['delta_C']-exact['delta_C']);arms[name]=q
    return dict(seed=int(m.cfg.seed),n_candidates=len(cands),exact=exact,arms=arms)


def selftest():
    assert ARMS['K8_m20']['sigma']==.2 and ARMS['K16_m30']['K']==16
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=412);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--max-candidates',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=40);ap.add_argument('--out',default='runs/operator_mismatch_learning/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps,a.max_candidates,a.eta,a.iterations)
        if r is None:continue
        rows.append(r);print('seed',seed,'exact',round(r['exact']['delta_C'],4),[(n,round(r['arms'][n]['delta_C'],4),round(r['arms'][n]['mean_map_corr'],4)) for n in ARMS],flush=True)
    summary=dict(bodies=len(rows),exact=dict(mean_delta_C=float(np.mean([r['exact']['delta_C'] for r in rows])),improved=int(np.sum([r['exact']['improved'] for r in rows]))),arms={})
    for n in ARMS:
        q=[r['arms'][n] for r in rows];md=float(np.mean([x['delta_C'] for x in q]));summary['arms'][n]=dict(mean_delta_C=md,improved=int(np.sum([x['improved'] for x in q])),group_gain_ratio=float(md/(summary['exact']['mean_delta_C']+1e-30)),
            mean_minus_exact=float(np.mean([x['minus_exact_delta'] for x in q])),mean_map_corr=float(np.mean([x['mean_map_corr'] for x in q])),mean_map_relative_l2=float(np.mean([x['mean_map_relative_l2'] for x in q])),mean_monotone_fraction=float(np.mean([x['monotone_fraction'] for x in q])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='operator_mismatch_learning_dev_v01',conditions=ARMS,summary=summary,rows=rows),indent=2))
    print('\nOPERATOR MISMATCH LEARNING DEV');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
