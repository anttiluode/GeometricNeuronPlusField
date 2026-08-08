"""Closed-loop compressed-gradient learning under device-inspired measurement errors.

Systematic controller errors are fixed for each absolute DFT bin across the whole
training run.  Local tap/readout noise is redrawn for every gradient measurement.
The port-selected frequency set itself may change as the structure changes.

Arms:
  exact       full discrete adjoint
  K8_ideal    8-bin boundary-selected ideal compressed gradient
  K8_mod      sigma_phi=.025 rad, sigma_amp=.025, tap=.010
  K8_high     sigma_phi=.050 rad, sigma_amp=.050, tap=.020
  K16_high    same high error levels with 16 bins

The amplitude/tap normalizations are explicit toy measurement-error models; the
numerical scales are device-inspired, not a claim to reproduce a specific chip.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from analog_frontier_learning import candidate_bonds,weights_from_rho,gradient_on_candidates,normalized_step,eval_C
from device_error_probe import order_complex
from transfer_decomposition_probe import safe_corr

ARMS={
 'K8_ideal':dict(K=8,phase=0.,amp=0.,tap=0.),
 'K8_mod':dict(K=8,phase=.025,amp=.025,tap=.010),
 'K8_high':dict(K=8,phase=.050,amp=.050,tap=.020),
 'K16_high':dict(K=16,phase=.050,amp=.050,tap=.020),
}


def weights_from_state(m,wh,wv,seqT,seqD):
    ET=ae.linear_forward(m,wh,wv,seqT,store=False);ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30;return ET,ED,2*ED/S**2,-2*ET/S**2


def compressed_gradient(m,wh,wv,seqT,seqD,pars,phase_err,amp_err,rng):
    ET,ED,aT,aD=weights_from_state(m,wh,wv,seqT,seqD)
    T=order_complex(m,wh,wv,seqT,aT);D=order_complex(m,wh,wv,seqD,aD)
    ZH=T['ZH']+D['ZH'];ZV=T['ZV']+D['ZV'];eh=T['exact'][0]+D['exact'][0];ev=T['exact'][1]+D['exact'][1]
    score=T['port_score']+D['port_score'];order=np.argsort(score)[::-1];kk=np.asarray(order[:int(pars['K'])],int)
    zh=ZH[kk].copy();zv=ZV[kk].copy()
    if pars['phase']>0:
        rot=np.exp(1j*phase_err[kk]);zh*=rot[:,None,None];zv*=rot[:,None,None]
    if pars['amp']>0:
        a=amp_err[kk];zh*=a[:,None,None];zv*=a[:,None,None]
    gh=np.real(zh);gv=np.real(zv)
    if pars['tap']>0:
        for j in range(len(kk)):
            sh=float(np.sqrt(np.mean(np.abs(zh[j])**2))+1e-30);sv=float(np.sqrt(np.mean(np.abs(zv[j])**2))+1e-30)
            gh[j]+=rng.normal(0.,float(pars['tap'])*sh,size=gh[j].shape)
            gv[j]+=rng.normal(0.,float(pars['tap'])*sv,size=gv[j].shape)
    ah=np.sum(gh,axis=0);av=np.sum(gv,axis=0)
    ex=np.concatenate([eh.ravel(),ev.ravel()]);ap=np.concatenate([ah.ravel(),av.ravel()])
    return dict(C=float((ET-ED)/(ET+ED+1e-30)),gh=ah,gv=av,map_corr=float(safe_corr(ex,ap)),map_rel=float(np.linalg.norm(ap-ex)/(np.linalg.norm(ex)+1e-30)),bins=[int(x) for x in kk])


def exact_gradient(m,wh,wv,seqT,seqD):
    ET,ED,aT,aD=weights_from_state(m,wh,wv,seqT,seqD)
    pT,vT,_=ae.linear_forward(m,wh,wv,seqT,store=True);pD,vD,_=ae.linear_forward(m,wh,wv,seqD,store=True)
    ghT,gvT=ae.adjoint_grad(m,wh,wv,pT,vT,aT);ghD,gvD=ae.adjoint_grad(m,wh,wv,pD,vD,aD)
    return dict(C=float((ET-ED)/(ET+ED+1e-30)),gh=ghT+ghD,gv=gvT+gvD,map_corr=1.,map_rel=0.)


def run_arm(m,bwh,bwv,cands,seqT,seqD,eta,iters,kbath,dk,name,pars=None):
    rho=np.zeros(len(cands),float);traj=[];corr=[];rels=[]
    Tlen=len(seqT);sysrng=np.random.default_rng(int(m.cfg.seed)*7919 + sum(map(ord,name))*101)
    if pars is None:
        phase_err=np.zeros(Tlen);amp_err=np.ones(Tlen)
    else:
        phase_err=sysrng.normal(0.,float(pars['phase']),size=Tlen) if pars['phase'] else np.zeros(Tlen)
        amp_err=1.+sysrng.normal(0.,float(pars['amp']),size=Tlen) if pars['amp'] else np.ones(Tlen)
    for it in range(int(iters)):
        wh,wv=weights_from_rho(bwh,bwv,cands,rho,kbath,dk)
        if pars is None:st=exact_gradient(m,wh,wv,seqT,seqD)
        else:
            rng=np.random.default_rng(int(m.cfg.seed)*1000003 + sum(map(ord,name))*1009 + it)
            st=compressed_gradient(m,wh,wv,seqT,seqD,pars,phase_err,amp_err,rng)
        if it==0:traj.append(float(st['C']))
        g=gradient_on_candidates(st['gh'],st['gv'],cands,dk);rho=np.clip(rho+normalized_step(g,eta),0.,1.)
        wh2,wv2=weights_from_rho(bwh,bwv,cands,rho,kbath,dk);traj.append(eval_C(m,wh2,wv2,seqT,seqD))
        corr.append(st['map_corr']);rels.append(st['map_rel'])
    tr=np.asarray(traj,float)
    return dict(start_C=float(tr[0]),final_C=float(tr[-1]),delta_C=float(tr[-1]-tr[0]),improved=bool(tr[-1]>tr[0]+1e-10),
                monotone_fraction=float(np.mean(np.diff(tr)>=-1e-8)),sum_rho=float(np.sum(rho)),rho=[float(x) for x in rho],
                mean_map_corr=float(np.mean(corr)),mean_map_relative_l2=float(np.mean(rels)),trajectory=[float(x) for x in tr])


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
    assert ARMS['K8_mod']['phase']==.025 and ARMS['K16_high']['K']==16
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=412);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--max-candidates',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=40);ap.add_argument('--out',default='runs/device_error_learning/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps,a.max_candidates,a.eta,a.iterations)
        if r is None:continue
        rows.append(r);print('seed',seed,'exact',round(r['exact']['delta_C'],4),[(n,round(r['arms'][n]['delta_C'],4),round(r['arms'][n]['mean_map_corr'],4)) for n in ARMS],flush=True)
    summary=dict(bodies=len(rows),arms={},exact=dict(mean_delta_C=float(np.mean([r['exact']['delta_C'] for r in rows])),improved=int(np.sum([r['exact']['improved'] for r in rows]))))
    for n in ARMS:
        q=[r['arms'][n] for r in rows];summary['arms'][n]=dict(mean_delta_C=float(np.mean([x['delta_C'] for x in q])),improved=int(np.sum([x['improved'] for x in q])),
            group_gain_ratio=float(np.mean([x['delta_C'] for x in q])/(summary['exact']['mean_delta_C']+1e-30)),mean_minus_exact=float(np.mean([x['minus_exact_delta'] for x in q])),
            mean_map_corr=float(np.mean([x['mean_map_corr'] for x in q])),mean_map_relative_l2=float(np.mean([x['mean_map_relative_l2'] for x in q])),mean_monotone_fraction=float(np.mean([x['monotone_fraction'] for x in q])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='device_error_learning_dev_v01',conditions=ARMS,summary=summary,rows=rows),indent=2))
    print('\nDEVICE ERROR LEARNING DEV');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
