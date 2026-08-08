"""Closed-loop analog frontier learning with spectrally compressed reciprocal gradients.

Development probe following SPECTRAL_CORRELATION_COMPRESSION_V01.md.

All arms use the same fixed candidate frontier bonds, eta, iteration count and
max-normalized local update rule.  They differ only in gradient readout:

  exact   full discrete adjoint
  K4      boundary-selected 4-bin spectral correlation
  K8      boundary-selected 8-bin spectral correlation
  K16     boundary-selected 16-bin spectral correlation

The boundary frequency ranking is recomputed from the current arm's external source
and soma-return spectra at every relinearization step.  It never inspects internal
per-bond gradient contributions.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from analog_frontier_learning import (candidate_bonds,weights_from_rho,gradient_on_candidates,
                                      normalized_step,eval_C)
from spectral_correlation_compression_probe import order_data
from transfer_decomposition_probe import safe_corr

KS=(4,8,16)


def gradient_state(m,wh,wv,seqT,seqD,K=None):
    # Energies determine the derivative weights for the current physical state.
    ET=ae.linear_forward(m,wh,wv,seqT,store=False)
    ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30; aT=2.0*ED/(S*S); aD=-2.0*ET/(S*S)
    T=order_data(m,wh,wv,seqT,aT); D=order_data(m,wh,wv,seqD,aD)
    eh=T['exact'][0]+D['exact'][0]; ev=T['exact'][1]+D['exact'][1]
    if K is None:
        return dict(C=float((ET-ED)/S),gh=eh,gv=ev,map_corr=1.0,map_rel=0.0,bins=None)
    GH=T['GH']+D['GH']; GV=T['GV']+D['GV']
    score=T['port_score']+D['port_score']; order=np.argsort(score)[::-1]
    kk=np.asarray(order[:min(int(K),len(order))],int)
    gh=np.sum(GH[kk],axis=0); gv=np.sum(GV[kk],axis=0)
    ex=np.concatenate([eh.ravel(),ev.ravel()]); ap=np.concatenate([gh.ravel(),gv.ravel()])
    rel=float(np.linalg.norm(ap-ex)/(np.linalg.norm(ex)+1e-30))
    return dict(C=float((ET-ED)/S),gh=gh,gv=gv,map_corr=float(safe_corr(ex,ap)),map_rel=rel,bins=[int(x) for x in kk])


def run_arm(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kbath,dk,K=None):
    rho=np.zeros(len(cands),float); traj=[]; mapcorr=[]; maprel=[]; binhist=[]
    for it in range(int(iters)):
        wh,wv=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        st=gradient_state(m,wh,wv,seqT,seqD,K)
        if it==0: traj.append(float(st['C']))
        g=gradient_on_candidates(st['gh'],st['gv'],cands,dk)
        rho=np.clip(rho+normalized_step(g,eta),0.,1.)
        wh2,wv2=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        traj.append(eval_C(m,wh2,wv2,seqT,seqD))
        mapcorr.append(float(st['map_corr'])); maprel.append(float(st['map_rel'])); binhist.append(st['bins'])
    tr=np.asarray(traj,float)
    return dict(start_C=float(tr[0]),final_C=float(tr[-1]),delta_C=float(tr[-1]-tr[0]),
                monotone_fraction=float(np.mean(np.diff(tr)>=-1e-8)),sum_rho=float(np.sum(rho)),
                rho=[float(x) for x in rho],trajectory=[float(x) for x in tr],
                mean_map_corr=float(np.mean(mapcorr)),mean_map_relative_l2=float(np.mean(maprel)),
                min_map_corr=float(np.min(mapcorr)),max_map_relative_l2=float(np.max(maprel)),bin_history=binhist)


def one(m,lag,steps,max_candidates,eta,iters):
    rng=np.random.default_rng(int(m.cfg.seed)+240240)
    cands=candidate_bonds(m,max_candidates,rng)
    if len(cands)<2:return None
    bwh,bwv=ae.bond_weights(m,m.body);kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-kb)
    seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    exact=run_arm(m,bwh,bwv,cands,seqT,seqD,eta,iters,kb,dk,None)
    arms={str(K):run_arm(m,bwh,bwv,cands,seqT,seqD,eta,iters,kb,dk,K) for K in KS}
    for K in KS:
        q=arms[str(K)];q['minus_exact_delta']=float(q['delta_C']-exact['delta_C'])
        q['gain_fraction_of_exact']=float(q['delta_C']/(exact['delta_C']+1e-30))
    return dict(seed=int(m.cfg.seed),cells=int(m.body.sum()),n_candidates=len(cands),
                candidates=[[list(p),list(q)] for p,q in cands],exact=exact,compressed=arms)


def selftest():
    g=np.array([2.,-1.,0.]);s=normalized_step(g,.01);assert abs(s[0]-.01)<1e-12
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--max-candidates',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=40);ap.add_argument('--out',default='runs/spectral_gradient_learning/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps,a.max_candidates,a.eta,a.iterations)
        if r is None:continue
        rows.append(r);print('seed',seed,'exact',round(r['exact']['delta_C'],4),[(K,round(r['compressed'][str(K)]['delta_C'],4),round(r['compressed'][str(K)]['mean_map_corr'],4)) for K in KS],flush=True)
    if not rows:raise SystemExit('No valid bodies')
    summary=dict(bodies=len(rows),eta=a.eta,iterations=a.iterations,max_candidates=a.max_candidates,
                 exact=dict(mean_delta_C=float(np.mean([r['exact']['delta_C'] for r in rows])),improved=int(np.sum([r['exact']['delta_C']>0 for r in rows])),mean_sum_rho=float(np.mean([r['exact']['sum_rho'] for r in rows]))),compressed={})
    for K in KS:
        q=[r['compressed'][str(K)] for r in rows]
        summary['compressed'][str(K)]=dict(mean_delta_C=float(np.mean([x['delta_C'] for x in q])),improved=int(np.sum([x['delta_C']>0 for x in q])),
            mean_gain_fraction_of_exact=float(np.mean([x['gain_fraction_of_exact'] for x in q])),median_gain_fraction_of_exact=float(np.median([x['gain_fraction_of_exact'] for x in q])),
            mean_minus_exact_delta=float(np.mean([x['minus_exact_delta'] for x in q])),beats_exact=int(np.sum([x['minus_exact_delta']>0 for x in q])),
            mean_map_corr=float(np.mean([x['mean_map_corr'] for x in q])),mean_map_relative_l2=float(np.mean([x['mean_map_relative_l2'] for x in q])),
            mean_monotone_fraction=float(np.mean([x['monotone_fraction'] for x in q])),mean_sum_rho=float(np.mean([x['sum_rho'] for x in q])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='spectral_gradient_learning_dev_v01',summary=summary,rows=rows),indent=2))
    print('\nSPECTRAL GRADIENT LEARNING DEV');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
