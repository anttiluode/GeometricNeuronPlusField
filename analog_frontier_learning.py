"""Joint analog frontier learning with relinearized soma adjoints.

See ANALOG_FRONTIER_LEARNING_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # patches ae.weighted_lap to the exact mature-boundary operator
from structural_interference_probe import event_candidates, n4


def contrast_adjoint_weights(m, wh, wv, seqT, seqD):
    pT,vT,ET=ae.linear_forward(m,wh,wv,seqT,store=True)
    pD,vD,ED=ae.linear_forward(m,wh,wv,seqD,store=True)
    S=ET+ED+1e-30
    aT=2.0*ED/(S*S)
    aD=-2.0*ET/(S*S)
    ghT,gvT=ae.adjoint_grad(m,wh,wv,pT,vT,aT)
    ghD,gvD=ae.adjoint_grad(m,wh,wv,pD,vD,aD)
    return dict(C=ae.contrast(ET,ED),ET=ET,ED=ED,gh=ghT+ghD,gv=gvT+gvD)


def candidate_bonds(m, max_candidates, rng):
    adds,_,_,_=event_candidates(m,int(max_candidates),rng)
    body=m.body.astype(bool)
    out=[]
    for p in adds:
        qs=[q for q in n4(*p,body.shape) if body[q]]
        if len(qs)==1:
            out.append((tuple(map(int,p)),tuple(map(int,qs[0]))))
    return out


def weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk):
    wh=base_wh.copy();wv=base_wv.copy()
    for (p,q),r in zip(cands,rho):
        ae.set_edge(wh,wv,p,q,float(kbath+float(r)*dk))
    return wh,wv


def gradient_on_candidates(gh,gv,cands,dk):
    return np.asarray([dk*ae.edge_lookup(gh,gv,p,q) for p,q in cands],float)


def normalized_step(g,eta):
    g=np.asarray(g,float)
    mx=float(np.max(np.abs(g))) if len(g) else 0.0
    if not np.isfinite(mx) or mx<=1e-30:
        return np.zeros_like(g)
    return float(eta)*g/mx


def eval_C(m,wh,wv,seqT,seqD):
    C,_,_=ae.linear_contrast(m,wh,wv,seqT,seqD)
    return float(C)


def run_arm_relinearized(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kbath,dk):
    rho=np.zeros(len(cands),float)
    traj=[]; grad_hist=[]
    initial_g=None
    initial_positive_later_negative=np.zeros(len(cands),bool)
    for it in range(int(iters)):
        wh,wv=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        st=contrast_adjoint_weights(m,wh,wv,seqT,seqD)
        if it==0: traj.append(float(st['C']))
        g=gradient_on_candidates(st['gh'],st['gv'],cands,dk)
        if initial_g is None: initial_g=g.copy()
        initial_positive_later_negative |= ((initial_g>0) & (g<0))
        grad_hist.append([float(x) for x in g])
        rho=np.clip(rho+normalized_step(g,eta),0.0,1.0)
        wh2,wv2=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        traj.append(eval_C(m,wh2,wv2,seqT,seqD))
    return dict(rho=rho,traj=np.asarray(traj,float),initial_g=np.asarray(initial_g,float),
                grad_hist=grad_hist,reversed_candidates=initial_positive_later_negative)


def run_arm_frozen(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kbath,dk,initial_g):
    rho=np.zeros(len(cands),float)
    C0=eval_C(m,base_wh,base_wv,seqT,seqD)
    traj=[C0]
    step=normalized_step(initial_g,eta)
    for _ in range(int(iters)):
        rho=np.clip(rho+step,0.0,1.0)
        wh,wv=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        traj.append(eval_C(m,wh,wv,seqT,seqD))
    return dict(rho=rho,traj=np.asarray(traj,float))


def run_arm_shuffled(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kbath,dk,rng):
    rho=np.zeros(len(cands),float)
    traj=[]
    for it in range(int(iters)):
        wh,wv=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        st=contrast_adjoint_weights(m,wh,wv,seqT,seqD)
        if it==0: traj.append(float(st['C']))
        g=gradient_on_candidates(st['gh'],st['gv'],cands,dk)
        gp=g[rng.permutation(len(g))] if len(g) else g
        rho=np.clip(rho+normalized_step(gp,eta),0.0,1.0)
        wh2,wv2=weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
        traj.append(eval_C(m,wh2,wv2,seqT,seqD))
    return dict(rho=rho,traj=np.asarray(traj,float))


def arm_summary(arm):
    rho=np.asarray(arm['rho'],float);tr=np.asarray(arm['traj'],float)
    return dict(
        start_C=float(tr[0]),final_C=float(tr[-1]),delta_C=float(tr[-1]-tr[0]),
        monotone_fraction=float(np.mean(np.diff(tr)>=-1e-8)) if len(tr)>1 else float('nan'),
        sum_rho=float(np.sum(rho)),
        active_fraction=float(np.mean(rho>1e-4)) if len(rho) else float('nan'),
        strongly_matured_fraction=float(np.mean(rho>.25)) if len(rho) else float('nan'),
        rho=[float(x) for x in rho],
        trajectory=[float(x) for x in tr],
    )


def body_probe(m,lag,steps,max_candidates,eta,iters):
    rng=np.random.default_rng(int(m.cfg.seed)+240240)
    cands=candidate_bonds(m,max_candidates,rng)
    if len(cands)<2:
        return None
    base_wh,base_wv=ae.bond_weights(m,m.body)
    kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)

    rel=run_arm_relinearized(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kb,dk)
    fro=run_arm_frozen(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kb,dk,rel['initial_g'])
    shu=run_arm_shuffled(m,base_wh,base_wv,cands,seqT,seqD,eta,iters,kb,dk,np.random.default_rng(int(m.cfg.seed)+777777))

    rels=arm_summary(rel);fros=arm_summary(fro);shus=arm_summary(shu)
    initcorr=float(np.corrcoef(rel['initial_g'],rel['rho'])[0,1]) if len(cands)>1 and np.std(rel['initial_g'])>1e-15 and np.std(rel['rho'])>1e-15 else float('nan')
    return dict(
        seed=int(m.cfg.seed),cells=int(m.body.sum()),n_candidates=len(cands),
        candidates=[[list(p),list(q)] for p,q in cands],
        initial_gradient=[float(x) for x in rel['initial_g']],
        initial_gradient_final_rho_corr=initcorr,
        initially_positive_later_negative_n=int(np.sum(rel['reversed_candidates'])),
        initially_positive_n=int(np.sum(rel['initial_g']>0)),
        relinearized=rels,frozen=fros,shuffled=shus,
        rel_minus_frozen=float(rels['delta_C']-fros['delta_C']),
        rel_minus_shuffled=float(rels['delta_C']-shus['delta_C']),
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=240)
    ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-candidates',type=int,default=8)
    ap.add_argument('--eta',type=float,default=.01)
    ap.add_argument('--iterations',type=int,default=40)
    ap.add_argument('--out',default='runs/analog_frontier_learning/analog_frontier_learning.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    g=np.array([2.,-1.,0.])
    s=normalized_step(g,.01)
    assert abs(s[0]-.01)<1e-12 and abs(s[1]+.005)<1e-12
    print('selftest ok')


def main():
    a=parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor

    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,a.lag,a.steps,a.max_candidates,a.eta,a.iterations)
        if r is None:continue
        rows.append(r)
        print(f"seed {seed}: n={r['n_candidates']} rel={r['relinearized']['delta_C']:+.4f} "
              f"frozen={r['frozen']['delta_C']:+.4f} shuffled={r['shuffled']['delta_C']:+.4f} "
              f"mono={r['relinearized']['monotone_fraction']:.2f}",flush=True)

    if not rows:raise SystemExit('No valid bodies')
    rel=np.asarray([r['relinearized']['delta_C'] for r in rows],float)
    df=np.asarray([r['rel_minus_frozen'] for r in rows],float)
    ds=np.asarray([r['rel_minus_shuffled'] for r in rows],float)
    mono=np.asarray([r['relinearized']['monotone_fraction'] for r in rows],float)
    summary=dict(
        bodies=len(rows),eta=a.eta,iterations=a.iterations,max_candidates=a.max_candidates,
        mean_delta_relinearized=float(np.mean(rel)),
        median_delta_relinearized=float(np.median(rel)),
        improved_relinearized_bodies=int(np.sum(rel>0)),
        mean_delta_frozen=float(np.mean([r['frozen']['delta_C'] for r in rows])),
        mean_delta_shuffled=float(np.mean([r['shuffled']['delta_C'] for r in rows])),
        mean_rel_minus_frozen=float(np.mean(df)),
        rel_beats_frozen_bodies=int(np.sum(df>0)),
        mean_rel_minus_shuffled=float(np.mean(ds)),
        rel_beats_shuffled_bodies=int(np.sum(ds>0)),
        median_relinearized_monotone_fraction=float(np.median(mono)),
        mean_relinearized_sum_rho=float(np.mean([r['relinearized']['sum_rho'] for r in rows])),
        mean_frozen_sum_rho=float(np.mean([r['frozen']['sum_rho'] for r in rows])),
        mean_shuffled_sum_rho=float(np.mean([r['shuffled']['sum_rho'] for r in rows])),
        mean_initial_gradient_final_rho_corr=float(np.nanmean([r['initial_gradient_final_rho_corr'] for r in rows])),
        total_initially_positive_later_negative=int(np.sum([r['initially_positive_later_negative_n'] for r in rows])),
        total_initially_positive=int(np.sum([r['initially_positive_n'] for r in rows])),
    )
    summary['D1_pass']=bool(summary['mean_delta_relinearized']>.005 and summary['improved_relinearized_bodies']>=10)
    summary['D2_pass']=bool(summary['mean_rel_minus_frozen']>.003 and summary['rel_beats_frozen_bodies']>=9)
    summary['D3_pass']=bool(summary['mean_rel_minus_shuffled']>.003 and summary['rel_beats_shuffled_bodies']>=9)
    summary['D4_pass']=bool(summary['median_relinearized_monotone_fraction']>.80)

    payload=dict(experiment='analog_frontier_learning_discovery_v01',prereg='ANALOG_FRONTIER_LEARNING_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,
                 max_candidates=a.max_candidates,eta=a.eta,iterations=a.iterations,
                 summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nANALOG FRONTIER LEARNING RECEIPT')
    for k,v in summary.items():print(f' {k}: {v}')

if __name__=='__main__':main()
