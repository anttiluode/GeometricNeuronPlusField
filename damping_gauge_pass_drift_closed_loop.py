"""Closed-loop conformal-echo learning with independent +/- pass drift.

Development only on reused bodies.

The static differential-pass probe found approximately:
    0.2% differential coupling drift -> mean gradient corr ~.94
    0.5% differential coupling drift -> mean gradient corr ~.80
under the rest of the combined hardware-error model.

This script asks whether those imperfect two-pass gradients still train the
actual physical device in loop.

Arms:
    physical_exact      digital exact gradient of actual lossy device
    echo_drift_002      two-pass echo with 0.2% differential +/- drift
    echo_drift_005      two-pass echo with 0.5% differential +/- drift
    shuffled_002        same 0.2%-drift echo credit, permuted among candidates

All echo arms also include:
    mean residual loss .005
    loss CV .20
    loss calibration +5%
    mirror alpha .95
    common reverse drift .02
    candidate gradient readout noise .05
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import analog_frontier_learning as afl
import adjoint_eligibility_probe as ae
import damping_gauge_combined_closed_loop as cc
import damping_gauge_loss_heterogeneity_probe as lh
import damping_gauge_pass_drift_probe as pd
from transfer_decomposition_probe import safe_corr


def noisy_candidate(gh,gv,cands,dk,noise,rng):
    g=afl.gradient_on_candidates(gh,gv,cands,dk)
    if noise>0 and len(g):
        g=g+rng.standard_normal(g.shape)*(float(noise)*float(np.sqrt(np.mean(g*g))))
    return g


def echo_two_pass_gradient(m,wh,wv,st,epsmap,eps_hat,a,diff,rng):
    gt=pd.two_pass_order(m,wh,wv,st['xT'],st['uT'],st['qT'],epsmap,eps_hat,a.mirror_alpha,a.common_drift,diff,rng)
    gd=pd.two_pass_order(m,wh,wv,st['xD'],st['uD'],st['qD'],epsmap,eps_hat,a.mirror_alpha,a.common_drift,diff,rng)
    return gt[0]+gd[0],gt[1]+gd[1]


def run_arm(m,base_wh,base_wv,cands,seqT,seqD,epsmap,a,mode,diff,seedoff):
    kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    eps_hat=float(a.mean_eps)*(1+float(a.loss_cal_error))
    rho=np.zeros(len(cands),float);hist=[]
    rng=np.random.default_rng(int(m.cfg.seed)*1_000_003+seedoff)
    for it in range(a.iterations):
        wh,wv=afl.weights_from_rho(base_wh,base_wv,cands,rho,kb,dk)
        st=cc.physical_state(m,wh,wv,seqT,seqD,epsmap)
        ghx,gvx=cc.exact_physical_gradient(m,wh,wv,st,epsmap)
        gx=afl.gradient_on_candidates(ghx,gvx,cands,dk)
        if mode=='physical_exact':
            g=gx.copy();corr=float('nan')
        else:
            gh,gv=echo_two_pass_gradient(m,wh,wv,st,epsmap,eps_hat,a,diff,rng)
            g=noisy_candidate(gh,gv,cands,dk,a.readout_noise,rng)
            corr=float(safe_corr(gx,g)) if len(g)>1 else float('nan')
            if mode=='shuffle': g=g[rng.permutation(len(g))]
        hist.append(dict(iteration=it,physical_C=float(st['C']),echo_corr=corr))
        rho=np.clip(rho+afl.normalized_step(g,a.eta),0,1)
    wh,wv=afl.weights_from_rho(base_wh,base_wv,cands,rho,kb,dk)
    fin=cc.physical_state(m,wh,wv,seqT,seqD,epsmap)
    corrs=[x['echo_corr'] for x in hist if np.isfinite(x['echo_corr'])]
    return dict(start_C=hist[0]['physical_C'],final_C=float(fin['C']),delta_C=float(fin['C']-hist[0]['physical_C']),
                mean_corr=float(np.mean(corrs)) if corrs else float('nan'),min_corr=float(np.min(corrs)) if corrs else float('nan'),rho=rho.tolist(),history=hist)


def body(m,a):
    rng=np.random.default_rng(int(m.cfg.seed)+240240);cands=afl.candidate_bonds(m,a.max_candidates,rng)
    if len(cands)<2:return None
    bwh,bwv=ae.bond_weights(m,m.body);seqT=ae.source_sequence(m,True,a.lag,a.steps);seqD=ae.source_sequence(m,False,a.lag,a.steps)
    epsmap=lh.make_epsmap(m,a.mean_eps,a.loss_cv,np.random.default_rng(int(m.cfg.seed)+8_818_181))
    arms={
      'physical_exact':run_arm(m,bwh,bwv,cands,seqT,seqD,epsmap,a,'physical_exact',0,100),
      'echo_002':run_arm(m,bwh,bwv,cands,seqT,seqD,epsmap,a,'echo',.002,200),
      'echo_005':run_arm(m,bwh,bwv,cands,seqT,seqD,epsmap,a,'echo',.005,300),
      'shuffle_002':run_arm(m,bwh,bwv,cands,seqT,seqD,epsmap,a,'shuffle',.002,400),
    }
    return dict(seed=int(m.cfg.seed),n_candidates=len(cands),arms=arms)


def summarize(rows):
    out={'bodies':len(rows),'arms':{}}
    for name in ('physical_exact','echo_002','echo_005','shuffle_002'):
        d=np.asarray([r['arms'][name]['delta_C'] for r in rows]);c=[r['arms'][name]['mean_corr'] for r in rows if np.isfinite(r['arms'][name]['mean_corr'])]
        out['arms'][name]=dict(mean_delta=float(np.mean(d)),median_delta=float(np.median(d)),positive=int(np.sum(d>0)),mean_corr=float(np.mean(c)) if c else float('nan'))
    ex=np.mean([r['arms']['physical_exact']['delta_C'] for r in rows])
    for name in ('echo_002','echo_005'):
        out[name+'_retained_fraction']=float(np.mean([r['arms'][name]['delta_C'] for r in rows])/(ex+1e-30))
    out['echo002_minus_shuffle']=dict(mean=float(np.mean([r['arms']['echo_002']['delta_C']-r['arms']['shuffle_002']['delta_C'] for r in rows])),positive=int(sum(r['arms']['echo_002']['delta_C']>r['arms']['shuffle_002']['delta_C'] for r in rows)))
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=472);ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--max-candidates',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=16)
    ap.add_argument('--mean-eps',type=float,default=.005);ap.add_argument('--loss-cv',type=float,default=.20);ap.add_argument('--loss-cal-error',type=float,default=.05);ap.add_argument('--mirror-alpha',type=float,default=.95);ap.add_argument('--common-drift',type=float,default=.02);ap.add_argument('--readout-noise',type=float,default=.05)
    ap.add_argument('--out',default='runs/damping_gauge_pass_drift_closed_loop/dev_472_477.json');a=ap.parse_args()
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));bb=m.bootstrap()
        if not bb.get('ok'):continue
        m.mature=True;r=body(m,a)
        if r is None:continue
        rows.append(r);A=r['arms'];print(f"seed {seed}: exact={A['physical_exact']['delta_C']:+.4f} d002={A['echo_002']['delta_C']:+.4f} d005={A['echo_005']['delta_C']:+.4f} shuffle={A['shuffle_002']['delta_C']:+.4f}",flush=True)
    out=dict(config=vars(a),rows=rows,summary=summarize(rows));q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2));print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':main()
