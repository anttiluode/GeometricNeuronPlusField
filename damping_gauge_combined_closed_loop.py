"""Closed-loop graded-bond learning with a combined hardware-error echo gradient.

Development-only on reused bodies.  This is the first test in the damping-gauge
branch where static gradient-map correlation is not the primary criterion.

We compare four learning arms on the same candidate frontier bonds:

    nominal_exact
        existing exact adjoint of the intended damped task

    physical_exact
        digital reference gradient of the *actual compiled physical device*
        including fixed spatial residual-loss disorder

    echo_combined
        memory-free physical-echo surrogate with all of the following present
        simultaneously:
            mean residual loss eps=.005
            cellwise loss CV=.20
            +5% scalar loss-calibration error
            5% terminal time-mirror error
            2% reverse-operator coupling drift
            5% local gradient-readout noise

    shuffled_echo
        same echo gradient values but randomly permuted among candidate bonds
        before every update (norm/histogram matched credit-placement control)

The physical objective is the contrast evaluated from the actual compiled
forward device.  The nominal objective is also recorded to ask whether training
the imperfect physical implementation transfers back to the intended task.

The goal is deliberately practical:

    does the corrupted constant-history echo gradient actually train the body
    for several steps, or does its attractive static correlation fail in loop?

No held-out seeds are used here and no novelty claim is attached.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path
import numpy as np

# Import analog_frontier_learning first: it imports adjoint_dose_probe, which
# patches ae.weighted_lap to the exact mature-boundary operator used by the
# graded frontier learner.
import analog_frontier_learning as afl
import adjoint_eligibility_probe as ae
import damping_gauge_reversal_probe as dg
import damping_gauge_residual_loss_v02 as dl
import damping_gauge_loss_heterogeneity_probe as lh
from transfer_decomposition_probe import safe_corr


def flat(h,v): return np.concatenate([np.ravel(h),np.ravel(v)])


def cosine(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))


def physical_energy(m,x,r_desired):
    soma=tuple(map(int,m.soma))
    T=len(x)-1
    wt=(float(r_desired)**(2*np.arange(1,T+1,dtype=float)))
    return float(np.sum(wt*np.abs(x[1:,soma[0],soma[1]])**2))


def perturb_weights(wh,wv,sigma,rng):
    if sigma<=0:return wh.copy(),wv.copy()
    ph=np.exp(float(sigma)*rng.standard_normal(wh.shape)-.5*float(sigma)**2)
    pv=np.exp(float(sigma)*rng.standard_normal(wv.shape)-.5*float(sigma)**2)
    return wh*ph,wv*pv


def damaged_compensated_echo(m,wh,wv,x,u,epsmap,eps_hat,mirror_alpha):
    """Same-device conformal echo with scalar calibration + terminal error."""
    ah=1.0-float(eps_hat); back=1.0-epsmap; T=len(x)-1
    y0=x[T].copy()
    ideal_y1=ah*x[T-1]
    # alpha=1 gives exact calibrated conformal terminal state; alpha=0 leaves
    # the second reverse state equal to terminal position (no momentum flip).
    y1=y0-float(mirror_alpha)*(y0-ideal_y1)
    out=[y0,y1]
    for j in range(T-1):
        src=(ah**(j+1))*u[T-1-j]
        nxt=lh.map_apply(m,wh,wv,out[j+1],epsmap)-back*out[j]+src
        out.append(nxt)
    return np.asarray(out)


def weighted_interference_gradient(m,y,b,eps_hat):
    *_,beta,_=dg.params(m); ah=1.0-float(eps_hat); T=len(y)-1
    yp=y[1:]+b[1:]; ym=y[1:]-b[1:]
    ph,pv=dg.edge_diffs(yp); mh,mv=dg.edge_diffs(ym)
    wt=(ah**(-np.arange(1,T+1,dtype=float))).reshape((T,1,1))
    ch=.25*np.sum(wt*(np.abs(ph)**2-np.abs(mh)**2),axis=0)
    cv=.25*np.sum(wt*(np.abs(pv)**2-np.abs(mv)**2),axis=0)
    return (-2.0*beta*ch).real,(-2.0*beta*cv).real


def compiled_sources(m,wh,wv,seq):
    # u depends only on intended uniform damping/source envelope, but using the
    # current operator here keeps the compiler path explicit.
    z,u=dg.gauge_forward(m,wh,wv,seq)
    *_,r,_,_=dg.params(m)
    return u,r


def physical_state(m,wh,wv,seqT,seqD,epsmap):
    uT,r=compiled_sources(m,wh,wv,seqT)
    uD,_=compiled_sources(m,wh,wv,seqD)
    xT=lh.forward_map(m,wh,wv,uT,epsmap)
    xD=lh.forward_map(m,wh,wv,uD,epsmap)
    ET=physical_energy(m,xT,r); ED=physical_energy(m,xD,r)
    S=ET+ED+1e-30
    C=float((ET-ED)/S)
    cT=2.0*ED/(S*S); cD=-2.0*ET/(S*S)
    qT=lh.objective_sources(m,xT,cT,r)
    qD=lh.objective_sources(m,xD,cD,r)
    return dict(C=C,ET=ET,ED=ED,r=r,uT=uT,uD=uD,xT=xT,xD=xD,qT=qT,qD=qD,cT=cT,cD=cD)


def exact_physical_gradient(m,wh,wv,st,epsmap):
    gT=lh.exact_gradient_map(m,wh,wv,st['xT'],st['qT'],epsmap)
    gD=lh.exact_gradient_map(m,wh,wv,st['xD'],st['qD'],epsmap)
    return gT[0]+gD[0],gT[1]+gD[1]


def echo_gradient(m,wh,wv,st,epsmap,eps_hat,mirror_alpha,operator_drift,rng):
    # One quasi-static reverse-operator mismatch is shared by the +/- phase
    # states.  Pass-to-pass differential drift is a separate future wall.
    rwh,rwv=perturb_weights(wh,wv,operator_drift,rng)
    out=[]
    for x,u,q in ((st['xT'],st['uT'],st['qT']),(st['xD'],st['uD'],st['qD'])):
        y=damaged_compensated_echo(m,rwh,rwv,x,u,epsmap,eps_hat,mirror_alpha)
        b=lh.physical_adjoint_map(m,rwh,rwv,q,epsmap)
        out.append(weighted_interference_gradient(m,y,b,eps_hat))
    return out[0][0]+out[1][0],out[0][1]+out[1][1]


def noisy_candidate_gradient(gh,gv,cands,dk,noise,rng):
    g=afl.gradient_on_candidates(gh,gv,cands,dk)
    if float(noise)>0 and len(g):
        rms=float(np.sqrt(np.mean(g*g)))
        g=g+rng.standard_normal(g.shape)*(float(noise)*rms)
    return g


def eval_nominal(m,wh,wv,seqT,seqD):
    return afl.eval_C(m,wh,wv,seqT,seqD)


def make_loss_map(m,mean_eps,cv,rng):
    return lh.make_epsmap(m,float(mean_eps),float(cv),rng)


def arm_step(m,base_wh,base_wv,cands,rho,seqT,seqD,epsmap,eta,kbath,dk,
             mode,eps_hat,mirror_alpha,operator_drift,readout_noise,rng):
    wh,wv=afl.weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
    pst=physical_state(m,wh,wv,seqT,seqD,epsmap)
    nominal=eval_nominal(m,wh,wv,seqT,seqD)

    gref_h,gref_v=exact_physical_gradient(m,wh,wv,pst,epsmap)
    gref=afl.gradient_on_candidates(gref_h,gref_v,cands,dk)

    if mode=='physical_exact':
        g=gref.copy(); gecho=None
    elif mode in ('echo_combined','shuffled_echo'):
        gh,gv=echo_gradient(m,wh,wv,pst,epsmap,eps_hat,mirror_alpha,operator_drift,rng)
        gecho=noisy_candidate_gradient(gh,gv,cands,dk,readout_noise,rng)
        g=gecho.copy()
        if mode=='shuffled_echo' and len(g):
            g=g[rng.permutation(len(g))]
    elif mode=='nominal_exact':
        st=afl.contrast_adjoint_weights(m,wh,wv,seqT,seqD)
        g=afl.gradient_on_candidates(st['gh'],st['gv'],cands,dk)
        gecho=None
    else: raise ValueError(mode)

    corr=float(safe_corr(gref,gecho)) if gecho is not None and len(g)>1 else float('nan')
    cos=cosine(gref,gecho) if gecho is not None else float('nan')
    rho2=np.clip(rho+afl.normalized_step(g,eta),0.0,1.0)
    return rho2,dict(physical_C=float(pst['C']),nominal_C=float(nominal),
                     exact_physical_grad_norm=float(np.linalg.norm(gref)),
                     used_grad_norm=float(np.linalg.norm(g)),
                     echo_vs_exact_corr=corr,echo_vs_exact_cosine=cos)


def final_eval(m,base_wh,base_wv,cands,rho,seqT,seqD,epsmap,kbath,dk):
    wh,wv=afl.weights_from_rho(base_wh,base_wv,cands,rho,kbath,dk)
    pst=physical_state(m,wh,wv,seqT,seqD,epsmap)
    return dict(physical_C=float(pst['C']),nominal_C=float(eval_nominal(m,wh,wv,seqT,seqD)),rho=[float(x) for x in rho])


def run_arm(m,base_wh,base_wv,cands,seqT,seqD,epsmap,eta,iters,kbath,dk,mode,
            eps_hat,mirror_alpha,operator_drift,readout_noise,seed_offset):
    rho=np.zeros(len(cands),float); hist=[]
    rng=np.random.default_rng(int(m.cfg.seed)*1_000_003+int(seed_offset))
    for it in range(int(iters)):
        rho,rec=arm_step(m,base_wh,base_wv,cands,rho,seqT,seqD,epsmap,eta,kbath,dk,
                         mode,eps_hat,mirror_alpha,operator_drift,readout_noise,rng)
        rec['iteration']=it;hist.append(rec)
    fin=final_eval(m,base_wh,base_wv,cands,rho,seqT,seqD,epsmap,kbath,dk)
    start=hist[0]
    corr=[x['echo_vs_exact_corr'] for x in hist if np.isfinite(x['echo_vs_exact_corr'])]
    cos=[x['echo_vs_exact_cosine'] for x in hist if np.isfinite(x['echo_vs_exact_cosine'])]
    return dict(mode=mode,start_physical_C=float(start['physical_C']),final_physical_C=float(fin['physical_C']),
                delta_physical_C=float(fin['physical_C']-start['physical_C']),
                start_nominal_C=float(start['nominal_C']),final_nominal_C=float(fin['nominal_C']),
                delta_nominal_C=float(fin['nominal_C']-start['nominal_C']),
                mean_echo_corr=float(np.mean(corr)) if corr else float('nan'),
                min_echo_corr=float(np.min(corr)) if corr else float('nan'),
                mean_echo_cosine=float(np.mean(cos)) if cos else float('nan'),
                rho=fin['rho'],history=hist)


def body_probe(m,a):
    rng=np.random.default_rng(int(m.cfg.seed)+240240)
    cands=afl.candidate_bonds(m,a.max_candidates,rng)
    if len(cands)<2:return None
    base_wh,base_wv=ae.bond_weights(m,m.body)
    kb=float(m.cfg.k_mature_bath);dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    seqT=ae.source_sequence(m,True,a.lag,a.steps);seqD=ae.source_sequence(m,False,a.lag,a.steps)

    # Fixed fabrication loss disorder for this physical body across all arms.
    epsmap=make_loss_map(m,a.mean_eps,a.loss_cv,np.random.default_rng(int(m.cfg.seed)+8_818_181))
    body=m.body.astype(bool); actual_mean=float(np.mean(epsmap[body]));actual_cv=float(np.std(epsmap[body])/(actual_mean+1e-30))
    eps_hat=float(a.mean_eps)*(1.0+float(a.loss_cal_error))

    arms={}
    specs=[('nominal_exact',101),('physical_exact',202),('echo_combined',303),('shuffled_echo',404)]
    for mode,off in specs:
        arms[mode]=run_arm(m,base_wh,base_wv,cands,seqT,seqD,epsmap,a.eta,a.iterations,kb,dk,mode,
                           eps_hat,a.mirror_alpha,a.operator_drift,a.readout_noise,off)
    return dict(seed=int(m.cfg.seed),cells=int(m.body.sum()),n_candidates=len(cands),
                actual_mean_eps=actual_mean,actual_loss_cv=actual_cv,eps_hat=eps_hat,arms=arms)


def summarize(rows):
    modes=['nominal_exact','physical_exact','echo_combined','shuffled_echo']
    out=dict(bodies=len(rows),arms={})
    for m in modes:
        dp=np.asarray([r['arms'][m]['delta_physical_C'] for r in rows],float)
        dn=np.asarray([r['arms'][m]['delta_nominal_C'] for r in rows],float)
        cs=[r['arms'][m]['mean_echo_corr'] for r in rows if np.isfinite(r['arms'][m]['mean_echo_corr'])]
        out['arms'][m]=dict(mean_delta_physical=float(np.mean(dp)),median_delta_physical=float(np.median(dp)),
                            improved_physical=int(np.sum(dp>0)),
                            mean_delta_nominal=float(np.mean(dn)),improved_nominal=int(np.sum(dn>0)),
                            mean_echo_corr=float(np.mean(cs)) if cs else float('nan'))
    eh=np.asarray([r['arms']['echo_combined']['delta_physical_C'] for r in rows],float)
    sh=np.asarray([r['arms']['shuffled_echo']['delta_physical_C'] for r in rows],float)
    ex=np.asarray([r['arms']['physical_exact']['delta_physical_C'] for r in rows],float)
    out['echo_minus_shuffle']=dict(mean=float(np.mean(eh-sh)),positive=int(np.sum(eh>sh)))
    out['echo_retained_fraction_of_exact']=float(np.mean(eh)/(np.mean(ex)+1e-30))
    return out


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472);ap.add_argument('--seeds',type=int,default=8)
    ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-candidates',type=int,default=8);ap.add_argument('--eta',type=float,default=.01);ap.add_argument('--iterations',type=int,default=24)
    ap.add_argument('--mean-eps',type=float,default=.005);ap.add_argument('--loss-cv',type=float,default=.20);ap.add_argument('--loss-cal-error',type=float,default=.05)
    ap.add_argument('--mirror-alpha',type=float,default=.95);ap.add_argument('--operator-drift',type=float,default=.02);ap.add_argument('--readout-noise',type=float,default=.05)
    ap.add_argument('--out',default='runs/damping_gauge_combined_closed_loop/dev_472_479.json')
    return ap.parse_args()


def main():
    a=parse_args();fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=body_probe(m,a)
        if r is None:continue
        rows.append(r)
        A=r['arms']
        print(f"seed {seed}: phys-exact={A['physical_exact']['delta_physical_C']:+.4f} echo={A['echo_combined']['delta_physical_C']:+.4f} shuffled={A['shuffled_echo']['delta_physical_C']:+.4f} corr={A['echo_combined']['mean_echo_corr']:.4f}",flush=True)
    if not rows:raise SystemExit('No valid bodies')
    out=dict(experiment='combined_echo_closed_loop_dev_v01',config=vars(a),rows=rows,summary=summarize(rows))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nCOMBINED ECHO CLOSED LOOP')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
