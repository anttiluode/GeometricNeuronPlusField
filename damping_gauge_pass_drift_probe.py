"""How sensitive is conformal echo gradient readout to +/- pass-to-pass drift?

Development-only on reused bodies.

The combined closed-loop result used one common reverse-operator mismatch for
both +/- interference states.  That is realistic for a quasi-static calibrated
mismatch, but optimistic if the two phase states require separate physical
trials and the mesh drifts between them.

This probe keeps the already hostile common conditions:

    mean residual loss .005
    spatial loss CV .20
    scalar loss calibration error +5%
    terminal mirror alpha .95
    common reverse operator mismatch sigma .02

Then it adds *independent differential coupling drift* to the + and - reverse
trials.  Because the gradient is a small difference of two branch energies,
pass-to-pass drift can leak the much larger self-energy terms into the result.

For each differential drift scale we compare the complete target+distractor
physical gradient against the exact digital gradient of the actual lossy
physical device.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import analog_frontier_learning as afl  # mature-boundary patch
import adjoint_eligibility_probe as ae
import damping_gauge_reversal_probe as dg
import damping_gauge_loss_heterogeneity_probe as lh
import damping_gauge_combined_closed_loop as cc
from transfer_decomposition_probe import safe_corr


def flat(h,v):return np.concatenate([np.ravel(h),np.ravel(v)])
def rel(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))
def cosine(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))


def perturb_from(wh,wv,sigma,rng):
    if sigma<=0:return wh.copy(),wv.copy()
    ph=np.exp(float(sigma)*rng.standard_normal(wh.shape)-.5*float(sigma)**2)
    pv=np.exp(float(sigma)*rng.standard_normal(wv.shape)-.5*float(sigma)**2)
    return wh*ph,wv*pv


def trial_energy(m,wh,wv,x,u,q,epsmap,eps_hat,mirror_alpha,sign):
    y=cc.damaged_compensated_echo(m,wh,wv,x,u,epsmap,eps_hat,mirror_alpha)
    b=lh.physical_adjoint_map(m,wh,wv,q,epsmap)
    field=y[1:]+float(sign)*b[1:]
    fh,fv=dg.edge_diffs(field)
    ah=1.0-float(eps_hat);T=len(field)
    wt=(ah**(-np.arange(1,T+1,dtype=float))).reshape((T,1,1))
    Eh=np.sum(wt*np.abs(fh)**2,axis=0)
    Ev=np.sum(wt*np.abs(fv)**2,axis=0)
    return Eh,Ev


def two_pass_order(m,base_wh,base_wv,x,u,q,epsmap,eps_hat,mirror_alpha,common_sigma,diff_sigma,rng):
    # Fabrication/calibration mismatch common to both phase states.
    cwh,cwv=perturb_from(base_wh,base_wv,common_sigma,rng)
    # Differential pass drift around that common reverse operator.
    pwh,pwv=perturb_from(cwh,cwv,diff_sigma,rng)
    mwh,mwv=perturb_from(cwh,cwv,diff_sigma,rng)
    Ep=trial_energy(m,pwh,pwv,x,u,q,epsmap,eps_hat,mirror_alpha,+1)
    Em=trial_energy(m,mwh,mwv,x,u,q,epsmap,eps_hat,mirror_alpha,-1)
    *_,beta,_=dg.params(m)
    gh=(-.5*beta*(Ep[0]-Em[0])).real
    gv=(-.5*beta*(Ep[1]-Em[1])).real
    return gh,gv


def one_body(m,a,diffs):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,a.lag,a.steps);seqD=ae.source_sequence(m,False,a.lag,a.steps)
    epsmap=lh.make_epsmap(m,a.mean_eps,a.loss_cv,np.random.default_rng(int(m.cfg.seed)+5_551_337))
    eps_hat=a.mean_eps*(1+a.loss_cal_error)
    st=cc.physical_state(m,wh,wv,seqT,seqD,epsmap)
    refh,refv=cc.exact_physical_gradient(m,wh,wv,st,epsmap);ref=flat(refh,refv)
    rows=[]
    for d in diffs:
        for rep in range(a.reps):
            rng=np.random.default_rng(int(m.cfg.seed)*1_000_003+rep*9127+int(round(d*1e8)))
            gt=two_pass_order(m,wh,wv,st['xT'],st['uT'],st['qT'],epsmap,eps_hat,a.mirror_alpha,a.common_drift,d,rng)
            gd=two_pass_order(m,wh,wv,st['xD'],st['uD'],st['qD'],epsmap,eps_hat,a.mirror_alpha,a.common_drift,d,rng)
            got=flat(gt[0]+gd[0],gt[1]+gd[1])
            rows.append(dict(seed=int(m.cfg.seed),differential_drift=float(d),rep=rep,
                             corr=float(safe_corr(ref,got)),cosine=cosine(ref,got),relative_l2=rel(ref,got),
                             norm_ratio=float(np.linalg.norm(got)/(np.linalg.norm(ref)+1e-30))))
    return rows


def summarize(rows):
    out={}
    for d in sorted(set(r['differential_drift'] for r in rows)):
        q=[r for r in rows if r['differential_drift']==d]
        out[str(d)]=dict(n=len(q),mean_corr=float(np.mean([x['corr'] for x in q])),min_corr=float(np.min([x['corr'] for x in q])),
                         mean_cosine=float(np.mean([x['cosine'] for x in q])),mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
                         median_norm_ratio=float(np.median([x['norm_ratio'] for x in q])))
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--mean-eps',type=float,default=.005);ap.add_argument('--loss-cv',type=float,default=.20);ap.add_argument('--loss-cal-error',type=float,default=.05)
    ap.add_argument('--mirror-alpha',type=float,default=.95);ap.add_argument('--common-drift',type=float,default=.02)
    ap.add_argument('--differential-drifts',default='0,0.00001,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01')
    ap.add_argument('--reps',type=int,default=4);ap.add_argument('--out',default='runs/damping_gauge_pass_drift/dev_472_475.json');a=ap.parse_args()
    diffs=[float(x) for x in a.differential_drifts.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;rows.extend(one_body(m,a,diffs));print('seed',seed,'done',flush=True)
    out=dict(config=vars(a),rows=rows,summary=summarize(rows));q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2));print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
