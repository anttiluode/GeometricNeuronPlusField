"""Development probe for TIME_DOMAIN_IN_SITU_WALL.md.

The same reciprocal medium generates the exact adjoint field, but can an O(1)-trial
TRIM-like *same-time* local interference measurement recover the dynamic internal
bond gradient without storing/reversing each local forward history?

Compare:
  exact/aligned      stored forward[n] x adjoint[n]
  causal+storedR     physically replayed adjoint[r] x stored forward[T-r] (positive control)
  causal+forward     replayed adjoint[r] x ordinary forward[r] (simple simultaneous trial)
  causal+revinput    replayed adjoint[r] x response to time-reversed INPUT waveform[r]
                     (naive attempt to avoid local forward memory)
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from reciprocal_adjoint_probe import retro_source_sequence,flat_pair,normalized_l2
from transfer_decomposition_probe import safe_corr


def raw_overlap(m,forward_samples,retro_samples):
    """Arrays length T of local field states to multiply samplewise."""
    dt=float(m.cfg.dt);stiff=float(m.cfg.stiffness)
    gh=np.zeros((m.body.shape[0],m.body.shape[1]-1),float)
    gv=np.zeros((m.body.shape[0]-1,m.body.shape[1]),float)
    for f,r in zip(forward_samples,retro_samples):
        dfh=f[:,1:]-f[:,:-1];drh=r[:,:-1]-r[:,1:]
        gh += 2*dt*stiff*np.real(np.conj(drh)*dfh)
        dfv=f[1:,:]-f[:-1,:];drv=r[:-1,:]-r[1:,:]
        gv += 2*dt*stiff*np.real(np.conj(drv)*dfv)
    return gh,gv


def order_fields(m,wh,wv,seq,weight):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    g=weight*np.asarray(p[1:,m.soma[0],m.soma[1]],np.complex128)
    retro_seq=retro_source_sequence(m,g,reverse=True)
    rp,rv,_=ae.linear_forward(m,wh,wv,retro_seq,store=True)
    T=len(seq)
    # exact algorithmic adjoint for this weighted energy contribution
    eh,ev=ae.adjoint_grad(m,wh,wv,p,v,weight)
    # Positive-control physical timeline: at retro state r=1..T use original forward state T-r.
    sh,sv=raw_overlap(m,[p[T-r] for r in range(1,T+1)],[rp[r] for r in range(1,T+1)])
    # Memory-free simple simultaneous pairing: ordinary forward timeline with retro timeline.
    qh,qv=raw_overlap(m,[p[r-1] for r in range(1,T+1)],[rp[r] for r in range(1,T+1)])
    # Naive reversed-input attempt: reverse external source sequence, but start medium from rest.
    revseq=list(seq[::-1]);pr,vr,_=ae.linear_forward(m,wh,wv,revseq,store=True)
    rh,rvv=raw_overlap(m,[pr[r-1] for r in range(1,T+1)],[rp[r] for r in range(1,T+1)])
    # How well did reversed input reproduce the actual reversed local forward state, globally?
    wanted=np.asarray([p[T-r] for r in range(1,T+1)])
    got=np.asarray([pr[r-1] for r in range(1,T+1)])
    hist_rel=float(np.linalg.norm(got-wanted)/(np.linalg.norm(wanted)+1e-30))
    return dict(exact=(eh,ev),storedR=(sh,sv),simultaneous=(qh,qv),revinput=(rh,rvv),revinput_history_rel=hist_rel,E=E)


def one(m,lag,steps):
    wh,wv=ae.bond_weights(m,m.body);seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    # First get unweighted energies for contrast objective weights.
    _,_,ET=ae.linear_forward(m,wh,wv,seqT,store=False);_,_,ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30;aT=2*ED/S**2;aD=-2*ET/S**2
    T=order_fields(m,wh,wv,seqT,aT);D=order_fields(m,wh,wv,seqD,aD)
    def add(name):return (T[name][0]+D[name][0],T[name][1]+D[name][1])
    exact=add('exact');methods={k:add(k) for k in ('storedR','simultaneous','revinput')};ex=flat_pair(*exact)
    out={}
    for k,pair in methods.items():
        z=flat_pair(*pair);out[k]=dict(corr=float(safe_corr(ex,z)),relative_l2=normalized_l2(ex,z))
    return dict(seed=int(m.cfg.seed),C=float((ET-ED)/S),methods=out,
                target_revinput_history_rel=float(T['revinput_history_rel']),distractor_revinput_history_rel=float(D['revinput_history_rel']))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--out',default='runs/temporal_trim_wall/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        x=np.arange(5);assert [x[5-r] for r in range(1,6)]==[4,3,2,1,0];print('selftest ok');return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps);rows.append(r);print('seed',seed,{k:(round(v['corr'],4),round(v['relative_l2'],4)) for k,v in r['methods'].items()},'revhist',round(r['target_revinput_history_rel'],3),flush=True)
    summary=dict(bodies=len(rows))
    for k in ('storedR','simultaneous','revinput'):
        summary[k]=dict(mean_corr=float(np.mean([r['methods'][k]['corr'] for r in rows])),mean_relative_l2=float(np.mean([r['methods'][k]['relative_l2'] for r in rows])),min_corr=float(np.min([r['methods'][k]['corr'] for r in rows])))
    summary['mean_revinput_history_rel']=float(np.mean([(r['target_revinput_history_rel']+r['distractor_revinput_history_rel'])/2 for r in rows]))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='temporal_trim_wall_dev_v01',summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
