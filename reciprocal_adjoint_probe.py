"""Same-medium time-reversed soma wave as the exact discrete adjoint credit field.

See RECIPROCAL_ADJOINT_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # patches ae.weighted_lap to exact mature-boundary operator
from structural_interference_probe import event_candidates, n4
from transfer_decomposition_probe import safe_corr


def retro_source_sequence(m, soma_wave, reverse=True):
    """Acceleration-source sequence s=g/dt placed only at the soma."""
    g=np.asarray(soma_wave,np.complex128)
    if reverse:
        g=g[::-1]
    out=[]
    dt=float(m.cfg.dt)
    for x in g:
        src=np.zeros(m.body.shape,np.complex128)
        src[m.soma]=x/dt
        out.append(src)
    return out


def physical_credit_history(m,wh,wv,g_states,reverse=True):
    """Return mu[n] aligned to forward transitions n=0..T-1.

    A reverse-wave run starts at chi_0=mu[T+1]=0.  After its first step,
    chi_1=mu[T].  Therefore transition n uses retro state T-n.
    """
    seq=retro_source_sequence(m,g_states,reverse=reverse)
    ps,vs,_=ae.linear_forward(m,wh,wv,seq,store=True)
    T=len(g_states)
    mu=np.asarray([ps[T-n] for n in range(T)],np.complex128)
    return mu,ps,vs


def overlap_gradient(m,forward_ps,credit_mu):
    """Bond gradient from local forward × retrograde field overlap."""
    dt=float(m.cfg.dt); stiff=float(m.cfg.stiffness)
    gh=np.zeros((m.body.shape[0],m.body.shape[1]-1),float)
    gv=np.zeros((m.body.shape[0]-1,m.body.shape[1]),float)
    T=len(credit_mu)
    for n in range(T):
        prev=forward_ps[n]
        mu=credit_mu[n]
        dpsi_h=prev[:,1:]-prev[:,:-1]
        dmu_h=mu[:,:-1]-mu[:,1:]
        gh += 2.0*dt*stiff*np.real(np.conj(dmu_h)*dpsi_h)
        dpsi_v=prev[1:,:]-prev[:-1,:]
        dmu_v=mu[:-1,:]-mu[1:,:]
        gv += 2.0*dt*stiff*np.real(np.conj(dmu_v)*dpsi_v)
    return gh,gv


def exact_and_reciprocal(m,lag,steps,reverse=True):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,lag,steps)
    seqD=ae.source_sequence(m,False,lag,steps)
    pT,vT,ET=ae.linear_forward(m,wh,wv,seqT,store=True)
    pD,vD,ED=ae.linear_forward(m,wh,wv,seqD,store=True)
    S=ET+ED+1e-30
    aT=2.0*ED/(S*S)
    aD=-2.0*ET/(S*S)

    ghT,gvT=ae.adjoint_grad(m,wh,wv,pT,vT,aT)
    ghD,gvD=ae.adjoint_grad(m,wh,wv,pD,vD,aD)
    exact_h=ghT+ghD; exact_v=gvT+gvD

    # Output derivative added at states 1..T in the exact adjoint.
    gT=aT*np.asarray(pT[1:,m.soma[0],m.soma[1]],np.complex128)
    gD=aD*np.asarray(pD[1:,m.soma[0],m.soma[1]],np.complex128)
    muT,_,_=physical_credit_history(m,wh,wv,gT,reverse=reverse)
    muD,_,_=physical_credit_history(m,wh,wv,gD,reverse=reverse)
    rhT,rvT=overlap_gradient(m,pT,muT)
    rhD,rvD=overlap_gradient(m,pD,muD)
    replay_h=rhT+rhD; replay_v=rvT+rvD

    return dict(wh=wh,wv=wv,seqT=seqT,seqD=seqD,pT=pT,pD=pD,
                exact_h=exact_h,exact_v=exact_v,replay_h=replay_h,replay_v=replay_v,
                C=ae.contrast(ET,ED),ET=ET,ED=ED)


def flat_pair(h,v):
    return np.concatenate([np.asarray(h,float).ravel(),np.asarray(v,float).ravel()])


def normalized_l2(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))


def addition_events(m,max_candidates,rng):
    adds,_,_,_=event_candidates(m,int(max_candidates),rng)
    body=m.body.astype(bool);events=[]
    for p in adds:
        p=tuple(map(int,p))
        qs=[q for q in n4(*p,body.shape) if body[q]]
        if len(qs)==1:
            events.append(dict(kind='add',cell=[p[0],p[1]]))
    return events


def norm_step(g,eta=.01):
    g=np.asarray(g,float)
    mx=float(np.max(np.abs(g))) if len(g) else 0.0
    if mx<=1e-30:return np.zeros_like(g)
    return np.clip(float(eta)*g/mx,0.0,1.0)


def body_probe(m,lag,steps,max_candidates):
    z=exact_and_reciprocal(m,lag,steps,reverse=True)
    zn=exact_and_reciprocal(m,lag,steps,reverse=False)
    ex=flat_pair(z['exact_h'],z['exact_v']); rp=flat_pair(z['replay_h'],z['replay_v'])
    nrp=flat_pair(zn['replay_h'],zn['replay_v'])

    rng=np.random.default_rng(int(m.cfg.seed)+264264)
    events=addition_events(m,max_candidates,rng)
    dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    gx=np.asarray([ae.event_score(z['exact_h'],z['exact_v'],m.body,e,dk) for e in events],float)
    gr=np.asarray([ae.event_score(z['replay_h'],z['replay_v'],m.body,e,dk) for e in events],float)
    scale=float(np.max(np.abs(gx))+1e-30) if len(gx) else 1.0
    max_norm_err=float(np.max(np.abs(gr-gx))/scale) if len(gx) else float('nan')
    rho_x=norm_step(gx,.01);rho_r=norm_step(gr,.01)

    return dict(
        seed=int(m.cfg.seed),cells=int(m.body.sum()),C=float(z['C']),
        bond_map_corr=float(safe_corr(ex,rp)),
        bond_map_relative_l2=normalized_l2(ex,rp),
        nonreversed_bond_map_corr=float(safe_corr(ex,nrp)),
        nonreversed_bond_map_relative_l2=normalized_l2(ex,nrp),
        frontier_n=len(events),
        frontier_corr=float(safe_corr(gx,gr)) if len(events)>=2 else float('nan'),
        frontier_max_normalized_abs_error=max_norm_err,
        one_step_max_abs_rho_difference=float(np.max(np.abs(rho_x-rho_r))) if len(events) else float('nan'),
        exact_frontier=[float(x) for x in gx],
        reciprocal_frontier=[float(x) for x in gr],
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=240)
    ap.add_argument('--seeds',type=int,default=2)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-candidates',type=int,default=8)
    ap.add_argument('--out',default='runs/reciprocal_adjoint/reciprocal_adjoint.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    g=np.array([-2.,1.,.5])
    r=norm_step(g,.01)
    assert abs(r[1]-.005)<1e-12 and r[0]==0
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
        r=body_probe(m,a.lag,a.steps,a.max_candidates);rows.append(r)
        print(f"seed {seed}: map r={r['bond_map_corr']:+.9f} rel={r['bond_map_relative_l2']:.2e} "
              f"front r={r['frontier_corr']:+.9f} step={r['one_step_max_abs_rho_difference']:.2e} "
              f"no-rev r={r['nonreversed_bond_map_corr']:+.3f}",flush=True)

    if not rows:raise SystemExit('No valid bodies')
    summary=dict(
        bodies=len(rows),
        mean_bond_map_corr=float(np.nanmean([r['bond_map_corr'] for r in rows])),
        min_bond_map_corr=float(np.nanmin([r['bond_map_corr'] for r in rows])),
        mean_bond_map_relative_l2=float(np.nanmean([r['bond_map_relative_l2'] for r in rows])),
        max_bond_map_relative_l2=float(np.nanmax([r['bond_map_relative_l2'] for r in rows])),
        pooled_frontier_corr=float(safe_corr([x for r in rows for x in r['exact_frontier']],
                                             [x for r in rows for x in r['reciprocal_frontier']])),
        max_frontier_normalized_abs_error=float(np.nanmax([r['frontier_max_normalized_abs_error'] for r in rows])),
        max_one_step_rho_difference=float(np.nanmax([r['one_step_max_abs_rho_difference'] for r in rows])),
        mean_nonreversed_bond_map_corr=float(np.nanmean([r['nonreversed_bond_map_corr'] for r in rows])),
        mean_nonreversed_bond_map_relative_l2=float(np.nanmean([r['nonreversed_bond_map_relative_l2'] for r in rows])),
    )
    summary['R1_pass']=bool(summary['mean_bond_map_corr']>.999999 and summary['mean_bond_map_relative_l2']<1e-8)
    summary['R2_pass']=bool(summary['pooled_frontier_corr']>.999999 and summary['max_frontier_normalized_abs_error']<1e-7)
    summary['R3_pass']=bool(summary['max_one_step_rho_difference']<1e-8)
    payload=dict(experiment='reciprocal_adjoint_identity_v01',prereg='RECIPROCAL_ADJOINT_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,
                 summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nRECIPROCAL ADJOINT RECEIPT')
    for k,v in summary.items():print(f' {k}: {v}')

if __name__=='__main__':main()
