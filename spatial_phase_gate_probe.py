"""Spatial phase gating of the local forward x return structural signal.

Motivation
----------
In-vivo CA1 voltage imaging reports a smooth intracellular theta phase gradient across
pyramidal-cell morphology (Liao et al., Nature Communications 2024).  That motivates a
strictly narrower model question: if plasticity is periodically gated, does a smooth
morphology-indexed phase field preserve/use the local structural signal differently
from the same phases shuffled across the arbor?

This probe does NOT inject a theta waveform at the soma and does NOT alter forward or
adjoint propagation.  It first computes the exact reciprocal local overlap density

    c_e(t) = 2 dt K Re[conj(dmu_e(t)) dpsi_e(t)]

for target+distractor.  The exact bond gradient is sum_t c_e(t).

Then a 50%-duty local gate is applied only to structural accumulation:

    grad_e(gated) = 2 sum_t m_e(t) c_e(t).

The factor 2 compensates the common 50% duty; correlation is scale-invariant anyway.

At period P=42 (chosen because global gating is known to be phase-sensitive), compare:

  global     one phase for all bonds
  smooth     phase = phi0 + alpha * graph_distance_from_soma
  shuffled   same smooth phase values, randomly permuted over real bonds
  random     independent uniform phase per real bond

`alpha` is in simulation frames of phase shift per graph edge.  No mapping to biological
Hz or micrometers is implied.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # mature-boundary operator patch
from reciprocal_adjoint_probe import physical_credit_history, flat_pair, normalized_l2
from transfer_decomposition_probe import safe_corr

PERIOD = 42
ALPHAS = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0)
N_SHUFFLES = 8
N_RANDOM = 8


def bond_mask(body):
    body = np.asarray(body, bool)
    return body[:, :-1] & body[:, 1:], body[:-1, :] & body[1:, :]


def graph_distances(body, soma):
    body = np.asarray(body, bool)
    H, W = body.shape
    d = np.full((H, W), np.inf, float)
    s = tuple(map(int, soma))
    d[s] = 0.0
    q = deque([s])
    while q:
        r, c = q.popleft()
        nd = d[r, c] + 1.0
        for rr, cc in ((r-1,c),(r+1,c),(r,c-1),(r,c+1)):
            if 0 <= rr < H and 0 <= cc < W and body[rr,cc] and nd < d[rr,cc]:
                d[rr,cc] = nd
                q.append((rr,cc))
    return d


def bond_distances(body, soma):
    d = graph_distances(body, soma)
    mh, mv = bond_mask(body)
    dh = 0.5 * (d[:, :-1] + d[:, 1:])
    dv = 0.5 * (d[:-1, :] + d[1:, :])
    dh[~mh] = 0.0
    dv[~mv] = 0.0
    return dh, dv, mh, mv


def overlap_time(m, forward_ps, credit_mu):
    z = np.asarray(forward_ps, np.complex128)
    mu = np.asarray(credit_mu, np.complex128)
    T = len(mu)
    dt = float(m.cfg.dt)
    stiff = float(m.cfg.stiffness)
    ch = np.zeros((T, m.body.shape[0], m.body.shape[1]-1), float)
    cv = np.zeros((T, m.body.shape[0]-1, m.body.shape[1]), float)
    for n in range(T):
        prev = z[n]
        mm = mu[n]
        dpsi_h = prev[:,1:] - prev[:,:-1]
        dmu_h = mm[:,:-1] - mm[:,1:]
        dpsi_v = prev[1:,:] - prev[:-1,:]
        dmu_v = mm[:-1,:] - mm[1:,:]
        ch[n] = 2.0*dt*stiff*np.real(np.conj(dmu_h)*dpsi_h)
        cv[n] = 2.0*dt*stiff*np.real(np.conj(dmu_v)*dpsi_v)
    return ch, cv


def exact_local_signal(m, lag, steps):
    wh, wv = ae.bond_weights(m, m.body)
    seqT = ae.source_sequence(m, True, lag, steps)
    seqD = ae.source_sequence(m, False, lag, steps)
    pT, vT, ET = ae.linear_forward(m, wh, wv, seqT, store=True)
    pD, vD, ED = ae.linear_forward(m, wh, wv, seqD, store=True)
    S = ET + ED + 1e-30
    aT = 2.0*ED/(S*S)
    aD = -2.0*ET/(S*S)
    exT_h, exT_v = ae.adjoint_grad(m, wh, wv, pT, vT, aT)
    exD_h, exD_v = ae.adjoint_grad(m, wh, wv, pD, vD, aD)
    gT = aT*np.asarray(pT[1:,m.soma[0],m.soma[1]], np.complex128)
    gD = aD*np.asarray(pD[1:,m.soma[0],m.soma[1]], np.complex128)
    muT, _, _ = physical_credit_history(m, wh, wv, gT, reverse=True)
    muD, _, _ = physical_credit_history(m, wh, wv, gD, reverse=True)
    cTh, cTv = overlap_time(m, pT, muT)
    cDh, cDv = overlap_time(m, pD, muD)
    ch, cv = cTh+cDh, cTv+cDv
    ex_h, ex_v = exT_h+exD_h, exT_v+exD_v
    # Positive control: local density must sum to the exact adjoint.
    rel = normalized_l2(flat_pair(ex_h,ex_v), flat_pair(np.sum(ch,axis=0),np.sum(cv,axis=0)))
    if rel > 1e-10:
        raise RuntimeError(f'local overlap density mismatch: {rel}')
    return dict(ch=ch, cv=cv, exact_h=ex_h, exact_v=ex_v, C=float((ET-ED)/S))


def gate_from_offsets(T, offsets_h, offsets_v, period=PERIOD):
    t = np.arange(T, dtype=float)[:,None,None]
    on = float(period)/2.0
    gh = np.mod(t + offsets_h[None,:,:], float(period)) < on
    gv = np.mod(t + offsets_v[None,:,:], float(period)) < on
    return gh, gv


def gated_map(ch, cv, offsets_h, offsets_v, period=PERIOD):
    gh, gv = gate_from_offsets(len(ch), offsets_h, offsets_v, period)
    return 2.0*np.sum(ch*gh, axis=0), 2.0*np.sum(cv*gv, axis=0)


def metrics(ex_h, ex_v, h, v):
    ex = flat_pair(ex_h,ex_v)
    ap = flat_pair(h,v)
    mx = float(np.max(np.abs(ex))+1e-30)
    mask = np.abs(ex) > .01*mx
    return dict(
        corr=float(safe_corr(ex,ap)),
        relative_l2=normalized_l2(ex,ap),
        strong_sign=float(np.mean(np.sign(ex[mask])==np.sign(ap[mask]))) if np.any(mask) else float('nan'),
    )


def phase_summary(q):
    return dict(
        mean_corr=float(np.mean([x['corr'] for x in q])),
        median_corr=float(np.median([x['corr'] for x in q])),
        min_corr=float(np.min([x['corr'] for x in q])),
        max_corr=float(np.max([x['corr'] for x in q])),
        mean_sign=float(np.mean([x['strong_sign'] for x in q])),
        mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
    )


def shuffle_offsets(oh, ov, mh, mv, rng):
    vals = np.concatenate([oh[mh], ov[mv]])
    vals = vals[rng.permutation(len(vals))]
    a = int(np.sum(mh))
    xh = np.zeros_like(oh)
    xv = np.zeros_like(ov)
    xh[mh] = vals[:a]
    xv[mv] = vals[a:]
    return xh, xv


def random_offsets(oh, ov, mh, mv, rng, period=PERIOD):
    xh = np.zeros_like(oh)
    xv = np.zeros_like(ov)
    xh[mh] = rng.uniform(0.0, float(period), size=int(np.sum(mh)))
    xv[mv] = rng.uniform(0.0, float(period), size=int(np.sum(mv)))
    return xh, xv


def scan_global_phase(z, base_h, base_v):
    q=[]
    for phi in range(PERIOD):
        h,v=gated_map(z['ch'],z['cv'],base_h+phi,base_v+phi)
        mm=metrics(z['exact_h'],z['exact_v'],h,v); mm['phi']=phi; q.append(mm)
    return phase_summary(q)


def one(m, lag, steps):
    z=exact_local_signal(m,lag,steps)
    dh,dv,mh,mv=bond_distances(m.body,m.soma)
    zeros_h=np.zeros_like(dh); zeros_v=np.zeros_like(dv)
    global_result=scan_global_phase(z,zeros_h,zeros_v)
    rng=np.random.default_rng(int(m.cfg.seed)+420024)
    alphas={}
    for alpha in ALPHAS:
        oh=float(alpha)*dh; ov=float(alpha)*dv
        smooth=scan_global_phase(z,oh,ov)
        sh=[]
        for _ in range(N_SHUFFLES):
            xh,xv=shuffle_offsets(oh,ov,mh,mv,rng)
            sh.append(scan_global_phase(z,xh,xv))
        rr=[]
        for _ in range(N_RANDOM):
            xh,xv=random_offsets(oh,ov,mh,mv,rng)
            rr.append(scan_global_phase(z,xh,xv))
        alphas[str(alpha)]=dict(
            smooth=smooth,
            shuffled=dict(
                mean_mean_corr=float(np.mean([x['mean_corr'] for x in sh])),
                mean_min_corr=float(np.mean([x['min_corr'] for x in sh])),
                best_mean_corr=float(np.max([x['mean_corr'] for x in sh])),
            ),
            random=dict(
                mean_mean_corr=float(np.mean([x['mean_corr'] for x in rr])),
                mean_min_corr=float(np.mean([x['min_corr'] for x in rr])),
                best_mean_corr=float(np.max([x['mean_corr'] for x in rr])),
            )
        )
    return dict(seed=int(m.cfg.seed),C=z['C'],global_gate=global_result,alphas=alphas)


def summarize(rows):
    out=dict(bodies=len(rows),global_gate={},alphas={})
    g=[r['global_gate'] for r in rows]
    out['global_gate']=dict(
        mean_corr=float(np.mean([x['mean_corr'] for x in g])),
        mean_min_corr=float(np.mean([x['min_corr'] for x in g])),
    )
    for a in map(str,ALPHAS):
        q=[r['alphas'][a] for r in rows]
        out['alphas'][a]=dict(
            smooth_mean_corr=float(np.mean([x['smooth']['mean_corr'] for x in q])),
            smooth_mean_min_corr=float(np.mean([x['smooth']['min_corr'] for x in q])),
            shuffled_mean_corr=float(np.mean([x['shuffled']['mean_mean_corr'] for x in q])),
            shuffled_mean_min_corr=float(np.mean([x['shuffled']['mean_min_corr'] for x in q])),
            random_mean_corr=float(np.mean([x['random']['mean_mean_corr'] for x in q])),
            smooth_minus_shuffled=float(np.mean([x['smooth']['mean_corr']-x['shuffled']['mean_mean_corr'] for x in q])),
            smooth_minus_global=float(np.mean([x['smooth']['mean_corr']-r['global_gate']['mean_corr'] for x,r in zip(q,rows)])),
        )
    return out


def selftest():
    body=np.zeros((5,5),bool); body[2,1:4]=1; body[1:4,2]=1
    dh,dv,mh,mv=bond_distances(body,(2,2))
    assert np.isfinite(dh[mh]).all() and np.isfinite(dv[mv]).all()
    oh=.5*dh; ov=.5*dv
    gh,gv=gate_from_offsets(210,oh,ov)
    assert gh.shape[0]==210 and gv.shape[0]==210
    assert np.allclose(np.mean(gh[:,mh],axis=0),.5)
    assert np.allclose(np.mean(gv[:,mv],axis=0),.5)
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=528)
    ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/spatial_phase_gate/dev.json')
    ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest: selftest(); return
    fa=Path(a.functional_arbors).resolve(); sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); b=m.bootstrap()
        if not b.get('ok'): continue
        m.mature=True; r=one(m,a.lag,a.steps); rows.append(r)
        best=max(r['alphas'].items(), key=lambda kv:kv[1]['smooth']['mean_corr'])
        print('seed',seed,'global',round(r['global_gate']['mean_corr'],3),
              'best-alpha',best[0],round(best[1]['smooth']['mean_corr'],3),
              'shuffle',round(best[1]['shuffled']['mean_mean_corr'],3),flush=True)
    if not rows: raise SystemExit('No valid bodies')
    s=summarize(rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(dict(experiment='spatial_phase_gate_dev_v01',summary=s,rows=rows),indent=2))
    print('\nSPATIAL PHASE GATE DEV')
    print(json.dumps(s,indent=2))


if __name__=='__main__': main()
