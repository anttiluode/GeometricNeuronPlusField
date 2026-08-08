"""Mechanism diagnostic for the morphology-indexed impedance candidate.

Development-only.  The synchrony probe found a candidate region where a smooth
soma->distal delayed-restorative gradient reduces location-dependent somatic
phase spread relative to uniform and shuffled same-material controls.

This script asks whether that improvement has the *HCN-like sign*:

  - distal inputs should acquire a phase advance at the soma relative to the
    no-gradient / uniform-material control;
  - the phase advance should increase with graph distance;
  - amplitude changes alone should not explain the phase correction;
  - local dendritic responses may remain phase-rich.

The model is still a quasi-active linear proxy, not a conductance-based HCN
current and not a biological claim.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import (
    weighted_laplacian_sparse, build_profiles, harmonic_matrix,
    injection_sites, circular_rms,
)


def safe_corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a = a[ok]; b = b[ok]
    if len(a) < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]; y = y[ok]
    if len(x) < 2 or np.var(x) < 1e-12:
        return float('nan')
    return float(np.polyfit(x, y, 1)[0])


def solve_transfer_vectors(m, L, gfield, omega, tau_h, sites):
    A = harmonic_matrix(m, L, gfield, omega, tau_h).tocsc()
    lu = splu(A)
    h, w = m.body.shape
    N = h * w
    idx = np.asarray([p[0] * w + p[1] for p in sites], dtype=int)
    B = np.zeros((N, len(sites)), np.complex128)
    B[idx, np.arange(len(sites))] = 1.0
    X = lu.solve(B)
    si = m.soma[0] * w + m.soma[1]
    hs = X[si, np.arange(len(sites))]
    hl = X[idx, np.arange(len(sites))]
    return hs, hl


def circ_delta(a, b):
    return np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b))))


def body_probe(m, g0, ratio, tau_h, omega, nshuffle):
    wh, wv = ae.bond_weights(m, m.body)
    L = weighted_laplacian_sparse(wh, wv)
    sites = injection_sites(m)
    distmap = m.graph_distance_from_soma()
    dist = np.asarray([distmap[p] for p in sites], float)
    rng = np.random.default_rng(int(m.cfg.seed) + 919_919)
    prof = build_profiles(m, g0, ratio, rng, nshuffle)

    hz, lz = solve_transfer_vectors(m, L, prof['zero'], omega, tau_h, sites)
    hu, lu = solve_transfer_vectors(m, L, prof['uniform'], omega, tau_h, sites)
    hs, ls = solve_transfer_vectors(m, L, prof['smooth'], omega, tau_h, sites)
    hr, lr = solve_transfer_vectors(m, L, prof['reverse'], omega, tau_h, sites)

    hsh = []; lsh = []
    for z in prof['shuffled']:
        a, b = solve_transfer_vectors(m, L, z, omega, tau_h, sites)
        hsh.append(a); lsh.append(b)
    hsh = np.asarray(hsh); lsh = np.asarray(lsh)

    # Phase advance of smooth relative to controls, wrapped into [-pi,pi].
    adv_zero = circ_delta(np.angle(hs), np.angle(hz))
    adv_uniform = circ_delta(np.angle(hs), np.angle(hu))
    # Circular mean of shuffled complex transfer then compare phase.
    hsh_mean = np.mean(hsh, axis=0)
    adv_shuffle = circ_delta(np.angle(hs), np.angle(hsh_mean))

    amp_ratio_zero = np.abs(hs) / (np.abs(hz) + 1e-30)
    amp_ratio_uniform = np.abs(hs) / (np.abs(hu) + 1e-30)

    # Distal quartile versus proximal quartile is less sensitive than a single
    # linear slope to branch-specific phase wrapping.
    q25, q75 = np.quantile(dist, [0.25, 0.75])
    prox = dist <= q25
    distal = dist >= q75

    def qsummary(x):
        return dict(
            mean=float(np.mean(x)), median=float(np.median(x)),
            prox_mean=float(np.mean(x[prox])), distal_mean=float(np.mean(x[distal])),
            distal_minus_prox=float(np.mean(x[distal]) - np.mean(x[prox])),
            positive_fraction=float(np.mean(x > 0)),
            corr_distance=safe_corr(dist, x), slope_per_edge=slope(dist, x),
        )

    return dict(
        seed=int(m.cfg.seed), cells=int(m.body.sum()), sites=len(sites),
        max_graph_distance=float(np.max(dist)),
        soma_phase_rms=dict(
            zero=circular_rms(np.angle(hz)), uniform=circular_rms(np.angle(hu)),
            smooth=circular_rms(np.angle(hs)), reverse=circular_rms(np.angle(hr)),
            shuffle=float(np.mean([circular_rms(np.angle(x)) for x in hsh])),
        ),
        local_phase_rms=dict(
            zero=circular_rms(np.angle(lz)), uniform=circular_rms(np.angle(lu)),
            smooth=circular_rms(np.angle(ls)), reverse=circular_rms(np.angle(lr)),
            shuffle=float(np.mean([circular_rms(np.angle(x)) for x in lsh])),
        ),
        advance_vs_zero=qsummary(adv_zero),
        advance_vs_uniform=qsummary(adv_uniform),
        advance_vs_shuffle=qsummary(adv_shuffle),
        amplitude_ratio_vs_zero=qsummary(np.log(np.maximum(amp_ratio_zero, 1e-30))),
        amplitude_ratio_vs_uniform=qsummary(np.log(np.maximum(amp_ratio_uniform, 1e-30))),
        per_site=[dict(
            cell=[int(p[0]), int(p[1])], distance=float(d),
            advance_zero=float(a0), advance_uniform=float(au), advance_shuffle=float(ash),
            abs_smooth=float(abs(h)), abs_zero=float(abs(h0)), abs_uniform=float(abs(h1)),
        ) for p,d,a0,au,ash,h,h0,h1 in zip(sites,dist,adv_zero,adv_uniform,adv_shuffle,hs,hz,hu)],
    )


def summarize(rows):
    def meanpath(*keys):
        z=[]
        for r in rows:
            q=r
            for k in keys:q=q[k]
            z.append(float(q))
        return float(np.mean(z))
    def vals(*keys):
        out=[]
        for r in rows:
            q=r
            for k in keys:q=q[k]
            out.append(float(q))
        return np.asarray(out,float)

    out=dict(
        bodies=len(rows),
        smooth_soma_phase_rms=meanpath('soma_phase_rms','smooth'),
        zero_soma_phase_rms=meanpath('soma_phase_rms','zero'),
        uniform_soma_phase_rms=meanpath('soma_phase_rms','uniform'),
        shuffle_soma_phase_rms=meanpath('soma_phase_rms','shuffle'),
        smooth_local_phase_rms=meanpath('local_phase_rms','smooth'),
        zero_local_phase_rms=meanpath('local_phase_rms','zero'),
        smooth_gain_vs_uniform=float(np.mean(vals('soma_phase_rms','uniform')-vals('soma_phase_rms','smooth'))),
        smooth_gain_vs_shuffle=float(np.mean(vals('soma_phase_rms','shuffle')-vals('soma_phase_rms','smooth'))),
        smooth_gain_vs_zero=float(np.mean(vals('soma_phase_rms','zero')-vals('soma_phase_rms','smooth'))),
        local_phase_retention=float(meanpath('local_phase_rms','smooth')/(meanpath('local_phase_rms','zero')+1e-30)),
        phase_advance_zero_mean=meanpath('advance_vs_zero','mean'),
        phase_advance_zero_distal_minus_prox=meanpath('advance_vs_zero','distal_minus_prox'),
        phase_advance_zero_corr_distance=meanpath('advance_vs_zero','corr_distance'),
        phase_advance_uniform_mean=meanpath('advance_vs_uniform','mean'),
        phase_advance_uniform_distal_minus_prox=meanpath('advance_vs_uniform','distal_minus_prox'),
        phase_advance_uniform_corr_distance=meanpath('advance_vs_uniform','corr_distance'),
        phase_advance_shuffle_distal_minus_prox=meanpath('advance_vs_shuffle','distal_minus_prox'),
        phase_advance_shuffle_corr_distance=meanpath('advance_vs_shuffle','corr_distance'),
        log_amp_ratio_zero_corr_distance=meanpath('amplitude_ratio_vs_zero','corr_distance'),
        log_amp_ratio_uniform_corr_distance=meanpath('amplitude_ratio_vs_uniform','corr_distance'),
    )
    return out


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=552)
    ap.add_argument('--seeds',type=int,default=8)
    ap.add_argument('--g0',type=float,default=.02)
    ap.add_argument('--ratio',type=float,default=7.0)
    ap.add_argument('--tau',type=float,default=2.0)
    ap.add_argument('--omega',type=float,default=.04)
    ap.add_argument('--nshuffle',type=int,default=8)
    ap.add_argument('--out',default='runs/hcn_inductive/dev.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    x=np.array([-.1,.1]);y=np.array([3.1,-3.1])
    d=circ_delta(x,y)
    assert np.all(np.abs(d)<=np.pi+1e-12)
    assert abs(safe_corr([1,2,3],[1,2,3])-1)<1e-12
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
        r=body_probe(m,a.g0,a.ratio,a.tau,a.omega,a.nshuffle);rows.append(r)
        print(f"seed {seed}: soma smooth={r['soma_phase_rms']['smooth']:.4f} "
              f"shuf={r['soma_phase_rms']['shuffle']:.4f} unif={r['soma_phase_rms']['uniform']:.4f} "
              f"adv0 distal-prox={r['advance_vs_zero']['distal_minus_prox']:+.4f} "
              f"r_dist={r['advance_vs_zero']['corr_distance']:+.3f}",flush=True)
    if not rows:raise SystemExit('No valid bodies')
    s=summarize(rows)
    payload=dict(experiment='hcn_inductive_phase_diagnostic_v01',development_only=True,
                 seed_start=a.seed_start,seeds_requested=a.seeds,g0=a.g0,ratio=a.ratio,
                 tau_h=a.tau,omega=a.omega,nshuffle=a.nshuffle,summary=s,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nHCN INDUCTIVE DIAGNOSTIC')
    print(json.dumps(s,indent=2))

if __name__=='__main__':main()
