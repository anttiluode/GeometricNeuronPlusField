"""Minimal morphology-indexed inductive-material probe.

This is NOT a conductance-based HCN model.  It tests only the linearized idea
motivated by Vaidya & Johnston (Nat Neurosci 2013, doi:10.1038/nn.3562): a
spatial gradient of delayed restorative/inductive membrane response may
compensate location-dependent transfer phase at a common somatic readout.

The frozen FunctionalArbor geometry is unchanged.  We add an auxiliary local
state z on occupied arbor cells:

    v[n+1]   = v[n] + dt*(K L psi[n] - gamma*v[n] - rho*psi[n]
                            - g_h(x)*z[n] + source[n])
    psi[n+1] = psi[n] + dt*v[n+1]
    z[n+1]   = z[n] + dt/tau_h*(psi[n] - z[n])

For harmonic drive exp(i*omega*n), z can be eliminated exactly, giving a sparse
frequency-domain transfer solve.  The primary comparison keeps the histogram
of g_h values fixed:

    smooth soma->distal gradient
    shuffled same values
    uniform same mean
    reversed same values

We inject the same harmonic current at many occupied locations and measure
phase both locally and at the soma.  The desired biological signature is:

    somatic phase spread decreases for the smooth morphology-indexed gradient
    while local dendritic phase spread remains substantial.

Development use only.  Freeze a candidate before held-out confirmation.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags, eye
from scipy.sparse.linalg import splu

import adjoint_eligibility_probe as ae


def parse_floats(s: str):
    return [float(x) for x in str(s).split(',') if str(x).strip()]


def weighted_laplacian_sparse(wh, wv):
    """Sparse matrix L matching ae.weighted_lap: negative semidefinite."""
    h = wh.shape[0]
    w = wh.shape[1] + 1
    N = h * w
    rows, cols, vals = [], [], []

    def add_edge(i, j, wt):
        wt = float(wt)
        rows.extend([i, i, j, j])
        cols.extend([i, j, j, i])
        vals.extend([-wt, wt, -wt, wt])

    for y in range(h):
        for x in range(w - 1):
            i = y * w + x
            add_edge(i, i + 1, wh[y, x])
    for y in range(h - 1):
        for x in range(w):
            i = y * w + x
            add_edge(i, i + w, wv[y, x])
    return coo_matrix((vals, (rows, cols)), shape=(N, N), dtype=np.complex128).tocsr()


def circular_rms(phases):
    p = np.asarray(phases, float)
    if p.size == 0:
        return float('nan')
    mu = np.angle(np.mean(np.exp(1j * p)))
    d = np.angle(np.exp(1j * (p - mu)))
    return float(np.sqrt(np.mean(d * d)))


def circ_gap(a, b):
    return float(abs(np.angle(np.exp(1j * (float(a) - float(b))))))


def build_profiles(m, g0, ratio, rng, nshuffle):
    body = m.body.astype(bool)
    dist = np.asarray(m.graph_distance_from_soma(), float)
    dmax = max(float(dist[body].max()), 1.0)
    dn = np.zeros_like(dist, float)
    dn[body] = np.clip(dist[body] / dmax, 0.0, 1.0)

    smooth = np.zeros_like(dist, float)
    smooth[body] = float(g0) * (1.0 + (float(ratio) - 1.0) * dn[body])
    vals = smooth[body].copy()

    uniform = np.zeros_like(dist, float)
    uniform[body] = float(np.mean(vals))

    # Same exact histogram, assigned in the opposite graph-distance rank order.
    reverse = np.zeros_like(dist, float)
    cells = [tuple(map(int, p)) for p in np.argwhere(body)]
    cells_sorted = sorted(cells, key=lambda p: (dist[p], p[0], p[1]))
    vals_sorted = np.sort(vals)[::-1]
    for p, v in zip(cells_sorted, vals_sorted):
        reverse[p] = float(v)

    shuffled = []
    for _ in range(int(nshuffle)):
        z = np.zeros_like(dist, float)
        vv = vals.copy()
        rng.shuffle(vv)
        z[body] = vv
        shuffled.append(z)

    zero = np.zeros_like(dist, float)
    return dict(zero=zero, uniform=uniform, smooth=smooth, reverse=reverse, shuffled=shuffled)


def harmonic_matrix(m, L, gfield, omega, tau_h):
    """Exact discrete-time harmonic operator for the linear auxiliary-state model."""
    c = m.cfg
    dt = float(c.dt)
    lam = np.exp(1j * float(omega))
    ah = (dt / float(tau_h)) / (lam - 1.0 + dt / float(tau_h))
    dyn = ((lam - 1.0 + dt * float(c.damping)) * (lam - 1.0) /
           (dt * dt * lam) + float(c.restoring))
    N = m.body.size
    gdiag = diags(np.asarray(gfield, float).ravel() * ah, 0, shape=(N, N), dtype=np.complex128)
    return (dyn * eye(N, dtype=np.complex128, format='csr')
            + gdiag
            - float(c.stiffness) * L)


def injection_sites(m, min_dist=2):
    body = m.body.astype(bool)
    d = m.graph_distance_from_soma()
    pts = [tuple(map(int, p)) for p in np.argwhere(body & (d >= int(min_dist)))]
    # Ensure the two task terminals are represented whenever available.
    for k in (0, 1):
        p = m.source_terminal(k)
        if p is not None and p not in pts:
            pts.append(tuple(map(int, p)))
    pts.sort(key=lambda p: (int(d[p]), p[0], p[1]))
    return pts


def transfer_metrics(m, L, gfield, omega, tau_h, sites):
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

    sp = np.angle(hs)
    lp = np.angle(hl)
    out = dict(
        soma_phase_rms=circular_rms(sp),
        local_phase_rms=circular_rms(lp),
        soma_amp_median=float(np.median(np.abs(hs))),
        soma_amp_cv=float(np.std(np.abs(hs)) / (np.mean(np.abs(hs)) + 1e-30)),
        local_amp_median=float(np.median(np.abs(hl))),
    )

    # Direct A-vs-B terminal phase gap as a familiar two-input secondary metric.
    col = {p: j for j, p in enumerate(sites)}
    ta, tb = m.source_terminal(0), m.source_terminal(1)
    if ta is not None and tb is not None and tuple(ta) in col and tuple(tb) in col:
        out['terminal_phase_gap'] = circ_gap(sp[col[tuple(ta)]], sp[col[tuple(tb)]])
    else:
        out['terminal_phase_gap'] = float('nan')
    return out


def one_body(m, g0, ratio, tau_h, omegas, nshuffle):
    wh, wv = ae.bond_weights(m, m.body)
    L = weighted_laplacian_sparse(wh, wv)
    sites = injection_sites(m)
    d = m.graph_distance_from_soma()
    rng = np.random.default_rng(int(m.cfg.seed) + 731_731 + int(round(g0 * 1e6)) + int(round(ratio * 100)) + int(round(tau_h * 10)))
    prof = build_profiles(m, g0, ratio, rng, nshuffle)

    byfreq = []
    for omega in omegas:
        row = dict(omega=float(omega))
        for name in ('zero', 'uniform', 'smooth', 'reverse'):
            row[name] = transfer_metrics(m, L, prof[name], omega, tau_h, sites)
        sh = [transfer_metrics(m, L, z, omega, tau_h, sites) for z in prof['shuffled']]
        row['shuffle'] = dict(
            soma_phase_rms=float(np.mean([q['soma_phase_rms'] for q in sh])),
            local_phase_rms=float(np.mean([q['local_phase_rms'] for q in sh])),
            terminal_phase_gap=float(np.nanmean([q['terminal_phase_gap'] for q in sh])),
            soma_amp_median=float(np.mean([q['soma_amp_median'] for q in sh])),
            soma_amp_cv=float(np.mean([q['soma_amp_cv'] for q in sh])),
        )
        byfreq.append(row)

    return dict(
        seed=int(m.cfg.seed),
        cells=int(m.body.sum()),
        sites=len(sites),
        max_graph_distance=int(np.max(d[m.body.astype(bool)])),
        g0=float(g0), ratio=float(ratio), tau_h=float(tau_h),
        byfreq=byfreq,
    )


def summarize(rows, omegas):
    def vals(profile, key):
        z = []
        for r in rows:
            for q in r['byfreq']:
                z.append(float(q[profile][key]))
        return np.asarray(z, float)

    z_s = vals('zero', 'soma_phase_rms')
    u_s = vals('uniform', 'soma_phase_rms')
    s_s = vals('smooth', 'soma_phase_rms')
    h_s = vals('shuffle', 'soma_phase_rms')
    r_s = vals('reverse', 'soma_phase_rms')
    z_l = vals('zero', 'local_phase_rms')
    s_l = vals('smooth', 'local_phase_rms')
    h_l = vals('shuffle', 'local_phase_rms')

    out = dict(
        bodies=len(rows),
        observations=int(len(s_s)),
        zero_soma_phase_rms=float(np.mean(z_s)),
        uniform_soma_phase_rms=float(np.mean(u_s)),
        smooth_soma_phase_rms=float(np.mean(s_s)),
        shuffle_soma_phase_rms=float(np.mean(h_s)),
        reverse_soma_phase_rms=float(np.mean(r_s)),
        smooth_gain_vs_zero=float(np.mean(z_s - s_s)),
        smooth_gain_vs_uniform=float(np.mean(u_s - s_s)),
        smooth_gain_vs_shuffle=float(np.mean(h_s - s_s)),
        smooth_gain_vs_reverse=float(np.mean(r_s - s_s)),
        zero_local_phase_rms=float(np.mean(z_l)),
        smooth_local_phase_rms=float(np.mean(s_l)),
        shuffle_local_phase_rms=float(np.mean(h_l)),
        local_phase_retention=float(np.mean(s_l) / (np.mean(z_l) + 1e-30)),
        smooth_terminal_gap=float(np.nanmean(vals('smooth', 'terminal_phase_gap'))),
        shuffle_terminal_gap=float(np.nanmean(vals('shuffle', 'terminal_phase_gap'))),
        uniform_terminal_gap=float(np.nanmean(vals('uniform', 'terminal_phase_gap'))),
    )
    # Positive score requires morphology-indexed placement to beat both controls;
    # local phase retention prevents a trivial globally flattened field from winning.
    retention_penalty = max(0.0, 0.70 - out['local_phase_retention'])
    out['development_score'] = float(
        out['smooth_gain_vs_shuffle'] + out['smooth_gain_vs_uniform']
        + 0.5 * out['smooth_gain_vs_zero'] - 2.0 * retention_penalty
    )

    # Frequency-resolved means are important: HCN-like compensation should be a
    # band effect, not an unexplained improvement at every omega.
    freq = {}
    for omega in omegas:
        zz = [q for r in rows for q in r['byfreq'] if abs(float(q['omega']) - float(omega)) < 1e-12]
        freq[str(float(omega))] = dict(
            zero=float(np.mean([q['zero']['soma_phase_rms'] for q in zz])),
            uniform=float(np.mean([q['uniform']['soma_phase_rms'] for q in zz])),
            smooth=float(np.mean([q['smooth']['soma_phase_rms'] for q in zz])),
            shuffle=float(np.mean([q['shuffle']['soma_phase_rms'] for q in zz])),
            reverse=float(np.mean([q['reverse']['soma_phase_rms'] for q in zz])),
            smooth_local=float(np.mean([q['smooth']['local_phase_rms'] for q in zz])),
            smooth_minus_shuffle=float(np.mean([q['smooth']['soma_phase_rms'] - q['shuffle']['soma_phase_rms'] for q in zz])),
        )
    out['frequency'] = freq
    return out


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=552)
    ap.add_argument('--seeds', type=int, default=4)
    ap.add_argument('--g0s', default='0.005,0.02')
    ap.add_argument('--ratios', default='2,7')
    ap.add_argument('--taus', default='4,12,32')
    ap.add_argument('--omegas', default='0.04,0.08,0.12,0.16,0.24')
    ap.add_argument('--nshuffle', type=int, default=3)
    ap.add_argument('--out', default='runs/hcn_impedance/dev.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    # Circular spread invariant to a global phase rotation.
    p = np.array([-.2, 0.0, .2])
    assert abs(circular_rms(p) - circular_rms(p + 1.3)) < 1e-12
    # Discrete auxiliary transfer tends to one at DC.
    dt, tau = .12, 12.0
    lam = np.exp(1j * 1e-9)
    ah = (dt/tau) / (lam - 1 + dt/tau)
    assert abs(ah - 1) < 1e-6
    print('selftest ok')


def main():
    a = parse_args()
    if a.selftest:
        selftest(); return
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    g0s = parse_floats(a.g0s)
    ratios = parse_floats(a.ratios)
    taus = parse_floats(a.taus)
    omegas = parse_floats(a.omegas)

    bodies = []
    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature = True
        bodies.append(m)
        print(f'body {seed}: cells={m.body.sum()} dmax={m.graph_distance_from_soma().max()}', flush=True)
    if not bodies:
        raise SystemExit('No valid bodies')

    candidates = []
    for g0 in g0s:
        for ratio in ratios:
            for tau in taus:
                rows = []
                for m in bodies:
                    rows.append(one_body(m, g0, ratio, tau, omegas, a.nshuffle))
                s = summarize(rows, omegas)
                rec = dict(g0=float(g0), ratio=float(ratio), tau_h=float(tau), summary=s, rows=rows)
                candidates.append(rec)
                print(
                    f"g0={g0:.4g} ratio={ratio:g} tau={tau:g}: "
                    f"smooth={s['smooth_soma_phase_rms']:.4f} shuffle={s['shuffle_soma_phase_rms']:.4f} "
                    f"uniform={s['uniform_soma_phase_rms']:.4f} local-ret={s['local_phase_retention']:.3f} "
                    f"score={s['development_score']:+.5f}", flush=True)

    candidates.sort(key=lambda r: r['summary']['development_score'], reverse=True)
    payload = dict(
        experiment='hcn_impedance_gradient_development_v01',
        biological_motivation='Vaidya & Johnston 2013 doi:10.1038/nn.3562',
        model_warning='minimal linear delayed-restorative impedance proxy; not conductance-based HCN',
        seed_start=a.seed_start, seeds_requested=a.seeds, bodies=len(bodies),
        omegas=omegas, nshuffle=a.nshuffle,
        candidates=candidates,
    )
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nTOP CANDIDATES')
    for r in candidates[:8]:
        s = r['summary']
        print(
            f"g0={r['g0']:.4g} ratio={r['ratio']:g} tau={r['tau_h']:g} "
            f"score={s['development_score']:+.6f} gain-shuf={s['smooth_gain_vs_shuffle']:+.6f} "
            f"gain-unif={s['smooth_gain_vs_uniform']:+.6f} local-ret={s['local_phase_retention']:.3f}"
        )


if __name__ == '__main__':
    main()
