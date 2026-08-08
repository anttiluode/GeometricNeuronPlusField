"""Mechanism audit for the single-run +/- lock-in gradient readout.

`polarization_identity_probe.py` found a sharp rate-dependent effect: fast balanced
0/pi toggling recovers the compressed K-bin adjoint gradient essentially exactly,
while slower toggling leaks badly.

The algebra says why.  With s(t) in {+1,-1},

    I(t) = |u(t) + s(t) v(t)|^2

and a signed power accumulator gives

    1/2 sum_t s(t) I(t)
      = sum_t Re(conj(u(t)) v(t))
        + 1/2 sum_t s(t) [|u(t)|^2 + |v(t)|^2].

The second term is the only error.  If u and v contain only selected spectral bins K,
then their self-energy spectra live on pairwise DIFFERENCE frequencies k_i-k_j.
The square-wave reference contributes its own harmonic set.  Leakage therefore requires
spectral overlap between those two sets (with amplitudes/phases determining its size).

This probe tests three levels of explanation on the existing development bodies:

1. exact frequency-domain Parseval reconstruction of the measured leakage;
2. actual weighted spectral-overlap magnitude versus measured leakage;
3. a structure-only collision score using only selected-bin differences and the
   reference spectrum, ignoring internal amplitudes.

If (1) is machine exact and (2)/(3) order the modulation rates correctly, the apparent
"fast rhythm" effect is an ordinary heterodyne/lock-in separation rule, not a new
learning principle.  That rule can then be preregistered on fresh bodies.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # exact mature-boundary patch
import polarization_identity_probe as pi
from reciprocal_adjoint_probe import flat_pair
from transfer_decomposition_probe import safe_corr


def reference(T: int, half_period: int):
    h = int(half_period)
    if T % (2 * h) != 0:
        raise ValueError(f'half_period {h} does not balance T={T}')
    s = np.where((np.arange(T) // h) % 2 == 0, 1.0, -1.0)
    assert abs(float(np.sum(s))) < 1e-12
    return s


def selected_time(U, V, kk):
    us = np.zeros_like(U)
    vs = np.zeros_like(V)
    us[kk] = U[kk]
    vs[kk] = V[kk]
    return np.fft.ifft(us, axis=0), np.fft.ifft(vs, axis=0)


def one_component(U, V, kk, td_fac, half_period):
    u, v = selected_time(U, V, kk)
    T = len(u)
    s = reference(T, half_period)
    ss = s.reshape((T,) + (1,) * (u.ndim - 1))
    e = np.abs(u) ** 2 + np.abs(v) ** 2

    # Direct time-domain leakage in gradient units.
    leak_td = 0.5 * td_fac * np.sum(ss * e, axis=0)

    # Frequency-domain identity. numpy convention gives
    # sum_t s[t] e[t] = (1/T) sum_k S[-k] E[k].
    S = np.fft.fft(s)
    E = np.fft.fft(e, axis=0)
    k = np.arange(T)
    leak_fd = 0.5 * td_fac * np.real(
        np.sum(S[(-k) % T].reshape((T,) + (1,) * (E.ndim - 1)) * E, axis=0) / T
    )

    # Actual weighted spectral overlap before signed cancellation.  This is an upper
    # bound/availability measure, not the leakage itself.
    shape = (T,) + (1,) * (E.ndim - 1)
    weighted_overlap = 0.5 * td_fac * np.sum(
        np.abs(S[(-k) % T]).reshape(shape) * np.abs(E), axis=0
    ) / T

    return leak_td, leak_fd, weighted_overlap, S, E


def difference_mask(T, kk):
    m = np.zeros(T, bool)
    kk = np.asarray(kk, int)
    for a in kk:
        for b in kk:
            m[(int(a) - int(b)) % T] = True
    return m


def structural_collision(T, kk, S, tol=1e-9):
    """Reference harmonic mass landing on the selected-bin difference set.

    Uses only selected frequency indices and the known modulation waveform.  It does
    not inspect local field amplitudes, so it is a deliberately weak predictor.
    """
    dm = difference_mask(T, kk)
    a = np.abs(S)
    active = a > float(tol) * (float(np.max(a)) + 1e-30)
    denom = float(np.sum(a[active]) + 1e-30)
    return dict(
        difference_bins=int(np.sum(dm)),
        reference_harmonics=int(np.sum(active)),
        colliding_harmonics=int(np.sum(active & dm)),
        reference_mass_on_difference_set=float(np.sum(a[dm]) / denom),
        exact_support_collision=bool(np.any(active & dm)),
    )


def combine_component(Td, Dd, component, kk, half_period):
    if component == 'h':
        keys = ('UH', 'VH')
    else:
        keys = ('UV', 'VV')
    lt, lf, wt, S, E = one_component(Td[keys[0]], Td[keys[1]], kk, Td['td_fac'], half_period)
    ld, lfd, wd, _, ED = one_component(Dd[keys[0]], Dd[keys[1]], kk, Dd['td_fac'], half_period)
    return lt + ld, lf + lfd, wt + wd, S, E + ED


def norm_ratio(xh, xv, refh, refv):
    x = flat_pair(xh, xv)
    r = flat_pair(refh, refv)
    return float(np.linalg.norm(x) / (np.linalg.norm(r) + 1e-30))


def one(m, lag, steps):
    wh, wv = ae.bond_weights(m, m.body)
    seqT = ae.source_sequence(m, True, lag, steps)
    seqD = ae.source_sequence(m, False, lag, steps)
    ET = ae.linear_forward(m, wh, wv, seqT, store=False)
    ED = ae.linear_forward(m, wh, wv, seqD, store=False)
    denom = ET + ED + 1e-30
    aT = 2.0 * ED / (denom * denom)
    aD = -2.0 * ET / (denom * denom)
    Td = pi.order_data(m, wh, wv, seqT, aT)
    Dd = pi.order_data(m, wh, wv, seqD, aD)
    order = np.argsort(Td['port_score'] + Dd['port_score'])[::-1]

    out = {}
    for K in pi.KS:
        kk = np.asarray(order[:K], int)
        gh, gv = pi.combine(pi.direct_map, Td, Dd, kk)
        kr = {}
        for h in pi.LOCKIN_HALF_PERIODS:
            lth, lfh, wth, S, _ = combine_component(Td, Dd, 'h', kk, h)
            ltv, lfv, wtv, _, _ = combine_component(Td, Dd, 'v', kk, h)
            fd_err = norm_ratio(lfh - lth, lfv - ltv, gh, gv)
            leak_ratio = norm_ratio(lth, ltv, gh, gv)
            overlap_ratio = norm_ratio(wth, wtv, gh, gv)
            sc = structural_collision(steps, kk, S)
            kr[str(h)] = dict(
                leakage_to_gradient_l2=leak_ratio,
                spectral_parseval_error=fd_err,
                weighted_overlap_to_gradient_l2=overlap_ratio,
                **sc,
            )
        out[str(K)] = dict(bins=[int(x) for x in kk], half_period=kr)
    return dict(seed=int(m.cfg.seed), C=float((ET - ED) / denom), results=out)


def summarize(rows):
    s = dict(bodies=len(rows), K={})
    for K in map(str, pi.KS):
        s['K'][K] = {}
        for h in map(str, pi.LOCKIN_HALF_PERIODS):
            q = [r['results'][K]['half_period'][h] for r in rows]
            s['K'][K][h] = dict(
                mean_leakage=float(np.mean([x['leakage_to_gradient_l2'] for x in q])),
                mean_weighted_overlap=float(np.mean([x['weighted_overlap_to_gradient_l2'] for x in q])),
                max_parseval_error=float(np.max([x['spectral_parseval_error'] for x in q])),
                mean_collision_mass=float(np.mean([x['reference_mass_on_difference_set'] for x in q])),
                bodies_with_support_collision=int(np.sum([x['exact_support_collision'] for x in q])),
                mean_colliding_harmonics=float(np.mean([x['colliding_harmonics'] for x in q])),
            )

        # Across all body/rate points, do the two predictors order leakage?
        leak = []
        weighted = []
        collision = []
        for r in rows:
            for h in map(str, pi.LOCKIN_HALF_PERIODS):
                z = r['results'][K]['half_period'][h]
                leak.append(z['leakage_to_gradient_l2'])
                weighted.append(z['weighted_overlap_to_gradient_l2'])
                collision.append(z['reference_mass_on_difference_set'])
        s['K'][K]['pooled_predictors'] = dict(
            weighted_overlap_corr=float(safe_corr(np.asarray(weighted), np.asarray(leak))),
            structural_collision_corr=float(safe_corr(np.asarray(collision), np.asarray(leak))),
        )
    return s


def selftest():
    rng = np.random.default_rng(8)
    T = 30
    U = rng.normal(size=(T, 2, 3)) + 1j * rng.normal(size=(T, 2, 3))
    V = rng.normal(size=(T, 2, 3)) + 1j * rng.normal(size=(T, 2, 3))
    kk = np.asarray([1, 4, 8, 12, 17, 21, 25, 28])
    td, fd, _, S, _ = one_component(U, V, kk, 0.31, 3)
    assert np.max(np.abs(td - fd)) < 1e-11
    sc = structural_collision(T, kk, S)
    assert sc['reference_harmonics'] > 0
    print('selftest ok')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=472)
    ap.add_argument('--seeds', type=int, default=8)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/lockin_overlap/dev_472_479.json')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return

    fa = Path(a.functional_arbors).resolve()
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    rows = []
    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature = True
        r = one(m, a.lag, a.steps)
        rows.append(r)
        print('seed', seed,
              'K8', [(h, round(r['results']['8']['half_period'][str(h)]['leakage_to_gradient_l2'], 4),
                       r['results']['8']['half_period'][str(h)]['colliding_harmonics'])
                      for h in pi.LOCKIN_HALF_PERIODS], flush=True)
    if not rows:
        raise SystemExit('No valid bodies')

    summary = summarize(rows)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(experiment='lockin_overlap_dev_v01', summary=summary, rows=rows), indent=2))
    print('\nLOCKIN OVERLAP DEV')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
