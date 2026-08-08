"""Physical-readout audit for the compressed transient adjoint gradient.

The compressed K-bin learner currently evaluates

    Re[conj(U_k) V_k]

numerically at every bond.  This probe asks how much of that multiplication can be
replaced by power/intensity measurements.

For each selected frequency bin the polarization identity gives

    (|U+V|^2 - |U-V|^2) / 4 = Re[conj(U)V].

Two stronger questions are tested as well:

1. Can all selected bins be measured with one GLOBAL 0/pi toggle rather than 2K
   independent bin toggles?  With a full-window power integral the answer should be
   yes algebraically: Parseval removes cross-frequency terms.

2. Can a single continuously running +/- phase modulation act as a lock-in reference,
   so that a signed power accumulator recovers the cross term without two separate
   measurements?  This is not guaranteed because the self-energy envelope leaks into
   the modulation channel.  Several exactly balanced square-wave periods are tested.

The second question is the deliberately speculative bridge to rhythmic gating.  It is
an engineering test only; it is NOT a model of theta/gamma, chandelier cells, GABA, or
biological backpropagation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # patch exact mature-boundary Laplacian
from reciprocal_adjoint_probe import flat_pair, normalized_l2, retro_source_sequence
from spectral_correlation_compression_probe import edge_series, port_spectrum_score
from transfer_decomposition_probe import safe_corr

KS = (8, 16)
PHASE_ERRORS = (0.025, 0.05, 0.10)
# T=210 in the registered task.  These half-periods all divide T/2=105 exactly,
# so the +/- reference has exactly zero DC component.
LOCKIN_HALF_PERIODS = (1, 3, 5, 7, 15, 21, 35)


def component_spectra(m, forward_states, retro_states):
    """Return local spectral factors U,V such that fac*Re(conj(U)V) is gradient."""
    Fh, Fv = edge_series(forward_states, False)
    Rh, Rv = edge_series(retro_states, True)
    T = len(Fh)
    k = np.arange(T)
    neg = (-k) % T

    FH = np.fft.fft(Fh, axis=0)
    FV = np.fft.fft(Fv, axis=0)
    RH = np.fft.fft(Rh, axis=0)
    RV = np.fft.fft(Rv, axis=0)

    # FFT(reverse(f))[k] = exp(+i 2pi k/T) FFT(f)[-k].
    # The apparently frequency-specific phase shift is therefore one deterministic
    # linear phase ramp produced by the reverse-time alignment, not K free phases.
    phase = np.exp(2j * np.pi * k / T)[:, None, None]
    VH = phase * FH[neg]
    VV = phase * FV[neg]
    UH = RH
    UV = RV

    fac = 2.0 * float(m.cfg.dt) * float(m.cfg.stiffness) / T
    td_fac = 2.0 * float(m.cfg.dt) * float(m.cfg.stiffness)
    return UH, VH, UV, VV, fac, td_fac


def order_data(m, wh, wv, seq, weight):
    p, v, E = ae.linear_forward(m, wh, wv, seq, store=True)
    g = weight * np.asarray(p[1:, m.soma[0], m.soma[1]], np.complex128)
    rseq = retro_source_sequence(m, g, reverse=True)
    rp, rv, _ = ae.linear_forward(m, wh, wv, rseq, store=True)
    UH, VH, UV, VV, fac, td_fac = component_spectra(m, p[:-1], rp[1:])
    eh, ev = ae.adjoint_grad(m, wh, wv, p, v, weight)
    return dict(
        UH=UH, VH=VH, UV=UV, VV=VV, fac=fac, td_fac=td_fac,
        exact=(eh, ev), port_score=port_spectrum_score(seq, rseq), E=E,
    )


def direct_map(U, V, kk, fac):
    return fac * np.sum(np.real(np.conj(U[kk]) * V[kk]), axis=0)


def binwise_polarization_map(U, V, kk, fac, phase_error=0.0):
    """Two intensity states per bin; phase_error perturbs the nominal pi state."""
    u = U[kk]
    v = V[kk]
    i0 = np.abs(u + v) ** 2
    rot = np.exp(1j * (np.pi + float(phase_error)))
    ipi = np.abs(u + rot * v) ** 2
    return 0.25 * fac * np.sum(i0 - ipi, axis=0)


def broadband_polarization_map(U, V, kk, td_fac, phase_error=0.0):
    """One global 0/pi state pair for all selected bins, integrated over time.

    The selected spectra are transformed back to a complex time waveform.  Parseval
    makes the full-window intensity difference equal the sum of the selected bin
    cross terms.  Therefore this is a two-state measurement independent of K, provided
    the K-bin local phasors can actually be retained/replayed and the detector can
    integrate the full aligned window.
    """
    us = np.zeros_like(U)
    vs = np.zeros_like(V)
    us[kk] = U[kk]
    vs[kk] = V[kk]
    u = np.fft.ifft(us, axis=0)
    v = np.fft.ifft(vs, axis=0)
    i0 = np.abs(u + v) ** 2
    rot = np.exp(1j * (np.pi + float(phase_error)))
    ipi = np.abs(u + rot * v) ** 2
    return 0.25 * td_fac * np.sum(i0 - ipi, axis=0)


def lockin_map(U, V, kk, td_fac, half_period):
    """Single-run +/- phase modulation with a signed power accumulator.

    If s(t)=+/-1 and has zero mean,

      s |u+s v|^2 = s(|u|^2+|v|^2) + 2 Re(conj(u)v).

    The desired cross term is exact; the only error is leakage of the self-energy
    envelope into the square-wave reference.  Faster toggling should help only if that
    envelope is slow relative to the modulation.
    """
    us = np.zeros_like(U)
    vs = np.zeros_like(V)
    us[kk] = U[kk]
    vs[kk] = V[kk]
    u = np.fft.ifft(us, axis=0)
    v = np.fft.ifft(vs, axis=0)
    T = len(u)
    h = int(half_period)
    if T % (2 * h) != 0:
        raise ValueError(f'half_period {h} does not exactly balance T={T}')
    s = np.where((np.arange(T) // h) % 2 == 0, 1.0, -1.0)
    if abs(float(np.sum(s))) > 1e-12:
        raise AssertionError('lock-in reference must be exactly balanced')
    ss = s.reshape((T,) + (1,) * (u.ndim - 1))
    measured = np.abs(u + ss * v) ** 2
    # target = td_fac * sum Re(conj(u)v)
    # signed-power estimate = td_fac/2 * sum s*I
    est = 0.5 * td_fac * np.sum(ss * measured, axis=0)
    leakage = 0.5 * td_fac * np.sum(ss * (np.abs(u) ** 2 + np.abs(v) ** 2), axis=0)
    return est, leakage


def metrics(ref_h, ref_v, test_h, test_v):
    ref = flat_pair(ref_h, ref_v)
    test = flat_pair(test_h, test_v)
    mx = float(np.max(np.abs(ref)) + 1e-30)
    mask = np.abs(ref) > 0.01 * mx
    return dict(
        corr=float(safe_corr(ref, test)),
        relative_l2=normalized_l2(ref, test),
        strong_sign_agreement=(
            float(np.mean(np.sign(ref[mask]) == np.sign(test[mask])))
            if np.any(mask) else float('nan')
        ),
    )


def combine(method, T, D, kk, **kwargs):
    gh = method(T['UH'], T['VH'], kk, T['fac'] if method in (direct_map, binwise_polarization_map) else T['td_fac'], **kwargs)
    gv = method(T['UV'], T['VV'], kk, T['fac'] if method in (direct_map, binwise_polarization_map) else T['td_fac'], **kwargs)
    gh += method(D['UH'], D['VH'], kk, D['fac'] if method in (direct_map, binwise_polarization_map) else D['td_fac'], **kwargs)
    gv += method(D['UV'], D['VV'], kk, D['fac'] if method in (direct_map, binwise_polarization_map) else D['td_fac'], **kwargs)
    return gh, gv


def combine_lockin(T, D, kk, half_period):
    th, tlh = lockin_map(T['UH'], T['VH'], kk, T['td_fac'], half_period)
    tv, tlv = lockin_map(T['UV'], T['VV'], kk, T['td_fac'], half_period)
    dh, dlh = lockin_map(D['UH'], D['VH'], kk, D['td_fac'], half_period)
    dv, dlv = lockin_map(D['UV'], D['VV'], kk, D['td_fac'], half_period)
    return th + dh, tv + dv, tlh + dlh, tlv + dlv


def one(m, lag, steps):
    wh, wv = ae.bond_weights(m, m.body)
    seqT = ae.source_sequence(m, True, lag, steps)
    seqD = ae.source_sequence(m, False, lag, steps)
    ET = ae.linear_forward(m, wh, wv, seqT, store=False)
    ED = ae.linear_forward(m, wh, wv, seqD, store=False)
    S = ET + ED + 1e-30
    aT = 2.0 * ED / (S * S)
    aD = -2.0 * ET / (S * S)
    T = order_data(m, wh, wv, seqT, aT)
    D = order_data(m, wh, wv, seqD, aD)
    exact_h = T['exact'][0] + D['exact'][0]
    exact_v = T['exact'][1] + D['exact'][1]
    order = np.argsort(T['port_score'] + D['port_score'])[::-1]

    kres = {}
    for K in KS:
        kk = np.asarray(order[:K], int)
        dh, dv = combine(direct_map, T, D, kk)
        bh, bv = combine(binwise_polarization_map, T, D, kk)
        ph, pv = combine(broadband_polarization_map, T, D, kk)
        q = dict(
            bins=[int(x) for x in kk],
            direct_vs_exact=metrics(exact_h, exact_v, dh, dv),
            binwise_vs_direct=metrics(dh, dv, bh, bv),
            broadband_two_state_vs_direct=metrics(dh, dv, ph, pv),
            phase_error={},
            lockin={},
        )
        for eps in PHASE_ERRORS:
            eh, ev = combine(broadband_polarization_map, T, D, kk, phase_error=eps)
            q['phase_error'][str(eps)] = metrics(dh, dv, eh, ev)
        for h in LOCKIN_HALF_PERIODS:
            lh, lv, leak_h, leak_v = combine_lockin(T, D, kk, h)
            mm = metrics(dh, dv, lh, lv)
            leak = normalized_l2(np.zeros_like(flat_pair(dh, dv)), flat_pair(leak_h, leak_v))
            # The normalized_l2 call above is not useful with a zero reference; report
            # leakage relative to desired compressed-gradient norm explicitly.
            desired = flat_pair(dh, dv)
            leakflat = flat_pair(leak_h, leak_v)
            mm['leakage_to_gradient_l2'] = float(np.linalg.norm(leakflat) / (np.linalg.norm(desired) + 1e-30))
            q['lockin'][str(h)] = mm
        kres[str(K)] = q

    return dict(seed=int(m.cfg.seed), C=float((ET - ED) / S), results=kres)


def selftest():
    rng = np.random.default_rng(123)
    T = 30
    U = rng.normal(size=(T, 3, 4)) + 1j * rng.normal(size=(T, 3, 4))
    V = rng.normal(size=(T, 3, 4)) + 1j * rng.normal(size=(T, 3, 4))
    kk = np.asarray([1, 3, 7, 11, 19, 22, 25, 29])
    fac = 0.37 / T
    td_fac = 0.37
    d = direct_map(U, V, kk, fac)
    b = binwise_polarization_map(U, V, kk, fac)
    p = broadband_polarization_map(U, V, kk, td_fac)
    assert np.max(np.abs(d - b)) < 1e-12
    assert np.max(np.abs(d - p)) < 1e-12
    # Exact identity of lock-in estimate and desired+leakage.
    l, leak = lockin_map(U, V, kk, td_fac, 1)
    assert np.max(np.abs(l - d - leak)) < 1e-12
    print('selftest ok')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=472)
    ap.add_argument('--seeds', type=int, default=4)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/polarization_identity/dev.json')
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
        print(
            'seed', seed,
            'K8 direct', round(r['results']['8']['direct_vs_exact']['corr'], 4),
            '2state rel', f"{r['results']['8']['broadband_two_state_vs_direct']['relative_l2']:.2e}",
            'lock1', round(r['results']['8']['lockin']['1']['corr'], 4),
            'lock7', round(r['results']['8']['lockin']['7']['corr'], 4),
            flush=True,
        )
    if not rows:
        raise SystemExit('No valid bodies')

    summary = dict(bodies=len(rows), steps=a.steps, K={})
    for K in map(str, KS):
        q = [r['results'][K] for r in rows]
        summary['K'][K] = dict(
            direct_vs_exact=dict(
                mean_corr=float(np.mean([x['direct_vs_exact']['corr'] for x in q])),
                mean_relative_l2=float(np.mean([x['direct_vs_exact']['relative_l2'] for x in q])),
            ),
            binwise_identity_max_relative_l2=float(np.max([x['binwise_vs_direct']['relative_l2'] for x in q])),
            broadband_two_state_max_relative_l2=float(np.max([x['broadband_two_state_vs_direct']['relative_l2'] for x in q])),
            phase_error={},
            lockin={},
        )
        for eps in PHASE_ERRORS:
            z = [x['phase_error'][str(eps)] for x in q]
            summary['K'][K]['phase_error'][str(eps)] = dict(
                mean_corr=float(np.mean([x['corr'] for x in z])),
                mean_relative_l2=float(np.mean([x['relative_l2'] for x in z])),
                mean_strong_sign_agreement=float(np.mean([x['strong_sign_agreement'] for x in z])),
            )
        for h in LOCKIN_HALF_PERIODS:
            z = [x['lockin'][str(h)] for x in q]
            summary['K'][K]['lockin'][str(h)] = dict(
                mean_corr=float(np.mean([x['corr'] for x in z])),
                mean_relative_l2=float(np.mean([x['relative_l2'] for x in z])),
                mean_strong_sign_agreement=float(np.mean([x['strong_sign_agreement'] for x in z])),
                mean_leakage_to_gradient_l2=float(np.mean([x['leakage_to_gradient_l2'] for x in z])),
            )

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(experiment='polarization_identity_dev_v01', summary=summary, rows=rows), indent=2))
    print('\nPOLARIZATION IDENTITY DEV')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
