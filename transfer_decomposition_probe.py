"""Single-source transfer decomposition of temporal-order selectivity.

For each frozen body/gain condition, record the two complex single-source
responses h_A(x,t), h_B(x,t), then synthesize the paired task either coherently
(with the cross-source complex term) or incoherently (power envelopes only).
Compare both with the actual paired simulation.

See TRANSFER_DECOMPOSITION_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


CONDS = {
    'baseline': (1.0, 1.0),
    'A_half': (0.5, 1.0),
    'A_double': (2.0, 1.0),
    'B_half': (1.0, 0.5),
    'B_double': (1.0, 2.0),
}


def addsrc(a, b):
    if isinstance(a, (float, int, np.floating)):
        return b
    if isinstance(b, (float, int, np.floating)):
        return a
    return a + b


def scaled_source(m, which, q, gain):
    src = m.pulse_source(which, q, False)
    if isinstance(src, (float, int, np.floating)):
        return src
    return src * float(gain)


def trace_single(m, which, gain, steps):
    m.reset_fast(True)
    out = np.zeros((int(steps),) + m.body.shape, np.complex64)
    for t in range(int(steps)):
        m.advance(scaled_source(m, which, t, gain), False, True, 'none')
        out[t] = m.psi
    return out


def trace_pair(m, gain_a, gain_b, lag, target, steps):
    m.reset_fast(True)
    out = np.zeros((int(steps),) + m.body.shape, np.complex64)
    first, second = ((0, 1) if target else (1, 0))
    gains = {0: float(gain_a), 1: float(gain_b)}
    for t in range(int(steps)):
        a = scaled_source(m, first, t, gains[first])
        b = scaled_source(m, second, t - int(lag), gains[second])
        m.advance(addsrc(a, b), False, True, 'none')
        out[t] = m.psi
    return out


def tshift(x, lag):
    x = np.asarray(x)
    out = np.zeros_like(x)
    lag = int(lag)
    if lag <= 0:
        out[:] = x
    elif lag < len(x):
        out[lag:] = x[:-lag]
    return out


def contrast_from_power(pT, pD):
    t = np.max(np.asarray(pT, float), axis=0)
    d = np.max(np.asarray(pD, float), axis=0)
    return (t - d) / (t + d + 1e-30)


def safe_corr(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) < 3 or np.std(a) < 1e-14 or np.std(b) < 1e-14:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w, l)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2 * tail)), w, l


def condition_decomposition(m, gain_a, gain_b, lag, steps):
    hA = trace_single(m, 0, gain_a, steps)
    hB = trace_single(m, 1, gain_b, steps)
    aT = trace_pair(m, gain_a, gain_b, lag, True, steps)
    aD = trace_pair(m, gain_a, gain_b, lag, False, steps)

    hA_l = tshift(hA, lag)
    hB_l = tshift(hB, lag)

    sT = hA + hB_l
    sD = hB + hA_l

    p_actual_T = np.abs(aT) ** 2
    p_actual_D = np.abs(aD) ** 2
    p_coh_T = np.abs(sT) ** 2
    p_coh_D = np.abs(sD) ** 2
    p_inc_T = np.abs(hA) ** 2 + np.abs(hB_l) ** 2
    p_inc_D = np.abs(hB) ** 2 + np.abs(hA_l) ** 2

    C_actual = contrast_from_power(p_actual_T, p_actual_D)
    C_coh = contrast_from_power(p_coh_T, p_coh_D)
    C_inc = contrast_from_power(p_inc_T, p_inc_D)

    cross_T = 2.0 * np.real(hA * np.conj(hB_l))
    cross_D = 2.0 * np.real(hB * np.conj(hA_l))
    it = np.argmax(p_actual_T, axis=0)
    id_ = np.argmax(p_actual_D, axis=0)
    xt = np.take_along_axis(cross_T, it[None, ...], axis=0)[0]
    xd = np.take_along_axis(cross_D, id_[None, ...], axis=0)[0]

    return C_actual, C_coh, C_inc, xt, xd


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=48)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/transfer_decomposition/transfer_decomposition.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    x = np.arange(5)
    y = tshift(x, 2)
    assert np.array_equal(y, [0,0,0,1,2])
    a = np.array([[[1.]], [[2.]]])
    b = np.array([[[1.]], [[1.]]])
    c = contrast_from_power(a, b)
    assert abs(float(c[0,0]) - 1/3) < 1e-12
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

    rows = []
    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature = True
        body = m.body.astype(bool)
        cells = np.argwhere(body)
        idx = tuple(cells.T)
        soma = tuple(m.soma)
        si = int(np.flatnonzero((cells[:,0] == soma[0]) & (cells[:,1] == soma[1]))[0])

        all_actual = []
        all_coh = []
        all_inc = []
        condition_rows = {}
        soma_err_coh = []
        soma_err_inc = []

        for name, (ga, gb) in CONDS.items():
            Ca, Cc, Ci, xt, xd = condition_decomposition(m, ga, gb, a.lag, a.steps)
            va = np.asarray(Ca[idx], float)
            vc = np.asarray(Cc[idx], float)
            vi = np.asarray(Ci[idx], float)
            all_actual.append(va); all_coh.append(vc); all_inc.append(vi)
            ec = np.abs(va - vc); ei = np.abs(va - vi)
            soma_err_coh.append(float(ec[si])); soma_err_inc.append(float(ei[si]))
            condition_rows[name] = dict(
                gain_A=ga, gain_B=gb,
                signed_corr_coherent=safe_corr(va, vc),
                signed_corr_incoherent=safe_corr(va, vi),
                abs_corr_coherent=safe_corr(np.abs(va), np.abs(vc)),
                abs_corr_incoherent=safe_corr(np.abs(va), np.abs(vi)),
                mae_coherent=float(np.mean(ec)),
                mae_incoherent=float(np.mean(ei)),
                body_absC_actual=float(np.mean(np.abs(va))),
                body_absC_coherent=float(np.mean(np.abs(vc))),
                body_absC_incoherent=float(np.mean(np.abs(vi))),
                soma_C_actual=float(va[si]),
                soma_C_coherent=float(vc[si]),
                soma_C_incoherent=float(vi[si]),
                soma_abs_error_coherent=float(ec[si]),
                soma_abs_error_incoherent=float(ei[si]),
                soma_absC_actual=float(abs(va[si])),
                soma_absC_coherent=float(abs(vc[si])),
                soma_absC_incoherent=float(abs(vi[si])),
                cross_at_actual_peaks_mean_abs=float(np.mean(0.5*(np.abs(xt[idx]) + np.abs(xd[idx])))),
                cross_at_actual_peaks_soma_abs=float(0.5*(abs(float(xt[soma])) + abs(float(xd[soma])))),
            )

        va = np.concatenate(all_actual)
        vc = np.concatenate(all_coh)
        vi = np.concatenate(all_inc)
        mae_c = float(np.mean(np.abs(va - vc)))
        mae_i = float(np.mean(np.abs(va - vi)))
        body_row = dict(
            seed=seed,
            cells=int(body.sum()),
            P1_signed_corr_coherent=safe_corr(va, vc),
            signed_corr_incoherent=safe_corr(va, vi),
            abs_corr_coherent=safe_corr(np.abs(va), np.abs(vc)),
            abs_corr_incoherent=safe_corr(np.abs(va), np.abs(vi)),
            mae_coherent=mae_c,
            mae_incoherent=mae_i,
            P2_incoherent_minus_coherent_mae=float(mae_i - mae_c),
            P3_soma_incoherent_minus_coherent_error=float(np.mean(soma_err_inc) - np.mean(soma_err_coh)),
            soma_mae_coherent=float(np.mean(soma_err_coh)),
            soma_mae_incoherent=float(np.mean(soma_err_inc)),
            conditions=condition_rows,
        )
        rows.append(body_row)
        print(f"seed {seed:2d}: coh r {body_row['P1_signed_corr_coherent']:+.4f} "
              f"MAE coh/inc {mae_c:.4f}/{mae_i:.4f} "
              f"Dbody {body_row['P2_incoherent_minus_coherent_mae']:+.4f} "
              f"Dsoma {body_row['P3_soma_incoherent_minus_coherent_error']:+.4f}", flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    p1 = np.asarray([r['P1_signed_corr_coherent'] for r in rows], float)
    p2 = np.asarray([r['P2_incoherent_minus_coherent_mae'] for r in rows], float)
    p3 = np.asarray([r['P3_soma_incoherent_minus_coherent_error'] for r in rows], float)
    p2p, p2w, p2l = sign_test_two_sided(p2)
    p3p, p3w, p3l = sign_test_two_sided(p3)

    summary = dict(
        bodies=len(rows),
        P1=dict(
            mean_signed_corr_coherent=float(np.mean(p1)),
            median_signed_corr_coherent=float(np.median(p1)),
            passed=bool(np.mean(p1) > .95),
        ),
        P2=dict(
            mean_incoherent_minus_coherent_mae=float(np.mean(p2)),
            median=float(np.median(p2)),
            positive_bodies=p2w,
            negative_bodies=p2l,
            sign_p=p2p,
            passed=bool(np.mean(p2) > 0 and np.isfinite(p2p) and p2p < .05 and np.mean(p1) > .95),
        ),
        P3=dict(
            mean_soma_incoherent_minus_coherent_error=float(np.mean(p3)),
            median=float(np.median(p3)),
            positive_bodies=p3w,
            negative_bodies=p3l,
            sign_p=p3p,
            passed=bool(np.mean(p3) > 0 and np.isfinite(p3p) and p3p < .05 and np.mean(p1) > .95),
        ),
        descriptive=dict(
            mean_signed_corr_incoherent=float(np.mean([r['signed_corr_incoherent'] for r in rows])),
            mean_abs_corr_coherent=float(np.mean([r['abs_corr_coherent'] for r in rows])),
            mean_abs_corr_incoherent=float(np.mean([r['abs_corr_incoherent'] for r in rows])),
            mean_mae_coherent=float(np.mean([r['mae_coherent'] for r in rows])),
            mean_mae_incoherent=float(np.mean([r['mae_incoherent'] for r in rows])),
            mean_soma_mae_coherent=float(np.mean([r['soma_mae_coherent'] for r in rows])),
            mean_soma_mae_incoherent=float(np.mean([r['soma_mae_incoherent'] for r in rows])),
            body_absC_actual=float(np.mean([q['body_absC_actual'] for r in rows for q in r['conditions'].values()])),
            body_absC_coherent=float(np.mean([q['body_absC_coherent'] for r in rows for q in r['conditions'].values()])),
            body_absC_incoherent=float(np.mean([q['body_absC_incoherent'] for r in rows for q in r['conditions'].values()])),
            soma_absC_actual=float(np.mean([q['soma_absC_actual'] for r in rows for q in r['conditions'].values()])),
            soma_absC_coherent=float(np.mean([q['soma_absC_coherent'] for r in rows for q in r['conditions'].values()])),
            soma_absC_incoherent=float(np.mean([q['soma_absC_incoherent'] for r in rows for q in r['conditions'].values()])),
        ),
    )

    payload = dict(
        experiment='transfer_decomposition_v01',
        prereg='TRANSFER_DECOMPOSITION_PREREG_V01.md',
        seed_start=a.seed_start,
        seeds_requested=a.seeds,
        lag=a.lag,
        steps=a.steps,
        conditions=CONDS,
        summary=summary,
        rows=rows,
    )
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nTRANSFER DECOMPOSITION RECEIPT')
    print(f" P1 coherent signed r {summary['P1']['mean_signed_corr_coherent']:+.5f} PASS={summary['P1']['passed']}")
    print(f" P2 Dbody {summary['P2']['mean_incoherent_minus_coherent_mae']:+.5f} "
          f"pos/neg {p2w}/{p2l} p={p2p:.5g} PASS={summary['P2']['passed']}")
    print(f" P3 Dsoma {summary['P3']['mean_soma_incoherent_minus_coherent_error']:+.5f} "
          f"pos/neg {p3w}/{p3l} p={p3p:.5g} PASS={summary['P3']['passed']}")
    print(' wrote', out)


if __name__ == '__main__':
    main()
