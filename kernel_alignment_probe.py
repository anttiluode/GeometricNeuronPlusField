"""Why is the soma point-power kernel aligned with the temporal-order task?

Frozen-body discovery probe.  See KERNEL_ALIGNMENT_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from graph_mode_probe import graph_laplacian_modes
from transfer_decomposition_probe import trace_single, tshift, safe_corr


def percentile(v, x):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) == 0 or not np.isfinite(x):
        return float('nan')
    return float(np.mean(v <= float(x)))


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w, l)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2.0 * tail)), w, l


def body_probe(m, lag=20, steps=210):
    hA = trace_single(m, 0, 1.0, steps)
    hB = trace_single(m, 1, 1.0, steps)

    coords, evals, evecs = graph_laplacian_modes(m.body)
    n = len(coords)
    ys = np.asarray([p[0] for p in coords], int)
    xs = np.asarray([p[1] for p in coords], int)
    cidx = {p: i for i, p in enumerate(coords)}
    soma = tuple(map(int, m.soma))
    si = cidx[soma]

    ZA = np.asarray(hA[:, ys, xs], np.complex128)
    ZB = np.asarray(hB[:, ys, xs], np.complex128)
    ZA_l = tshift(ZA, lag)
    ZB_l = tshift(ZB, lag)

    # Full graph-modal source histories.  evecs is real orthogonal.
    QA = ZA @ evecs
    QB = ZB @ evecs
    QA_l = tshift(QA, lag)
    QB_l = tshift(QB, lag)

    # Historical coherent pair task reconstructed from the independently measured
    # single-source fields.
    psiT = ZA + ZB_l
    psiD = ZB + ZA_l
    pT = np.abs(psiT) ** 2
    pD = np.abs(psiD) ** 2
    peakT = np.max(pT, axis=0)
    peakD = np.max(pD, axis=0)
    C = (peakT - peakD) / (peakT + peakD + 1e-30)
    absC = np.abs(C)

    # Existing amplitude-balance opportunity map.
    PA = np.max(np.abs(ZA) ** 2, axis=0)
    PB = np.max(np.abs(ZB) ** 2, axis=0)
    balance = np.minimum(PA, PB) / (np.maximum(PA, PB) + 1e-30)

    # Direct coordinate-space order-sensitive cross-source difference at each time.
    crossT = 2.0 * np.real(ZA * np.conj(ZB_l))
    crossD = 2.0 * np.real(ZB * np.conj(ZA_l))
    crossdiff = crossT - crossD

    # D(t) is a global Hermitian task-interaction operator in modal coordinates.
    # Its Frobenius norm removes global instantaneous interaction magnitude so the
    # remaining score asks how well K_x=v_x v_x^T is oriented toward D(t).
    dnorm = np.zeros(int(steps), float)
    for t in range(int(steps)):
        a = QA[t]; b = QB_l[t]; c = QB[t]; d = QA_l[t]
        D = (np.outer(a, np.conj(b)) + np.outer(b, np.conj(a))
             - np.outer(c, np.conj(d)) - np.outer(d, np.conj(c)))
        dnorm[t] = float(np.linalg.norm(D, ord='fro'))

    mx = float(np.max(dnorm)) if len(dnorm) else 0.0
    valid = dnorm > max(1e-30, mx * 1e-10)
    if not np.any(valid):
        alignment = np.zeros(n, float)
    else:
        z = crossdiff[valid] / dnorm[valid, None]
        alignment = np.sqrt(np.mean(z * z, axis=0))

    cross_rms = np.sqrt(np.mean(crossdiff * crossdiff, axis=0))
    combined = balance * alignment

    # Modal reconstruction control for the coordinate-space cross difference.
    recA = QA @ evecs.T
    recB = QB @ evecs.T
    recA_l = tshift(recA, lag)
    recB_l = tshift(recB, lag)
    rec_crossdiff = (2.0 * np.real(recA * np.conj(recB_l))
                     - 2.0 * np.real(recB * np.conj(recA_l)))
    modal_cross_mae = float(np.mean(np.abs(rec_crossdiff - crossdiff)))
    cross_scale = float(np.mean(np.abs(crossdiff)) + 1e-30)
    modal_cross_rel = float(modal_cross_mae / cross_scale)

    rA = safe_corr(alignment, absC)
    rB = safe_corr(balance, absC)
    rQ = safe_corr(combined, absC)
    rU = safe_corr(cross_rms, absC)
    rAB = safe_corr(alignment, balance)

    return dict(
        cells=int(n),
        soma=[int(soma[0]), int(soma[1])],
        valid_operator_times=int(np.sum(valid)),
        modal_cross_reconstruction_mae=modal_cross_mae,
        modal_cross_reconstruction_relative=modal_cross_rel,
        r_alignment_absC=float(rA),
        r_balance_absC=float(rB),
        r_combined_absC=float(rQ),
        r_crossrms_absC=float(rU),
        r_alignment_balance=float(rAB),
        delta_combined_minus_balance=float(rQ - rB),
        soma_absC=float(absC[si]),
        soma_C=float(C[si]),
        soma_alignment=float(alignment[si]),
        soma_balance=float(balance[si]),
        soma_combined=float(combined[si]),
        soma_crossrms=float(cross_rms[si]),
        soma_absC_percentile=percentile(absC, absC[si]),
        soma_alignment_percentile=percentile(alignment, alignment[si]),
        soma_balance_percentile=percentile(balance, balance[si]),
        soma_combined_percentile=percentile(combined, combined[si]),
        descriptive=dict(
            alignment_mean=float(np.mean(alignment)),
            alignment_median=float(np.median(alignment)),
            balance_mean=float(np.mean(balance)),
            absC_mean=float(np.mean(absC)),
            combined_mean=float(np.mean(combined)),
            max_dnorm=mx,
        ),
    )


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=108)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/kernel_alignment_discovery/kernel_alignment.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    assert abs(percentile([1,2,3,4], 3) - .75) < 1e-12
    p, w, l = sign_test_two_sided([1,1,-1])
    assert w == 2 and l == 1 and 0 <= p <= 1
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
        r = body_probe(m, a.lag, a.steps)
        r['seed'] = int(seed)
        rows.append(r)
        print(
            f"seed {seed}: rA={r['r_alignment_absC']:+.3f} "
            f"rB={r['r_balance_absC']:+.3f} rQ={r['r_combined_absC']:+.3f} "
            f"dQ-B={r['delta_combined_minus_balance']:+.3f} "
            f"soma pct A/B/Q/C={r['soma_alignment_percentile']:.2f}/"
            f"{r['soma_balance_percentile']:.2f}/{r['soma_combined_percentile']:.2f}/"
            f"{r['soma_absC_percentile']:.2f}", flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    def vals(k):
        return np.asarray([r[k] for r in rows], float)

    rA = vals('r_alignment_absC')
    rB = vals('r_balance_absC')
    rQ = vals('r_combined_absC')
    d = vals('delta_combined_minus_balance')
    somaA = vals('soma_alignment_percentile')
    dp, dw, dl = sign_test_two_sided(d)

    summary = dict(
        bodies=len(rows),
        mean_r_alignment_absC=float(np.nanmean(rA)),
        median_r_alignment_absC=float(np.nanmedian(rA)),
        alignment_positive_bodies=int(np.sum(rA > 0)),
        mean_r_balance_absC=float(np.nanmean(rB)),
        mean_r_combined_absC=float(np.nanmean(rQ)),
        mean_delta_combined_minus_balance=float(np.nanmean(d)),
        median_delta_combined_minus_balance=float(np.nanmedian(d)),
        combined_improved_bodies=int(np.sum(d > 0)),
        combined_worse_bodies=int(np.sum(d < 0)),
        combined_improvement_sign_p=float(dp),
        median_soma_alignment_percentile=float(np.nanmedian(somaA)),
        mean_soma_alignment_percentile=float(np.nanmean(somaA)),
        soma_alignment_above_median_bodies=int(np.sum(somaA > .5)),
        mean_soma_balance_percentile=float(np.nanmean(vals('soma_balance_percentile'))),
        mean_soma_combined_percentile=float(np.nanmean(vals('soma_combined_percentile'))),
        mean_soma_absC_percentile=float(np.nanmean(vals('soma_absC_percentile'))),
        mean_r_alignment_balance=float(np.nanmean(vals('r_alignment_balance'))),
        mean_r_crossrms_absC=float(np.nanmean(vals('r_crossrms_absC'))),
        modal_cross_reconstruction_relative_mean=float(np.nanmean(vals('modal_cross_reconstruction_relative'))),
    )
    summary['D1_pass'] = bool(summary['mean_r_alignment_absC'] > .50 and summary['alignment_positive_bodies'] >= 9)
    summary['D2_pass'] = bool(summary['median_soma_alignment_percentile'] > .75 and summary['soma_alignment_above_median_bodies'] >= 9)
    summary['D3_pass'] = bool(summary['mean_delta_combined_minus_balance'] > 0 and summary['combined_improved_bodies'] >= 8)

    payload = dict(
        experiment='kernel_alignment_discovery_v01',
        prereg='KERNEL_ALIGNMENT_DISCOVERY_PREREG_V01.md',
        seed_start=a.seed_start,
        seeds_requested=a.seeds,
        lag=a.lag,
        steps=a.steps,
        summary=summary,
        rows=rows,
    )
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nKERNEL ALIGNMENT DISCOVERY RECEIPT')
    for k, v in summary.items():
        print(f' {k}: {v}')


if __name__ == '__main__':
    main()
