"""Gain-only causal test of the amplitude-balance explanation.

Fresh frozen bodies, fixed source locations, fixed lag and field physics.
Only one source pulse amplitude is multiplied by 0.5 or 2.0.
See SOURCE_GAIN_BALANCE_PREREG_V01.md.
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


def peak_single(m, which, gain, steps):
    m.reset_fast(True)
    acc = np.zeros(m.body.shape, float)
    for t in range(int(steps)):
        m.advance(scaled_source(m, which, t, gain), False, True, 'none')
        acc = np.maximum(acc, np.abs(m.psi) ** 2)
    return acc


def peak_order(m, gain_a, gain_b, lag, target, steps):
    m.reset_fast(True)
    acc = np.zeros(m.body.shape, float)
    first, second = ((0, 1) if target else (1, 0))
    gains = {0: float(gain_a), 1: float(gain_b)}
    for t in range(int(steps)):
        a = scaled_source(m, first, t, gains[first])
        b = scaled_source(m, second, t - int(lag), gains[second])
        m.advance(addsrc(a, b), False, True, 'none')
        acc = np.maximum(acc, np.abs(m.psi) ** 2)
    return acc


def safe_corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) < 3 or np.std(a) < 1e-14 or np.std(b) < 1e-14:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def rankdata(x):
    """Average ranks, enough for a tiny Spearman receipt without scipy."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind='mergesort')
    ranks = np.empty(len(x), float)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        r = 0.5 * (i + j - 1)
        ranks[order[i:j]] = r
        i = j
    return ranks


def spearman(a, b):
    return safe_corr(rankdata(a), rankdata(b))


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w, l)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2 * tail)), w, l


def top_overlap(a, b, frac=.10):
    a = np.asarray(a, float); b = np.asarray(b, float)
    k = max(1, int(math.ceil(frac * len(a))))
    ia = set(np.argsort(a, kind='mergesort')[-k:].tolist())
    ib = set(np.argsort(b, kind='mergesort')[-k:].tolist())
    return float(len(ia & ib) / k)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=24)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/source_gain_balance/source_gain_balance.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    assert abs(safe_corr([0,1,2],[0,2,4]) - 1.0) < 1e-12
    assert abs(spearman([3,1,2],[9,1,4]) - 1.0) < 1e-12
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

        recs = {}
        maps = {}
        for name, (ga, gb) in CONDS.items():
            pA = peak_single(m, 0, ga, a.steps)
            pB = peak_single(m, 1, gb, a.steps)
            T = peak_order(m, ga, gb, a.lag, True, a.steps)
            D = peak_order(m, ga, gb, a.lag, False, a.steps)
            absC_map = np.abs((T - D) / (T + D + 1e-30))
            bal_map = np.minimum(pA, pB) / (np.maximum(pA, pB) + 1e-30)
            cv = np.asarray(absC_map[idx], float)
            bv = np.asarray(bal_map[idx], float)
            recs[name] = dict(
                gain_A=ga, gain_B=gb,
                soma_balance=float(bv[si]),
                soma_selectivity=float(cv[si]),
                soma_balance_percentile=float(np.mean(bv < bv[si])),
                soma_selectivity_percentile=float(np.mean(cv < cv[si])),
                corr_balance_selectivity=safe_corr(bv, cv),
                top_decile_overlap=top_overlap(bv, cv),
                best_balance_selectivity=float(cv[int(np.argmax(bv))]),
                best_selectivity=float(np.max(cv)),
            )
            maps[name] = (bv, cv)

        base_b, base_c = maps['baseline']
        delta_soma_b = []
        delta_soma_c = []
        delta_map_rs = []
        for name in CONDS:
            if name == 'baseline':
                recs[name]['delta_map_corr'] = float('nan')
                continue
            bv, cv = maps[name]
            recs[name]['delta_map_corr'] = safe_corr(bv - base_b, cv - base_c)
            delta_map_rs.append(recs[name]['delta_map_corr'])
            delta_soma_b.append(recs[name]['soma_balance'] - recs['baseline']['soma_balance'])
            delta_soma_c.append(recs[name]['soma_selectivity'] - recs['baseline']['soma_selectivity'])

        soma_bal = [recs[n]['soma_balance'] for n in CONDS]
        soma_sel = [recs[n]['soma_selectivity'] for n in CONDS]
        body = dict(
            seed=seed,
            cells=int(m.body.sum()),
            P1_soma_delta_corr=safe_corr(delta_soma_b, delta_soma_c),
            soma_absolute_pearson=safe_corr(soma_bal, soma_sel),
            soma_absolute_spearman=spearman(soma_bal, soma_sel),
            max_balance_condition=max(CONDS, key=lambda n: recs[n]['soma_balance']),
            max_selectivity_condition=max(CONDS, key=lambda n: recs[n]['soma_selectivity']),
            P2_mean_delta_map_corr=float(np.nanmean(delta_map_rs)),
            conditions=recs,
        )
        rows.append(body)
        print(f"seed {seed:2d}: soma delta r {body['P1_soma_delta_corr']:+.3f}, "
              f"map delta r {body['P2_mean_delta_map_corr']:+.3f}, "
              f"max b/sel {body['max_balance_condition']}/{body['max_selectivity_condition']}", flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    p1 = np.asarray([r['P1_soma_delta_corr'] for r in rows], float)
    p2 = np.asarray([r['P2_mean_delta_map_corr'] for r in rows], float)
    p1p, p1w, p1l = sign_test_two_sided(p1)
    p2p, p2w, p2l = sign_test_two_sided(p2)
    same_max = int(sum(r['max_balance_condition'] == r['max_selectivity_condition'] for r in rows))
    summary = dict(
        bodies=len(rows),
        P1=dict(
            mean_soma_delta_corr=float(np.nanmean(p1)),
            median_soma_delta_corr=float(np.nanmedian(p1)),
            positive_bodies=p1w,
            negative_bodies=p1l,
            sign_p=p1p,
            passed=bool(np.nanmean(p1) > 0 and np.isfinite(p1p) and p1p < .05),
        ),
        P2=dict(
            mean_body_delta_map_corr=float(np.nanmean(p2)),
            median_body_delta_map_corr=float(np.nanmedian(p2)),
            positive_bodies=p2w,
            negative_bodies=p2l,
            sign_p=p2p,
            passed=bool(np.nanmean(p2) > 0 and np.isfinite(p2p) and p2p < .05),
        ),
        descriptive=dict(
            mean_soma_absolute_pearson=float(np.nanmean([r['soma_absolute_pearson'] for r in rows])),
            mean_soma_absolute_spearman=float(np.nanmean([r['soma_absolute_spearman'] for r in rows])),
            same_max_condition=same_max,
            mean_cell_corr_balance_selectivity=float(np.mean([q['corr_balance_selectivity'] for r in rows for q in r['conditions'].values()])),
            mean_top_decile_overlap=float(np.mean([q['top_decile_overlap'] for r in rows for q in r['conditions'].values()])),
        ),
    )

    payload = dict(
        experiment='source_gain_balance_v01',
        prereg='SOURCE_GAIN_BALANCE_PREREG_V01.md',
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

    print('\nSOURCE-GAIN BALANCE RECEIPT')
    print(f" P1 soma delta r {summary['P1']['mean_soma_delta_corr']:+.4f} "
          f"pos/neg {p1w}/{p1l} p={p1p:.5g} PASS={summary['P1']['passed']}")
    print(f" P2 map delta r {summary['P2']['mean_body_delta_map_corr']:+.4f} "
          f"pos/neg {p2w}/{p2l} p={p2p:.5g} PASS={summary['P2']['passed']}")
    print(f" same max balance/selectivity condition {same_max}/{len(rows)}")
    print(' wrote', out)


if __name__ == '__main__':
    main()
