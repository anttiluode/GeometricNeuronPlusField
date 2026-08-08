"""Hostile controls for the fast 50%-duty return-gating result.

The held-out coarse-return test showed that a P=14 periodic 50%-duty gate preserved the
exact bond-gradient direction at ~.999 correlation for every tested phase offset,
while slower P=42 gating became phase-sensitive.

Before interpreting that as a rhythmic mechanism, compare equal-duty masks:

  periodic   deterministic square-wave masks at several periods
  random     exactly half the samples selected uniformly without replacement
  block      one contiguous circular half-window
  comb       evenly spaced isolated samples (same 50% duty: alternating samples)

All masks are applied to the exact soma return waveform separately for target and
distractor and then L2-dose matched exactly as in biological_return_code_probe.py.

The core question is whether *rhythm* is special or whether the gradient is simply
redundant under distributed temporal subsampling.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from biological_return_code_probe import build_exact, gradient_from_codes, map_metrics, l2_match

PERIODS = (2, 6, 10, 14, 30, 42, 70)
N_RANDOM = 24
N_BLOCK = 24


def apply_mask(g, mask):
    g = np.asarray(g, np.complex128)
    mask = np.asarray(mask, bool)
    return l2_match(g * mask, g)


def evaluate_mask(m, z, mask):
    gt = apply_mask(z['gT'], mask)
    gd = apply_mask(z['gD'], mask)
    h, v = gradient_from_codes(m, z['wh'], z['wv'], z['pT'], z['pD'], gt, gd)
    mm = map_metrics(z['exact_h'], z['exact_v'], h, v)
    mm['duty'] = float(np.mean(mask))
    return mm


def periodic_mask(T, period, offset):
    p = int(period)
    on = p // 2
    t = np.arange(T)
    return ((t + int(offset)) % p) < on


def block_mask(T, offset):
    n = T // 2
    idx = (np.arange(n) + int(offset)) % T
    m = np.zeros(T, bool)
    m[idx] = True
    return m


def random_mask(T, rng):
    n = T // 2
    idx = rng.choice(T, size=n, replace=False)
    m = np.zeros(T, bool)
    m[idx] = True
    return m


def summarize(q):
    return dict(
        n=len(q),
        mean_corr=float(np.mean([x['corr'] for x in q])),
        median_corr=float(np.median([x['corr'] for x in q])),
        min_corr=float(np.min([x['corr'] for x in q])),
        max_corr=float(np.max([x['corr'] for x in q])),
        mean_sign=float(np.mean([x['strong_sign_agreement'] for x in q])),
        mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
    )


def one(m, lag, steps):
    z = build_exact(m, lag, steps)
    T = len(z['gT'])
    if T % 2:
        raise ValueError('control requires even T')

    periodic = {}
    for P in PERIODS:
        if P % 2:
            continue
        q = [evaluate_mask(m, z, periodic_mask(T, P, off)) for off in range(P)]
        periodic[str(P)] = summarize(q)

    rng = np.random.default_rng(int(m.cfg.seed) + 20260808)
    random_q = [evaluate_mask(m, z, random_mask(T, rng)) for _ in range(N_RANDOM)]

    offsets = np.linspace(0, T-1, N_BLOCK, dtype=int)
    block_q = [evaluate_mask(m, z, block_mask(T, int(off))) for off in offsets]

    # Alternating samples are exactly P=2 periodic gating, but retain an explicit
    # named control because this is the extreme distributed-subsampling limit.
    comb_q = [evaluate_mask(m, z, periodic_mask(T, 2, off)) for off in (0, 1)]

    return dict(
        seed=int(m.cfg.seed),
        C=float(z['C']),
        periodic=periodic,
        random=summarize(random_q),
        block=summarize(block_q),
        comb=summarize(comb_q),
    )


def group_summary(rows):
    out = dict(bodies=len(rows), periodic={})
    for P in map(str, PERIODS):
        q = [r['periodic'][P] for r in rows]
        out['periodic'][P] = dict(
            mean_of_mean_corr=float(np.mean([x['mean_corr'] for x in q])),
            mean_of_min_corr=float(np.mean([x['min_corr'] for x in q])),
            worst_body_min_corr=float(np.min([x['min_corr'] for x in q])),
            mean_sign=float(np.mean([x['mean_sign'] for x in q])),
        )
    for name in ('random', 'block', 'comb'):
        q = [r[name] for r in rows]
        out[name] = dict(
            mean_of_mean_corr=float(np.mean([x['mean_corr'] for x in q])),
            mean_of_min_corr=float(np.mean([x['min_corr'] for x in q])),
            worst_body_min_corr=float(np.min([x['min_corr'] for x in q])),
            mean_sign=float(np.mean([x['mean_sign'] for x in q])),
        )
    return out


def selftest():
    T = 210
    for P in PERIODS:
        for off in (0, 1, P-1):
            m = periodic_mask(T, P, off)
            # Some periods do not divide T, so duty can differ from .5 by at most
            # one partial-period edge; P=2/6/10/14/30/42/70 all divide 210 here.
            assert abs(np.mean(m) - .5) < 1e-12
    for off in (0, 17, 209):
        assert int(np.sum(block_mask(T, off))) == T//2
    rng = np.random.default_rng(1)
    assert int(np.sum(random_mask(T, rng))) == T//2
    print('selftest ok')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=510)
    ap.add_argument('--seeds', type=int, default=6)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/return_gate_controls/dev.json')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest(); return

    fa = Path(a.functional_arbors).resolve()
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    rows = []
    for seed in range(a.seed_start, a.seed_start+a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        b = m.bootstrap()
        if not b.get('ok'):
            continue
        m.mature = True
        r = one(m, a.lag, a.steps)
        rows.append(r)
        print(
            'seed', seed,
            'P2', round(r['periodic']['2']['min_corr'], 4),
            'P14', round(r['periodic']['14']['min_corr'], 4),
            'rand', round(r['random']['mean_corr'], 4),
            'rand-min', round(r['random']['min_corr'], 4),
            'block', round(r['block']['mean_corr'], 4),
            'block-min', round(r['block']['min_corr'], 4),
            flush=True,
        )
    if not rows:
        raise SystemExit('No valid bodies')
    s = group_summary(rows)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(experiment='return_gate_controls_dev_v01', summary=s, rows=rows), indent=2))
    print('\nRETURN GATE CONTROLS DEV')
    print(json.dumps(s, indent=2))


if __name__ == '__main__':
    main()
