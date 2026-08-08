"""Evaluate BIOLOGICAL_RETURN_CONFIRM_PREREG_V01.md from a JSON receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--receipt', required=True)
    ap.add_argument('--out', default='runs/biological_return_confirm/verdict.json')
    a = ap.parse_args()

    p = json.loads(Path(a.receipt).read_text())
    rows = p['rows']
    s = p['summary']
    c = s['codes']
    g = s['gates']

    R0 = c['exact']['mean_corr'] > .999999 and c['exact']['mean_relative_l2'] < 1e-10
    R1 = c['real_wave']['mean_corr'] > .995 and c['real_wave']['mean_strong_sign'] > .99

    env_body = [r['codes']['envelope_signed']['corr'] for r in rows]
    R2 = (c['envelope_signed']['mean_corr'] > .80 and
          c['envelope_signed']['median_corr'] > .80 and
          sum(x > .75 for x in env_body) >= 9)

    R3 = (g['14']['mean_median_corr'] > .99 and
          g['14']['mean_worst_corr'] > .98 and
          g['14']['mean_median_sign'] > .96)

    fast_slow_median = g['14']['mean_median_corr'] - g['42']['mean_median_corr']
    fast_slow_worst = g['14']['mean_worst_corr'] - g['42']['mean_worst_corr']
    R4 = fast_slow_median > .03 and fast_slow_worst > .15

    sp32_body = [r['codes']['sparse_phase_32']['corr'] for r in rows]
    R5 = c['sparse_phase_32']['mean_corr'] > .85 and sum(x > .80 for x in sp32_body) >= 8

    signed_minus_positive = c['sparse_signed_32']['mean_corr'] - c['sparse_positive_32']['mean_corr']
    R6 = signed_minus_positive > .50 and abs(c['sparse_positive_32']['mean_corr']) < .25

    env_minus_phase = c['envelope_signed']['mean_corr'] - c['phase_only']['mean_corr']
    R7 = env_minus_phase > .08

    verdict = dict(
        prereg='BIOLOGICAL_RETURN_CONFIRM_PREREG_V01.md',
        bodies=len(rows),
        criteria={
            'R0_exact_positive_control': dict(pass_=bool(R0), corr=c['exact']['mean_corr'], rel_l2=c['exact']['mean_relative_l2']),
            'R1_real_wave': dict(pass_=bool(R1), corr=c['real_wave']['mean_corr'], strong_sign=c['real_wave']['mean_strong_sign'], rel_l2=c['real_wave']['mean_relative_l2']),
            'R2_envelope_signed': dict(pass_=bool(R2), mean_corr=c['envelope_signed']['mean_corr'], median_corr=c['envelope_signed']['median_corr'], bodies_gt_0p75=int(sum(x > .75 for x in env_body))),
            'R3_fast_gate': dict(pass_=bool(R3), median_corr=g['14']['mean_median_corr'], worst_corr=g['14']['mean_worst_corr'], median_sign=g['14']['mean_median_sign']),
            'R4_scale_separation': dict(pass_=bool(R4), median_advantage=float(fast_slow_median), worst_advantage=float(fast_slow_worst), slow_median=g['42']['mean_median_corr'], slow_worst=g['42']['mean_worst_corr']),
            'R5_sparse_phase_32': dict(pass_=bool(R5), mean_corr=c['sparse_phase_32']['mean_corr'], bodies_gt_0p8=int(sum(x > .80 for x in sp32_body))),
            'R6_consequence_sign': dict(pass_=bool(R6), signed_corr=c['sparse_signed_32']['mean_corr'], positive_corr=c['sparse_positive_32']['mean_corr'], difference=float(signed_minus_positive)),
            'R7_envelope_beats_phase_only': dict(pass_=bool(R7), envelope_corr=c['envelope_signed']['mean_corr'], phase_only_corr=c['phase_only']['mean_corr'], difference=float(env_minus_phase)),
        },
        descriptive=dict(
            delayed_sparse8={k: v for k, v in c.items() if k.startswith('delayed_sparse8_')},
            sparse_curve={k: v for k, v in c.items() if k.startswith('sparse_phase_') or k.startswith('sparse_signed_') or k.startswith('sparse_positive_')},
            gates=g,
        )
    )
    verdict['passed'] = sum(int(x['pass_']) for x in verdict['criteria'].values())
    verdict['total'] = len(verdict['criteria'])

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))
    if verdict['passed'] != verdict['total']:
        raise SystemExit(f"registered confirmation failed: {verdict['passed']}/{verdict['total']}")


if __name__ == '__main__':
    main()
