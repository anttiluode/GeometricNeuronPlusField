"""Evaluate RETURN_GATE_CONTROLS_CONFIRM_PREREG_V01.md from a JSON receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--receipt', required=True)
    ap.add_argument('--out', default='runs/return_gate_controls_confirm/verdict.json')
    a = ap.parse_args()

    p = json.loads(Path(a.receipt).read_text())
    s = p['summary']
    q = s['periodic']
    rnd = s['random']
    blk = s['block']
    comb = s['comb']

    G0 = (q['6']['mean_of_mean_corr'] > .999 and
          q['6']['mean_of_min_corr'] > .999 and
          q['6']['worst_body_min_corr'] > .995)

    G1 = (q['14']['mean_of_mean_corr'] > .995 and
          q['14']['mean_of_min_corr'] > .990 and
          q['6']['mean_of_mean_corr'] - q['14']['mean_of_mean_corr'] > .0005)

    p14_rand_mean = q['14']['mean_of_mean_corr'] - rnd['mean_of_mean_corr']
    p14_rand_min = q['14']['mean_of_min_corr'] - rnd['mean_of_min_corr']
    G2 = p14_rand_mean > .008 and p14_rand_min > .05

    rand_block_mean = rnd['mean_of_mean_corr'] - blk['mean_of_mean_corr']
    rand_block_min = rnd['mean_of_min_corr'] - blk['mean_of_min_corr']
    G3 = rand_block_mean > .15 and rand_block_min > .30

    means = [q[str(P)]['mean_of_mean_corr'] for P in (6,14,30,42,70)]
    ordering = all(means[i] > means[i+1] for i in range(len(means)-1))
    p14_p42 = q['14']['mean_of_mean_corr'] - q['42']['mean_of_mean_corr']
    p42_p70 = q['42']['mean_of_mean_corr'] - q['70']['mean_of_mean_corr']
    G4 = ordering and p14_p42 > .10 and p42_p70 > .05

    comb_rand_mean = comb['mean_of_mean_corr'] - rnd['mean_of_mean_corr']
    comb_rand_min = comb['mean_of_min_corr'] - rnd['mean_of_min_corr']
    G5 = comb_rand_mean > .01 and comb_rand_min > .05

    verdict = {
        'prereg': 'RETURN_GATE_CONTROLS_CONFIRM_PREREG_V01.md',
        'bodies': s['bodies'],
        'criteria': {
            'G0_fast_regular_transparent': {
                'pass_': G0, 'P6': q['6']},
            'G1_P14_robust_but_below_fast': {
                'pass_': G1, 'P14': q['14'],
                'P6_minus_P14_mean': q['6']['mean_of_mean_corr'] - q['14']['mean_of_mean_corr']},
            'G2_P14_beats_random': {
                'pass_': G2,
                'mean_advantage': p14_rand_mean,
                'min_advantage': p14_rand_min},
            'G3_random_beats_block': {
                'pass_': G3,
                'mean_advantage': rand_block_mean,
                'min_advantage': rand_block_min},
            'G4_period_orders_degradation': {
                'pass_': G4,
                'means_P6_P14_P30_P42_P70': means,
                'P14_minus_P42': p14_p42,
                'P42_minus_P70': p42_p70},
            'G5_comb_beats_random': {
                'pass_': G5,
                'mean_advantage': comb_rand_mean,
                'min_advantage': comb_rand_min},
        },
        'summary': s,
    }
    verdict['passed'] = sum(int(v['pass_']) for v in verdict['criteria'].values())
    verdict['total'] = len(verdict['criteria'])

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))
    if verdict['passed'] != verdict['total']:
        raise SystemExit(f"registered confirmation failed: {verdict['passed']}/{verdict['total']}")


if __name__ == '__main__':
    main()
