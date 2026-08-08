"""Evaluate LOCKIN_OVERLAP_CONFIRM_PREREG_V01.md on frozen JSON receipts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def corr(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or np.std(a) < 1e-30 or np.std(b) < 1e-30:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--overlap', required=True)
    ap.add_argument('--polarization', required=True)
    ap.add_argument('--out', default='runs/lockin_overlap_confirm/verdict.json')
    a = ap.parse_args()

    ov = json.loads(Path(a.overlap).read_text())
    po = json.loads(Path(a.polarization).read_text())

    rows = ov['rows']
    all_points = []
    for r in rows:
        for K, kres in r['results'].items():
            for h, z in kres['half_period'].items():
                all_points.append(dict(seed=r['seed'], K=int(K), h=int(h), **z))

    max_parseval = max(z['spectral_parseval_error'] for z in all_points)
    c0 = max_parseval < 1e-10

    free = [z for z in all_points if not z['exact_support_collision']]
    hit = [z for z in all_points if z['exact_support_collision']]
    max_free = max([z['leakage_to_gradient_l2'] for z in free], default=float('inf'))
    med_free = float(np.median([z['leakage_to_gradient_l2'] for z in free])) if free else float('inf')
    med_hit = float(np.median([z['leakage_to_gradient_l2'] for z in hit])) if hit else 0.0
    c1 = len(free) >= 24 and max_free < 1e-10
    ratio = med_hit / (med_free + 1e-30)
    c2 = len(hit) >= 24 and med_hit > 0.01 and ratio > 1e6

    overlap_corr = corr(
        [z['weighted_overlap_to_gradient_l2'] for z in all_points],
        [z['leakage_to_gradient_l2'] for z in all_points],
    )
    c3 = overlap_corr > 0.90

    ps = po['summary']['K']
    two_state_max = max(ps[str(K)]['broadband_two_state_max_relative_l2'] for K in (8, 16))
    c4 = two_state_max < 1e-10

    phase = {}
    c5 = True
    for K in (8, 16):
        z = ps[str(K)]['phase_error']['0.1']
        phase[str(K)] = z
        c5 = c5 and z['mean_corr'] > 0.999 and z['mean_relative_l2'] < 0.01 and z['mean_strong_sign_agreement'] > 0.99

    verdict = dict(
        prereg='LOCKIN_OVERLAP_CONFIRM_PREREG_V01.md',
        bodies=len(rows),
        criteria={
            'C0_parseval_identity': dict(pass_=bool(c0), max_normalized_error=float(max_parseval)),
            'C1_collision_free_exact': dict(pass_=bool(c1), points=len(free), max_leakage=float(max_free), median_leakage=float(med_free)),
            'C2_collision_failure_boundary': dict(pass_=bool(c2), points=len(hit), median_leakage=float(med_hit), median_ratio=float(ratio)),
            'C3_weighted_overlap_predicts': dict(pass_=bool(c3), pooled_corr=float(overlap_corr)),
            'C4_broadband_two_state_identity': dict(pass_=bool(c4), max_relative_l2=float(two_state_max)),
            'C5_phase_error_0p1': dict(pass_=bool(c5), K=phase),
        },
    )
    verdict['passed'] = int(sum(int(x['pass_']) for x in verdict['criteria'].values()))
    verdict['total'] = len(verdict['criteria'])

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))
    if verdict['passed'] != verdict['total']:
        raise SystemExit(f"registered confirmation failed: {verdict['passed']}/{verdict['total']}")


if __name__ == '__main__':
    main()
