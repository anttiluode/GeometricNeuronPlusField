"""Held-out confirmation for GeometricNeuronPlusField discovery v0.1.

The choices in this file are frozen from seeds 0-11 before seeds 12-23 are run.
See DISCOVERY_V01.md.  This script intentionally does NOT search for a better band.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from graph_mode_probe import graph_laplacian_modes, trace_order, contrast_from_traces


BAND = (18, 19, 20)
MODE0_LIMIT = 0.005
BAND_GAIN_MIN_POSITIVE = 9
BAND_GAIN_SIGN_P_LIMIT = 0.05
BAND_ABSC_MIN = 0.05


def exact_sign_test(values):
    x = np.asarray(values, float)
    x = x[np.abs(x) > 1e-12]
    n = len(x)
    if n == 0:
        return 1.0, 0, 0
    pos = int((x > 0).sum())
    k = min(pos, n - pos)
    tail = sum(math.comb(n, i) for i in range(k + 1))
    return min(1.0, 2.0 * tail / (2 ** n)), pos, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='../FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=12)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--modes', type=int, default=24)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--probe-steps', type=int, default=150)
    ap.add_argument('--outdir', default='runs/confirm_graph_band')
    a = ap.parse_args()

    if a.modes <= max(BAND):
        raise SystemExit(f'--modes must be > {max(BAND)} for registered band {BAND}')

    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            print(f'seed {seed}: bootstrap FAILED')
            continue
        m.mature = True
        coords, evals, evecs = graph_laplacian_modes(m.body)
        nmodes = min(int(a.modes), len(coords))
        if nmodes <= max(BAND):
            print(f'seed {seed}: only {nmodes} modes; skip')
            continue

        tg = trace_order(m, coords, evecs, nmodes, int(a.lag), True, int(a.probe_steps))
        ds = trace_order(m, coords, evecs, nmodes, int(a.lag), False, int(a.probe_steps))
        point_c = float(contrast_from_traces(tg['point'][:, None], ds['point'][:, None])[0])
        coh = np.asarray(contrast_from_traces(tg['coherent'], ds['coherent']), float)
        incoh = np.asarray(contrast_from_traces(tg['incoherent'], ds['incoherent']), float)

        band = np.asarray(BAND, int)
        band_gain = float(np.mean(np.abs(coh[band]) - np.abs(incoh[band])))
        band_abs_c = float(np.mean(np.abs(coh[band])))
        mode0_abs_c = float(abs(coh[0]))
        band_lambdas = [float(evals[i]) for i in BAND]

        rows.append(dict(
            seed=int(seed),
            boot=boot,
            point_absC=float(abs(point_c)),
            mode0_absC=mode0_abs_c,
            band_absC=band_abs_c,
            band_coherence_gain=band_gain,
            band_eigenvalues=band_lambdas,
            per_mode=[dict(index=i,
                           eigenvalue=float(evals[i]),
                           coherent_absC=float(abs(coh[i])),
                           incoherent_absC=float(abs(incoh[i])),
                           gain=float(abs(coh[i]) - abs(incoh[i]))) for i in BAND],
        ))
        print(f'seed {seed:2d}  mode0 |C|={mode0_abs_c:.6f}  '
              f'band |C|={band_abs_c:.4f}  band gain={band_gain:+.4f}')

    if not rows:
        raise SystemExit('No completed held-out bodies.')

    mode0 = np.asarray([r['mode0_absC'] for r in rows])
    gain = np.asarray([r['band_coherence_gain'] for r in rows])
    band_c = np.asarray([r['band_absC'] for r in rows])
    point = np.asarray([r['point_absC'] for r in rows])
    p, pos, n = exact_sign_test(gain)

    tests = dict(
        common_mode_blind=dict(
            criterion=f'mean mode0 |C| < {MODE0_LIMIT}',
            value=float(mode0.mean()),
            pass_=bool(mode0.mean() < MODE0_LIMIT),
        ),
        band_coherence_advantage=dict(
            criterion=f'>= {BAND_GAIN_MIN_POSITIVE}/{len(rows)} positive and sign p < {BAND_GAIN_SIGN_P_LIMIT}',
            mean_gain=float(gain.mean()),
            positive=int(pos),
            n=int(n),
            sign_test_p=float(p),
            pass_=bool(pos >= BAND_GAIN_MIN_POSITIVE and p < BAND_GAIN_SIGN_P_LIMIT),
        ),
        band_informative=dict(
            criterion=f'mean coherent band |C| > {BAND_ABSC_MIN}',
            value=float(band_c.mean()),
            pass_=bool(band_c.mean() > BAND_ABSC_MIN),
        ),
    )
    all_pass = all(t['pass_'] for t in tests.values())

    payload = dict(
        experiment='heldout_graph_band_v01',
        discovery_seeds='0-11',
        heldout_seed_start=int(a.seed_start),
        heldout_seeds_requested=int(a.seeds),
        heldout_seeds_completed=len(rows),
        registered_band=list(BAND),
        lag=int(a.lag),
        no_claim_band_beats_point=True,
        summary=dict(
            mode0_absC_mean=float(mode0.mean()),
            band_absC_mean=float(band_c.mean()),
            band_gain_mean=float(gain.mean()),
            point_absC_mean=float(point.mean()),
        ),
        tests=tests,
        all_registered_tests_pass=bool(all_pass),
        rows=rows,
    )

    with open(outdir / 'confirmation_results.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2)

    with open(outdir / 'CONFIRMATION.txt', 'w', encoding='utf-8') as f:
        f.write('GeometricNeuronPlusField held-out graph-band confirmation\n')
        f.write(f'seeds {a.seed_start}-{a.seed_start + a.seeds - 1}; completed {len(rows)}\n')
        f.write(f'mode0 mean |C|: {mode0.mean():.8f}\n')
        f.write(f'band 18-20 mean |C|: {band_c.mean():.6f}\n')
        f.write(f'band mean coherence gain: {gain.mean():+.6f}\n')
        f.write(f'band gain positive: {pos}/{n}; sign p={p:.8f}\n')
        f.write(f'point mean |C|: {point.mean():.6f}\n')
        f.write(f'all registered tests pass: {all_pass}\n')
        for name, test in tests.items():
            f.write(f'{name}: {test}\n')

    print('\nHELD-OUT CONFIRMATION')
    print(f'  mode0 mean |C|          {mode0.mean():.8f}  pass={tests["common_mode_blind"]["pass_"]}')
    print(f'  band 18-20 mean |C|     {band_c.mean():.6f}  pass={tests["band_informative"]["pass_"]}')
    print(f'  band coherence gain     {gain.mean():+.6f}')
    print(f'  gain positive / sign p  {pos}/{n}  p={p:.8f}  pass={tests["band_coherence_advantage"]["pass_"]}')
    print(f'  point mean |C|          {point.mean():.6f}  (context only; no superiority claim)')
    print(f'  ALL REGISTERED PASS     {all_pass}')


if __name__ == '__main__':
    main()
