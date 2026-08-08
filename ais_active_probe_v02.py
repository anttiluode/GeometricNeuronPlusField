"""AIS active-boundary v0.2: timing precision under exact per-frequency rate matching.

The active membrane and all upstream physics are imported unchanged from
`ais_active_probe.py` (v0.1).  This file changes only the controls/estimators:

* memoryless and linearized controls get an oracle top-k rate match separately
  at each frequency;
* timing is interpreted only when the active boundary emits >= 4 events;
* vector strength is accompanied by spike-count-unbiased PPC and cycle coverage;
* A->B/B->A controls are rate-matched over the task pair only.

See AIS_ACTIVE_PREREG_V02.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from ais_active_probe import (
    AISConfig, ActiveAIS, task_trace, modulated_trace, normalize_traces,
    linear_score, in_window, vector_strength, first_or_nan, centroid_or_nan,
    safe_absdiff, count_contrast,
)


def topk_events(score: np.ndarray, mask: np.ndarray, k: int) -> np.ndarray:
    """Exact event-count match using the k largest allowed scalar samples.

    This is an intentionally strong oracle control.  Stable sorting provides a
    deterministic tie break; event times are frame centers.
    """
    score = np.asarray(score, float)
    ids = np.flatnonzero(np.asarray(mask, bool))
    k = int(max(0, min(int(k), len(ids))))
    if k == 0:
        return np.asarray([], float)
    order = np.argsort(-score[ids], kind='mergesort')[:k]
    return np.sort(ids[order].astype(float) + 0.5)


def pooled_topk(score_a: np.ndarray, score_b: np.ndarray, k: int):
    """Exact total-rate match over a pair of traces using one pooled ranking."""
    a = np.asarray(score_a, float); b = np.asarray(score_b, float)
    vals = np.concatenate([a, b])
    k = int(max(0, min(int(k), len(vals))))
    if k == 0:
        return np.asarray([], float), np.asarray([], float)
    ids = np.argsort(-vals, kind='mergesort')[:k]
    ia = np.sort(ids[ids < len(a)].astype(float) + 0.5)
    ib0 = ids[ids >= len(a)] - len(a)
    ib = np.sort(ib0.astype(float) + 0.5)
    return ia, ib


def ppc(times: np.ndarray, freq: float) -> float:
    """Pairwise phase consistency; unbiased by event count for N>=2."""
    t = np.asarray(times, float)
    n = len(t)
    if n < 2:
        return float('nan')
    z = np.exp(2j * math.pi * freq * t).sum()
    return float((abs(z) ** 2 - n) / (n * (n - 1)))


def phase_metrics(times: np.ndarray, freq: float, start: float, end: float):
    z = in_window(times, start, end)
    n = len(z)
    cycles = max((end - start) * freq, 1e-12)
    vs = vector_strength(z, freq)
    return dict(
        events=int(n),
        vector_strength=float(vs),
        ppc=ppc(z, freq),
        events_per_cycle=float(n / cycles),
        cycle_coverage=float(min(1.0, n / cycles)),
        coverage_score=float(vs * min(1.0, n / cycles)),
    )


def sign_test_two_sided(wins: int, losses: int) -> float:
    n = int(wins + losses)
    if n <= 0:
        return float('nan')
    k = min(int(wins), int(losses))
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2.0 * tail))


def task_metrics(a: np.ndarray, b: np.ndarray):
    return dict(
        events_T=int(len(a)), events_D=int(len(b)),
        count_contrast=count_contrast(a, b),
        first_latency_delta=safe_absdiff(first_or_nan(a), first_or_nan(b)),
        centroid_delta=safe_absdiff(centroid_or_nan(a), centroid_or_nan(b)),
    )


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='../FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=0)
    ap.add_argument('--seeds', type=int, default=24)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--task-steps', type=int, default=180)
    ap.add_argument('--freqs', default='0.00625,0.0125,0.025,0.05,0.0833333,0.125')
    ap.add_argument('--freq-steps', type=int, default=640)
    ap.add_argument('--burn', type=int, default=160)
    ap.add_argument('--min-events', type=int, default=4)
    ap.add_argument('--upper-min-freq', type=float, default=0.05)
    ap.add_argument('--out', default='runs/ais_active_v02/ais_active_v02.json')
    return ap.parse_args()


def main():
    a = parse_args()
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    freqs = [float(x) for x in a.freqs.split(',') if x.strip()]
    cfg = AISConfig()  # MUST remain identical to v0.1
    ais = ActiveAIS(cfg)
    kernel = ais.linear_kernel()
    rows = []

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature = True

        raw = {
            'task_T': task_trace(m, a.lag, True, a.task_steps),
            'task_D': task_trace(m, a.lag, False, a.task_steps),
        }
        for j, f in enumerate(freqs):
            raw[f'f{j}'] = modulated_trace(m, f, a.freq_steps, source=0)
        x, scale = normalize_traces(raw)

        active = {}
        linear = {}
        for key, tr in x.items():
            _, sp = ais.run(tr)
            active[key] = sp
            linear[key] = linear_score(tr, kernel, cfg.input_gain)

        freq_rows = []
        for j, f in enumerate(freqs):
            key = f'f{j}'; end = len(x[key])
            az = in_window(active[key], a.burn, end)
            k = len(az)
            mask = np.arange(end) >= a.burn
            mz = topk_events(x[key], mask, k)
            lz = topk_events(linear[key], mask, k)
            # Exact rate matching is a hard validity condition.
            if len(mz) != k or len(lz) != k:
                raise RuntimeError('per-frequency rate match failed')
            freq_rows.append(dict(
                freq=f,
                valid_timing=bool(k >= a.min_events),
                active=phase_metrics(active[key], f, a.burn, end),
                memoryless=phase_metrics(mz, f, a.burn, end),
                linearized=phase_metrics(lz, f, a.burn, end),
            ))

        # Task controls: one pooled top-k over T+D, matched to total active task events.
        at = active['task_T']; ad = active['task_D']; kt = len(at) + len(ad)
        mt, md = pooled_topk(x['task_T'], x['task_D'], kt)
        lt, ld = pooled_topk(linear['task_T'], linear['task_D'], kt)
        if len(mt) + len(md) != kt or len(lt) + len(ld) != kt:
            raise RuntimeError('task-pair rate match failed')
        task = dict(
            active=task_metrics(at, ad),
            memoryless=task_metrics(mt, md),
            linearized=task_metrics(lt, ld),
            matched_total_events=int(kt),
        )

        rows.append(dict(seed=seed, cells=int(m.body.sum()), scale=scale,
                         frequency=freq_rows, task=task))
        valid_upper = [q for q in freq_rows if q['freq'] >= a.upper_min_freq and q['valid_timing']]
        if valid_upper:
            da = np.mean([q['active']['ppc'] - q['linearized']['ppc'] for q in valid_upper])
            print(f'seed {seed:2d}: valid upper {len(valid_upper)}  mean active-linear PPC {da:+.3f}', flush=True)
        else:
            print(f'seed {seed:2d}: valid upper 0', flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    # Per-frequency receipts.
    freq_summary = []
    for j, f in enumerate(freqs):
        qs = [r['frequency'][j] for r in rows]
        valid = [q for q in qs if q['valid_timing']]
        item = dict(freq=f, bodies=len(qs), valid_pairs=len(valid),
                    mean_active_events=float(np.mean([q['active']['events'] for q in qs])))
        for name in ('memoryless', 'linearized', 'active'):
            item[name] = dict(
                vs_mean_valid=float(np.mean([q[name]['vector_strength'] for q in valid])) if valid else float('nan'),
                ppc_mean_valid=float(np.mean([q[name]['ppc'] for q in valid])) if valid else float('nan'),
                coverage_score_mean=float(np.mean([q[name]['coverage_score'] for q in qs])),
            )
        if valid:
            dlin = np.asarray([q['active']['ppc'] - q['linearized']['ppc'] for q in valid], float)
            dmem = np.asarray([q['active']['ppc'] - q['memoryless']['ppc'] for q in valid], float)
            item['active_minus_linear_ppc_mean'] = float(dlin.mean())
            item['active_minus_memoryless_ppc_mean'] = float(dmem.mean())
            item['active_beats_linear'] = int(np.sum(dlin > 0))
            item['active_beats_memoryless'] = int(np.sum(dmem > 0))
        freq_summary.append(item)

    # Registered upper-band test, pooled over all valid body/frequency pairs.
    upper = []
    for r in rows:
        for q in r['frequency']:
            if q['freq'] >= a.upper_min_freq and q['valid_timing']:
                upper.append((r['seed'], q))
    dlin = np.asarray([q['active']['ppc'] - q['linearized']['ppc'] for _, q in upper], float)
    dmem = np.asarray([q['active']['ppc'] - q['memoryless']['ppc'] for _, q in upper], float)
    wlin = int(np.sum(dlin > 0)); llin = int(np.sum(dlin < 0))
    wmem = int(np.sum(dmem > 0)); lmem = int(np.sum(dmem < 0))
    upper_summary = dict(
        valid_pairs=len(upper),
        unique_bodies=len(set(seed for seed, _ in upper)),
        active_ppc_mean=float(np.mean([q['active']['ppc'] for _, q in upper])) if upper else float('nan'),
        linearized_ppc_mean=float(np.mean([q['linearized']['ppc'] for _, q in upper])) if upper else float('nan'),
        memoryless_ppc_mean=float(np.mean([q['memoryless']['ppc'] for _, q in upper])) if upper else float('nan'),
        active_minus_linear_ppc_mean=float(dlin.mean()) if len(dlin) else float('nan'),
        active_minus_memoryless_ppc_mean=float(dmem.mean()) if len(dmem) else float('nan'),
        active_beats_linear=wlin, active_loses_linear=llin,
        active_beats_linear_sign_p=sign_test_two_sided(wlin, llin),
        active_beats_memoryless=wmem, active_loses_memoryless=lmem,
        active_beats_memoryless_sign_p=sign_test_two_sided(wmem, lmem),
    )

    task_summary = {}
    for name in ('memoryless', 'linearized', 'active'):
        cc = np.asarray([abs(r['task'][name]['count_contrast']) for r in rows], float)
        fd = np.asarray([r['task'][name]['first_latency_delta'] for r in rows], float)
        cd = np.asarray([r['task'][name]['centroid_delta'] for r in rows], float)
        task_summary[name] = dict(
            abs_count_contrast_mean=float(np.nanmean(cc)),
            first_latency_delta_mean=float(np.nanmean(fd)) if np.isfinite(fd).any() else float('nan'),
            first_latency_valid=int(np.isfinite(fd).sum()),
            centroid_delta_mean=float(np.nanmean(cd)) if np.isfinite(cd).any() else float('nan'),
            centroid_valid=int(np.isfinite(cd).sum()),
        )

    payload = dict(
        experiment='ais_active_boundary_v02',
        mechanism_frozen_from='ais_active_boundary_v01',
        registered=dict(min_events=a.min_events, upper_min_freq=a.upper_min_freq),
        config=dict(lag=a.lag, task_steps=a.task_steps, freqs=freqs,
                    freq_steps=a.freq_steps, burn=a.burn),
        summary=dict(frequency=freq_summary, upper_band=upper_summary, task=task_summary),
        rows=rows,
    )
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nAIS ACTIVE v0.2 — PER-FREQUENCY RATE-MATCHED PRECISION')
    print(f'  upper-band valid pairs / bodies  {upper_summary["valid_pairs"]} / {upper_summary["unique_bodies"]}')
    print(f'  PPC active       {upper_summary["active_ppc_mean"]:+.4f}')
    print(f'  PPC linearized   {upper_summary["linearized_ppc_mean"]:+.4f}')
    print(f'  PPC memoryless   {upper_summary["memoryless_ppc_mean"]:+.4f}')
    print(f'  active-linear    {upper_summary["active_minus_linear_ppc_mean"]:+.4f}  '
          f'w/l {wlin}/{llin}  sign p={upper_summary["active_beats_linear_sign_p"]:.5g}')
    print(f'  active-memoryless {upper_summary["active_minus_memoryless_ppc_mean"]:+.4f}  '
          f'w/l {wmem}/{lmem}  sign p={upper_summary["active_beats_memoryless_sign_p"]:.5g}')
    print(f'  wrote {out}')


if __name__ == '__main__':
    main()
