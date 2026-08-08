"""Active-boundary test downstream of the frozen FunctionalArbor soma.

This is deliberately a *boundary* experiment, not another growth experiment.
The arbor is bootstrapped once and frozen.  Its ordinary soma power trace then
feeds three event encoders:

  1. memoryless: the raw soma drive, thresholded at a globally rate-matched level;
  2. linearized: the small-signal impulse response of the SAME active membrane,
     convolved with the drive and thresholded at a globally rate-matched level;
  3. active: a compact Hodgkin-Huxley-like membrane with Na activation/inactivation
     and K activation state.  No derivative/high-pass term is inserted by hand.

The comparison asks whether active state buys anything beyond (a) instantaneous
amplitude selection and (b) the membrane's own passive/small-signal temporal
filter.  The main readout is phase locking to envelope modulation of the frozen
arbor field, with an A->B/B->A task receipt as a secondary check.

Important: one field frame is mapped to 1 ms *only as an explicit simulation
scale for the boundary ODE*.  Reported envelope frequencies are kept in
cycles/frame; no biological calibration is claimed.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Frozen-arbor drives


def reset_fast(m):
    try:
        m.reset_fast(clear_traces=True)
    except TypeError:
        m.reset_fast(True)


def addsrc(a, b):
    if isinstance(a, (float, int, np.floating)):
        return b
    if isinstance(b, (float, int, np.floating)):
        return a
    return a + b


def task_trace(m, lag: int, target: bool, steps: int) -> np.ndarray:
    """Ordinary soma power for A->B or B->A."""
    reset_fast(m)
    first, second = (0, 1) if target else (1, 0)
    out = np.zeros(steps, float)
    for t in range(steps):
        a = m.pulse_source(first, t, False)
        b = m.pulse_source(second, t - lag, False)
        out[t] = m.advance(addsrc(a, b), False, True, 'none')
    return out


def modulated_trace(m, freq: float, steps: int, source: int = 0,
                    floor: float = 0.12, depth: float = 0.88) -> np.ndarray:
    """Continuously drive one source with a sinusoidally modulated carrier.

    The body and carrier physics are unchanged.  Only the source envelope is
    extended from a short pulse to a positive periodic modulation so the
    downstream boundary can be tested over a frequency sweep.
    """
    reset_fast(m)
    p = m.source_terminal(source)
    if p is None:
        return np.zeros(steps, float)
    out = np.zeros(steps, float)
    c = m.cfg
    for t in range(steps):
        # Starts at the envelope floor to avoid a gratuitous onset transient.
        cyc = 0.5 * (1.0 - math.cos(2.0 * math.pi * freq * t))
        env = floor + depth * cyc
        src = np.zeros_like(m.psi)
        src[p] = c.source_amp * env * np.exp(1j * c.carrier_omega * t)
        out[t] = m.advance(src, False, True, 'none')
    return out


# ---------------------------------------------------------------------------
# Active boundary


@dataclass
class AISConfig:
    # Classic compact HH conductances.  The model is used here as a stateful
    # event boundary, not as a fitted reconstruction of a particular AIS.
    C: float = 1.0
    gNa: float = 120.0
    gK: float = 36.0
    gL: float = 0.3
    ENa: float = 50.0
    EK: float = -77.0
    EL: float = -54.4
    V0: float = -65.0
    input_gain: float = 15.0
    frame_ms: float = 1.0
    ode_dt_ms: float = 0.025
    spike_threshold_mv: float = 0.0
    clip_drive: float = 2.5


def _vtrap(x: float, y: float) -> float:
    z = x / y
    if abs(z) < 1e-6:
        return y * (1.0 + 0.5 * z)
    return x / (1.0 - math.exp(-z))


def _rates(V: float):
    am = 0.1 * _vtrap(V + 40.0, 10.0)
    bm = 4.0 * math.exp(-(V + 65.0) / 18.0)
    ah = 0.07 * math.exp(-(V + 65.0) / 20.0)
    bh = 1.0 / (1.0 + math.exp(-(V + 35.0) / 10.0))
    an = 0.01 * _vtrap(V + 55.0, 10.0)
    bn = 0.125 * math.exp(-(V + 65.0) / 80.0)
    return am, bm, ah, bh, an, bn


def _steady_gates(V: float):
    am, bm, ah, bh, an, bn = _rates(V)
    return am / (am + bm), ah / (ah + bh), an / (an + bn)


class ActiveAIS:
    """Small stateful output compartment with Na/K gating state."""
    def __init__(self, cfg: AISConfig | None = None):
        self.cfg = cfg or AISConfig()

    def run_current(self, current_frames: np.ndarray, detect_spikes: bool = True):
        c = self.cfg
        V = float(c.V0)
        m, h, n = _steady_gates(V)
        sub = max(1, int(round(c.frame_ms / c.ode_dt_ms)))
        dt = c.frame_ms / sub
        volts = np.zeros(len(current_frames), float)
        spikes = []
        tm = 0.0
        prev = V
        for i, I in enumerate(np.asarray(current_frames, float)):
            for _ in range(sub):
                am, bm, ah, bh, an, bn = _rates(V)
                dm = am * (1.0 - m) - bm * m
                dh = ah * (1.0 - h) - bh * h
                dn = an * (1.0 - n) - bn * n
                INa = c.gNa * (m ** 3) * h * (V - c.ENa)
                IK = c.gK * (n ** 4) * (V - c.EK)
                IL = c.gL * (V - c.EL)
                dV = (float(I) - INa - IK - IL) / c.C
                prev = V
                V += dt * dV
                m = float(np.clip(m + dt * dm, 0.0, 1.0))
                h = float(np.clip(h + dt * dh, 0.0, 1.0))
                n = float(np.clip(n + dt * dn, 0.0, 1.0))
                tm += dt
                if detect_spikes and prev < c.spike_threshold_mv <= V:
                    spikes.append(tm / c.frame_ms)  # field-frame coordinates
            volts[i] = V
        return volts, np.asarray(spikes, float)

    def run(self, normalized_drive: np.ndarray):
        x = np.clip(np.asarray(normalized_drive, float), 0.0, self.cfg.clip_drive)
        return self.run_current(self.cfg.input_gain * x, True)

    def linear_kernel(self, frames: int = 180, eps_current: float = 0.02):
        """Finite-difference small-signal impulse response about quiescent rest.

        This is the control Claude requested: the filter is not chosen by hand.
        It is measured from this membrane's own state equations below spike
        threshold, then used as an LTI filter before a rate-matched threshold.
        """
        z = np.zeros(frames, float)
        base, _ = self.run_current(z, False)
        imp = z.copy(); imp[0] = eps_current
        pert, _ = self.run_current(imp, False)
        return (pert - base) / eps_current


# ---------------------------------------------------------------------------
# Controls and metrics


def normalize_traces(traces: dict[str, np.ndarray], q: float = 0.95):
    vals = np.concatenate([v[np.isfinite(v)] for v in traces.values() if len(v)])
    pos = vals[vals > 0]
    scale = float(np.quantile(pos, q)) if len(pos) else 1.0
    scale = max(scale, 1e-12)
    return {k: np.asarray(v, float) / scale for k, v in traces.items()}, scale


def linear_score(x: np.ndarray, kernel: np.ndarray, input_gain: float):
    cur = input_gain * np.asarray(x, float)
    return np.convolve(cur, kernel, mode='full')[:len(cur)]


def threshold_for_total(scores: list[np.ndarray], masks: list[np.ndarray], target: int):
    pool = np.concatenate([s[m] for s, m in zip(scores, masks)])
    if target <= 0 or len(pool) == 0:
        return math.inf
    if target >= len(pool):
        return -math.inf
    ss = np.sort(pool)[::-1]
    a = float(ss[target - 1]); b = float(ss[target])
    if a == b:
        return a - max(1e-12, abs(a) * 1e-12)
    return 0.5 * (a + b)


def frame_events(score: np.ndarray, threshold: float):
    if not np.isfinite(threshold):
        return np.asarray([], float) if threshold > 0 else np.arange(len(score), dtype=float) + 0.5
    return np.flatnonzero(np.asarray(score) > threshold).astype(float) + 0.5


def in_window(times: np.ndarray, start: float, end: float):
    t = np.asarray(times, float)
    return t[(t >= start) & (t < end)]


def vector_strength(times: np.ndarray, freq: float):
    t = np.asarray(times, float)
    if len(t) == 0:
        return 0.0
    return float(abs(np.mean(np.exp(2j * math.pi * freq * t))))


def event_rate(times: np.ndarray, start: float, end: float):
    return float(100.0 * len(in_window(times, start, end)) / max(end - start, 1e-12))


def first_or_nan(times: np.ndarray):
    return float(times[0]) if len(times) else float('nan')


def centroid_or_nan(times: np.ndarray):
    return float(np.mean(times)) if len(times) else float('nan')


def safe_absdiff(a, b):
    return float(abs(a - b)) if np.isfinite(a) and np.isfinite(b) else float('nan')


def count_contrast(a: np.ndarray, b: np.ndarray):
    na, nb = len(a), len(b)
    return float((na - nb) / (na + nb + 1e-12))


def selftest():
    ais = ActiveAIS()
    x = np.zeros(300, float)
    x[40:90] = 1.0
    x[160:210] = 1.0
    v, sp = ais.run(x)
    k = ais.linear_kernel(120)
    assert np.isfinite(v).all() and np.isfinite(k).all()
    assert len(sp) > 0, 'active boundary produced no events in selftest'
    assert np.max(np.abs(k)) > 1e-6, 'linearization is degenerate'
    print(f'selftest ok: {len(sp)} active events, kernel peak {np.max(np.abs(k)):.4g}')


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
    ap.add_argument('--out', default='runs/ais_active/ais_active.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def main():
    a = parse_args()
    if a.selftest:
        selftest(); return
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    freqs = [float(x) for x in a.freqs.split(',') if x.strip()]
    ais_cfg = AISConfig()
    ais = ActiveAIS(ais_cfg)
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
        active_v = {}
        linear = {}
        for key, tr in x.items():
            vv, ss = ais.run(tr)
            active_v[key] = vv
            active[key] = ss
            linear[key] = linear_score(tr, kernel, ais_cfg.input_gain)

        # One threshold per control and body, fitted only to TOTAL active event
        # count across the battery.  It is not retuned per frequency/order.
        keys = ['task_T', 'task_D'] + [f'f{j}' for j in range(len(freqs))]
        masks = []
        target_events = 0
        scores_raw = []
        scores_lin = []
        for key in keys:
            if key.startswith('f'):
                mask = np.arange(len(x[key])) >= a.burn
                target_events += len(in_window(active[key], a.burn, len(x[key])))
            else:
                mask = np.ones(len(x[key]), bool)
                target_events += len(active[key])
            masks.append(mask)
            scores_raw.append(x[key])
            scores_lin.append(linear[key])
        thr_raw = threshold_for_total(scores_raw, masks, target_events)
        thr_lin = threshold_for_total(scores_lin, masks, target_events)
        mem = {k: frame_events(x[k], thr_raw) for k in keys}
        lin = {k: frame_events(linear[k], thr_lin) for k in keys}

        freq_rows = []
        for j, f in enumerate(freqs):
            key = f'f{j}'; end = len(x[key])
            enc = {}
            for name, ev in [('memoryless', mem[key]), ('linearized', lin[key]), ('active', active[key])]:
                z = in_window(ev, a.burn, end)
                enc[name] = dict(events=int(len(z)),
                                 rate_per_100_frames=event_rate(ev, a.burn, end),
                                 vector_strength=vector_strength(z, f))
            freq_rows.append(dict(freq=f, **enc))

        task = {}
        for name, evT, evD in [('memoryless', mem['task_T'], mem['task_D']),
                               ('linearized', lin['task_T'], lin['task_D']),
                               ('active', active['task_T'], active['task_D'])]:
            task[name] = dict(
                events_T=int(len(evT)), events_D=int(len(evD)),
                count_contrast=count_contrast(evT, evD),
                first_latency_delta=safe_absdiff(first_or_nan(evT), first_or_nan(evD)),
                centroid_delta=safe_absdiff(centroid_or_nan(evT), centroid_or_nan(evD)),
            )

        rows.append(dict(seed=seed, cells=int(m.body.sum()), scale=scale,
                         thresholds=dict(memoryless=thr_raw, linearized=thr_lin),
                         active_total_events=int(target_events),
                         frequency=freq_rows, task=task))
        hf = freq_rows[-2:]
        print(f'seed {seed:2d}: active events {target_events:3d}  '
              f'HF VS mem/lin/active '
              f'{np.mean([q["memoryless"]["vector_strength"] for q in hf]):.3f}/'
              f'{np.mean([q["linearized"]["vector_strength"] for q in hf]):.3f}/'
              f'{np.mean([q["active"]["vector_strength"] for q in hf]):.3f}', flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    summary_freq = []
    for j, f in enumerate(freqs):
        item = dict(freq=f)
        for name in ('memoryless', 'linearized', 'active'):
            vs = np.asarray([r['frequency'][j][name]['vector_strength'] for r in rows], float)
            rr = np.asarray([r['frequency'][j][name]['rate_per_100_frames'] for r in rows], float)
            item[name] = dict(vs_mean=float(vs.mean()), vs_median=float(np.median(vs)),
                              rate_mean=float(rr.mean()))
        av = np.asarray([r['frequency'][j]['active']['vector_strength'] for r in rows])
        lv = np.asarray([r['frequency'][j]['linearized']['vector_strength'] for r in rows])
        mv = np.asarray([r['frequency'][j]['memoryless']['vector_strength'] for r in rows])
        item['active_minus_linear_mean'] = float(np.mean(av-lv))
        item['active_minus_memoryless_mean'] = float(np.mean(av-mv))
        item['active_beats_linear_count'] = int(np.sum(av > lv))
        item['active_beats_memoryless_count'] = int(np.sum(av > mv))
        summary_freq.append(item)

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
        )

    # Pre-registered high-frequency score: mean of the two highest sweep freqs.
    hf_idx = list(range(max(0, len(freqs)-2), len(freqs)))
    per_body_hf = []
    for r in rows:
        d = {}
        for name in ('memoryless', 'linearized', 'active'):
            d[name] = float(np.mean([r['frequency'][j][name]['vector_strength'] for j in hf_idx]))
        per_body_hf.append(d)
    ah = np.asarray([q['active'] for q in per_body_hf])
    lh = np.asarray([q['linearized'] for q in per_body_hf])
    mh = np.asarray([q['memoryless'] for q in per_body_hf])
    high_frequency = dict(
        active_mean=float(ah.mean()), linearized_mean=float(lh.mean()), memoryless_mean=float(mh.mean()),
        active_minus_linear_mean=float(np.mean(ah-lh)),
        active_minus_memoryless_mean=float(np.mean(ah-mh)),
        active_beats_linear_count=int(np.sum(ah>lh)),
        active_beats_memoryless_count=int(np.sum(ah>mh)),
        bodies=len(rows),
    )

    payload = dict(experiment='ais_active_boundary_v01',
                   interpretation='simulation units; no biological calibration',
                   config=dict(ais=asdict(ais_cfg), lag=a.lag, task_steps=a.task_steps,
                               freqs=freqs, freq_steps=a.freq_steps, burn=a.burn),
                   linear_kernel=kernel.tolist(), summary=dict(frequency=summary_freq,
                   high_frequency=high_frequency, task=task_summary), rows=rows)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nAIS ACTIVE-BOUNDARY RECEIPT')
    print(f'  bodies {len(rows)}')
    print(f'  high-frequency VS memoryless {high_frequency["memoryless_mean"]:.4f}')
    print(f'  high-frequency VS linearized {high_frequency["linearized_mean"]:.4f}')
    print(f'  high-frequency VS active     {high_frequency["active_mean"]:.4f}')
    print(f'  active - linearized          {high_frequency["active_minus_linear_mean"]:+.4f} '
          f'({high_frequency["active_beats_linear_count"]}/{len(rows)} bodies)')
    print(f'  active - memoryless          {high_frequency["active_minus_memoryless_mean"]:+.4f} '
          f'({high_frequency["active_beats_memoryless_count"]}/{len(rows)} bodies)')
    print(f'  wrote {out}')


if __name__ == '__main__':
    main()
