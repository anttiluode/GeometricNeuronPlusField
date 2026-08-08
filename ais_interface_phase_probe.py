"""Final phase/interface test for the AIS-like boundary.

Upstream anatomy, carrier physics, HH parameters, gain, frequency battery and
rate-matched controls are frozen from AIS active v0.2.  The only experimental
factor is the scalar handed from the frozen soma field to the boundary:

    power      |psi_soma|^2
    magnitude  |psi_soma|
    real       Re psi_soma

The real feed is signed and is therefore injected as signed current.  Power and
magnitude remain nonnegative.  Every feed is normalized by its own body-level
95th percentile absolute scale across the registered frequency battery.  The
same clipping magnitude and input gain are then used for all feeds.

See AIS_FINAL_PHASE_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from ais_active_probe import AISConfig, ActiveAIS
from ais_active_probe_v02 import topk_events, ppc, sign_test_two_sided


FEEDS = ("power", "magnitude", "real")


def reset_fast(m):
    try:
        m.reset_fast(clear_traces=True)
    except TypeError:
        m.reset_fast(True)


def modulated_complex_trace(m, freq: float, steps: int, source: int = 0,
                            floor: float = 0.12, depth: float = 0.88) -> np.ndarray:
    """Return the actual complex soma field, not the historical power readout."""
    reset_fast(m)
    p = m.source_terminal(source)
    if p is None:
        return np.zeros(steps, np.complex128)
    out = np.zeros(steps, np.complex128)
    c = m.cfg
    for t in range(steps):
        cyc = 0.5 * (1.0 - math.cos(2.0 * math.pi * freq * t))
        env = floor + depth * cyc
        src = np.zeros_like(m.psi)
        src[p] = c.source_amp * env * np.exp(1j * c.carrier_omega * t)
        m.advance(src, False, True, "none")
        out[t] = complex(m.psi[m.soma])
    return out


def scalar_feed(z: np.ndarray, kind: str) -> np.ndarray:
    z = np.asarray(z)
    if kind == "power":
        return np.abs(z) ** 2
    if kind == "magnitude":
        return np.abs(z)
    if kind == "real":
        return np.real(z)
    raise ValueError(kind)


def normalize_feed(traces: list[np.ndarray], q: float = 0.95):
    vals = np.concatenate([np.abs(np.asarray(x, float)[np.isfinite(x)]) for x in traces])
    vals = vals[vals > 0]
    scale = float(np.quantile(vals, q)) if len(vals) else 1.0
    scale = max(scale, 1e-12)
    return [np.asarray(x, float) / scale for x in traces], scale


def current_from_normalized(x: np.ndarray, kind: str, cfg: AISConfig) -> np.ndarray:
    x = np.asarray(x, float)
    if kind == "real":
        y = np.clip(x, -cfg.clip_drive, cfg.clip_drive)
    else:
        y = np.clip(x, 0.0, cfg.clip_drive)
    return cfg.input_gain * y


def linear_score_from_current(current: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    return np.convolve(np.asarray(current, float), np.asarray(kernel, float), mode="full")[:len(current)]


def in_window(times: np.ndarray, start: float, end: float) -> np.ndarray:
    t = np.asarray(times, float)
    return t[(t >= start) & (t < end)]


def vector_strength(times: np.ndarray, freq: float) -> float:
    t = np.asarray(times, float)
    if not len(t):
        return float("nan")
    return float(abs(np.mean(np.exp(2j * math.pi * freq * t))))


def timing_metrics(events: np.ndarray, freq: float, carrier_freq: float,
                   start: float, end: float) -> dict:
    z = in_window(events, start, end)
    return {
        "events": int(len(z)),
        "ppc_envelope": ppc(z, freq),
        "vs_envelope": vector_strength(z, freq),
        "ppc_carrier": ppc(z, carrier_freq),
        "vs_carrier": vector_strength(z, carrier_freq),
    }


def safe_wilcoxon(d: np.ndarray) -> dict:
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if len(d) == 0 or np.allclose(d, 0):
        return {"n": int(len(d)), "stat": float("nan"), "p_two_sided": float("nan"),
                "p_active_greater": float("nan")}
    two = wilcoxon(d, alternative="two-sided", zero_method="wilcox")
    greater = wilcoxon(d, alternative="greater", zero_method="wilcox")
    return {"n": int(len(d)), "stat": float(two.statistic), "p_two_sided": float(two.pvalue),
            "p_active_greater": float(greater.pvalue)}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--functional-arbors", default="../FunctionalArbors")
    ap.add_argument("--seed-start", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=24)
    ap.add_argument("--freqs", default="0.00625,0.0125,0.025,0.05,0.0833333,0.125")
    ap.add_argument("--freq-steps", type=int, default=640)
    ap.add_argument("--burn", type=int, default=160)
    ap.add_argument("--min-events", type=int, default=4)
    ap.add_argument("--upper-min-freq", type=float, default=0.05)
    ap.add_argument("--min-valid-bodies", type=int, default=8)
    ap.add_argument("--out", default="runs/ais_final_phase/ais_final_phase.json")
    ap.add_argument("--selftest", action="store_true")
    return ap.parse_args()


def selftest():
    cfg = AISConfig()
    ais = ActiveAIS(cfg)
    kernel = ais.linear_kernel()
    t = np.arange(300.0)
    z = (0.6 + 0.4 * np.cos(2 * np.pi * 0.05 * t)) * np.exp(1j * 0.16 * t)
    for kind in FEEDS:
        raw = scalar_feed(z, kind)
        (x,), scale = normalize_feed([raw])
        cur = current_from_normalized(x, kind, cfg)
        lin = linear_score_from_current(cur, kernel)
        _, sp = ais.run_current(cur, True)
        mask = np.arange(len(x)) >= 50
        k = len(in_window(sp, 50, len(x)))
        ev = topk_events(cur, mask, k)
        assert len(ev) == k
        assert np.isfinite(cur).all() and np.isfinite(lin).all() and scale > 0
    print("selftest ok")


def main():
    a = parse_args()
    if a.selftest:
        selftest()
        return
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f"FunctionalArbors not found at {fa}")
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    freqs = [float(x) for x in a.freqs.split(",") if x.strip()]
    cfg = AISConfig()  # frozen
    ais = ActiveAIS(cfg)
    kernel = ais.linear_kernel()
    rows = []

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get("ok"):
            continue
        m.mature = True
        carrier_freq = float(m.cfg.carrier_omega / (2.0 * math.pi))

        complex_traces = [modulated_complex_trace(m, f, a.freq_steps, source=0) for f in freqs]
        feed_rows = {}

        for kind in FEEDS:
            raw = [scalar_feed(z, kind) for z in complex_traces]
            xlist, scale = normalize_feed(raw)
            freq_rows = []
            for f, x in zip(freqs, xlist):
                cur = current_from_normalized(x, kind, cfg)
                _, active = ais.run_current(cur, True)
                lin_score = linear_score_from_current(cur, kernel)

                az = in_window(active, a.burn, len(x))
                k = len(az)
                mask = np.arange(len(x)) >= a.burn
                memoryless = topk_events(cur, mask, k)
                linearized = topk_events(lin_score, mask, k)
                if len(memoryless) != k or len(linearized) != k:
                    raise RuntimeError("exact event-count match failed")

                freq_rows.append({
                    "freq": float(f),
                    "valid_timing": bool(k >= a.min_events),
                    "active": timing_metrics(active, f, carrier_freq, a.burn, len(x)),
                    "memoryless": timing_metrics(memoryless, f, carrier_freq, a.burn, len(x)),
                    "linearized": timing_metrics(linearized, f, carrier_freq, a.burn, len(x)),
                })
            feed_rows[kind] = {"scale": float(scale), "frequency": freq_rows}

        rows.append({
            "seed": int(seed), "cells": int(m.body.sum()),
            "carrier_freq_cycles_per_frame": carrier_freq,
            "feeds": feed_rows,
        })
        re_valid = [q for q in feed_rows["real"]["frequency"]
                    if q["freq"] >= a.upper_min_freq and q["valid_timing"]]
        delta = np.mean([q["active"]["ppc_envelope"] - q["linearized"]["ppc_envelope"]
                         for q in re_valid]) if re_valid else float("nan")
        print(f"seed {seed:2d}: Re valid upper {len(re_valid)} mean active-linear PPC {delta:+.3f}",
              flush=True)

    if not rows:
        raise SystemExit("No valid bodies")

    feed_summary = {}
    for kind in FEEDS:
        per_freq = []
        body_deltas = []
        for j, f in enumerate(freqs):
            qs = [r["feeds"][kind]["frequency"][j] for r in rows]
            valid = [q for q in qs if q["valid_timing"]]
            rec = {
                "freq": float(f),
                "bodies": len(qs),
                "valid_pairs": len(valid),
                "mean_active_events": float(np.mean([q["active"]["events"] for q in qs])),
            }
            for name in ("memoryless", "linearized", "active"):
                rec[name] = {
                    "ppc_envelope_mean_valid": float(np.nanmean([q[name]["ppc_envelope"] for q in valid]))
                        if valid else float("nan"),
                    "ppc_carrier_mean_valid": float(np.nanmean([q[name]["ppc_carrier"] for q in valid]))
                        if valid else float("nan"),
                }
            if valid:
                d = np.asarray([q["active"]["ppc_envelope"] - q["linearized"]["ppc_envelope"]
                                for q in valid], float)
                rec["active_minus_linear_ppc_mean"] = float(np.nanmean(d))
            per_freq.append(rec)

        for r in rows:
            valid = [q for q in r["feeds"][kind]["frequency"]
                     if q["freq"] >= a.upper_min_freq and q["valid_timing"]]
            if valid:
                body_deltas.append(float(np.mean(
                    [q["active"]["ppc_envelope"] - q["linearized"]["ppc_envelope"] for q in valid]
                )))

        d = np.asarray(body_deltas, float)
        wins = int(np.sum(d > 0)); losses = int(np.sum(d < 0))
        feed_summary[kind] = {
            "frequency": per_freq,
            "upper_body_level": {
                "valid_bodies": int(len(d)),
                "mean_active_minus_linear_ppc": float(np.mean(d)) if len(d) else float("nan"),
                "median_active_minus_linear_ppc": float(np.median(d)) if len(d) else float("nan"),
                "active_beats_linear": wins,
                "active_loses_linear": losses,
                "sign_p_two_sided": sign_test_two_sided(wins, losses),
                "wilcoxon": safe_wilcoxon(d),
            },
        }

    primary = feed_summary["real"]["upper_body_level"]
    phase_earns_role = bool(
        primary["valid_bodies"] >= a.min_valid_bodies
        and primary["median_active_minus_linear_ppc"] > 0
        and np.isfinite(primary["wilcoxon"]["p_active_greater"])
        and primary["wilcoxon"]["p_active_greater"] < 0.05
    )

    payload = {
        "experiment": "ais_final_phase_interface_v01",
        "frozen_from": "ais_active_boundary_v02",
        "registered": {
            "feeds": list(FEEDS),
            "min_events": a.min_events,
            "upper_min_freq": a.upper_min_freq,
            "min_valid_bodies": a.min_valid_bodies,
            "primary": "Re(psi) active-vs-own-linearization, body-level mean upper-band envelope PPC",
            "success_rule": ">=8 valid bodies, median delta > 0, one-sided Wilcoxon(active>linear) p<0.05",
        },
        "config": {
            "freqs": freqs, "freq_steps": a.freq_steps, "burn": a.burn,
            "ais": cfg.__dict__,
        },
        "summary": {
            "feeds": feed_summary,
            "phase_earns_active_eventization_role": phase_earns_role,
        },
        "rows": rows,
    }
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("\nAIS FINAL PHASE / INTERFACE RECEIPT")
    for kind in FEEDS:
        q = feed_summary[kind]["upper_body_level"]
        print(f" {kind:9s}: valid bodies {q['valid_bodies']:2d}  "
              f"mean dPPC {q['mean_active_minus_linear_ppc']:+.4f}  "
              f"median {q['median_active_minus_linear_ppc']:+.4f}  "
              f"w/l {q['active_beats_linear']}/{q['active_loses_linear']}  "
              f"Wilcoxon greater p={q['wilcoxon']['p_active_greater']:.5g}")
    print(" phase earns active eventization role:", phase_earns_role)
    print(" wrote", out)


if __name__ == "__main__":
    main()
