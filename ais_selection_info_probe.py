"""Does nonlinear AIS-like eventization buy frequency selection per spike?

The full active boundary is compared with a memoryless detector and the same
membrane's measured small-signal linearization.  Unlike the v0.2 timing test,
controls are matched to the active arm's TOTAL event count over the whole
frequency battery, not separately at each frequency.  This leaves frequency
allocation free to differ while holding the event budget fixed.

Primary statistic: KL[p(f|spike) || Uniform] in bits/spike for Re(psi_soma).
See AIS_SELECTION_INFO_PREREG_V01.md.
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
from ais_active_probe_v02 import sign_test_two_sided
from ais_interface_phase_probe import (
    FEEDS, modulated_complex_trace, scalar_feed, normalize_feed,
    current_from_normalized, linear_score_from_current,
)


def in_window(times, start, end):
    t = np.asarray(times, float)
    return t[(t >= start) & (t < end)]


def global_topk(scores: list[np.ndarray], start: int, k: int) -> list[np.ndarray]:
    """Exactly k frame-center events across equal-length condition traces."""
    vals = []
    refs = []
    for j, score in enumerate(scores):
        s = np.asarray(score, float)
        ids = np.arange(start, len(s), dtype=int)
        vals.append(s[ids])
        refs.extend((j, int(i)) for i in ids)
    if not vals:
        return [np.asarray([], float) for _ in scores]
    pool = np.concatenate(vals)
    k = max(0, min(int(k), len(pool)))
    out = [[] for _ in scores]
    if k == 0:
        return [np.asarray([], float) for _ in scores]
    # Stable sort gives deterministic tie handling.
    order = np.argsort(-pool, kind="mergesort")[:k]
    # refs follows the same concatenation order as pool.
    for idx in order:
        j, frame = refs[int(idx)]
        out[j].append(float(frame) + 0.5)
    return [np.asarray(sorted(x), float) for x in out]


def entropy_bits(p):
    p = np.asarray(p, float)
    z = p[p > 0]
    return float(-np.sum(z * np.log2(z))) if len(z) else 0.0


def bits_per_spike(counts):
    c = np.asarray(counts, float)
    total = float(c.sum())
    if total <= 0:
        return float("nan")
    p = c / total
    return float(math.log2(len(c)) - entropy_bits(p))


def binary_event_mi(counts, frames_per_condition):
    """I(F;Y) with F uniform over conditions, Y=event/no-event per frame."""
    c = np.asarray(counts, float)
    nfreq = len(c)
    T = float(frames_per_condition)
    if T <= 0 or np.any(c < 0) or np.any(c > T):
        return float("nan")
    q = c / T  # P(Y=1|F=f), at most one control/HH spike per field frame here.
    py1 = float(np.mean(q)); py0 = 1.0 - py1
    mi = 0.0
    for qi in q:
        if qi > 0 and py1 > 0:
            mi += (1.0 / nfreq) * qi * math.log2(qi / py1)
        if qi < 1 and py0 > 0:
            mi += (1.0 / nfreq) * (1.0 - qi) * math.log2((1.0 - qi) / py0)
    return float(mi)


def encoder_metrics(event_lists, start, end):
    counts = [int(len(in_window(e, start, end))) for e in event_lists]
    return {
        "counts": counts,
        "total": int(sum(counts)),
        "allocation": (np.asarray(counts, float) / sum(counts)).tolist() if sum(counts) else [float("nan")] * len(counts),
        "bits_per_spike": bits_per_spike(counts),
        "binary_event_mi_bits_per_frame": binary_event_mi(counts, end - start),
    }


def safe_wilcoxon(d, alternative="greater"):
    x = np.asarray(d, float)
    x = x[np.isfinite(x)]
    if len(x) == 0 or np.allclose(x, 0):
        return {"n": int(len(x)), "stat": float("nan"), "p": float("nan")}
    q = wilcoxon(x, alternative=alternative, zero_method="wilcox")
    return {"n": int(len(x)), "stat": float(q.statistic), "p": float(q.pvalue)}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--functional-arbors", default="../FunctionalArbors")
    ap.add_argument("--seed-start", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=24)
    ap.add_argument("--freqs", default="0.00625,0.0125,0.025,0.05,0.0833333,0.125")
    ap.add_argument("--freq-steps", type=int, default=640)
    ap.add_argument("--burn", type=int, default=160)
    ap.add_argument("--min-total-events", type=int, default=12)
    ap.add_argument("--min-valid-bodies", type=int, default=12)
    ap.add_argument("--out", default="runs/ais_selection_info/ais_selection_info.json")
    ap.add_argument("--selftest", action="store_true")
    return ap.parse_args()


def selftest():
    scores = [np.arange(20, dtype=float), np.arange(20, dtype=float)[::-1]]
    ev = global_topk(scores, 5, 7)
    assert sum(len(x) for x in ev) == 7
    assert abs(bits_per_spike([10, 10, 10]) - 0.0) < 1e-12
    assert abs(bits_per_spike([30, 0, 0]) - math.log2(3)) < 1e-12
    assert binary_event_mi([10, 10, 10], 100) < 1e-12
    print("selftest ok")


def main():
    a = parse_args()
    if a.selftest:
        selftest(); return
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f"FunctionalArbors not found at {fa}")
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    freqs = [float(x) for x in a.freqs.split(",") if x.strip()]
    cfg = AISConfig()
    ais = ActiveAIS(cfg)
    kernel = ais.linear_kernel()
    rows = []

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get("ok"):
            continue
        m.mature = True
        complex_traces = [modulated_complex_trace(m, f, a.freq_steps, source=0) for f in freqs]
        feed_rec = {}

        for kind in FEEDS:
            raw = [scalar_feed(z, kind) for z in complex_traces]
            xlist, scale = normalize_feed(raw)
            currents = [current_from_normalized(x, kind, cfg) for x in xlist]
            linear_scores = [linear_score_from_current(cur, kernel) for cur in currents]

            active_events = []
            K = 0
            for cur in currents:
                _, sp = ais.run_current(cur, True)
                active_events.append(sp)
                K += len(in_window(sp, a.burn, a.freq_steps))

            memoryless_events = global_topk(currents, a.burn, K)
            linearized_events = global_topk(linear_scores, a.burn, K)
            if sum(len(x) for x in memoryless_events) != K or sum(len(x) for x in linearized_events) != K:
                raise RuntimeError("global total-rate match failed")

            active = encoder_metrics(active_events, a.burn, a.freq_steps)
            memoryless = encoder_metrics(memoryless_events, a.burn, a.freq_steps)
            linearized = encoder_metrics(linearized_events, a.burn, a.freq_steps)
            if active["total"] != memoryless["total"] or active["total"] != linearized["total"]:
                raise RuntimeError("matched total counts disagree")

            feed_rec[kind] = {
                "scale": float(scale),
                "valid": bool(K >= a.min_total_events),
                "active": active,
                "memoryless": memoryless,
                "linearized": linearized,
            }

        rows.append({"seed": int(seed), "cells": int(m.body.sum()), "feeds": feed_rec})
        q = feed_rec["real"]
        if q["valid"]:
            d = q["active"]["bits_per_spike"] - q["linearized"]["bits_per_spike"]
            print(f"seed {seed:2d}: Re K={q['active']['total']:3d} dI_spike={d:+.4f} bits", flush=True)
        else:
            print(f"seed {seed:2d}: Re underexposed K={q['active']['total']}", flush=True)

    if not rows:
        raise SystemExit("No valid bodies")

    summary = {}
    for kind in FEEDS:
        valid = [r for r in rows if r["feeds"][kind]["valid"]]
        dlin = np.asarray([
            r["feeds"][kind]["active"]["bits_per_spike"] - r["feeds"][kind]["linearized"]["bits_per_spike"]
            for r in valid
        ], float)
        dmem = np.asarray([
            r["feeds"][kind]["active"]["bits_per_spike"] - r["feeds"][kind]["memoryless"]["bits_per_spike"]
            for r in valid
        ], float)
        dmi_lin = np.asarray([
            r["feeds"][kind]["active"]["binary_event_mi_bits_per_frame"] - r["feeds"][kind]["linearized"]["binary_event_mi_bits_per_frame"]
            for r in valid
        ], float)
        wl = (int(np.sum(dlin > 0)), int(np.sum(dlin < 0)))
        wm = (int(np.sum(dmem > 0)), int(np.sum(dmem < 0)))
        enc_mean = {}
        for enc in ("active", "linearized", "memoryless"):
            enc_mean[enc] = {
                "bits_per_spike_mean": float(np.mean([r["feeds"][kind][enc]["bits_per_spike"] for r in valid])) if valid else float("nan"),
                "binary_event_mi_mean": float(np.mean([r["feeds"][kind][enc]["binary_event_mi_bits_per_frame"] for r in valid])) if valid else float("nan"),
                "total_events_mean": float(np.mean([r["feeds"][kind][enc]["total"] for r in valid])) if valid else float("nan"),
            }
        summary[kind] = {
            "valid_bodies": len(valid),
            "encoder_means": enc_mean,
            "active_minus_linear_bits_per_spike_mean": float(np.mean(dlin)) if len(dlin) else float("nan"),
            "active_minus_linear_bits_per_spike_median": float(np.median(dlin)) if len(dlin) else float("nan"),
            "active_beats_linear": wl[0], "active_loses_linear": wl[1],
            "active_vs_linear_sign_p_two_sided": sign_test_two_sided(*wl),
            "active_gt_linear_wilcoxon": safe_wilcoxon(dlin, "greater"),
            "active_minus_memoryless_bits_per_spike_mean": float(np.mean(dmem)) if len(dmem) else float("nan"),
            "active_minus_memoryless_bits_per_spike_median": float(np.median(dmem)) if len(dmem) else float("nan"),
            "active_beats_memoryless": wm[0], "active_loses_memoryless": wm[1],
            "active_vs_memoryless_sign_p_two_sided": sign_test_two_sided(*wm),
            "active_gt_memoryless_wilcoxon": safe_wilcoxon(dmem, "greater"),
            "active_minus_linear_binary_mi_mean": float(np.mean(dmi_lin)) if len(dmi_lin) else float("nan"),
        }

    primary = summary["real"]
    earns = bool(
        primary["valid_bodies"] >= a.min_valid_bodies
        and primary["active_minus_linear_bits_per_spike_median"] > 0
        and np.isfinite(primary["active_gt_linear_wilcoxon"]["p"])
        and primary["active_gt_linear_wilcoxon"]["p"] < 0.05
    )

    payload = {
        "experiment": "ais_selection_information_v01",
        "registered": {
            "primary_feed": "real",
            "min_total_events": a.min_total_events,
            "min_valid_bodies": a.min_valid_bodies,
            "metric": "KL[p(f|spike)||Uniform(6)] bits/spike",
            "success_rule": ">=12 valid Re bodies, median active-linear >0, one-sided Wilcoxon active>linear p<.05",
        },
        "config": {"freqs": freqs, "freq_steps": a.freq_steps, "burn": a.burn, "ais": cfg.__dict__},
        "summary": {"feeds": summary, "active_earns_frequency_selection_role": earns},
        "rows": rows,
    }
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("\nAIS SELECTION-INFORMATION RECEIPT")
    for kind in FEEDS:
        q = summary[kind]
        print(f" {kind:9s}: valid {q['valid_bodies']:2d}  "
              f"I active/linear {q['encoder_means']['active']['bits_per_spike_mean']:.4f}/"
              f"{q['encoder_means']['linearized']['bits_per_spike_mean']:.4f}  "
              f"d={q['active_minus_linear_bits_per_spike_mean']:+.4f}  "
              f"w/l {q['active_beats_linear']}/{q['active_loses_linear']}  "
              f"p_gt={q['active_gt_linear_wilcoxon']['p']:.5g}")
    print(" active earns frequency-selection role:", earns)
    print(" wrote", out)


if __name__ == "__main__":
    main()
