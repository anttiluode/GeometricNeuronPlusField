"""Does n-gate kinetic speed set a refractory floor or mainly excitability?

This is a follow-up to AIS_KINETICS_V01.  Upstream bodies, power interface,
normalization, HH parameters, gain, frequencies, burn-in and the n_scale values
are frozen.  The only new readout is the interspike-interval distribution.

Primary test: at the already dominant/exposed f=0.025 cycles/frame, compare the
minimum ISI for n_slow (0.5x kinetics) and n_fast (2x kinetics) on bodies where
both emit at least two spikes after burn-in.

See AIS_N_ISI_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from ais_active_probe import AISConfig, task_trace, modulated_trace, normalize_traces
from ais_kinetics_probe import KineticAIS


VARIANTS = {"n_slow": 0.5, "full": 1.0, "n_fast": 2.0}


def in_window(times, start, end):
    t = np.asarray(times, float)
    return t[(t >= start) & (t < end)]


def isi_metrics(times):
    t = np.asarray(times, float)
    d = np.diff(t)
    if len(d) == 0:
        return {
            "events": int(len(t)), "isi_count": 0,
            "min_isi": float("nan"), "median_isi": float("nan"),
            "p10_isi": float("nan"), "mean_isi": float("nan"),
        }
    return {
        "events": int(len(t)), "isi_count": int(len(d)),
        "min_isi": float(np.min(d)),
        "median_isi": float(np.median(d)),
        "p10_isi": float(np.quantile(d, 0.10)),
        "mean_isi": float(np.mean(d)),
    }


def sign_two_sided(w, l):
    n = int(w + l)
    if n == 0:
        return float("nan")
    k = min(int(w), int(l))
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2 * tail))


def safe_wilcoxon(arr, alternative):
    x = np.asarray(arr, float)
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
    ap.add_argument("--lag", type=int, default=20)
    ap.add_argument("--task-steps", type=int, default=180)
    ap.add_argument("--freqs", default="0.00625,0.0125,0.025,0.05,0.0833333,0.125")
    ap.add_argument("--primary-freq", type=float, default=0.025)
    ap.add_argument("--freq-steps", type=int, default=640)
    ap.add_argument("--burn", type=int, default=160)
    ap.add_argument("--out", default="runs/ais_n_isi/ais_n_isi.json")
    return ap.parse_args()


def main():
    a = parse_args()
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f"FunctionalArbors not found at {fa}")
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    freqs = [float(x) for x in a.freqs.split(",") if x.strip()]
    j_primary = int(np.argmin(np.abs(np.asarray(freqs) - a.primary_freq)))
    if abs(freqs[j_primary] - a.primary_freq) > 1e-9:
        raise SystemExit("primary frequency must be one of --freqs")

    ais = KineticAIS(AISConfig())
    rows = []

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get("ok"):
            continue
        m.mature = True

        # Keep the exact old normalization battery, including the two task traces.
        raw = {
            "task_T": task_trace(m, a.lag, True, a.task_steps),
            "task_D": task_trace(m, a.lag, False, a.task_steps),
        }
        for j, f in enumerate(freqs):
            raw[f"f{j}"] = modulated_trace(m, f, a.freq_steps, source=0)
        x, scale = normalize_traces(raw)

        variants = {}
        for name, n_scale in VARIANTS.items():
            frows = []
            all_isi = []
            total_events = 0
            for j, f in enumerate(freqs):
                _, sp = ais.run_scaled(x[f"f{j}"], h_scale=1.0, n_scale=n_scale)
                z = in_window(sp, a.burn, a.freq_steps)
                met = isi_metrics(z)
                met["freq"] = float(f)
                frows.append(met)
                total_events += len(z)
                if len(z) >= 2:
                    all_isi.extend(np.diff(z).tolist())
            all_isi = np.asarray(all_isi, float)
            variants[name] = {
                "n_scale": float(n_scale),
                "total_events": int(total_events),
                "frequency": frows,
                "pooled_within_condition_isi_count": int(len(all_isi)),
                "pooled_min_isi": float(np.min(all_isi)) if len(all_isi) else float("nan"),
                "pooled_median_isi": float(np.median(all_isi)) if len(all_isi) else float("nan"),
            }

        rows.append({
            "seed": int(seed), "cells": int(m.body.sum()), "scale": float(scale),
            "variants": variants,
        })

        ps = variants["n_slow"]["frequency"][j_primary]["min_isi"]
        pf = variants["n_fast"]["frequency"][j_primary]["min_isi"]
        print(f"seed {seed:2d}: primary min ISI slow/full/fast "
              f"{ps!s:>8} / {variants['full']['frequency'][j_primary]['min_isi']!s:>8} / {pf!s:>8}",
              flush=True)

    if not rows:
        raise SystemExit("No valid bodies")

    # Primary paired slow-vs-fast test at f=0.025.
    log_ratio_fast_slow = []
    ordering = 0
    valid_ordering = 0
    primary_rows = []
    for r in rows:
        s = r["variants"]["n_slow"]["frequency"][j_primary]
        u = r["variants"]["full"]["frequency"][j_primary]
        f = r["variants"]["n_fast"]["frequency"][j_primary]
        if np.isfinite(s["min_isi"]) and np.isfinite(f["min_isi"]) and s["min_isi"] > 0 and f["min_isi"] > 0:
            lr = float(np.log2(f["min_isi"] / s["min_isi"]))
            log_ratio_fast_slow.append(lr)
            primary_rows.append({
                "seed": r["seed"], "slow_min_isi": s["min_isi"],
                "full_min_isi": u["min_isi"], "fast_min_isi": f["min_isi"],
                "log2_fast_over_slow": lr,
            })
            if np.isfinite(u["min_isi"]):
                valid_ordering += 1
                if s["min_isi"] > u["min_isi"] > f["min_isi"]:
                    ordering += 1

    arr = np.asarray(log_ratio_fast_slow, float)
    neg = int(np.sum(arr < 0))
    pos = int(np.sum(arr > 0))
    primary = {
        "freq": float(a.primary_freq),
        "valid_slow_fast_bodies": int(len(arr)),
        "mean_log2_fast_over_slow_min_isi": float(np.mean(arr)) if len(arr) else float("nan"),
        "median_log2_fast_over_slow_min_isi": float(np.median(arr)) if len(arr) else float("nan"),
        "fast_shorter": neg, "fast_longer": pos,
        "sign_p_two_sided": sign_two_sided(neg, pos),
        "wilcoxon_fast_shorter": safe_wilcoxon(arr, "less"),
        "strict_slow_gt_full_gt_fast": int(ordering),
        "strict_ordering_valid": int(valid_ordering),
        "pairs": primary_rows,
    }

    summaries = {}
    for name in VARIANTS:
        total = np.asarray([r["variants"][name]["total_events"] for r in rows], float)
        pmin = np.asarray([r["variants"][name]["frequency"][j_primary]["min_isi"] for r in rows], float)
        pmed = np.asarray([r["variants"][name]["frequency"][j_primary]["median_isi"] for r in rows], float)
        summaries[name] = {
            "emitting_bodies_anywhere": int(np.sum(total > 0)),
            "mean_total_events": float(np.mean(total)),
            "primary_min_isi_valid": int(np.isfinite(pmin).sum()),
            "primary_min_isi_mean": float(np.nanmean(pmin)) if np.isfinite(pmin).any() else float("nan"),
            "primary_median_isi_mean": float(np.nanmean(pmed)) if np.isfinite(pmed).any() else float("nan"),
        }

    payload = {
        "experiment": "ais_n_kinetics_isi_v01",
        "frozen_from": "ais_kinetics_v01",
        "registered": {
            "primary_freq": a.primary_freq,
            "n_scales": VARIANTS,
            "prediction_refractory": "n_slow minISI > full > n_fast; primary log2(fast/slow) < 0",
            "interpretation": {
                "supported": "minimum ISI changes in the predicted kinetic direction",
                "not_supported": "minimum ISI stable while event counts change strongly -> excitability/threshold effect",
            },
        },
        "config": {
            "freqs": freqs, "freq_steps": a.freq_steps, "burn": a.burn,
            "lag": a.lag, "task_steps": a.task_steps,
        },
        "summary": {"variants": summaries, "primary": primary},
        "rows": rows,
    }
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("\nAIS n-KINETICS ISI RECEIPT")
    print(f" valid slow/fast primary pairs {primary['valid_slow_fast_bodies']}")
    print(f" mean log2(fast/slow minISI) {primary['mean_log2_fast_over_slow_min_isi']:+.4f}")
    print(f" fast shorter/longer {primary['fast_shorter']}/{primary['fast_longer']}  "
          f"sign p={primary['sign_p_two_sided']:.5g}")
    print(f" Wilcoxon fast-shorter p={primary['wilcoxon_fast_shorter']['p']:.5g}")
    print(f" strict slow>full>fast {primary['strict_slow_gt_full_gt_fast']}/{primary['strict_ordering_valid']}")
    print(" wrote", out)


if __name__ == "__main__":
    main()
