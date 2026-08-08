"""Mechanism probe: can graph geometry predict the measured modal computation?

After graph-mode discovery/confirmation, the next question is WHY a higher spectral
band sees temporal order.  This script tests a deliberately simple reduction:

    q_n'' + damping q_n' + (restoring + stiffness*K*lambda_n) q_n
        = phi_n(A) s_A(t) + phi_n(B) s_B(t)

No fitted parameters are introduced.  All coefficients come from the existing
FunctionalArbor configuration and each frozen body's graph spectrum.

If the reduced scalar oscillators predict the measured modal A/B contrast, then
"geometry + field" has a concrete meaning: anatomy defines modal coordinates,
eigenvalue defines their time scale, and terminal geometry defines source coupling.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from graph_mode_probe import graph_laplacian_modes, trace_order, contrast_from_traces


def pulse_scalar(cfg, q):
    if not (0 <= q < cfg.pulse_frames):
        return 0.0j
    env = math.sin(math.pi * (q + 1) / (cfg.pulse_frames + 1)) ** 2
    phase = np.exp(1j * cfg.carrier_omega * q)
    return complex(cfg.source_amp * env * phase)


def reduced_modal_trace(cfg, evals, phi_a, phi_b, lag, target, steps):
    """Vectorized isolated-body linear modal predictor, no fitted parameters."""
    evals = np.asarray(evals, float)
    phi_a = np.asarray(phi_a, float)
    phi_b = np.asarray(phi_b, float)
    q = np.zeros(len(evals), np.complex128)
    v = np.zeros_like(q)
    out = np.zeros((steps, len(evals)), np.float64)

    for t in range(steps):
        if target:
            sa = pulse_scalar(cfg, t)
            sb = pulse_scalar(cfg, t - lag)
        else:
            sb = pulse_scalar(cfg, t)
            sa = pulse_scalar(cfg, t - lag)
        source = phi_a * sa + phi_b * sb
        v += cfg.dt * (
            -cfg.stiffness * cfg.k_arbor * evals * q
            - cfg.damping * v
            - cfg.restoring * q
            + source
        )
        q += cfg.dt * v
        out[t] = np.abs(q) ** 2
    return out


def corr(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 3 or x.std() < 1e-12 or y.std() < 1e-12:
        return float('nan')
    return float(np.corrcoef(x, y)[0, 1])


def rank_overlap(a, b, k=10):
    k = min(k, len(a), len(b))
    ia = set(np.argsort(np.abs(a))[-k:])
    ib = set(np.argsort(np.abs(b))[-k:])
    return len(ia & ib), k


def path_sign_vector(m, coords):
    pa = set(m.path(0) or [])
    pb = set(m.path(1) or [])
    x = np.zeros(len(coords), np.float64)
    for i, p in enumerate(coords):
        if p in pa and p not in pb:
            x[i] = 1.0
        elif p in pb and p not in pa:
            x[i] = -1.0
    n = np.linalg.norm(x)
    return x / n if n > 0 else x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='../FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=0)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--probe-steps', type=int, default=150)
    ap.add_argument('--outdir', default='runs/modal_mechanism')
    a = ap.parse_args()

    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []

    pooled_actual = []
    pooled_pred = []
    pooled_lambda = []
    pooled_source_diff = []
    pooled_source_interaction = []
    pooled_path_alignment = []

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            print(f'seed {seed}: bootstrap FAILED')
            continue
        m.mature = True
        coords, evals, evecs = graph_laplacian_modes(m.body)
        nmodes = len(coords)

        # Actual full spatial simulation, projected onto every body graph mode.
        tg = trace_order(m, coords, evecs, nmodes, int(a.lag), True, int(a.probe_steps))
        ds = trace_order(m, coords, evecs, nmodes, int(a.lag), False, int(a.probe_steps))
        actual = np.asarray(contrast_from_traces(tg['coherent'], ds['coherent']), float)

        # Geometry-derived terminal coupling.
        ta = tuple(m.source_terminal(0))
        tb = tuple(m.source_terminal(1))
        idx = {p: i for i, p in enumerate(coords)}
        ia, ib = idx[ta], idx[tb]
        phi_a = evecs[ia, :]
        phi_b = evecs[ib, :]

        # Parameter-free reduced modal prediction.
        pt = reduced_modal_trace(m.cfg, evals, phi_a, phi_b, int(a.lag), True, int(a.probe_steps))
        pd = reduced_modal_trace(m.cfg, evals, phi_a, phi_b, int(a.lag), False, int(a.probe_steps))
        predicted = np.asarray(contrast_from_traces(pt, pd), float)

        source_diff = np.abs(phi_a - phi_b)
        source_interaction = np.abs(phi_a * phi_b)
        path_sign = path_sign_vector(m, coords)
        path_alignment = np.abs(evecs.T @ path_sign)

        # Exclude the constant mode from correlations: both models predict it is
        # essentially order-blind, so including it can inflate agreement trivially.
        sl = slice(1, None)
        signed_r = corr(actual[sl], predicted[sl])
        abs_r = corr(np.abs(actual[sl]), np.abs(predicted[sl]))
        lambda_r = corr(np.asarray(evals[sl]), np.abs(actual[sl]))
        diff_r = corr(source_diff[sl], np.abs(actual[sl]))
        interaction_r = corr(source_interaction[sl], np.abs(actual[sl]))
        path_r = corr(path_alignment[sl], np.abs(actual[sl]))
        overlap, k = rank_overlap(actual[sl], predicted[sl], 10)

        # rank_overlap indices are relative after dropping mode 0; count is what matters.
        row = dict(
            seed=int(seed),
            boot=boot,
            terminal_A=list(ta),
            terminal_B=list(tb),
            measured_mode0_absC=float(abs(actual[0])),
            predicted_mode0_absC=float(abs(predicted[0])),
            signed_predicted_vs_measured_corr=signed_r,
            abs_predicted_vs_measured_corr=abs_r,
            absC_vs_lambda_corr=lambda_r,
            absC_vs_source_difference_corr=diff_r,
            absC_vs_source_interaction_corr=interaction_r,
            absC_vs_path_alignment_corr=path_r,
            top10_overlap=int(overlap),
            top10_k=int(k),
            modes=[dict(
                index=int(i),
                eigenvalue=float(evals[i]),
                measured_contrast=float(actual[i]),
                predicted_contrast=float(predicted[i]),
                source_A_coupling=float(phi_a[i]),
                source_B_coupling=float(phi_b[i]),
                source_difference=float(source_diff[i]),
                source_interaction=float(source_interaction[i]),
                path_alignment=float(path_alignment[i]),
            ) for i in range(nmodes)],
        )
        rows.append(row)

        pooled_actual.extend(actual[1:].tolist())
        pooled_pred.extend(predicted[1:].tolist())
        pooled_lambda.extend(evals[1:].tolist())
        pooled_source_diff.extend(source_diff[1:].tolist())
        pooled_source_interaction.extend(source_interaction[1:].tolist())
        pooled_path_alignment.extend(path_alignment[1:].tolist())

        print(
            f'seed {seed:2d}  signed r={signed_r:+.3f}  abs r={abs_r:+.3f}  '
            f'top10 overlap={overlap}/{k}  source-diff r={diff_r:+.3f}'
        )

    if not rows:
        raise SystemExit('No completed bodies.')

    def finite_mean(key):
        x = np.asarray([r[key] for r in rows], float)
        return float(np.nanmean(x)), float(np.nanmedian(x))

    signed_mean, signed_median = finite_mean('signed_predicted_vs_measured_corr')
    abs_mean, abs_median = finite_mean('abs_predicted_vs_measured_corr')
    diff_mean, diff_median = finite_mean('absC_vs_source_difference_corr')
    interact_mean, interact_median = finite_mean('absC_vs_source_interaction_corr')
    path_mean, path_median = finite_mean('absC_vs_path_alignment_corr')

    summary = dict(
        seeds_completed=len(rows),
        signed_prediction_corr_mean=signed_mean,
        signed_prediction_corr_median=signed_median,
        abs_prediction_corr_mean=abs_mean,
        abs_prediction_corr_median=abs_median,
        top10_overlap_mean=float(np.mean([r['top10_overlap'] for r in rows])),
        source_difference_corr_median=diff_median,
        source_interaction_corr_median=interact_median,
        path_alignment_corr_median=path_median,
        pooled_signed_prediction_corr=corr(pooled_actual, pooled_pred),
        pooled_abs_prediction_corr=corr(np.abs(pooled_actual), np.abs(pooled_pred)),
        pooled_absC_vs_lambda_corr=corr(pooled_lambda, np.abs(pooled_actual)),
        pooled_absC_vs_source_difference_corr=corr(pooled_source_diff, np.abs(pooled_actual)),
        pooled_absC_vs_source_interaction_corr=corr(pooled_source_interaction, np.abs(pooled_actual)),
        pooled_absC_vs_path_alignment_corr=corr(pooled_path_alignment, np.abs(pooled_actual)),
    )

    payload = dict(
        experiment='modal_mechanism_probe_v01',
        model='parameter-free isolated-body graph oscillator predictor',
        equation='q_n_ddot + damping*q_n_dot + (restoring + stiffness*K*lambda_n)*q_n = phi_n(A)*s_A + phi_n(B)*s_B',
        seed_start=int(a.seed_start),
        seeds_requested=int(a.seeds),
        lag=int(a.lag),
        probe_steps=int(a.probe_steps),
        summary=summary,
        rows=rows,
    )
    with open(outdir / 'modal_mechanism_results.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2)

    # Compact scatter plot: every non-constant mode from every body.
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(pooled_pred, pooled_actual, alpha=.28, s=14)
    lim = max(np.max(np.abs(pooled_pred)), np.max(np.abs(pooled_actual)), 1e-6)
    ax.plot([-lim, lim], [-lim, lim], linestyle='--', linewidth=1)
    ax.set_xlabel('reduced graph-oscillator predicted contrast')
    ax.set_ylabel('measured full-field modal contrast')
    ax.set_title('Does graph geometry predict the modal computation?')
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(outdir / 'predicted_vs_measured.png', dpi=180)
    plt.close(fig)

    print('\nMODAL MECHANISM RECEIPT')
    for k, v in summary.items():
        print(f'  {k:42s} {v:+.5f}' if isinstance(v, float) else f'  {k:42s} {v}')
    print(f'  wrote {outdir / "modal_mechanism_results.json"}')


if __name__ == '__main__':
    main()
