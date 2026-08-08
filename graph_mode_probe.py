"""GeometricNeuronPlusField — graph-mode microscope.

This experiment deliberately freezes morphology and changes only the observation
basis.  Instead of hand-selecting point/cross/ring soma taps, it diagonalizes the
unweighted graph Laplacian of the FunctionalArbor body and reads the SAME complex
field through the body's own spatial modes.

Two questions are kept separate:

1) Where does A/B temporal-order selectivity live in graph-mode space?
2) Can a cycle-level modal observable settle while the instantaneous field keeps
   moving under repeated drive?

The script imports the existing v0.9 FunctionalArbor implementation from a sibling
checkout.  It does not copy or alter its wave physics, growth physics, or readout.
No growth, credit, or learning occurs here.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


def graph_laplacian_modes(body: np.ndarray):
    """Return (coords, eigenvalues, eigenvectors) for the 4-neighbour body graph.

    eigenvectors[:, n] is phi_n evaluated in `coords` order.  np.linalg.eigh gives
    an orthonormal real basis.  For a connected body, mode 0 is the constant mode.
    """
    b = np.asarray(body, bool)
    coords = [tuple(map(int, p)) for p in np.argwhere(b)]
    if not coords:
        raise ValueError("empty body")
    idx = {p: i for i, p in enumerate(coords)}
    n = len(coords)
    A = np.zeros((n, n), np.float64)
    for i, (y, x) in enumerate(coords):
        for q in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
            j = idx.get(q)
            if j is not None:
                A[i, j] = 1.0
    A = np.maximum(A, A.T)
    degree = A.sum(axis=1)
    L = np.diag(degree) - A
    evals, evecs = np.linalg.eigh(L)
    order = np.argsort(evals)
    return coords, evals[order], evecs[:, order]


def _reset_fast(m):
    # FunctionalArbors versions have used both positional and keyword spelling.
    try:
        m.reset_fast(clear_traces=True)
    except TypeError:
        m.reset_fast(True)


def _add_sources(*parts):
    src = 0.0
    for p in parts:
        if isinstance(p, (float, int, np.floating)):
            if float(p) == 0.0:
                continue
        if isinstance(src, (float, int, np.floating)):
            src = p
        else:
            src = src + p
    return src


def project_field(psi, coords, evecs, soma, nmodes):
    """Read one field state in point, coherent-mode and incoherent-mode bases."""
    z = np.asarray([psi[p] for p in coords], np.complex128)
    V = evecs[:, :nmodes]

    # Graph-mode coherent projection: signed/phase-aware spatial filter.
    a = V.T @ z
    coherent = np.abs(a) ** 2

    # Matched phase-destroyed control.  Since each phi_n has L2 norm 1,
    # phi_n^2 is a non-negative unit-mass spatial aperture.
    incoherent = (V * V).T @ (np.abs(z) ** 2)

    point = float(np.abs(psi[soma]) ** 2)
    total_power = float(np.vdot(z, z).real)
    return point, coherent.real, incoherent.real, total_power, z


def trace_order(m, coords, evecs, nmodes, lag, target, steps):
    """One A->B or B->A probe; every readout watches the identical psi(t)."""
    _reset_fast(m)
    first, second = (0, 1) if target else (1, 0)
    point = np.zeros(steps, np.float64)
    coh = np.zeros((steps, nmodes), np.float64)
    incoh = np.zeros((steps, nmodes), np.float64)
    total = np.zeros(steps, np.float64)

    for t in range(steps):
        a = m.pulse_source(first, t, False)
        b = m.pulse_source(second, t - int(lag), False)
        src = _add_sources(a, b)
        m.advance(src, False, True, 'none')
        p, c, ic, tp, _ = project_field(m.psi, coords, evecs, m.soma, nmodes)
        point[t] = p
        coh[t] = c
        incoh[t] = ic
        total[t] = tp
    return dict(point=point, coherent=coh, incoherent=incoh, total_power=total)


def contrast_from_traces(tg, ds):
    """Same peak contrast used by the soma tap test, vectorized over modes."""
    a = np.max(tg, axis=0)
    b = np.max(ds, axis=0)
    return (a - b) / (a + b + 1e-12)


def sweep_landscape(m, coords, evals, evecs, nmodes, lags, steps):
    point_C = []
    coherent_C = []
    incoherent_C = []
    for lag in lags:
        tg = trace_order(m, coords, evecs, nmodes, lag, True, steps)
        ds = trace_order(m, coords, evecs, nmodes, lag, False, steps)
        point_C.append(float(contrast_from_traces(tg['point'][:, None], ds['point'][:, None])[0]))
        coherent_C.append(contrast_from_traces(tg['coherent'], ds['coherent']))
        incoherent_C.append(contrast_from_traces(tg['incoherent'], ds['incoherent']))

    point_C = np.asarray(point_C)
    coherent_C = np.asarray(coherent_C)
    incoherent_C = np.asarray(incoherent_C)
    return dict(
        lags=np.asarray(lags, int),
        point=point_C,
        coherent=coherent_C,
        incoherent=incoherent_C,
        eigenvalues=np.asarray(evals[:nmodes], float),
    )


def _scheduled_source(m, t, events):
    parts = []
    for which, t0 in events:
        p = m.pulse_source(which, t - t0, False)
        if isinstance(p, (float, int, np.floating)) and float(p) == 0.0:
            continue
        parts.append(p)
    return _add_sources(*parts) if parts else 0.0


def live_field_probe(m, coords, evecs, nmodes, lag, cycles, period, burn_cycles,
                     field_live_threshold, settled_cv_threshold, energy_threshold):
    """Repeated A->B drive: look for stable cycle observables of a moving field.

    "Settled" here is deliberately operational rather than metaphysical: after
    burn-in, the *cycle-mean modal power* has low coefficient of variation while
    normalized instantaneous field motion remains nonzero.
    """
    _reset_fast(m)
    total_steps = int(cycles * period)
    events = []
    for c in range(cycles):
        t0 = c * period
        events.append((0, t0))
        events.append((1, t0 + int(lag)))

    coh = np.zeros((total_steps, nmodes), np.float64)
    incoh = np.zeros((total_steps, nmodes), np.float64)
    total_power = np.zeros(total_steps, np.float64)
    point = np.zeros(total_steps, np.float64)
    field_motion = np.zeros(total_steps, np.float64)
    prev = None

    for t in range(total_steps):
        src = _scheduled_source(m, t, events)
        m.advance(src, False, True, 'none')
        p, c, ic, tp, z = project_field(m.psi, coords, evecs, m.soma, nmodes)
        point[t] = p
        coh[t] = c
        incoh[t] = ic
        total_power[t] = tp
        if prev is not None:
            field_motion[t] = np.linalg.norm(z - prev) / (np.linalg.norm(z) + 1e-12)
        prev = z

    cycle_coh = np.asarray([
        coh[c * period:(c + 1) * period].mean(axis=0) for c in range(cycles)
    ])
    cycle_incoh = np.asarray([
        incoh[c * period:(c + 1) * period].mean(axis=0) for c in range(cycles)
    ])
    cycle_point = np.asarray([
        point[c * period:(c + 1) * period].mean() for c in range(cycles)
    ])

    start_cycle = min(max(int(burn_cycles), 0), max(cycles - 1, 0))
    start = start_cycle * period
    post = cycle_coh[start_cycle:]
    post_ic = cycle_incoh[start_cycle:]
    post_point = cycle_point[start_cycle:]

    mean_coh = post.mean(axis=0)
    cv_coh = post.std(axis=0) / (np.abs(mean_coh) + 1e-12)
    mean_incoh = post_ic.mean(axis=0)
    cv_incoh = post_ic.std(axis=0) / (np.abs(mean_incoh) + 1e-12)
    point_cv = float(post_point.std() / (abs(post_point.mean()) + 1e-12))

    mean_total = float(total_power[start:].mean()) + 1e-12
    # Parseval makes this a principled modal energy fraction for the coherent basis.
    energy_fraction = mean_coh / mean_total

    # Instantaneous output motion is recorded too.  The primary "settled" test is
    # cycle-level because a settled measurement may summarize a still-oscillating field.
    dcoh = np.abs(np.diff(coh[start:], axis=0))
    inst_motion = np.median(dcoh, axis=0) / (np.median(np.abs(coh[start:]), axis=0) + 1e-12)
    live_motion = float(np.median(field_motion[max(start, 1):]))

    settled = ((live_motion >= field_live_threshold) &
               (cv_coh <= settled_cv_threshold) &
               (energy_fraction >= energy_threshold))

    return dict(
        field_motion_median=live_motion,
        field_live=bool(live_motion >= field_live_threshold),
        cycle_mean_coherent=mean_coh,
        cycle_cv_coherent=cv_coh,
        cycle_mean_incoherent=mean_incoh,
        cycle_cv_incoherent=cv_incoh,
        instantaneous_motion_coherent=inst_motion,
        modal_energy_fraction=energy_fraction,
        settled_energy_nontrivial=settled,
        point_cycle_cv=point_cv,
        thresholds=dict(field_live=float(field_live_threshold),
                        settled_cv=float(settled_cv_threshold),
                        energy_fraction=float(energy_threshold)),
    )


def _f(x):
    return float(x) if np.isfinite(x) else None


def seed_receipt(seed, boot, landscape, live, lags, target_lag, selectivity_threshold):
    L = list(lags).index(int(target_lag))
    out_modes = []
    nmodes = landscape['coherent'].shape[1]
    for n in range(nmodes):
        c_curve = landscape['coherent'][:, n]
        ic_curve = landscape['incoherent'][:, n]
        best = int(lags[int(np.argmax(np.abs(c_curve)))])
        c = float(c_curve[L])
        ic = float(ic_curve[L])
        settled = bool(live['settled_energy_nontrivial'][n])
        informative = abs(c) >= selectivity_threshold
        out_modes.append(dict(
            index=n,
            eigenvalue=float(landscape['eigenvalues'][n]),
            contrast=float(c),
            abs_contrast=float(abs(c)),
            incoherent_contrast=float(ic),
            coherence_gain_absC=float(abs(c) - abs(ic)),
            argmax_absC_lag=best,
            live_cycle_cv=float(live['cycle_cv_coherent'][n]),
            live_incoherent_cycle_cv=float(live['cycle_cv_incoherent'][n]),
            instantaneous_readout_motion=float(live['instantaneous_motion_coherent'][n]),
            energy_fraction=float(live['modal_energy_fraction'][n]),
            settled_energy_nontrivial=settled,
            settled_and_informative=bool(settled and informative),
        ))
    return dict(
        seed=int(seed),
        boot=boot,
        point_landscape={str(int(l)): float(v) for l, v in zip(lags, landscape['point'])},
        field_motion_median=float(live['field_motion_median']),
        field_live=bool(live['field_live']),
        point_cycle_cv=float(live['point_cycle_cv']),
        modes=out_modes,
    )


def summarize(rows, nmodes):
    summary = []
    for n in range(nmodes):
        modes = [r['modes'][n] for r in rows if len(r['modes']) > n]
        if not modes:
            continue
        vals = lambda key: np.asarray([m[key] for m in modes], float)
        summary.append(dict(
            index=n,
            eigenvalue_mean=float(vals('eigenvalue').mean()),
            absC_mean=float(vals('abs_contrast').mean()),
            absC_median=float(np.median(vals('abs_contrast'))),
            coherence_gain_mean=float(vals('coherence_gain_absC').mean()),
            coherence_gain_positive=int((vals('coherence_gain_absC') > 0).sum()),
            n=len(modes),
            energy_fraction_mean=float(vals('energy_fraction').mean()),
            live_cycle_cv_median=float(np.median(vals('live_cycle_cv'))),
            settled_informative_count=int(sum(bool(m['settled_and_informative']) for m in modes)),
        ))
    return summary


def plot_summary(summary, rows, outdir):
    import matplotlib.pyplot as plt

    idx = np.asarray([s['index'] for s in summary], int)
    absC = np.asarray([s['absC_mean'] for s in summary], float)
    gain = np.asarray([s['coherence_gain_mean'] for s in summary], float)
    energy = np.asarray([s['energy_fraction_mean'] for s in summary], float)
    cv = np.asarray([s['live_cycle_cv_median'] for s in summary], float)
    count = np.asarray([s['settled_informative_count'] for s in summary], int)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    axes[0].plot(idx, absC, marker='o', label='mean |C| coherent mode')
    axes[0].plot(idx, np.maximum(absC - gain, 0), marker='.', label='approx mean |C| incoherent')
    axes[0].set_ylabel('temporal selectivity')
    axes[0].legend()
    axes[0].grid(alpha=.25)

    axes[1].plot(idx, gain, marker='o', label='coherent - incoherent |C|')
    axes[1].axhline(0, linewidth=1)
    axes[1].plot(idx, energy, marker='.', label='modal energy fraction')
    axes[1].set_ylabel('gain / energy')
    axes[1].legend()
    axes[1].grid(alpha=.25)

    axes[2].plot(idx, cv, marker='o', label='median cycle CV')
    axes[2].plot(idx, count / max(len(rows), 1), marker='s', label='fraction settled + informative')
    axes[2].set_xlabel('graph Laplacian mode index')
    axes[2].set_ylabel('live-field readout')
    axes[2].legend()
    axes[2].grid(alpha=.25)

    fig.suptitle('GeometricNeuronPlusField — graph-mode microscope')
    fig.tight_layout()
    fig.savefig(outdir / 'graph_mode_summary.png', dpi=180)
    plt.close(fig)


def plot_mode_maps(body, coords, evals, evecs, outdir, count=8):
    import matplotlib.pyplot as plt

    count = min(int(count), evecs.shape[1])
    cols = 4
    rows = int(math.ceil(count / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(12, 3 * rows))
    axes = np.asarray(axes).reshape(-1)
    for n in range(count):
        img = np.full(body.shape, np.nan, np.float64)
        for value, p in zip(evecs[:, n], coords):
            img[p] = value
        axes[n].imshow(img)
        axes[n].set_title(f'mode {n}  lambda={evals[n]:.4g}')
        axes[n].axis('off')
    for ax in axes[count:]:
        ax.axis('off')
    fig.suptitle('First body graph modes (seed 0)')
    fig.tight_layout()
    fig.savefig(outdir / 'mode_maps_seed0.png', dpi=180)
    plt.close(fig)


def selftest():
    body = np.zeros((7, 7), np.uint8)
    body[3, 1:6] = 1
    body[2, 3] = 1
    body[1, 3] = 1
    coords, evals, evecs = graph_laplacian_modes(body)
    assert len(coords) == int(body.sum())
    assert evals[0] > -1e-10
    assert np.max(np.abs(evecs.T @ evecs - np.eye(len(coords)))) < 1e-10
    # Connected graph: first eigenvector is constant up to arbitrary sign.
    assert np.std(np.abs(evecs[:, 0])) < 1e-10

    rng = np.random.default_rng(123)
    z = rng.normal(size=len(coords)) + 1j * rng.normal(size=len(coords))
    a = evecs.T @ z
    assert abs(np.sum(np.abs(a) ** 2) - np.sum(np.abs(z) ** 2)) < 1e-9
    print('SELFTEST PASS: Laplacian basis, constant mode, orthogonality, Parseval')


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='../FunctionalArbors')
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--modes', type=int, default=24)
    ap.add_argument('--lags', default='0,4,8,12,16,20,24,28,32,40')
    ap.add_argument('--target-lag', type=int, default=20)
    ap.add_argument('--probe-steps', type=int, default=150)
    ap.add_argument('--cycles', type=int, default=12)
    ap.add_argument('--period', type=int, default=64)
    ap.add_argument('--burn-cycles', type=int, default=4)
    ap.add_argument('--field-live-threshold', type=float, default=0.01)
    ap.add_argument('--settled-cv-threshold', type=float, default=0.05)
    ap.add_argument('--energy-threshold', type=float, default=0.005)
    ap.add_argument('--selectivity-threshold', type=float, default=0.05)
    ap.add_argument('--outdir', default='runs/graph_modes')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def main():
    a = parse_args()
    if a.selftest:
        selftest()
        return

    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(
            f'FunctionalArbors not found at {fa}. Clone it beside this repo or pass '
            '--functional-arbors PATH.'
        )
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    lags = [int(x) for x in a.lags.split(',') if x.strip()]
    if a.target_lag not in lags:
        raise SystemExit('--target-lag must be present in --lags')

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    first_map = None

    for seed in range(a.seeds):
        cfg = V09Config(seed=seed)
        m = CausalEligibilityArbor(cfg)
        boot = m.bootstrap()
        if not boot.get('ok'):
            print(f'seed {seed}: bootstrap FAILED')
            continue
        m.mature = True

        coords, evals, evecs = graph_laplacian_modes(m.body)
        nmodes = min(int(a.modes), len(coords))
        if first_map is None:
            first_map = (m.body.copy(), coords, evals.copy(), evecs.copy())

        landscape = sweep_landscape(m, coords, evals, evecs, nmodes,
                                    lags, int(a.probe_steps))
        live = live_field_probe(
            m, coords, evecs, nmodes, int(a.target_lag), int(a.cycles), int(a.period),
            int(a.burn_cycles), float(a.field_live_threshold),
            float(a.settled_cv_threshold), float(a.energy_threshold)
        )
        row = seed_receipt(seed, boot, landscape, live, lags, a.target_lag,
                           float(a.selectivity_threshold))
        rows.append(row)

        target_i = lags.index(a.target_lag)
        mode_abs = np.abs(landscape['coherent'][target_i])
        best = int(np.argmax(mode_abs))
        settled_n = sum(mo['settled_and_informative'] for mo in row['modes'])
        print(
            f'seed {seed:2d}  cells={len(coords):2d}  '
            f'field_motion={row["field_motion_median"]:.4f}  '
            f'best_mode={best:2d} |C|={mode_abs[best]:.4f}  '
            f'settled+informative={settled_n}'
        )

    if not rows:
        raise SystemExit('No successful bootstrap bodies.')

    common_modes = min(len(r['modes']) for r in rows)
    summary = summarize(rows, common_modes)

    # Explicit constant-mode diagnostic.  Mode 0 should correspond to the smooth
    # spatial common mode on every connected body.
    c0 = np.asarray([r['modes'][0]['abs_contrast'] for r in rows], float)
    c0_ic = np.asarray([abs(r['modes'][0]['incoherent_contrast']) for r in rows], float)
    point = np.asarray([abs(r['point_landscape'][str(a.target_lag)]) for r in rows], float)
    live_motion = np.asarray([r['field_motion_median'] for r in rows], float)

    payload = dict(
        experiment='graph_mode_probe_v01',
        source_model='FunctionalArbors/v09_causal_eligibility frozen bootstrap body',
        seeds_requested=int(a.seeds),
        seeds_completed=len(rows),
        target_lag=int(a.target_lag),
        lags=lags,
        modes_requested=int(a.modes),
        live_protocol=dict(cycles=int(a.cycles), period=int(a.period), burn_cycles=int(a.burn_cycles)),
        thresholds=dict(field_live=float(a.field_live_threshold),
                        settled_cv=float(a.settled_cv_threshold),
                        energy_fraction=float(a.energy_threshold),
                        selectivity=float(a.selectivity_threshold)),
        diagnostics=dict(
            constant_mode_absC_mean=float(c0.mean()),
            constant_mode_incoherent_absC_mean=float(c0_ic.mean()),
            point_absC_mean=float(point.mean()),
            field_motion_median_across_seeds=float(np.median(live_motion)),
            field_live_seed_count=int(sum(r['field_live'] for r in rows)),
        ),
        summary=summary,
        rows=rows,
    )

    with open(outdir / 'graph_mode_results.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2)

    plot_summary(summary, rows, outdir)
    if first_map is not None:
        plot_mode_maps(*first_map, outdir=outdir, count=min(8, common_modes))

    print('\nGRAPH-MODE RECEIPT')
    print(f'  completed seeds            {len(rows)}')
    print(f'  mean |C| point @ lag       {point.mean():.5f}')
    print(f'  mean |C| constant mode     {c0.mean():.5f}')
    print(f'  mean |C| mode0 incoherent  {c0_ic.mean():.5f}')
    print(f'  median live field motion   {np.median(live_motion):.5f}')
    print(f'  live-field seeds           {sum(r["field_live"] for r in rows)}/{len(rows)}')
    ranked = sorted(summary, key=lambda s: (s['settled_informative_count'], s['absC_mean']), reverse=True)
    print('  top modes by settled+informative count:')
    for s in ranked[:5]:
        print(f'    mode {s["index"]:2d}: count={s["settled_informative_count"]:2d}/{len(rows)} '
              f'mean|C|={s["absC_mean"]:.4f} gain={s["coherence_gain_mean"]:+.4f} '
              f'energy={s["energy_fraction_mean"]:.4f} cv={s["live_cycle_cv_median"]:.4f}')
    print(f'  wrote {outdir / "graph_mode_results.json"}')


if __name__ == '__main__':
    main()
