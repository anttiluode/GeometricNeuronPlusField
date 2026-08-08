"""Discover which graph-mode pairs create the soma order-sensitive cross term.

See MODE_PAIR_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from graph_mode_probe import graph_laplacian_modes
from transfer_decomposition_probe import trace_single, tshift


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w, l)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2 * tail)), w, l


def concentration(absvals):
    a = np.asarray(absvals, float)
    a = a[np.isfinite(a) & (a >= 0)]
    total = float(a.sum())
    n = len(a)
    if n == 0 or total <= 0:
        return dict(total_abs=0.0, k50=0, k80=0, frac50=0.0, frac80=0.0,
                    participation=0.0, participation_fraction=0.0)
    s = np.sort(a)[::-1]
    cs = np.cumsum(s) / total
    k50 = int(np.searchsorted(cs, .5) + 1)
    k80 = int(np.searchsorted(cs, .8) + 1)
    pr = float(total * total / (np.sum(a * a) + 1e-30))
    return dict(total_abs=total, k50=k50, k80=k80,
                frac50=float(k50 / n), frac80=float(k80 / n),
                participation=pr, participation_fraction=float(pr / n))


def band_label(i):
    i = int(i)
    if i == 0: return 'm0'
    if i <= 5: return 'm1-5'
    if i <= 11: return 'm6-11'
    if i <= 17: return 'm12-17'
    if i <= 20: return 'm18-20'
    if i <= 31: return 'm21-31'
    return 'm32+'


BANDS = ['m0', 'm1-5', 'm6-11', 'm12-17', 'm18-20', 'm21-31', 'm32+']


def body_probe(m, lag, steps, top_pairs=30, top_modes=12):
    hA = trace_single(m, 0, 1.0, steps)
    hB = trace_single(m, 1, 1.0, steps)

    coords, evals, evecs = graph_laplacian_modes(m.body)
    n = len(coords)
    cidx = {p: i for i, p in enumerate(coords)}
    soma = tuple(map(int, m.soma))
    if soma not in cidx:
        raise RuntimeError('soma missing from graph coordinates')
    si = cidx[soma]

    ys = np.asarray([p[0] for p in coords], int)
    xs = np.asarray([p[1] for p in coords], int)
    ZA = np.asarray(hA[:, ys, xs], np.complex128)
    ZB = np.asarray(hB[:, ys, xs], np.complex128)

    # For real orthonormal eigenvectors, q_n=<phi_n,z> and z = sum q_n phi_n.
    QA = ZA @ evecs
    QB = ZB @ evecs
    phi_s = np.asarray(evecs[si, :], float)
    uA = QA * phi_s[None, :]
    uB = QB * phi_s[None, :]

    hsA = np.asarray(hA[:, soma[0], soma[1]], np.complex128)
    hsB = np.asarray(hB[:, soma[0], soma[1]], np.complex128)
    recA = uA.sum(axis=1)
    recB = uB.sum(axis=1)
    source_rec_mae = float(0.5 * (np.mean(np.abs(recA - hsA)) + np.mean(np.abs(recB - hsB))))
    source_scale = float(0.5 * (np.mean(np.abs(hsA)) + np.mean(np.abs(hsB))) + 1e-30)
    source_rec_rel = float(source_rec_mae / source_scale)

    uA_l = tshift(uA, lag)
    uB_l = tshift(uB, lag)
    hsA_l = tshift(hsA, lag)
    hsB_l = tshift(hsB, lag)

    psiT = recA + uB_l.sum(axis=1)
    psiD = recB + uA_l.sum(axis=1)
    pT = np.abs(psiT) ** 2
    pD = np.abs(psiD) ** 2
    tT = int(np.argmax(pT)); tD = int(np.argmax(pD))
    soma_C = float((pT[tT] - pD[tD]) / (pT[tT] + pD[tD] + 1e-30))

    MT = 2.0 * np.real(np.outer(uA[tT], np.conj(uB_l[tT])))
    MD = 2.0 * np.real(np.outer(uB[tD], np.conj(uA_l[tD])))
    M = MT - MD

    direct_cross_delta = float(
        2.0 * np.real(hsA[tT] * np.conj(hsB_l[tT]))
        - 2.0 * np.real(hsB[tD] * np.conj(hsA_l[tD]))
    )
    pair_sum = float(M.sum())
    pair_rec_abs_error = float(abs(pair_sum - direct_cross_delta))
    pair_rec_rel_error = float(pair_rec_abs_error / (abs(direct_cross_delta) + 1e-30))

    # Unordered pair contributions preserve the signed sum while avoiding double counting.
    pair_rows = []
    absvals = []
    mode_involve = np.zeros(n, float)
    diag_abs = 0.0
    for i in range(n):
        v = float(M[i, i])
        av = abs(v)
        pair_rows.append((i, i, v, av))
        absvals.append(av)
        mode_involve[i] += av
        diag_abs += av
        for j in range(i + 1, n):
            v = float(M[i, j] + M[j, i])
            av = abs(v)
            pair_rows.append((i, j, v, av))
            absvals.append(av)
            mode_involve[i] += 0.5 * av
            mode_involve[j] += 0.5 * av

    conc = concentration(absvals)
    total_abs = conc['total_abs'] + 1e-30
    cancellation = float(abs(pair_sum) / total_abs)
    diagonal_fraction = float(diag_abs / total_abs)
    mode_frac = mode_involve / (mode_involve.sum() + 1e-30)

    task_indices = [i for i in (18, 19, 20) if i < n]
    task_frac = float(mode_frac[task_indices].sum()) if task_indices else 0.0
    expected_task_frac = float(len(task_indices) / n) if n else 0.0
    task_enrichment = float(task_frac / expected_task_frac) if expected_task_frac > 0 else float('nan')

    pair_rows_sorted = sorted(pair_rows, key=lambda r: r[3], reverse=True)
    strongest_pairs = [dict(
        i=int(i), j=int(j), signed=float(v), abs=float(av), abs_fraction=float(av / total_abs),
        lambda_i=float(evals[i]), lambda_j=float(evals[j]),
        band_i=band_label(i), band_j=band_label(j),
    ) for i, j, v, av in pair_rows_sorted[:int(top_pairs)]]

    mode_order = np.argsort(mode_frac)[::-1]
    strongest_modes = [dict(
        index=int(i), involvement_fraction=float(mode_frac[i]), eigenvalue=float(evals[i]), band=band_label(i)
    ) for i in mode_order[:int(top_modes)]]

    # Ordered spectral-block absolute interaction matrix, descriptive only.
    block = {bi: {bj: 0.0 for bj in BANDS} for bi in BANDS}
    ordered_abs_total = float(np.abs(M).sum()) + 1e-30
    for i in range(n):
        bi = band_label(i)
        for j in range(n):
            bj = band_label(j)
            block[bi][bj] += float(abs(M[i, j]) / ordered_abs_total)

    return dict(
        cells=n,
        soma=[int(soma[0]), int(soma[1])],
        peak_t_target=tT,
        peak_t_distractor=tD,
        soma_contrast=soma_C,
        source_modal_reconstruction_mae=source_rec_mae,
        source_modal_reconstruction_rel=source_rec_rel,
        cross_delta_direct=direct_cross_delta,
        cross_delta_pair_sum=pair_sum,
        pair_reconstruction_abs_error=pair_rec_abs_error,
        pair_reconstruction_rel_error=pair_rec_rel_error,
        pair_count=len(pair_rows),
        concentration=conc,
        cancellation_ratio=cancellation,
        diagonal_abs_fraction=diagonal_fraction,
        task_band_18_20_involvement_fraction=task_frac,
        task_band_dimensional_fraction=expected_task_frac,
        task_band_enrichment=task_enrichment,
        mode_involvement_fraction=[float(x) for x in mode_frac],
        strongest_modes=strongest_modes,
        strongest_pairs=strongest_pairs,
        spectral_block_abs_fraction=block,
    )


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=72)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/mode_pair_discovery/mode_pair_discovery.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    a = np.array([1., 2., 3., 4.])
    c = concentration(a)
    assert c['k50'] == 2
    assert c['k80'] == 3
    assert band_label(18) == 'm18-20'
    print('selftest ok')


def main():
    a = parse_args()
    if a.selftest:
        selftest(); return
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    rows = []
    mode_sum = defaultdict(float); mode_n = defaultdict(int)
    pair_sum = defaultdict(float); pair_n = defaultdict(int)

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature = True
        r = body_probe(m, a.lag, a.steps)
        r['seed'] = int(seed)
        rows.append(r)

        for i, x in enumerate(r['mode_involvement_fraction']):
            mode_sum[i] += float(x); mode_n[i] += 1
        for q in r['strongest_pairs']:
            k = (int(q['i']), int(q['j']))
            pair_sum[k] += float(q['abs_fraction']); pair_n[k] += 1

        c = r['concentration']
        print(f"seed {seed}: N={r['cells']} C={r['soma_contrast']:+.3f} "
              f"k50={c['k50']}/{r['pair_count']} ({c['frac50']:.3f}) "
              f"k80={c['k80']} ({c['frac80']:.3f}) PRfrac={c['participation_fraction']:.3f} "
              f"cancel={r['cancellation_ratio']:.3f} diag={r['diagonal_abs_fraction']:.3f} "
              f"band18-20 enrich={r['task_band_enrichment']:.2f}", flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    def arr(path):
        out = []
        for r in rows:
            x = r
            for k in path:
                x = x[k]
            out.append(float(x))
        return np.asarray(out, float)

    mean_mode = sorted(((i, mode_sum[i] / mode_n[i], mode_n[i]) for i in mode_sum),
                       key=lambda x: x[1], reverse=True)
    # Pair aggregation uses the per-body strongest-pair list; missing there means < top-pairs cutoff.
    mean_pair = sorted(((k[0], k[1], pair_sum[k] / len(rows), pair_n[k]) for k in pair_sum),
                       key=lambda x: x[2], reverse=True)

    enrich = arr(['task_band_enrichment'])
    ep, ew, el = sign_test_two_sided(enrich - 1.0)
    summary = dict(
        bodies=len(rows),
        seed_start=a.seed_start,
        mean_cells=float(np.mean([r['cells'] for r in rows])),
        source_reconstruction_rel_mean=float(arr(['source_modal_reconstruction_rel']).mean()),
        pair_reconstruction_rel_mean=float(arr(['pair_reconstruction_rel_error']).mean()),
        mean_pair_fraction_50=float(arr(['concentration','frac50']).mean()),
        median_pair_fraction_50=float(np.median(arr(['concentration','frac50']))),
        mean_pair_fraction_80=float(arr(['concentration','frac80']).mean()),
        median_pair_fraction_80=float(np.median(arr(['concentration','frac80']))),
        mean_effective_pair_fraction=float(arr(['concentration','participation_fraction']).mean()),
        median_effective_pair_fraction=float(np.median(arr(['concentration','participation_fraction']))),
        mean_cancellation_ratio=float(arr(['cancellation_ratio']).mean()),
        median_cancellation_ratio=float(np.median(arr(['cancellation_ratio']))),
        mean_diagonal_abs_fraction=float(arr(['diagonal_abs_fraction']).mean()),
        mean_task_band_18_20_involvement=float(arr(['task_band_18_20_involvement_fraction']).mean()),
        mean_task_band_18_20_enrichment=float(enrich.mean()),
        median_task_band_18_20_enrichment=float(np.median(enrich)),
        task_band_enrichment_gt1_bodies=ew,
        task_band_enrichment_lt1_bodies=el,
        task_band_enrichment_sign_p=ep,
        discovery_top_modes=[dict(index=int(i), mean_involvement_fraction=float(v), bodies=int(nn))
                             for i, v, nn in mean_mode[:16]],
        discovery_top_pairs=[dict(i=int(i), j=int(j), mean_abs_fraction_over_all_bodies=float(v),
                                  appeared_in_top_list_bodies=int(nn))
                             for i, j, v, nn in mean_pair[:30]],
    )

    payload = dict(
        experiment='mode_pair_discovery_v01',
        prereg='MODE_PAIR_DISCOVERY_PREREG_V01.md',
        lag=a.lag,
        steps=a.steps,
        summary=summary,
        rows=rows,
    )
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nMODE-PAIR DISCOVERY RECEIPT')
    print(f" bodies={len(rows)}")
    print(f" mean k50 fraction={summary['mean_pair_fraction_50']:.4f}")
    print(f" mean k80 fraction={summary['mean_pair_fraction_80']:.4f}")
    print(f" mean effective-pair fraction={summary['mean_effective_pair_fraction']:.4f}")
    print(f" mean cancellation={summary['mean_cancellation_ratio']:.4f}")
    print(f" mean diagonal fraction={summary['mean_diagonal_abs_fraction']:.4f}")
    print(f" band18-20 involvement={summary['mean_task_band_18_20_involvement']:.4f} "
          f"enrichment={summary['mean_task_band_18_20_enrichment']:.3f} "
          f"gt/lt1={ew}/{el} p={ep:.5g}")
    print(' top modes:', [(q['index'], round(q['mean_involvement_fraction'],4)) for q in summary['discovery_top_modes'][:10]])
    print(' top pairs:', [(q['i'], q['j'], round(q['mean_abs_fraction_over_all_bodies'],4)) for q in summary['discovery_top_pairs'][:10]])
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
