"""Is a one-cell structural change a LOCAL event in the arbor's modal basis?

This is the locality audit prompted by the graph-mode reduction.  It asks whether a
single local anatomical change has a local effect in the geometry-defined normal
coordinates, or whether it perturbs a substantial fraction of the global spectrum.

The calculation is pure graph linear algebra on frozen FunctionalArbor bodies:
no wave simulation, growth, credit, or fitted constants.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np


def n4(y, x, h, w):
    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        yy, xx = y + dy, x + dx
        if 0 <= yy < h and 0 <= xx < w:
            yield yy, xx


def laplacian(cells, shape):
    idx = {tuple(c): i for i, c in enumerate(cells)}
    n = len(cells)
    L = np.zeros((n, n), np.float64)
    for c, i in idx.items():
        for r in n4(*c, *shape):
            j = idx.get(r)
            if j is not None:
                L[i, i] += 1.0
                L[i, j] -= 1.0
    return L, idx


def gdist(cells, shape, start):
    s = set(map(tuple, cells))
    d = {tuple(start): 0}
    q = deque([tuple(start)])
    while q:
        p = q.popleft()
        for r in n4(*p, *shape):
            if r in s and r not in d:
                d[r] = d[p] + 1
                q.append(r)
    return d


def connected(cells, shape):
    return bool(cells) and len(gdist(cells, shape, cells[0])) == len(cells)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='../FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=0)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--perts', type=int, default=25)
    ap.add_argument('--band', default='18,19,20')
    ap.add_argument('--out', default='runs/modal_locality/modal_locality.json')
    return ap.parse_args()


def main():
    a = parse_args()
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    band = [int(x) for x in a.band.split(',') if x.strip()]
    rows = []

    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        shape = m.body.shape
        cells = [tuple(map(int, c)) for c in np.argwhere(m.body > 0)]
        L0, idx0 = laplacian(cells, shape)
        w0, V0 = np.linalg.eigh(L0)
        gap = float(np.median(np.diff(w0)))
        dsoma = gdist(cells, shape, m.soma)

        occ = set(cells)
        cand = sorted({r for c in cells for r in n4(*c, *shape)
                       if r not in occ and 0 < r[0] < shape[0]-1 and 0 < r[1] < shape[1]-1})
        rng = np.random.default_rng(seed + 90210)
        if len(cand) > a.perts:
            cand = [cand[i] for i in rng.choice(len(cand), a.perts, replace=False)]

        recs = []
        for newc in cand:
            cells1 = cells + [newc]
            if not connected(cells1, shape):
                continue
            L1, idx1 = laplacian(cells1, shape)
            w1, V1 = np.linalg.eigh(L1)

            P = np.zeros((len(cells1), len(cells)), np.float64)
            for c, i in idx0.items():
                P[idx1[c], i] = 1.0
            V0e = P @ V0
            O = np.abs(V1.T @ V0e)
            best = O.max(axis=0)
            scrambled = int((best < 0.90).sum())
            degen = int((np.diff(w0) < 0.02).sum())
            k = min(len(w0), len(w1))
            dlam = float(np.linalg.norm(w1[:k] - w0[:k]) / (gap + 1e-12))

            signs = np.sign(np.einsum('ij,ij->j', V1[:, :len(w0)], V0e))
            signs[signs == 0] = 1.0
            dphi = V1[:, :len(w0)] * signs - V0e
            mass = (dphi ** 2).sum(axis=1)
            mass /= mass.sum() + 1e-30
            dn = gdist(cells1, shape, newc)
            loc = {}
            for R in (1, 2, 3, 5):
                inball = [i for c, i in idx1.items() if dn.get(c, 999) <= R]
                loc[str(R)] = [float(mass[inball].sum()), float(len(inball) / len(cells1))]

            bm = [n for n in band if n < len(w0)]
            bandloss = float(np.mean([1.0 - best[n] for n in bm])) if bm else float('nan')
            near = [dsoma.get(r, 999) for r in n4(*newc, *shape)]
            recs.append(dict(
                cell=list(map(int, newc)), scrambled=scrambled, n_modes=len(w0),
                degen=degen, dlam_over_gap=dlam, band_identity_loss=bandloss,
                mean_best_overlap=float(best.mean()), dist_soma=int(min(near) + 1), loc=loc,
            ))

        if recs:
            rows.append(dict(seed=seed, n_cells=len(cells), n_perts=len(recs), recs=recs))
            print(f'seed {seed:2d}: {len(recs)} perturbations')

    allr = [r for row in rows for r in row['recs']]
    if not allr:
        raise SystemExit('No valid perturbations')
    N = allr[0]['n_modes']
    sc = np.asarray([r['scrambled'] for r in allr], float)
    bl = np.asarray([r['band_identity_loss'] for r in allr], float)
    dl = np.asarray([r['dlam_over_gap'] for r in allr], float)
    ds = np.asarray([r['dist_soma'] for r in allr], float)

    locality = {}
    for R in (1, 2, 3, 5):
        mm = np.asarray([r['loc'][str(R)][0] for r in allr])
        ff = np.asarray([r['loc'][str(R)][1] for r in allr])
        locality[str(R)] = dict(dphi_mass_mean=float(mm.mean()), cell_share_mean=float(ff.mean()),
                                enrichment=float(mm.mean() / ff.mean()))

    summary = dict(
        bodies=len(rows), perturbations=len(allr), modes=N,
        scrambled_mean=float(sc.mean()), scrambled_median=float(np.median(sc)),
        scrambled_min=float(sc.min()), scrambled_max=float(sc.max()),
        spectrum_fraction_scrambled=float(sc.mean() / N),
        mean_best_overlap=float(np.mean([r['mean_best_overlap'] for r in allr])),
        task_band_identity_loss=float(bl.mean()),
        dlam_over_gap_mean=float(dl.mean()),
        near_degenerate_pairs_mean=float(np.mean([r['degen'] for r in allr])),
        corr_distance_soma_scrambled=float(np.corrcoef(ds, sc)[0, 1]),
        corr_distance_soma_band_loss=float(np.corrcoef(ds, bl)[0, 1]),
        locality=locality,
    )
    payload = dict(experiment='modal_locality_v01', band=band, summary=summary, rows=rows)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nMODAL LOCALITY RECEIPT')
    print(f'  one cell scrambles         {summary["scrambled_mean"]:.2f}/{N} modes')
    print(f'  spectrum fraction          {summary["spectrum_fraction_scrambled"]:.3f}')
    print(f'  task-band identity loss    {summary["task_band_identity_loss"]:.4f}')
    print(f'  ||d lambda|| / gap         {summary["dlam_over_gap_mean"]:.2f}')
    print(f'  corr(distance,scramble)    {summary["corr_distance_soma_scrambled"]:+.3f}')
    for R, v in locality.items():
        print(f'  locality R={R}: enrichment {v["enrichment"]:.3f}x')
    print(f'  wrote {out}')


if __name__ == '__main__':
    main()
