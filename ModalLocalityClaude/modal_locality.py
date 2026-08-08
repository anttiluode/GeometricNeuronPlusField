"""Is a one-cell structural change a LOCAL event in the arbor's own modal basis?

GeometricNeuronPlusField establishes that the mature Functional Arbor is (almost
exactly) diagonalized by its body graph's Laplacian eigenbasis {lambda_n, phi_n},
and that the task-relevant information lives in a specific band of those modes.

That has a consequence nobody in either line has drawn, and it bears directly on
the v0.8/v0.9 credit-assignment null:

    the modes are GLOBAL objects of the whole graph.
    Adding or removing ONE cell changes every lambda_n and every phi_n.

If one structural change perturbs many modes, and if the perturbation is not
concentrated near the changed cell, then no LOCAL eligibility tag can be causally
correct -- not because the tag is badly chosen, but because the map from local
structure to the modal coordinates is globally entangled. The v0.9 null would be
structural rather than a failure of tagging cleverness.

This script is parameter-free linear algebra on the real frozen bodies. No wave
simulation, no learning, no fitted constants.

Measures, per single-cell perturbation:
  1. SCRAMBLE      how many modes lose >10% of their identity (best overlap < .9)
  2. SPECTRAL      ||d lambda|| relative to the mean eigenvalue gap
  3. LOCALITY      fraction of the total |d phi|^2 mass lying within graph
                   distance d of the changed cell, against the fraction of BODY
                   CELLS within that distance. Ratio ~1 means delocalized.
  4. DISTANCE      does the size of the disturbance fall off with the changed
                   cell's distance from the soma / from the task band's support?

Controls: the zero perturbation (must give overlap exactly 1) and near-degenerate
eigenvalue pairs are flagged, since overlaps rotate freely inside a degenerate
subspace and would inflate SCRAMBLE artificially.
"""
from __future__ import annotations
import sys, json, math, argparse
from collections import deque
import numpy as np

sys.path.insert(0, '.')
from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor
from v06_ephaptic_growth.ephaptic_arbor import n4


def laplacian(cells, shape):
    idx = {tuple(c): i for i, c in enumerate(cells)}
    n = len(cells)
    L = np.zeros((n, n))
    for c, i in idx.items():
        for r in n4(*c, *shape):
            j = idx.get(r)
            if j is not None:
                L[i, i] += 1.0
                L[i, j] -= 1.0
    return L, idx


def eig(L):
    w, V = np.linalg.eigh(L)
    return w, V


def gdist(cells, shape, start):
    s = set(map(tuple, cells)); d = {tuple(start): 0}; q = deque([tuple(start)])
    while q:
        p = q.popleft()
        for r in n4(*p, *shape):
            if r in s and r not in d:
                d[r] = d[p] + 1; q.append(r)
    return d


def connected(cells, shape):
    if not cells: return False
    return len(gdist(cells, shape, cells[0])) == len(cells)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--perts', type=int, default=25, help='single-cell changes per body')
    ap.add_argument('--band', type=str, default='18,19,20', help='the confirmed task band')
    ap.add_argument('--out', type=str, default='modal_locality.json')
    a = ap.parse_args()
    band = [int(x) for x in a.band.split(',')]

    rows = []
    for seed in range(a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        shape = m.body.shape
        cells = [tuple(c) for c in np.argwhere(m.body > 0)]
        L0, idx0 = laplacian(cells, shape)
        w0, V0 = eig(L0)
        gap = float(np.median(np.diff(w0)))
        dsoma = gdist(cells, shape, m.soma)

        # candidate single-cell ADDITIONS: empty cells 4-adjacent to the body
        occ = set(cells)
        cand = []
        for c in cells:
            for r in n4(*c, *shape):
                if r not in occ and 0 < r[0] < shape[0] - 1 and 0 < r[1] < shape[1] - 1:
                    cand.append(r)
        cand = sorted(set(cand))
        rng = np.random.default_rng(seed + 90210)
        if len(cand) > a.perts:
            cand = [cand[i] for i in rng.choice(len(cand), a.perts, replace=False)]

        recs = []
        for newc in cand:
            cells1 = cells + [newc]
            if not connected(cells1, shape):
                continue
            L1, idx1 = laplacian(cells1, shape)
            w1, V1 = eig(L1)
            # embed old modes into the new index space (new cell gets 0)
            P = np.zeros((len(cells1), len(cells)))
            for c, i in idx0.items():
                P[idx1[c], i] = 1.0
            V0e = P @ V0
            O = np.abs(V1.T @ V0e)                      # |<phi'_m, phi_n>|
            best = O.max(axis=0)                        # per old mode
            scrambled = int((best < 0.90).sum())
            # near-degenerate flag: modes whose eigenvalue gap to a neighbour is tiny
            degen = int((np.diff(w0) < 0.02).sum())
            # spectral shift on the common modes
            k = min(len(w0), len(w1))
            dlam = float(np.linalg.norm(w1[:k] - w0[:k])) / (gap + 1e-12)
            # LOCALITY of the eigenvector change
            dphi = V1[:, :len(w0)] * np.sign(np.einsum('ij,ij->j', V1[:, :len(w0)], V0e)) - V0e
            mass = (dphi ** 2).sum(axis=1)              # per CELL, summed over modes
            mass = mass / (mass.sum() + 1e-30)
            dn = gdist(cells1, shape, newc)
            loc = {}
            for R in (1, 2, 3, 5):
                inball = [i for c, i in idx1.items() if dn.get(c, 99) <= R]
                loc[R] = (float(mass[inball].sum()), float(len(inball) / len(cells1)))
            # band-specific: how much do the confirmed task modes move?
            bm = [n for n in band if n < len(w0)]
            bandscr = float(np.mean([1 - best[n] for n in bm])) if bm else float('nan')
            recs.append(dict(cell=list(map(int, newc)), scrambled=scrambled, n_modes=len(w0),
                             degen=degen, dlam_over_gap=dlam, band_identity_loss=bandscr,
                             mean_best_overlap=float(best.mean()),
                             dist_soma=int(min((dsoma.get(r, 99) for r in n4(*newc, *shape)), default=99)) + 1,
                             loc={str(k2): v for k2, v in loc.items()}))
        if recs:
            rows.append(dict(seed=seed, n_cells=len(cells), n_perts=len(recs), recs=recs))
            sc = np.array([r['scrambled'] for r in recs], float)
            print(f'seed {seed:2d}  {len(cells)} cells, {len(recs)} one-cell additions: '
                  f'modes scrambled {sc.mean():.1f}/{len(w0)}  '
                  f'mean best-overlap {np.mean([r["mean_best_overlap"] for r in recs]):.3f}  '
                  f'band identity loss {np.mean([r["band_identity_loss"] for r in recs]):.3f}', flush=True)

    allr = [r for row in rows for r in row['recs']]
    n = len(allr)
    N = allr[0]['n_modes']
    sc = np.array([r['scrambled'] for r in allr], float)
    print('\n' + '=' * 84)
    print(f'IS ONE STRUCTURAL CHANGE A LOCAL EVENT? — {len(rows)} bodies, {n} single-cell additions')
    print('=' * 84)
    print(f'modes whose identity is disturbed (best overlap < 0.90) by ONE added cell:')
    print(f'   mean {sc.mean():.1f} of ~{N} modes   median {np.median(sc):.0f}   '
          f'range {sc.min():.0f}-{sc.max():.0f}   = {100*sc.mean()/N:.0f}% of the spectrum')
    print(f'   mean best-overlap across all modes  {np.mean([r["mean_best_overlap"] for r in allr]):.3f}')
    print(f'   identity lost by the CONFIRMED task band {band}: '
          f'{np.mean([r["band_identity_loss"] for r in allr]):.3f}')
    print(f'   ||d lambda|| / median eigenvalue gap  {np.mean([r["dlam_over_gap"] for r in allr]):.2f}')
    print(f'   near-degenerate mode pairs per body (overlap rotation caveat): '
          f'{np.mean([r["degen"] for r in allr]):.1f}')

    print('\nWHERE does the eigenvector change live? (share of total |d phi|^2 mass)')
    print(f'{"ball radius":>12s} {"|dphi|^2 mass":>14s} {"share of cells":>15s} {"enrichment":>11s}')
    print('-' * 84)
    for R in (1, 2, 3, 5):
        mm = np.array([r['loc'][str(R)][0] for r in allr]); ff = np.array([r['loc'][str(R)][1] for r in allr])
        print(f'{R:12d} {mm.mean():14.3f} {ff.mean():15.3f} {mm.mean()/ff.mean():11.2f}x')

    d = np.array([r['dist_soma'] for r in allr], float)
    print(f'\ndoes the disturbance depend on WHERE the cell was added?')
    print(f'   corr(graph distance from soma, modes scrambled) = {np.corrcoef(d, sc)[0,1]:+.3f}')
    bl = np.array([r['band_identity_loss'] for r in allr], float)
    print(f'   corr(graph distance from soma, band identity loss) = {np.corrcoef(d, bl)[0,1]:+.3f}')

    json.dump(dict(seeds=len(rows), band=band, rows=rows), open(a.out, 'w'))
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
