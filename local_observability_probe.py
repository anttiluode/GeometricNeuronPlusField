"""Can a small local soma aperture observe the confirmed global task-mode band?

The graph eigenbasis is an omniscient microscope. A physical soma/AIS interface is
local. This probe restricts confirmed modes 18-20 to graph-distance balls around
the soma and measures whether those three global coordinates remain independently
visible and well-conditioned.

Two controls are kept separate:

1. equal-size RANDOM CELL SETS -- distributed apertures with the same tap count;
2. graph-radius balls around EVERY OTHER BODY CELL -- the important locality-
   matched control, asking whether the soma is special among equally local apertures.

No wave simulation, learning, fitted readout, or active AIS model is used here.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np


def n4(y, x, h, w):
    for dy, dx in ((-1,0),(1,0),(0,-1),(0,1)):
        yy, xx = y+dy, x+dx
        if 0 <= yy < h and 0 <= xx < w:
            yield yy, xx


def graph_modes(body):
    coords = [tuple(map(int,p)) for p in np.argwhere(body > 0)]
    idx = {p:i for i,p in enumerate(coords)}
    A = np.zeros((len(coords),len(coords)), float)
    for p,i in idx.items():
        for q in n4(*p,*body.shape):
            j = idx.get(q)
            if j is not None:
                A[i,j] = 1.0
    A = np.maximum(A,A.T)
    L = np.diag(A.sum(1)) - A
    w,V = np.linalg.eigh(L)
    o = np.argsort(w)
    return coords, idx, w[o], V[:,o]


def distances(body, start):
    b = body.astype(bool)
    d = {tuple(start):0}
    q = deque([tuple(start)])
    while q:
        p = q.popleft()
        for r in n4(*p,*b.shape):
            if b[r] and r not in d:
                d[r] = d[p] + 1
                q.append(r)
    return d


def aperture_metrics(Vband, ids, total_cells):
    ids = np.asarray(ids, int)
    if ids.size == 0:
        return dict(n=0, cell_fraction=0.0, capture=0.0, smin=0.0, smax=0.0, isotropy=0.0, rank=0)
    M = Vband[ids, :]
    s = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape) * np.finfo(float).eps * (s[0] if len(s) else 1.0)
    rank = int((s > tol).sum())
    capture = float(np.sum(M*M) / Vband.shape[1])
    smax = float(s[0]) if len(s) else 0.0
    smin = float(s[-1]) if len(s) >= Vband.shape[1] else 0.0
    iso = float(smin / (smax + 1e-15))
    return dict(n=int(ids.size), cell_fraction=float(ids.size/total_cells),
                capture=capture, smin=smin, smax=smax, isotropy=iso, rank=rank)


def percentile(x, arr):
    arr = np.asarray(arr,float)
    return float((np.sum(arr < x) + 0.5*np.sum(arr == x)) / len(arr)) if len(arr) else float('nan')


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='../FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=0)
    ap.add_argument('--seeds', type=int, default=24)
    ap.add_argument('--band', default='18,19,20')
    ap.add_argument('--radii', default='0,1,2,3,4,5,6')
    ap.add_argument('--random-controls', type=int, default=250)
    ap.add_argument('--out', default='runs/local_observability/local_observability.json')
    return ap.parse_args()


def main():
    a = parse_args()
    fa = Path(a.functional_arbors).resolve()
    if not fa.exists():
        raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    band = [int(x) for x in a.band.split(',') if x.strip()]
    radii = [int(x) for x in a.radii.split(',') if x.strip()]
    rows=[]

    for seed in range(a.seed_start, a.seed_start+a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        coords, idx, w, V = graph_modes(m.body)
        if max(band) >= V.shape[1]:
            continue
        Vband = V[:,band]
        all_dist = {p: distances(m.body,p) for p in coords}
        ds = all_dist[tuple(m.soma)]
        rng = np.random.default_rng(seed+1701)
        per=[]

        for R in radii:
            soma_cells = [p for p in coords if ds.get(p,999) <= R]
            ids = [idx[p] for p in soma_cells]
            met = aperture_metrics(Vband,ids,len(coords))

            # Distributed equal-count controls.
            random_sets=[]; k=len(ids)
            if k>0:
                for _ in range(a.random_controls):
                    rid = rng.choice(len(coords), k, replace=False)
                    random_sets.append(aperture_metrics(Vband,rid,len(coords)))

            # Locality-matched controls: same graph radius around every possible center.
            local_balls=[]
            for center in coords:
                dc = all_dist[center]
                cid=[idx[p] for p in coords if dc.get(p,999) <= R]
                local_balls.append(aperture_metrics(Vband,cid,len(coords)))

            for key in ('capture','smin','isotropy'):
                met[f'{key}_random_mean'] = float(np.mean([q[key] for q in random_sets])) if random_sets else 0.0
                met[f'{key}_random_percentile'] = percentile(met[key],[q[key] for q in random_sets])
                met[f'{key}_localball_mean'] = float(np.mean([q[key] for q in local_balls])) if local_balls else 0.0
                met[f'{key}_localball_percentile'] = percentile(met[key],[q[key] for q in local_balls])
            met['localball_n_mean'] = float(np.mean([q['n'] for q in local_balls]))
            met['localball_n_median'] = float(np.median([q['n'] for q in local_balls]))
            per.append(dict(radius=R, **met))

        rows.append(dict(seed=seed, cells=len(coords), soma=list(map(int,m.soma)), radii=per))
        print(f'seed {seed:2d}: ' + ' '.join(
            f'R{q["radius"]}:n{q["n"]}/local-smin%{q["smin_localball_percentile"]:.2f}' for q in per))

    if not rows:
        raise SystemExit('No bodies')
    summary=[]
    for R in radii:
        rr=[next(q for q in row['radii'] if q['radius']==R) for row in rows]
        def mean(k): return float(np.mean([q[k] for q in rr]))
        def med(k): return float(np.median([q[k] for q in rr]))
        summary.append(dict(
            radius=R,
            n_mean=mean('n'),
            cell_fraction_mean=mean('cell_fraction'),
            capture_mean=mean('capture'),
            capture_random_mean=mean('capture_random_mean'),
            capture_random_enrichment=float(mean('capture')/(mean('capture_random_mean')+1e-15)),
            capture_random_percentile_median=med('capture_random_percentile'),
            capture_localball_mean=mean('capture_localball_mean'),
            capture_localball_percentile_median=med('capture_localball_percentile'),
            smin_mean=mean('smin'),
            smin_random_mean=mean('smin_random_mean'),
            smin_random_percentile_median=med('smin_random_percentile'),
            smin_localball_mean=mean('smin_localball_mean'),
            smin_localball_percentile_median=med('smin_localball_percentile'),
            isotropy_mean=mean('isotropy'),
            isotropy_random_mean=mean('isotropy_random_mean'),
            isotropy_random_percentile_median=med('isotropy_random_percentile'),
            isotropy_localball_mean=mean('isotropy_localball_mean'),
            isotropy_localball_percentile_median=med('isotropy_localball_percentile'),
            full_rank_count=int(sum(q['rank']==len(band) for q in rr)),
            bodies=len(rr),
        ))

    payload=dict(experiment='local_observability_v02', band=band, radii=radii,
                 random_controls=a.random_controls, summary=summary, rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')

    print('\nLOCAL OBSERVABILITY RECEIPT v0.2')
    for s in summary:
        print(f'  R={s["radius"]}: n={s["n_mean"]:.1f} '
              f'smin pct random={s["smin_random_percentile_median"]:.3f} '
              f'local-ball={s["smin_localball_percentile_median"]:.3f} '
              f'iso pct local-ball={s["isotropy_localball_percentile_median"]:.3f} '
              f'full-rank={s["full_rank_count"]}/{s["bodies"]}')
    print(f'  wrote {out}')


if __name__=='__main__':
    main()
