"""Can a small local soma aperture observe the confirmed global task-mode band?

The graph eigenbasis is an omniscient microscope.  A physical soma/AIS interface is
local.  This probe asks a geometry-only observability question on frozen bodies:
restrict modes 18-20 to graph-distance balls around the soma and measure whether
those three global coordinates remain visible and well-conditioned locally.

Controls are equal-size random cell sets on the SAME body.  No wave simulation,
learning, fitted readout, or active AIS model is used here.  A positive result would
say the body itself brings the task band into a locally observable mixture near the
soma.  A null says passive geometry alone does not privilege the soma; an active
boundary would then have to supply the missing selectivity.
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
        return dict(n=0, capture=0.0, smin=0.0, smax=0.0, isotropy=0.0, rank=0)
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
        ds = distances(m.body, m.soma)
        rng = np.random.default_rng(seed+1701)
        per=[]
        for R in radii:
            cells = [p for p in coords if ds.get(p,999) <= R]
            ids = [idx[p] for p in cells]
            met = aperture_metrics(Vband,ids,len(coords))
            controls=[]
            k=len(ids)
            if k>0:
                for _ in range(a.random_controls):
                    rid = rng.choice(len(coords), k, replace=False)
                    controls.append(aperture_metrics(Vband,rid,len(coords)))
            met['capture_random_mean'] = float(np.mean([q['capture'] for q in controls])) if controls else 0.0
            met['smin_random_mean'] = float(np.mean([q['smin'] for q in controls])) if controls else 0.0
            met['isotropy_random_mean'] = float(np.mean([q['isotropy'] for q in controls])) if controls else 0.0
            met['capture_percentile'] = percentile(met['capture'], [q['capture'] for q in controls])
            met['smin_percentile'] = percentile(met['smin'], [q['smin'] for q in controls])
            met['isotropy_percentile'] = percentile(met['isotropy'], [q['isotropy'] for q in controls])
            per.append(dict(radius=R, **met))
        rows.append(dict(seed=seed, cells=len(coords), soma=list(map(int,m.soma)), radii=per))
        print(f'seed {seed:2d}: ' + ' '.join(
            f'R{q["radius"]}:n{q["n"]}/smin%{q["smin_percentile"]:.2f}' for q in per))

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
            capture_enrichment=float(mean('capture')/(mean('capture_random_mean')+1e-15)),
            smin_mean=mean('smin'),
            smin_random_mean=mean('smin_random_mean'),
            smin_percentile_median=med('smin_percentile'),
            isotropy_mean=mean('isotropy'),
            isotropy_random_mean=mean('isotropy_random_mean'),
            isotropy_percentile_median=med('isotropy_percentile'),
            full_rank_count=int(sum(q['rank']==len(band) for q in rr)),
            bodies=len(rr),
        ))

    payload=dict(experiment='local_observability_v01', band=band, radii=radii,
                 random_controls=a.random_controls, summary=summary, rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')

    print('\nLOCAL OBSERVABILITY RECEIPT')
    for s in summary:
        print(f'  R={s["radius"]}: n={s["n_mean"]:.1f} capture x{s["capture_enrichment"]:.2f} '
              f'smin pct={s["smin_percentile_median"]:.2f} iso pct={s["isotropy_percentile_median"]:.2f} '
              f'full-rank={s["full_rank_count"]}/{s["bodies"]}')
    print(f'  wrote {out}')


if __name__=='__main__':
    main()
