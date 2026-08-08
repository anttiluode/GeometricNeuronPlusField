"""Causal source-move test of the soma amplitude-balance explanation.

Freeze a FunctionalArbor body, move one source inward along its existing
source-to-soma path, and ask whether the *change* in the amplitude-balance map
predicts the *change* in the temporal-order selectivity map.

See SOURCE_MOVE_BALANCE_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import deque
from pathlib import Path

import numpy as np


FRACTIONS = (0.75, 0.50, 0.25)


def n4_local(y, x, h, w):
    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        yy, xx = y + dy, x + dx
        if 0 <= yy < h and 0 <= xx < w:
            yield yy, xx


def gdist(body, start):
    b = np.asarray(body, bool)
    d = np.full(b.shape, -1, np.int32)
    start = tuple(start)
    if not b[start]:
        return d
    d[start] = 0
    q = deque([start])
    while q:
        p = q.popleft()
        for r in n4_local(*p, *b.shape):
            if b[r] and d[r] < 0:
                d[r] = d[p] + 1
                q.append(r)
    return d


def shortest_path(body, start, goal):
    """Occupied-cell shortest path from start to goal, inclusive."""
    b = np.asarray(body, bool)
    start, goal = tuple(start), tuple(goal)
    if not b[start] or not b[goal]:
        return None
    q = deque([start])
    prev = {start: None}
    while q:
        p = q.popleft()
        if p == goal:
            out = []
            while p is not None:
                out.append(p)
                p = prev[p]
            return out[::-1]
        for r in n4_local(*p, *b.shape):
            if b[r] and r not in prev:
                prev[r] = p
                q.append(r)
    return None


def pulse_at(m, cell, q):
    """Historical mature pulse shape injected at an arbitrary occupied cell."""
    c = m.cfg
    if not (0 <= q < c.pulse_frames):
        return 0.0
    env = math.sin(math.pi * (q + 1) / (c.pulse_frames + 1)) ** 2
    phase = np.exp(1j * c.carrier_omega * q)
    src = np.zeros_like(m.psi)
    src[tuple(cell)] = c.source_amp * env * phase
    return src


def addsrc(a, b):
    if isinstance(a, (float, int, np.floating)):
        return b
    if isinstance(b, (float, int, np.floating)):
        return a
    return a + b


def peak_single(m, cell, steps):
    m.reset_fast(True)
    acc = np.zeros(m.body.shape, float)
    for t in range(int(steps)):
        m.advance(pulse_at(m, cell, t), False, True, 'none')
        acc = np.maximum(acc, np.abs(m.psi) ** 2)
    return acc


def peak_order(m, cell_a, cell_b, lag, target, steps):
    m.reset_fast(True)
    acc = np.zeros(m.body.shape, float)
    first, second = ((cell_a, cell_b) if target else (cell_b, cell_a))
    for t in range(int(steps)):
        a = pulse_at(m, first, t)
        b = pulse_at(m, second, t - int(lag))
        m.advance(addsrc(a, b), False, True, 'none')
        acc = np.maximum(acc, np.abs(m.psi) ** 2)
    return acc


def safe_corr(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) < 3 or np.std(a) < 1e-14 or np.std(b) < 1e-14:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w, l)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2 * tail)), w, l


def wilcoxon_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (np.abs(z) > 1e-15)]
    if len(z) == 0:
        return float('nan')
    try:
        from scipy.stats import wilcoxon
        return float(wilcoxon(z, alternative='two-sided').pvalue)
    except Exception:
        return float('nan')


def top_overlap(a, b, frac=0.10):
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = len(a); k = max(1, int(math.ceil(frac * n)))
    ia = set(np.argsort(a, kind='mergesort')[-k:].tolist())
    ib = set(np.argsort(b, kind='mergesort')[-k:].tolist())
    return float(len(ia & ib) / k)


def position_on_path_from_soma(body, soma, terminal, frac):
    p = shortest_path(body, soma, terminal)
    if p is None or len(p) < 2:
        return None
    idx = int(round(float(frac) * (len(p) - 1)))
    idx = max(1, min(len(p) - 1, idx))
    return tuple(p[idx])


def condition_maps(m, cell_a, cell_b, lag, steps):
    pA = peak_single(m, cell_a, steps)
    pB = peak_single(m, cell_b, steps)
    T = peak_order(m, cell_a, cell_b, lag, True, steps)
    D = peak_order(m, cell_a, cell_b, lag, False, steps)
    C = (T - D) / (T + D + 1e-30)
    bal = np.minimum(pA, pB) / (np.maximum(pA, pB) + 1e-30)
    return pA, pB, np.abs(C), bal


def analyze_condition(body, soma, cells, absC, bal):
    idx = tuple(cells.T)
    cv = np.asarray(absC[idx], float)
    bv = np.asarray(bal[idx], float)
    si = int(np.flatnonzero((cells[:, 0] == soma[0]) & (cells[:, 1] == soma[1]))[0])
    ib = int(np.argmax(bv)); ic = int(np.argmax(cv))
    best_bal = tuple(map(int, cells[ib])); best_sel = tuple(map(int, cells[ic]))
    ds = gdist(body, soma)
    db = gdist(body, best_bal)
    return dict(
        corr_balance_selectivity=safe_corr(bv, cv),
        soma_balance=float(bv[si]),
        soma_selectivity=float(cv[si]),
        soma_balance_percentile=float(np.mean(bv < bv[si])),
        soma_selectivity_percentile=float(np.mean(cv < cv[si])),
        best_balance_cell=list(best_bal),
        best_selectivity_cell=list(best_sel),
        best_balance=float(bv[ib]),
        best_balance_selectivity=float(cv[ib]),
        best_selectivity=float(cv[ic]),
        best_balance_dist_soma=int(ds[best_bal]),
        best_selectivity_dist_soma=int(ds[best_sel]),
        best_to_best_distance=int(db[best_sel]),
        top_decile_overlap=top_overlap(bv, cv, 0.10),
    ), cv, bv


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=0)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/source_move_balance/source_move_balance.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    b = np.zeros((5, 5), np.uint8)
    b[2, :] = 1
    p = shortest_path(b, (2, 0), (2, 4))
    assert p == [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]
    assert position_on_path_from_soma(b, (2, 0), (2, 4), .5) == (2, 2)
    assert abs(safe_corr([0, 1, 2], [0, 2, 4]) - 1) < 1e-12
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
    for seed in range(a.seed_start, a.seed_start + a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature = True
        body = m.body.astype(bool)
        soma = tuple(m.soma)
        orig_a = tuple(m.source_terminal(0))
        orig_b = tuple(m.source_terminal(1))
        cells = np.argwhere(body)

        conditions = [('original', orig_a, orig_b)]
        seen = {(orig_a, orig_b)}
        for which, terminal in [('B', orig_b), ('A', orig_a)]:
            for frac in FRACTIONS:
                moved = position_on_path_from_soma(body, soma, terminal, frac)
                if moved is None:
                    continue
                ca, cb = (orig_a, moved) if which == 'B' else (moved, orig_b)
                if (ca, cb) in seen:
                    continue
                seen.add((ca, cb))
                conditions.append((f'{which}_{frac:.2f}', ca, cb))

        crows = []
        original_cv = original_bv = None
        for name, ca, cb in conditions:
            _, _, absC, bal = condition_maps(m, ca, cb, a.lag, a.steps)
            rec, cv, bv = analyze_condition(body, soma, cells, absC, bal)
            rec.update(name=name, source_A=list(ca), source_B=list(cb),
                       moved=bool(name != 'original'))
            if name == 'original':
                original_cv, original_bv = cv.copy(), bv.copy()
                rec['relocation_corr'] = float('nan')
                rec['displaced'] = bool(rec['best_balance_dist_soma'] >= 3)
                rec['balance_vs_soma_selectivity_delta'] = float(rec['best_balance_selectivity'] - rec['soma_selectivity'])
            else:
                rec['relocation_corr'] = safe_corr(bv - original_bv, cv - original_cv)
                rec['displaced'] = bool(rec['best_balance_dist_soma'] >= 3)
                rec['balance_vs_soma_selectivity_delta'] = float(rec['best_balance_selectivity'] - rec['soma_selectivity'])
            crows.append(rec)

        moved = [q for q in crows if q['moved']]
        reloc = np.asarray([q['relocation_corr'] for q in moved], float)
        displaced = [q for q in moved if q['displaced']]
        row = dict(
            seed=seed,
            cells=int(body.sum()),
            soma=list(soma),
            original_A=list(orig_a),
            original_B=list(orig_b),
            moved_conditions=len(moved),
            mean_relocation_corr=float(np.nanmean(reloc)) if np.isfinite(reloc).any() else float('nan'),
            positive_relocation_conditions=int(np.sum(reloc > 0)),
            displaced_conditions=len(displaced),
            mean_displaced_balance_vs_soma_delta=float(np.mean([q['balance_vs_soma_selectivity_delta'] for q in displaced])) if displaced else float('nan'),
            conditions=crows,
        )
        rows.append(row)
        print(f"seed {seed:2d}: reloc r {row['mean_relocation_corr']:+.3f} "
              f"({row['positive_relocation_conditions']}/{row['moved_conditions']} positive), "
              f"displaced {row['displaced_conditions']} "
              f"delta {row['mean_displaced_balance_vs_soma_delta']:+.3f}", flush=True)

    if not rows:
        raise SystemExit('No valid bodies')

    body_reloc = np.asarray([r['mean_relocation_corr'] for r in rows], float)
    p_reloc, w_reloc, l_reloc = sign_test_two_sided(body_reloc)
    body_disp = np.asarray([r['mean_displaced_balance_vs_soma_delta'] for r in rows], float)
    valid_disp = body_disp[np.isfinite(body_disp)]
    p_disp_sign, w_disp, l_disp = sign_test_two_sided(valid_disp)
    p_disp_wil = wilcoxon_two_sided(valid_disp)

    all_moved = [q for r in rows for q in r['conditions'] if q['moved']]
    all_displaced = [q for q in all_moved if q['displaced']]
    summary = dict(
        bodies=len(rows),
        moved_conditions=len(all_moved),
        displaced_conditions=len(all_displaced),
        P1=dict(
            body_mean_relocation_corr=float(np.nanmean(body_reloc)),
            body_median_relocation_corr=float(np.nanmedian(body_reloc)),
            positive_bodies=w_reloc,
            negative_bodies=l_reloc,
            sign_p=p_reloc,
            passed=bool(np.nanmean(body_reloc) > 0 and np.isfinite(p_reloc) and p_reloc < .05),
        ),
        P2=dict(
            valid_bodies=int(len(valid_disp)),
            body_mean_balance_vs_soma_delta=float(np.mean(valid_disp)) if len(valid_disp) else float('nan'),
            body_median_balance_vs_soma_delta=float(np.median(valid_disp)) if len(valid_disp) else float('nan'),
            positive_bodies=w_disp,
            negative_bodies=l_disp,
            sign_p=p_disp_sign,
            wilcoxon_two_sided_p=p_disp_wil,
            passed=bool(len(valid_disp) > 0 and np.mean(valid_disp) > 0 and
                        ((np.isfinite(p_disp_wil) and p_disp_wil < .05) or
                         (np.isfinite(p_disp_sign) and p_disp_sign < .05))),
        ),
        descriptive=dict(
            moved_corr_balance_selectivity_mean=float(np.mean([q['corr_balance_selectivity'] for q in all_moved])),
            moved_top_decile_overlap_mean=float(np.mean([q['top_decile_overlap'] for q in all_moved])),
            moved_best_to_best_distance_mean=float(np.mean([q['best_to_best_distance'] for q in all_moved])),
            moved_soma_selectivity_percentile_mean=float(np.mean([q['soma_selectivity_percentile'] for q in all_moved])),
            moved_soma_balance_percentile_mean=float(np.mean([q['soma_balance_percentile'] for q in all_moved])),
        ),
    )

    payload = dict(
        experiment='source_move_balance_v01',
        prereg='SOURCE_MOVE_BALANCE_PREREG_V01.md',
        seed_start=a.seed_start,
        seeds_requested=a.seeds,
        lag=a.lag,
        steps=a.steps,
        fractions=list(FRACTIONS),
        summary=summary,
        rows=rows,
    )
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print('\nSOURCE-MOVE BALANCE RECEIPT')
    print(f" bodies {summary['bodies']}  moved conditions {summary['moved_conditions']}")
    print(f" P1 relocation r mean {summary['P1']['body_mean_relocation_corr']:+.4f} "
          f"pos/neg {w_reloc}/{l_reloc} sign p={p_reloc:.5g} PASS={summary['P1']['passed']}")
    print(f" P2 displaced valid bodies {summary['P2']['valid_bodies']} "
          f"delta {summary['P2']['body_mean_balance_vs_soma_delta']:+.4f} "
          f"pos/neg {w_disp}/{l_disp} sign p={p_disp_sign:.5g} "
          f"wilcoxon p={p_disp_wil:.5g} PASS={summary['P2']['passed']}")
    print(f' wrote {out}')


if __name__ == '__main__':
    main()
