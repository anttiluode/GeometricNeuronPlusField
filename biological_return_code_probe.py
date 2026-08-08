"""How much of the exact local adjoint survives biologically coarse soma return codes?

This is a hostile bridge test, not a biological model.

The exact reciprocal result injects the full complex soma derivative waveform back into
an otherwise reciprocal arbor.  A biological back-propagating consequence event is
much coarser: action potentials are stereotyped, sparse, delayed, actively propagated,
and gated by inhibition/rhythms.

This probe freezes the forward field and changes ONLY the soma return code.  Every
variant is L2 dose-matched to the exact return separately for target and distractor so
map changes cannot be explained by total return energy alone.

Variants:
  exact              full complex derivative waveform (positive control)
  phase_only         full timing + phase, constant amplitude
  envelope_signed    amplitude envelope + task sign, no carrier phase
  real_wave          real part only
  sparse_phase_N     N fixed-amplitude events, event times from |g| peaks, local phase kept
  sparse_signed_N    N fixed-amplitude events, task sign only, carrier phase discarded
  sparse_positive_N  N fixed positive events, even task sign discarded
  delayed_sparse8_d  sparse_phase_8 shifted by d frames with zero fill
  gate_P             exact waveform retained only in 50%-duty periodic windows; every
                     possible phase offset is evaluated and summarized best/median/worst

No theta/gamma frequency claim is made: simulation frame rate is not biologically
calibrated.  The gate tests only whether periodic temporal multiplexing can discard
large fractions of the return waveform without destroying structural direction.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # exact mature-boundary operator patch
from reciprocal_adjoint_probe import physical_credit_history, overlap_gradient, flat_pair, normalized_l2
from transfer_decomposition_probe import safe_corr

SPARSE_NS = (1, 2, 4, 8, 16, 32)
DELAYS = (1, 2, 4, 8, 16)
GATE_PERIODS = (14, 30, 42)


def l2_match(x, ref):
    x = np.asarray(x, np.complex128)
    ref = np.asarray(ref, np.complex128)
    nx = float(np.linalg.norm(x))
    nr = float(np.linalg.norm(ref))
    if nx <= 1e-30 or nr <= 1e-30:
        return np.zeros_like(x)
    return x * (nr / nx)


def top_peaks(g, n, min_sep=2):
    a = np.abs(np.asarray(g))
    order = np.argsort(a)[::-1]
    chosen = []
    for i in order:
        i = int(i)
        if a[i] <= 0:
            break
        if all(abs(i-j) >= min_sep for j in chosen):
            chosen.append(i)
            if len(chosen) >= int(n):
                break
    return np.asarray(sorted(chosen), int)


def sparse_code(g, n, mode, task_sign):
    g = np.asarray(g, np.complex128)
    idx = top_peaks(g, n)
    out = np.zeros_like(g)
    if len(idx) == 0:
        return out
    if mode == 'phase':
        ph = g[idx] / (np.abs(g[idx]) + 1e-30)
        out[idx] = ph
    elif mode == 'signed':
        out[idx] = float(np.sign(task_sign) or 1.0)
    elif mode == 'positive':
        out[idx] = 1.0
    else:
        raise ValueError(mode)
    return l2_match(out, g)


def shift_zero(g, d):
    g = np.asarray(g, np.complex128)
    out = np.zeros_like(g)
    d = int(d)
    if d == 0:
        return g.copy()
    if d > 0 and d < len(g):
        out[d:] = g[:-d]
    elif d < 0 and -d < len(g):
        out[:d] = g[-d:]
    return l2_match(out, g)


def gate_code(g, period, offset):
    g = np.asarray(g, np.complex128)
    p = int(period)
    on = p // 2
    t = np.arange(len(g))
    mask = ((t + int(offset)) % p) < on
    return l2_match(g * mask, g), float(np.mean(mask))


def map_metrics(ex_h, ex_v, ap_h, ap_v):
    ex = flat_pair(ex_h, ex_v)
    ap = flat_pair(ap_h, ap_v)
    mx = float(np.max(np.abs(ex)) + 1e-30)
    mask = np.abs(ex) > 0.01 * mx
    return dict(
        corr=float(safe_corr(ex, ap)),
        relative_l2=normalized_l2(ex, ap),
        strong_sign_agreement=float(np.mean(np.sign(ex[mask]) == np.sign(ap[mask]))) if np.any(mask) else float('nan'),
    )


def gradient_from_codes(m, wh, wv, pT, pD, gT, gD):
    muT, _, _ = physical_credit_history(m, wh, wv, gT, reverse=True)
    muD, _, _ = physical_credit_history(m, wh, wv, gD, reverse=True)
    hT, vT = overlap_gradient(m, pT, muT)
    hD, vD = overlap_gradient(m, pD, muD)
    return hT + hD, vT + vD


def build_exact(m, lag, steps):
    wh, wv = ae.bond_weights(m, m.body)
    seqT = ae.source_sequence(m, True, lag, steps)
    seqD = ae.source_sequence(m, False, lag, steps)
    pT, vT, ET = ae.linear_forward(m, wh, wv, seqT, store=True)
    pD, vD, ED = ae.linear_forward(m, wh, wv, seqD, store=True)
    S = ET + ED + 1e-30
    aT = 2.0 * ED / (S*S)
    aD = -2.0 * ET / (S*S)
    ehT, evT = ae.adjoint_grad(m, wh, wv, pT, vT, aT)
    ehD, evD = ae.adjoint_grad(m, wh, wv, pD, vD, aD)
    ex_h, ex_v = ehT + ehD, evT + evD
    gT = aT * np.asarray(pT[1:, m.soma[0], m.soma[1]], np.complex128)
    gD = aD * np.asarray(pD[1:, m.soma[0], m.soma[1]], np.complex128)
    return dict(wh=wh, wv=wv, pT=pT, pD=pD, ET=ET, ED=ED, C=float((ET-ED)/S),
                aT=float(aT), aD=float(aD), gT=gT, gD=gD, exact_h=ex_h, exact_v=ex_v)


def eval_pair(m, z, gT, gD):
    h, v = gradient_from_codes(m, z['wh'], z['wv'], z['pT'], z['pD'], gT, gD)
    return map_metrics(z['exact_h'], z['exact_v'], h, v)


def one(m, lag, steps):
    z = build_exact(m, lag, steps)
    gT, gD = z['gT'], z['gD']
    rows = {}

    rows['exact'] = eval_pair(m, z, gT, gD)

    # Keep temporal phase everywhere but erase amplitude modulation.
    pt = gT / (np.abs(gT) + 1e-30)
    pd = gD / (np.abs(gD) + 1e-30)
    rows['phase_only'] = eval_pair(m, z, l2_match(pt, gT), l2_match(pd, gD))

    # Keep only the magnitude envelope and one task-level sign.
    et = np.sign(z['aT']) * np.abs(gT)
    ed = np.sign(z['aD']) * np.abs(gD)
    rows['envelope_signed'] = eval_pair(m, z, l2_match(et, gT), l2_match(ed, gD))

    rows['real_wave'] = eval_pair(m, z, l2_match(np.real(gT), gT), l2_match(np.real(gD), gD))

    for N in SPARSE_NS:
        st = sparse_code(gT, N, 'phase', z['aT'])
        sd = sparse_code(gD, N, 'phase', z['aD'])
        rows[f'sparse_phase_{N}'] = eval_pair(m, z, st, sd)

        st = sparse_code(gT, N, 'signed', z['aT'])
        sd = sparse_code(gD, N, 'signed', z['aD'])
        rows[f'sparse_signed_{N}'] = eval_pair(m, z, st, sd)

        st = sparse_code(gT, N, 'positive', z['aT'])
        sd = sparse_code(gD, N, 'positive', z['aD'])
        rows[f'sparse_positive_{N}'] = eval_pair(m, z, st, sd)

    baseT = sparse_code(gT, 8, 'phase', z['aT'])
    baseD = sparse_code(gD, 8, 'phase', z['aD'])
    for d in DELAYS:
        rows[f'delayed_sparse8_{d}'] = eval_pair(m, z, shift_zero(baseT, d), shift_zero(baseD, d))

    gates = {}
    for P in GATE_PERIODS:
        q = []
        for off in range(P):
            gt, dutyT = gate_code(gT, P, off)
            gd, dutyD = gate_code(gD, P, off)
            mm = eval_pair(m, z, gt, gd)
            mm['offset'] = int(off)
            mm['duty'] = 0.5 * (dutyT + dutyD)
            q.append(mm)
        gates[str(P)] = dict(
            best=max(q, key=lambda x: x['corr']),
            worst=min(q, key=lambda x: x['corr']),
            median_corr=float(np.median([x['corr'] for x in q])),
            mean_corr=float(np.mean([x['corr'] for x in q])),
            median_sign=float(np.median([x['strong_sign_agreement'] for x in q])),
            all=q,
        )

    return dict(seed=int(m.cfg.seed), cells=int(m.body.sum()), C=z['C'], codes=rows, gates=gates)


def summary(rows):
    names = list(rows[0]['codes'])
    out = dict(bodies=len(rows), codes={}, gates={})
    for name in names:
        q = [r['codes'][name] for r in rows]
        out['codes'][name] = dict(
            mean_corr=float(np.mean([x['corr'] for x in q])),
            median_corr=float(np.median([x['corr'] for x in q])),
            mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
            mean_strong_sign=float(np.mean([x['strong_sign_agreement'] for x in q])),
            bodies_corr_gt_0p8=int(np.sum([x['corr'] > .8 for x in q])),
        )
    for P in map(str, GATE_PERIODS):
        q = [r['gates'][P] for r in rows]
        out['gates'][P] = dict(
            mean_best_corr=float(np.mean([x['best']['corr'] for x in q])),
            mean_median_corr=float(np.mean([x['median_corr'] for x in q])),
            mean_worst_corr=float(np.mean([x['worst']['corr'] for x in q])),
            mean_median_sign=float(np.mean([x['median_sign'] for x in q])),
        )
    return out


def selftest():
    g = np.asarray([1+1j, 0, 2-1j, 0, .2+.4j], np.complex128)
    for mode in ('phase', 'signed', 'positive'):
        z = sparse_code(g, 2, mode, -1)
        assert abs(np.linalg.norm(z)-np.linalg.norm(g)) < 1e-12
        assert np.count_nonzero(z) <= 2
    for d in (1, 2):
        z = shift_zero(g, d)
        assert abs(np.linalg.norm(z)-np.linalg.norm(g)) < 1e-12
    print('selftest ok')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=492)
    ap.add_argument('--seeds', type=int, default=4)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/biological_return_code/dev.json')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest(); return

    fa = Path(a.functional_arbors).resolve()
    sys.path.insert(0, str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor

    rows = []
    for seed in range(a.seed_start, a.seed_start+a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        b = m.bootstrap()
        if not b.get('ok'):
            continue
        m.mature = True
        r = one(m, a.lag, a.steps)
        rows.append(r)
        print('seed', seed,
              'phase', round(r['codes']['phase_only']['corr'], 3),
              'env', round(r['codes']['envelope_signed']['corr'], 3),
              'sp8', round(r['codes']['sparse_phase_8']['corr'], 3),
              'ss8', round(r['codes']['sparse_signed_8']['corr'], 3),
              'delay4', round(r['codes']['delayed_sparse8_4']['corr'], 3),
              flush=True)
    if not rows:
        raise SystemExit('No valid bodies')
    s = summary(rows)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(experiment='biological_return_code_dev_v01', summary=s, rows=rows), indent=2))
    print('\nBIOLOGICAL RETURN CODE DEV')
    print(json.dumps(s, indent=2))


if __name__ == '__main__':
    main()
