"""Counterfactual one-cell structural events in interference coordinates.

See STRUCTURAL_INTERFERENCE_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import deque
from pathlib import Path

import numpy as np

from transfer_decomposition_probe import safe_corr


def n4(y, x, shape):
    h, w = shape
    for dy, dx in ((-1,0),(1,0),(0,-1),(0,1)):
        yy, xx = y + dy, x + dx
        if 0 <= yy < h and 0 <= xx < w:
            yield yy, xx


def tshift1(x, lag):
    x = np.asarray(x)
    out = np.zeros_like(x)
    lag = int(lag)
    if lag <= 0:
        out[:] = x
    elif lag < len(x):
        out[lag:] = x[:-lag]
    return out


def trace_soma(m, which, steps):
    m.reset_fast(True)
    out = np.zeros(int(steps), np.complex128)
    for t in range(int(steps)):
        src = m.pulse_source(int(which), t, False)
        m.advance(src, False, True, 'none')
        out[t] = complex(m.psi[m.soma])
    return out


def f_interference(V, rp, rm):
    return float(V * (rp-rm) / (2.0 + V*(rp+rm) + 1e-30))


def soma_metrics(m, lag, steps):
    A = trace_soma(m, 0, steps)
    B = trace_soma(m, 1, steps)
    Al = tshift1(A, lag)
    Bl = tshift1(B, lag)

    pT = np.abs(A + Bl) ** 2
    pD = np.abs(B + Al) ** 2
    mxT = float(np.max(pT)); mxD = float(np.max(pD))
    Cpeak = float((mxT-mxD)/(mxT+mxD+1e-30))

    EA = float(np.sum(np.abs(A)**2)); EB = float(np.sum(np.abs(B)**2))
    rootE = math.sqrt(max(EA*EB, 0.0))
    V = float(2.0*rootE/(EA+EB+1e-30))

    lag = int(lag)
    if lag <= 0:
        Anow, Aprev = A, A
        Bnow, Bprev = B, B
    else:
        Anow, Aprev = A[lag:], A[:-lag]
        Bnow, Bprev = B[lag:], B[:-lag]
    Rplus = np.sum(Anow * np.conj(Bprev))
    Rminus = np.sum(Bnow * np.conj(Aprev))
    rho_plus = Rplus/(rootE+1e-30)
    rho_minus = Rminus/(rootE+1e-30)
    rp = float(np.real(rho_plus)); rm = float(np.real(rho_minus))
    Cint = f_interference(V, rp, rm)

    return dict(
        Cpeak=Cpeak,
        Cint=Cint,
        V=V,
        rp=rp,
        rm=rm,
        abs_delta_rho=float(abs(rp-rm)),
        energy_A=EA,
        energy_B=EB,
    )


def connected_from_soma(m, body):
    b = np.asarray(body, bool)
    seen = np.zeros_like(b, bool)
    soma = tuple(map(int, m.soma))
    if not b[soma]: return seen
    q = deque([soma]); seen[soma] = True
    while q:
        p = q.popleft()
        for r in n4(*p, b.shape):
            if b[r] and not seen[r]:
                seen[r] = True; q.append(r)
    return seen


def graph_distances(m, body):
    b = np.asarray(body, bool)
    soma = tuple(map(int, m.soma))
    dist = {soma: 0}
    q = deque([soma])
    while q:
        p = q.popleft()
        for r in n4(*p, b.shape):
            if b[r] and r not in dist:
                dist[r] = dist[p] + 1; q.append(r)
    return dist


def event_candidates(m, max_each, rng):
    body = m.body.astype(bool)
    shape = body.shape
    term0 = m.source_terminal(0, body)
    term1 = m.source_terminal(1, body)
    target_union = np.any(m.target_masks, axis=0)

    additions = []
    for yy, xx in np.argwhere((~body) & m.window & (~target_union)):
        p = (int(yy), int(xx))
        occn = sum(bool(body[q]) for q in n4(*p, shape))
        if occn != 1:
            continue
        test = body.copy(); test[p] = True
        if m.source_terminal(0, test) != term0 or m.source_terminal(1, test) != term1:
            continue
        additions.append(p)

    deletions = []
    for yy, xx in np.argwhere(body & (~target_union)):
        p = (int(yy), int(xx))
        if p == tuple(map(int, m.soma)) or bool(m.protect[p]):
            continue
        test = body.copy(); test[p] = False
        if int(connected_from_soma(m, test).sum()) != int(test.sum()):
            continue
        if not m.both_connected(test):
            continue
        if m.source_terminal(0, test) != term0 or m.source_terminal(1, test) != term1:
            continue
        deletions.append(p)

    def choose(seq):
        seq = list(seq)
        if len(seq) <= int(max_each):
            return seq
        ii = rng.choice(len(seq), int(max_each), replace=False)
        return [seq[int(i)] for i in sorted(ii)]

    return choose(additions), choose(deletions), term0, term1


def shapley(base, after):
    V0, V1 = float(base['V']), float(after['V'])
    R0 = (float(base['rp']), float(base['rm']))
    R1 = (float(after['rp']), float(after['rm']))
    f00 = f_interference(V0, *R0)
    f10 = f_interference(V1, *R0)
    f01 = f_interference(V0, *R1)
    f11 = f_interference(V1, *R1)
    cV = 0.5*((f10-f00) + (f11-f01))
    cR = 0.5*((f01-f00) + (f11-f10))
    return float(cV), float(cR), float((cV+cR) - (f11-f00))


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w,l)
    tail = sum(math.comb(n,i) for i in range(k+1))/(2**n)
    return float(min(1.0, 2.0*tail)), w, l


def body_probe(m, lag, steps, max_each=6):
    base = soma_metrics(m, lag, steps)
    rng = np.random.default_rng(int(m.cfg.seed) + 156_156)
    adds, dels, term0, term1 = event_candidates(m, max_each, rng)
    dist = graph_distances(m, m.body)

    events = []
    for kind, seq in (('add', adds), ('delete', dels)):
        for p in seq:
            z = m.copy()
            z.body = m.body.copy()
            z.body[p] = 1 if kind == 'add' else 0
            z.mature = True
            # Exact source terminals are part of the registered event definition.
            if z.source_terminal(0) != term0 or z.source_terminal(1) != term1:
                continue
            after = soma_metrics(z, lag, steps)
            cV, cR, sherr = shapley(base, after)
            if kind == 'add':
                nd = [dist[q] for q in n4(*p, m.body.shape) if q in dist]
                ds = (min(nd)+1) if nd else -1
            else:
                ds = int(dist.get(p, -1))
            events.append(dict(
                kind=kind,
                cell=[int(p[0]),int(p[1])],
                dist_soma=int(ds),
                dCpeak=float(after['Cpeak']-base['Cpeak']),
                dCint=float(after['Cint']-base['Cint']),
                dV=float(after['V']-base['V']),
                d_abs_delta_rho=float(after['abs_delta_rho']-base['abs_delta_rho']),
                drp=float(after['rp']-base['rp']),
                drm=float(after['rm']-base['rm']),
                contrib_V=cV,
                contrib_R=cR,
                shapley_reconstruction_error=sherr,
                after=after,
            ))

    dcpeak = np.asarray([e['dCpeak'] for e in events], float)
    dcint = np.asarray([e['dCint'] for e in events], float)
    dv = np.asarray([e['dV'] for e in events], float)
    rI = safe_corr(dcint, dcpeak)
    rV = safe_corr(dv, dcpeak)
    return dict(
        seed=int(m.cfg.seed),
        cells=int(m.body.sum()),
        base=base,
        source_terminals=[list(map(int,term0)),list(map(int,term1))],
        candidate_counts=dict(add=len(adds), delete=len(dels)),
        events=events,
        n_events=len(events),
        r_dCint_dCpeak=float(rI),
        r_dV_dCpeak=float(rV),
        delta_r_full_minus_visibility=float(rI-rV) if np.isfinite(rI) and np.isfinite(rV) else float('nan'),
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=156)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--max-each', type=int, default=6)
    ap.add_argument('--out', default='runs/structural_interference_discovery/structural_interference.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    b=dict(V=.8,rp=.2,rm=-.1,Cint=f_interference(.8,.2,-.1))
    a=dict(V=.7,rp=.25,rm=-.08,Cint=f_interference(.7,.25,-.08))
    cV,cR,e=shapley(b,a)
    assert abs(cV+cR-(a['Cint']-b['Cint'])) < 1e-12
    assert abs(e)<1e-12
    print('selftest ok')


def main():
    a=parse_args()
    if a.selftest:
        selftest(); return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists(): raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor

    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed))
        boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True
        r=body_probe(m,a.lag,a.steps,a.max_each); rows.append(r)
        print(f"seed {seed}: events={r['n_events']} rInt={r['r_dCint_dCpeak']:+.3f} "
              f"rV={r['r_dV_dCpeak']:+.3f} d={r['delta_r_full_minus_visibility']:+.3f}",flush=True)

    if not rows: raise SystemExit('No valid bodies')
    valid=[r for r in rows if np.isfinite(r['r_dCint_dCpeak'])]
    rI=np.asarray([r['r_dCint_dCpeak'] for r in valid],float)
    deltas=np.asarray([r['delta_r_full_minus_visibility'] for r in valid if np.isfinite(r['delta_r_full_minus_visibility'])],float)
    dp,dw,dl=sign_test_two_sided(deltas)
    events=[e for r in rows for e in r['events']]
    nontriv=[e for e in events if abs(e['dCpeak'])>1e-5]
    signmatch=[np.sign(e['dCint'])==np.sign(e['dCpeak']) for e in nontriv]
    cV=np.asarray([abs(e['contrib_V']) for e in events],float)
    cR=np.asarray([abs(e['contrib_R']) for e in events],float)
    sh=np.asarray([abs(e['shapley_reconstruction_error']) for e in events],float)

    kind_summary={}
    for kind in ('add','delete'):
        es=[e for e in events if e['kind']==kind]
        if es:
            kind_summary[kind]=dict(
                n=len(es),
                mean_abs_dCpeak=float(np.mean([abs(e['dCpeak']) for e in es])),
                mean_abs_dCint=float(np.mean([abs(e['dCint']) for e in es])),
                mean_abs_contrib_V=float(np.mean([abs(e['contrib_V']) for e in es])),
                mean_abs_contrib_R=float(np.mean([abs(e['contrib_R']) for e in es])),
            )

    summary=dict(
        bodies=len(rows),
        valid_correlation_bodies=len(valid),
        total_events=len(events),
        mean_r_dCint_dCpeak=float(np.nanmean(rI)),
        median_r_dCint_dCpeak=float(np.nanmedian(rI)),
        positive_r_bodies=int(np.sum(rI>0)),
        pooled_nontrivial_events=len(nontriv),
        pooled_sign_agreement=float(np.mean(signmatch)) if signmatch else float('nan'),
        mean_r_improvement_full_over_visibility=float(np.nanmean(deltas)) if len(deltas) else float('nan'),
        improved_bodies=int(np.sum(deltas>0)),
        worse_bodies=int(np.sum(deltas<0)),
        improvement_sign_p=float(dp),
        mean_abs_contrib_V=float(np.mean(cV)),
        mean_abs_contrib_R=float(np.mean(cR)),
        compatibility_to_visibility_abs_ratio=float(np.sum(cR)/(np.sum(cV)+1e-30)),
        mean_abs_shapley_error=float(np.mean(sh)),
        max_abs_shapley_error=float(np.max(sh)),
        kind_summary=kind_summary,
    )
    summary['D1_pass']=bool(summary['mean_r_dCint_dCpeak']>.50 and summary['positive_r_bodies']>=10)
    summary['D2_pass']=bool(np.isfinite(summary['pooled_sign_agreement']) and summary['pooled_sign_agreement']>.65)
    summary['D3_pass']=bool(summary['mean_r_improvement_full_over_visibility']>.10 and summary['improved_bodies']>=9)

    payload=dict(experiment='structural_interference_discovery_v01',
                 prereg='STRUCTURAL_INTERFERENCE_DISCOVERY_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,
                 max_each=a.max_each,summary=summary,rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nSTRUCTURAL INTERFERENCE DISCOVERY RECEIPT')
    for k,v in summary.items(): print(f' {k}: {v}')

if __name__=='__main__':
    main()
