"""Explicit visibility x lagged-complex-compatibility factorization.

See INTERFERENCE_FACTORIZATION_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from transfer_decomposition_probe import trace_single, tshift, safe_corr


def percentile(v, x):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) == 0 or not np.isfinite(x):
        return float('nan')
    return float(np.mean(v <= float(x)))


def sign_test_two_sided(vals):
    z = np.asarray(vals, float)
    z = z[np.isfinite(z) & (z != 0)]
    if len(z) == 0:
        return float('nan'), 0, 0
    w = int(np.sum(z > 0)); l = int(np.sum(z < 0)); n = w + l
    k = min(w, l)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return float(min(1.0, 2.0 * tail)), w, l


def body_probe(m, lag=20, steps=210):
    hA = trace_single(m, 0, 1.0, steps)
    hB = trace_single(m, 1, 1.0, steps)
    body = m.body.astype(bool)
    cells = np.argwhere(body)
    idx = tuple(cells.T)
    soma = tuple(map(int, m.soma))
    si = int(np.flatnonzero((cells[:,0] == soma[0]) & (cells[:,1] == soma[1]))[0])

    A = np.asarray(hA[:, idx[0], idx[1]], np.complex128)
    B = np.asarray(hB[:, idx[0], idx[1]], np.complex128)
    Al = tshift(A, lag)
    Bl = tshift(B, lag)

    # Historical peak-power order contrast.
    psiT = A + Bl
    psiD = B + Al
    pT = np.abs(psiT) ** 2
    pD = np.abs(psiD) ** 2
    peakT = np.max(pT, axis=0)
    peakD = np.max(pD, axis=0)
    Cpeak = (peakT - peakD) / (peakT + peakD + 1e-30)
    absCpeak = np.abs(Cpeak)

    # Peak amplitude-balance variable retained as a control.
    PApeak = np.max(np.abs(A) ** 2, axis=0)
    PBpeak = np.max(np.abs(B) ** 2, axis=0)
    peak_balance = np.minimum(PApeak, PBpeak) / (np.maximum(PApeak, PBpeak) + 1e-30)

    # Energy visibility.  Zero-padding conceptually preserves the full single-source
    # energies when one trace is delayed; the overlap sums below are the only cross terms.
    EA = np.sum(np.abs(A) ** 2, axis=0)
    EB = np.sum(np.abs(B) ** 2, axis=0)
    rootE = np.sqrt(EA * EB)
    visibility = 2.0 * rootE / (EA + EB + 1e-30)

    lag = int(lag)
    if lag <= 0:
        A_now, A_prev = A, A
        B_now, B_prev = B, B
    else:
        A_now, A_prev = A[lag:], A[:-lag]
        B_now, B_prev = B[lag:], B[:-lag]

    Rplus = np.sum(A_now * np.conj(B_prev), axis=0)
    Rminus = np.sum(B_now * np.conj(A_prev), axis=0)
    rho_plus = Rplus / (rootE + 1e-30)
    rho_minus = Rminus / (rootE + 1e-30)
    rp = np.real(rho_plus)
    rm = np.real(rho_minus)
    delta_rho = rp - rm

    denom = 2.0 + visibility * (rp + rm)
    Cint = visibility * delta_rho / (denom + 1e-30)
    absCint = np.abs(Cint)
    numerator = visibility * np.abs(delta_rho)

    # Explicit zero-padded integrated-energy construction verifies the identity.
    L = int(steps + max(lag, 0))
    AT = np.zeros((L, len(cells)), np.complex128)
    BT = np.zeros_like(AT)
    AD = np.zeros_like(AT)
    BD = np.zeros_like(AT)
    AT[:steps] = A
    BD[:steps] = B
    if lag <= 0:
        BT[:steps] = B
        AD[:steps] = A
    else:
        BT[lag:lag+steps] = B
        AD[lag:lag+steps] = A
    ET = np.sum(np.abs(AT + BT) ** 2, axis=0)
    ED = np.sum(np.abs(BD + AD) ** 2, axis=0)
    Cint_direct = (ET - ED) / (ET + ED + 1e-30)
    identity_mae = float(np.mean(np.abs(Cint_direct - Cint)))
    identity_rel = float(identity_mae / (np.mean(np.abs(Cint_direct)) + 1e-30))

    r_signed = safe_corr(Cint, Cpeak)
    r_abs = safe_corr(absCint, absCpeak)
    rV = safe_corr(visibility, absCpeak)
    rN = safe_corr(numerator, absCpeak)
    rVB = safe_corr(visibility, peak_balance)

    return dict(
        cells=int(len(cells)),
        soma=[int(soma[0]), int(soma[1])],
        integrated_identity_mae=identity_mae,
        integrated_identity_relative=identity_rel,
        r_signed_Cint_Cpeak=float(r_signed),
        r_abs_Cint_Cpeak=float(r_abs),
        r_visibility_absCpeak=float(rV),
        r_num_absCpeak=float(rN),
        delta_num_minus_visibility=float(rN-rV),
        r_visibility_peakbalance=float(rVB),
        soma_Cpeak=float(Cpeak[si]),
        soma_Cint=float(Cint[si]),
        soma_visibility=float(visibility[si]),
        soma_abs_delta_rho=float(abs(delta_rho[si])),
        soma_numerator=float(numerator[si]),
        soma_peak_balance=float(peak_balance[si]),
        soma_Cpeak_percentile=percentile(absCpeak, absCpeak[si]),
        soma_Cint_percentile=percentile(absCint, absCint[si]),
        soma_visibility_percentile=percentile(visibility, visibility[si]),
        soma_delta_rho_percentile=percentile(np.abs(delta_rho), abs(delta_rho[si])),
        soma_numerator_percentile=percentile(numerator, numerator[si]),
        descriptive=dict(
            mean_visibility=float(np.mean(visibility)),
            mean_abs_delta_rho=float(np.mean(np.abs(delta_rho))),
            mean_abs_Cint=float(np.mean(absCint)),
            mean_abs_Cpeak=float(np.mean(absCpeak)),
            mean_numerator=float(np.mean(numerator)),
        ),
    )


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--functional-arbors', default='FunctionalArbors')
    ap.add_argument('--seed-start', type=int, default=132)
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--steps', type=int, default=210)
    ap.add_argument('--out', default='runs/interference_factorization_discovery/interference_factorization.json')
    ap.add_argument('--selftest', action='store_true')
    return ap.parse_args()


def selftest():
    x = np.array([1+0j, 2+0j, 0+0j])[:,None]
    y = np.array([0+0j, 1+0j, 1+0j])[:,None]
    EA = np.sum(abs(x)**2,axis=0); EB=np.sum(abs(y)**2,axis=0)
    V=2*np.sqrt(EA*EB)/(EA+EB)
    assert np.all((V>=0)&(V<=1+1e-12))
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

    rows=[]
    for seed in range(a.seed_start, a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed))
        boot=m.bootstrap()
        if not boot.get('ok'):
            continue
        m.mature=True
        r=body_probe(m,a.lag,a.steps); r['seed']=int(seed); rows.append(r)
        print(f"seed {seed}: rSigned={r['r_signed_Cint_Cpeak']:+.3f} "
              f"rAbs={r['r_abs_Cint_Cpeak']:+.3f} rV={r['r_visibility_absCpeak']:+.3f} "
              f"rN={r['r_num_absCpeak']:+.3f} d={r['delta_num_minus_visibility']:+.3f} "
              f"soma pct V/drho/N/C={r['soma_visibility_percentile']:.2f}/"
              f"{r['soma_delta_rho_percentile']:.2f}/{r['soma_numerator_percentile']:.2f}/"
              f"{r['soma_Cpeak_percentile']:.2f}", flush=True)

    if not rows:
        raise SystemExit('No valid bodies')
    def vals(k): return np.asarray([r[k] for r in rows],float)
    rs=vals('r_signed_Cint_Cpeak'); ra=vals('r_abs_Cint_Cpeak')
    d=vals('delta_num_minus_visibility'); dp,dw,dl=sign_test_two_sided(d)
    sV=vals('soma_visibility_percentile')
    summary=dict(
        bodies=len(rows),
        mean_r_signed_Cint_Cpeak=float(np.nanmean(rs)),
        signed_positive_bodies=int(np.sum(rs>0)),
        mean_r_abs_Cint_Cpeak=float(np.nanmean(ra)),
        mean_r_visibility_absCpeak=float(np.nanmean(vals('r_visibility_absCpeak'))),
        mean_r_num_absCpeak=float(np.nanmean(vals('r_num_absCpeak'))),
        mean_delta_num_minus_visibility=float(np.nanmean(d)),
        improved_bodies=int(np.sum(d>0)),
        worse_bodies=int(np.sum(d<0)),
        improvement_sign_p=float(dp),
        median_soma_visibility_percentile=float(np.nanmedian(sV)),
        mean_soma_visibility_percentile=float(np.nanmean(sV)),
        mean_soma_delta_rho_percentile=float(np.nanmean(vals('soma_delta_rho_percentile'))),
        mean_soma_numerator_percentile=float(np.nanmean(vals('soma_numerator_percentile'))),
        mean_soma_Cpeak_percentile=float(np.nanmean(vals('soma_Cpeak_percentile'))),
        mean_r_visibility_peakbalance=float(np.nanmean(vals('r_visibility_peakbalance'))),
        integrated_identity_relative_mean=float(np.nanmean(vals('integrated_identity_relative'))),
    )
    summary['D1_pass']=bool(summary['mean_r_signed_Cint_Cpeak']>.65 and summary['signed_positive_bodies']>=10)
    summary['D2_pass']=bool(summary['mean_r_abs_Cint_Cpeak']>.60)
    summary['D3_pass']=bool(summary['mean_delta_num_minus_visibility']>.05 and summary['improved_bodies']>=9)
    summary['D4_pass']=bool(summary['median_soma_visibility_percentile']>.85)

    payload=dict(experiment='interference_factorization_discovery_v01',
                 prereg='INTERFERENCE_FACTORIZATION_DISCOVERY_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,
                 summary=summary,rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nINTERFERENCE FACTORIZATION DISCOVERY RECEIPT')
    for k,v in summary.items(): print(f' {k}: {v}')

if __name__=='__main__':
    main()
