"""Exploratory test of the constant graph mode as a quadratic reference/mixer.

This follows the confirmed mode-pair result. It is exploratory on reused bodies;
any positive quantitative prediction must be frozen on fresh seeds later.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from graph_mode_probe import graph_laplacian_modes
from transfer_decomposition_probe import trace_single, tshift


def contrast(pT,pD):
    a=float(np.max(np.asarray(pT,float))); b=float(np.max(np.asarray(pD,float)))
    return float((a-b)/(a+b+1e-30))


def safe_corr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)<3 or np.std(a)<1e-14 or np.std(b)<1e-14: return float('nan')
    return float(np.corrcoef(a,b)[0,1])


def one(m,lag,steps):
    hA=trace_single(m,0,1.0,steps); hB=trace_single(m,1,1.0,steps)
    coords,evals,V=graph_laplacian_modes(m.body)
    idx={p:i for i,p in enumerate(coords)}; soma=tuple(map(int,m.soma)); si=idx[soma]
    ys=np.array([p[0] for p in coords]); xs=np.array([p[1] for p in coords])
    QA=np.asarray(hA[:,ys,xs],complex)@V; QB=np.asarray(hB[:,ys,xs],complex)@V
    phis=V[si]
    uA=QA*phis[None,:]; uB=QB*phis[None,:]
    UT=uA+tshift(uB,lag); UD=uB+tshift(uA,lag)
    cT=UT[:,0]; cD=UD[:,0]
    rT=UT[:,1:].sum(axis=1); rD=UD[:,1:].sum(axis=1)
    pfullT=np.abs(cT+rT)**2; pfullD=np.abs(cD+rD)**2
    pcommonT=np.abs(cT)**2; pcommonD=np.abs(cD)**2
    presidT=np.abs(rT)**2; presidD=np.abs(rD)**2
    pmixT=pcommonT+2*np.real(cT*np.conj(rT))
    pmixD=pcommonD+2*np.real(cD*np.conj(rD))
    Cfull=contrast(pfullT,pfullD); C0=contrast(pcommonT,pcommonD)
    Cres=contrast(presidT,presidD); Cmix=contrast(pmixT,pmixD)
    # amplitude dominance measured before summing sources: modal soma amplitudes
    source_common=np.concatenate([np.abs(uA[:,0]),np.abs(uB[:,0])])
    source_resid=np.concatenate([np.abs(uA[:,1:].sum(axis=1)),np.abs(uB[:,1:].sum(axis=1))])
    return dict(
        cells=len(coords), C_full=Cfull, C_common_only=C0, C_residual_only=Cres,
        C_first_order_common_mixer=Cmix,
        abs_error_mixer=float(abs(Cmix-Cfull)),
        absC_full=float(abs(Cfull)), absC_common=float(abs(C0)),
        absC_residual=float(abs(Cres)), absC_mixer=float(abs(Cmix)),
        sign_mixer_matches=bool(np.sign(Cmix)==np.sign(Cfull) or Cfull==0),
        common_to_residual_mean_amplitude_ratio=float(np.mean(source_common)/(np.mean(source_resid)+1e-30)),
        common_source_amplitude_mean=float(np.mean(source_common)),
        residual_source_amplitude_mean=float(np.mean(source_resid)),
        residual_power_fraction_at_full_target_peak=float(
            presidT[int(np.argmax(pfullT))]/(pfullT[int(np.argmax(pfullT))]+1e-30)),
        residual_power_fraction_at_full_distractor_peak=float(
            presidD[int(np.argmax(pfullD))]/(pfullD[int(np.argmax(pfullD))]+1e-30)),
    )


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=72); ap.add_argument('--seeds',type=int,default=24)
    ap.add_argument('--lag',type=int,default=20); ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/common_mode_mixer/common_mode_mixer.json'); ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest:
        assert abs(contrast([1,3],[1,1])-.5)<1e-12; print('selftest ok'); return
    fa=Path(a.functional_arbors).resolve(); sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True; r=one(m,a.lag,a.steps); r['seed']=seed; rows.append(r)
        print(f"seed {seed}: C full/common/resid/mix {r['C_full']:+.3f}/{r['C_common_only']:+.3f}/{r['C_residual_only']:+.3f}/{r['C_first_order_common_mixer']:+.3f} err={r['abs_error_mixer']:.3f}",flush=True)
    def v(k): return np.array([r[k] for r in rows],float)
    summary=dict(
        bodies=len(rows),
        mean_C_full=float(v('C_full').mean()), mean_absC_full=float(v('absC_full').mean()),
        mean_absC_common=float(v('absC_common').mean()), mean_absC_residual=float(v('absC_residual').mean()),
        mean_absC_mixer=float(v('absC_mixer').mean()),
        mixer_signed_corr_full=safe_corr(v('C_first_order_common_mixer'),v('C_full')),
        residual_signed_corr_full=safe_corr(v('C_residual_only'),v('C_full')),
        mixer_mean_abs_error=float(v('abs_error_mixer').mean()),
        mixer_median_abs_error=float(np.median(v('abs_error_mixer'))),
        mixer_sign_matches=int(sum(r['sign_mixer_matches'] for r in rows)),
        mean_common_to_residual_amplitude_ratio=float(v('common_to_residual_mean_amplitude_ratio').mean()),
        mean_residual_power_fraction_target_peak=float(v('residual_power_fraction_at_full_target_peak').mean()),
        mean_residual_power_fraction_distractor_peak=float(v('residual_power_fraction_at_full_distractor_peak').mean()),
    )
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(dict(experiment='common_mode_mixer_exploratory',seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,summary=summary,rows=rows),indent=2),encoding='utf-8')
    print('\nCOMMON-MODE MIXER EXPLORATORY RECEIPT')
    print(json.dumps(summary,indent=2)); print('wrote',out)

if __name__=='__main__': main()
