"""Dense response curves for one structural conductance degree of freedom.

See BOND_RESPONSE_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # patches ae.weighted_lap to exact mature-boundary version
from structural_interference_probe import body_probe as nonlinear_event_probe
from transfer_decomposition_probe import safe_corr

ALPHAS = (
    0.0,
    0.005, 0.01, 0.02, 0.03, 0.05, 0.075,
    0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60,
    0.75, 1.00,
)


def response_curve(m, base, seqT, seqD, event, dk):
    vals=[]
    for alpha in ALPHAS:
        if alpha == 0:
            c=float(base['C'])
        else:
            wh,wv=ad.apply_event_fraction(base['wh'],base['wv'],m.body,event,dk,alpha)
            c,_,_=ae.linear_contrast(m,wh,wv,seqT,seqD)
        vals.append(float(c-base['C']))
    return np.asarray(vals,float)


def slope_sign_reversal(vals):
    vals=np.asarray(vals,float)
    a=np.asarray(ALPHAS,float)
    s=np.diff(vals)/np.diff(a)
    tol=max(1e-12, float(np.max(np.abs(s)))*1e-8 if len(s) else 1e-12)
    sg=np.sign(np.where(np.abs(s)>tol,s,0.0)).astype(int)
    nz=sg[sg!=0]
    if len(nz)<2:
        return False
    return bool(np.any(nz[1:] != nz[:-1]))


def summarize_favored(events, kind):
    es=[e for e in events if e['kind']==kind and e['base_directional_derivative']>0]
    if not es:
        return dict(n=0)
    interior=[e for e in es if 0 < e['alpha_best'] < 1]
    full_negative=[e for e in es if e['dC_lin_full']<0]
    regrets=np.asarray([e['regret_binary'] for e in es],float)
    interior_alphas=np.asarray([e['alpha_best'] for e in interior],float)
    return dict(
        n=len(es),
        interior_best_n=len(interior),
        interior_best_fraction=float(len(interior)/len(es)),
        full_negative_n=len(full_negative),
        full_negative_fraction=float(len(full_negative)/len(es)),
        mean_regret_binary=float(np.mean(regrets)),
        median_regret_binary=float(np.median(regrets)),
        regret_positive_n=int(np.sum(regrets>1e-5)),
        regret_positive_fraction=float(np.mean(regrets>1e-5)),
        median_alpha_best_interior=float(np.median(interior_alphas)) if len(interior_alphas) else float('nan'),
        mean_alpha_best=float(np.mean([e['alpha_best'] for e in es])),
        slope_reversal_n=int(np.sum([e['slope_sign_reversal'] for e in es])),
        slope_reversal_fraction=float(np.mean([e['slope_sign_reversal'] for e in es])),
        mean_best_gain=float(np.mean([e['best_gain'] for e in es])),
        mean_full_gain=float(np.mean([e['dC_lin_full'] for e in es])),
        pooled_corr_full_with_dCint=float(safe_corr([e['dC_lin_full'] for e in es],[e['dC_int'] for e in es])),
        pooled_corr_full_with_dCpeak=float(safe_corr([e['dC_lin_full'] for e in es],[e['dC_peak'] for e in es])),
    )


def body_probe(m,lag,steps,max_each):
    observed=nonlinear_event_probe(m,lag,steps,max_each)
    seqT=ae.source_sequence(m,True,lag,steps)
    seqD=ae.source_sequence(m,False,lag,steps)
    base=ae.contrast_adjoint(m,m.body,seqT,seqD)
    dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)

    events=[]
    for e in observed['events']:
        deriv=float(ae.event_score(base['gh'],base['gv'],m.body,e,dk))
        curve=response_curve(m,base,seqT,seqD,e,dk)
        ib=int(np.argmax(curve))
        best=float(curve[ib]); full=float(curve[-1])
        events.append(dict(
            kind=e['kind'],cell=e['cell'],dist_soma=e['dist_soma'],
            base_directional_derivative=deriv,
            alphas=list(ALPHAS),
            dC_lin_curve=[float(x) for x in curve],
            alpha_best=float(ALPHAS[ib]),
            best_gain=best,
            dC_lin_full=full,
            regret_binary=float(best-full),
            slope_sign_reversal=slope_sign_reversal(curve),
            dC_int=float(e['dCint']),
            dC_peak=float(e['dCpeak']),
        ))

    return dict(
        seed=int(m.cfg.seed),cells=int(m.body.sum()),n_events=len(events),
        base_C_lin=float(base['C']),base_C_int=float(observed['base']['Cint']),base_C_peak=float(observed['base']['Cpeak']),
        favored_add=summarize_favored(events,'add'),
        favored_delete=summarize_favored(events,'delete'),
        events=events,
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=216)
    ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-each',type=int,default=6)
    ap.add_argument('--out',default='runs/bond_response_discovery/bond_response.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    y=np.array([0,.1,.2,.15,.1])
    assert slope_sign_reversal(y)
    y=np.array([0,.1,.2,.3,.4])
    assert not slope_sign_reversal(y)
    print('selftest ok')


def main():
    a=parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor

    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,a.lag,a.steps,a.max_each);rows.append(r)
        fa0=r['favored_add'];fd0=r['favored_delete']
        print(f"seed {seed}: fav add={fa0.get('n',0)} int={fa0.get('interior_best_fraction',float('nan')):.2f} "
              f"full-={fa0.get('full_negative_fraction',float('nan')):.2f} regret={fa0.get('mean_regret_binary',float('nan')):.4f}; "
              f"fav del={fd0.get('n',0)}",flush=True)

    events=[e for r in rows for e in r['events']]
    add=summarize_favored(events,'add');dele=summarize_favored(events,'delete')
    summary=dict(bodies=len(rows),total_events=len(events),alphas=list(ALPHAS),
                 favored_add=add,favored_delete=dele)
    summary['D1_pass']=bool(add.get('n',0)>0 and add['interior_best_fraction']>.70)
    summary['D2_pass']=bool(add.get('n',0)>0 and add['full_negative_fraction']>.40)
    summary['D3_pass']=bool(add.get('n',0)>0 and add['mean_regret_binary']>.005 and add['regret_positive_fraction']>=.75)
    summary['D4_pass']=bool(add.get('n',0)>0 and np.isfinite(add['median_alpha_best_interior']) and add['median_alpha_best_interior']>=.03)

    payload=dict(experiment='bond_response_discovery_v01',prereg='BOND_RESPONSE_DISCOVERY_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,max_each=a.max_each,
                 summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nBOND RESPONSE DISCOVERY RECEIPT')
    for k,v in summary.items():print(f' {k}: {v}')

if __name__=='__main__':main()
