from __future__ import annotations

import argparse, json, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
import bond_response_probe as br
from structural_interference_probe import body_probe as event_probe


def summarize(events):
    if not events:
        return {'n': 0}
    interior=[e for e in events if 0 < e['alpha_best'] < 1]
    a0=[e for e in events if e['alpha_best'] == 0]
    a1=[e for e in events if e['alpha_best'] == 1]
    ia=np.asarray([e['alpha_best'] for e in interior],float)
    regret=np.asarray([e['regret_binary'] for e in events],float)
    return dict(
        n=len(events),
        interior_n=len(interior), interior_fraction=float(len(interior)/len(events)),
        alpha0_n=len(a0), alpha0_fraction=float(len(a0)/len(events)),
        alpha1_n=len(a1), alpha1_fraction=float(len(a1)/len(events)),
        median_alpha_best_interior=float(np.median(ia)) if len(ia) else float('nan'),
        mean_binary_regret=float(np.mean(regret)),
        regret_positive_fraction=float(np.mean(regret>1e-5)),
        slope_reversal_fraction=float(np.mean([e['slope_sign_reversal'] for e in events])),
        mean_base_derivative=float(np.mean([e['base_directional_derivative'] for e in events])),
    )


def body_control(m, lag, steps, max_each):
    observed=event_probe(m,lag,steps,max_each)
    seqT=ae.source_sequence(m,True,lag,steps)
    seqD=ae.source_sequence(m,False,lag,steps)
    base=ae.contrast_adjoint(m,m.body,seqT,seqD)
    dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    events=[]
    for e in observed['events']:
        if e['kind'] != 'add':
            continue
        deriv=float(ae.event_score(base['gh'],base['gv'],m.body,e,dk))
        curve=br.response_curve(m,base,seqT,seqD,e,dk)
        ib=int(np.argmax(curve))
        events.append(dict(
            cell=e['cell'], dist_soma=e['dist_soma'],
            base_directional_derivative=deriv,
            alpha_best=float(br.ALPHAS[ib]),
            best_gain=float(curve[ib]), dC_lin_full=float(curve[-1]),
            regret_binary=float(curve[ib]-curve[-1]),
            slope_sign_reversal=br.slope_sign_reversal(curve),
            alphas=list(br.ALPHAS), dC_lin_curve=[float(x) for x in curve],
        ))
    pos=[e for e in events if e['base_directional_derivative']>0]
    neg=[e for e in events if e['base_directional_derivative']<0]
    return dict(seed=int(m.cfg.seed), random_add=summarize(events), positive_derivative=summarize(pos),
                negative_derivative=summarize(neg), events=events)


def selftest():
    e=[dict(alpha_best=.1,regret_binary=.2,slope_sign_reversal=True,base_directional_derivative=1),
       dict(alpha_best=0.,regret_binary=0.,slope_sign_reversal=False,base_directional_derivative=-1)]
    s=summarize(e)
    assert s['interior_fraction']==.5 and s['alpha0_fraction']==.5
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=276)
    ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-each',type=int,default=6)
    ap.add_argument('--out',default='runs/random_frontier_control/random_frontier.json')
    ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest: selftest(); return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists(): raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True
        r=body_control(m,a.lag,a.steps,a.max_each); rows.append(r)
        s=r['random_add']
        print(f"seed {seed}: n={s['n']} interior={s.get('interior_fraction',float('nan')):.2f} a0={s.get('alpha0_fraction',float('nan')):.2f}",flush=True)
    events=[e for r in rows for e in r['events']]
    pos=[e for e in events if e['base_directional_derivative']>0]
    neg=[e for e in events if e['base_directional_derivative']<0]
    rnd=summarize(events); fav=summarize(pos); bad=summarize(neg)
    if rnd.get('n',0):
        if rnd['interior_fraction']>=.70: verdict='BROAD_GRADED'
        elif rnd['interior_fraction']<=.55 and rnd['alpha0_fraction']>=.30: verdict='SELECTION_DOMINATED'
        else: verdict='MIXED'
    else: verdict='NO_DATA'
    summary=dict(bodies=len(rows),random=rnd,positive_derivative=fav,negative_derivative=bad,verdict=verdict,
                 favored_minus_random_interior=float(fav.get('interior_fraction',float('nan'))-rnd.get('interior_fraction',float('nan'))) if fav.get('n',0) and rnd.get('n',0) else float('nan'))
    payload=dict(experiment='random_frontier_control_v01',prereg='RANDOM_FRONTIER_CONTROL_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,max_each=a.max_each,
                 summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nRANDOM FRONTIER CONTROL RECEIPT')
    for k,v in summary.items(): print(f' {k}: {v}')

if __name__=='__main__': main()
