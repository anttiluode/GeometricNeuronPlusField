"""Evaluate the frozen HCN_QUASIACTIVE_CONFIRM_PREREG_V01 criteria."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

FREQS=(0.03,0.04)


def load_candidate(payload):
    cs=payload.get('candidates',[])
    if len(cs)!=1:
        raise SystemExit(f'Expected exactly one frozen candidate, found {len(cs)}')
    c=cs[0]
    wanted=(0.005,10.0,2.0,0.5)
    got=(float(c['g0']),float(c['ratio']),float(c['tau_h']),float(c['mu_ratio']))
    if any(abs(a-b)>1e-12 for a,b in zip(got,wanted)):
        raise SystemExit(f'Frozen candidate mismatch: got {got}, expected {wanted}')
    return c


def observations(c):
    out=[]
    for r in c['rows']:
        seed=int(r['seed'])
        for q in r['byfreq']:
            om=float(q['omega'])
            if not any(abs(om-x)<1e-12 for x in FREQS):
                continue
            sp=q['soma_phase_rms']; lp=q['local_phase_rms']; amp=q['soma_amp_median']
            out.append(dict(
                seed=seed,omega=om,
                gain_shuffle=float(sp['shuffle']-sp['smooth']),
                gain_uniform=float(sp['uniform']-sp['smooth']),
                gain_reverse=float(sp['reverse']-sp['smooth']),
                lead=float(q['advance_vs_zero']['distal_minus_prox']),
                local_retention=float(lp['smooth']/(lp['zero']+1e-30)),
                amp_ratio_shuffle=float(amp['smooth']/(amp['shuffle']+1e-30)),
                amp_ratio_uniform=float(amp['smooth']/(amp['uniform']+1e-30)),
            ))
    return out


def eval_receipt(payload):
    c=load_candidate(payload); obs=observations(c)
    if not obs: raise SystemExit('No frozen-frequency observations found')
    n=len(obs); frac=.75
    def arr(k):return np.asarray([z[k] for z in obs],float)
    gs=arr('gain_shuffle');gu=arr('gain_uniform');gr=arr('gain_reverse')
    lead=arr('lead');lr=arr('local_retention');ars=arr('amp_ratio_shuffle');aru=arr('amp_ratio_uniform')

    criteria={}
    criteria['Q0_placement_advantage']={
        'pass_':bool(np.mean(gs)>.05 and np.mean(gu)>.05),
        'mean_gain_shuffle':float(np.mean(gs)),'mean_gain_uniform':float(np.mean(gu))}
    criteria['Q1_robust_observations']={
        'pass_':bool(np.mean(gs>0)>=frac and np.mean(gu>0)>=frac),
        'gain_shuffle_positive':int(np.sum(gs>0)),'gain_uniform_positive':int(np.sum(gu>0)),
        'observations':n,'fraction_shuffle':float(np.mean(gs>0)),'fraction_uniform':float(np.mean(gu>0))}
    criteria['Q2_hcn_distance_sign']={
        'pass_':bool(np.mean(lead)>.30 and np.mean(lead>0)>=frac),
        'mean_lead':float(np.mean(lead)),'positive':int(np.sum(lead>0)),'fraction_positive':float(np.mean(lead>0))}
    criteria['Q3_local_phase_retention']={
        'pass_':bool(np.mean(lr)>.70),'mean_local_retention':float(np.mean(lr)),
        'median_local_retention':float(np.median(lr))}
    criteria['Q4_amplitude_control']={
        'pass_':bool(.5<np.median(ars)<2.0 and .5<np.median(aru)<2.0),
        'median_amp_ratio_shuffle':float(np.median(ars)),
        'median_amp_ratio_uniform':float(np.median(aru))}

    per_freq={}
    q5=True
    for om in FREQS:
        z=[x for x in obs if abs(x['omega']-om)<1e-12]
        f=dict(
            observations=len(z),
            mean_gain_shuffle=float(np.mean([x['gain_shuffle'] for x in z])),
            mean_gain_uniform=float(np.mean([x['gain_uniform'] for x in z])),
            mean_lead=float(np.mean([x['lead'] for x in z])),
        )
        f['pass_']=bool(f['mean_gain_shuffle']>.02 and f['mean_gain_uniform']>.02 and f['mean_lead']>.10)
        q5=q5 and f['pass_'];per_freq[str(om)]=f
    criteria['Q5_both_frequencies']={'pass_':bool(q5),'frequency':per_freq}
    criteria['Q6_gradient_direction']={
        'pass_':bool(np.mean(gr)>.05),'mean_reverse_minus_smooth':float(np.mean(gr))}

    passed=sum(int(v['pass_']) for v in criteria.values())
    summary=dict(
        bodies=len(set(z['seed'] for z in obs)),observations=n,
        candidate=dict(g0=.005,ratio=10.0,tau_h=2.0,mu_ratio=.5,omegas=list(FREQS)),
        criteria=criteria,passed=passed,total=len(criteria),all_pass=bool(passed==len(criteria)),
        pooled=dict(
            mean_gain_shuffle=float(np.mean(gs)),mean_gain_uniform=float(np.mean(gu)),
            mean_gain_reverse=float(np.mean(gr)),mean_lead=float(np.mean(lead)),
            mean_local_retention=float(np.mean(lr)),median_amp_ratio_shuffle=float(np.median(ars)),
            median_amp_ratio_uniform=float(np.median(aru))),
        observations_detail=obs,
    )
    return summary


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--receipt',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    payload=json.loads(Path(a.receipt).read_text(encoding='utf-8'))
    s=eval_receipt(payload)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(s,indent=2),encoding='utf-8')
    print(json.dumps(s,indent=2))
    if not s['all_pass']:
        raise SystemExit('Frozen quasi-active HCN confirmation failed one or more preregistered criteria')

if __name__=='__main__':main()
