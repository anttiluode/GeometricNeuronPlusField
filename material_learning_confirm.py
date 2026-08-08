"""Evaluate MATERIAL_LEARNING_CONFIRM_PREREG_V01 on one frozen receipt."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np

FREQS=(0.03,0.04)


def evaluate(p):
    rows=p['rows'];n=len(rows)
    if n<1:raise SystemExit('No bodies in receipt')
    frac=.75
    def cval(r,name,key):return float(r['controls'][name][key])
    gu=np.array([cval(r,'learned','coherence2')-cval(r,'uniform','coherence2') for r in rows])
    gs=np.array([cval(r,'learned','coherence2')-cval(r,'shuffle_learned','coherence2') for r in rows])
    gh=np.array([cval(r,'learned','coherence2')-cval(r,'hand','coherence2') for r in rows])
    rho=np.array([float(r['morphology']['spearman_distance']) for r in rows])
    ar=np.array([cval(r,'learned','mean_median_amp')/(cval(r,'uniform','mean_median_amp')+1e-30) for r in rows])
    phasegain=np.array([cval(r,'uniform','mean_phase_rms')-cval(r,'learned','mean_phase_rms') for r in rows])
    graderr=max(float(r['gradient_check']['max_relative_error']) for r in rows)
    crit={}
    crit['L0_gradient_audit']={'pass_':bool(graderr<1e-5),'max_relative_error':graderr}
    crit['L1_learns_objective']={'pass_':bool(gu.mean()>.05 and np.mean(gu>0)>=frac),'mean_gain':float(gu.mean()),'positive':int(np.sum(gu>0)),'fraction_positive':float(np.mean(gu>0))}
    crit['L2_beats_shuffled']={'pass_':bool(gs.mean()>.08 and np.mean(gs>0)>=frac),'mean_gain':float(gs.mean()),'positive':int(np.sum(gs>0)),'fraction_positive':float(np.mean(gs>0))}
    crit['L3_beats_hand']={'pass_':bool(gh.mean()>.025 and np.mean(gh>0)>=frac),'mean_gain':float(gh.mean()),'positive':int(np.sum(gh>0)),'fraction_positive':float(np.mean(gh>0))}
    crit['L4_distance_emerges']={'pass_':bool(rho.mean()>.40 and np.mean(rho>0)>=frac),'mean_spearman':float(rho.mean()),'positive':int(np.sum(rho>0)),'fraction_positive':float(np.mean(rho>0)),'min_spearman':float(np.min(rho))}
    crit['L5_amplitude_sane']={'pass_':bool(.40<np.median(ar)<1.50 and np.mean((ar>.30)&(ar<2.0))>=frac),'median_ratio':float(np.median(ar)),'mean_ratio':float(np.mean(ar)),'in_range':int(np.sum((ar>.30)&(ar<2.0))),'fraction_in_range':float(np.mean((ar>.30)&(ar<2.0)))}
    per={};ok=True
    for om in FREQS:
        d=[]
        for r in rows:
            lu={float(x['omega']):float(x['coherence2']) for x in r['controls']['learned']['frequency']}
            uu={float(x['omega']):float(x['coherence2']) for x in r['controls']['uniform']['frequency']}
            d.append(lu[om]-uu[om])
        m=float(np.mean(d));per[str(om)]={'mean_gain':m,'pass_':bool(m>.02)};ok=ok and m>.02
    crit['L6_both_frequencies']={'pass_':bool(ok),'frequency':per}
    crit['L7_phase_rms_improves']={'pass_':bool(phasegain.mean()>.08),'mean_phase_rms_gain':float(phasegain.mean()),'positive':int(np.sum(phasegain>0))}
    passed=sum(int(x['pass_']) for x in crit.values())
    return dict(prereg='MATERIAL_LEARNING_CONFIRM_PREREG_V01.md',bodies=n,criteria=crit,passed=passed,total=len(crit),all_pass=bool(passed==len(crit)),pooled=dict(mean_gain_uniform=float(gu.mean()),mean_gain_shuffle=float(gs.mean()),mean_gain_hand=float(gh.mean()),mean_spearman=float(rho.mean()),median_amp_ratio=float(np.median(ar)),mean_phase_rms_gain=float(phasegain.mean())))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--receipt',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    p=json.loads(Path(a.receipt).read_text(encoding='utf-8'));s=evaluate(p);Path(a.out).parent.mkdir(parents=True,exist_ok=True);Path(a.out).write_text(json.dumps(s,indent=2),encoding='utf-8');print(json.dumps(s,indent=2))
    if not s['all_pass']:raise SystemExit('Held-out material-learning confirmation failed one or more frozen criteria')
if __name__=='__main__':main()
