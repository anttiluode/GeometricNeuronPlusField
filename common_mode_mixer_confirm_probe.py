"""Held-out confirmation of the constant-mode reference/mixer hypothesis.
See COMMON_MODE_MIXER_CONFIRM_PREREG_V01.md.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
from common_mode_mixer_probe import one,safe_corr


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=96); ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20); ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/common_mode_mixer_confirm/common_mode_mixer_confirm.json')
    ap.add_argument('--selftest',action='store_true'); a=ap.parse_args()
    if a.selftest:
        print('selftest ok'); return
    fa=Path(a.functional_arbors).resolve(); sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True; r=one(m,a.lag,a.steps); r['seed']=seed; rows.append(r)
        print(f"seed {seed}: full={r['C_full']:+.3f} ref={r['C_first_order_common_mixer']:+.3f} residual={r['C_residual_only']:+.3f} common={r['C_common_only']:+.5f}",flush=True)
    if not rows: raise SystemExit('No valid bodies')
    full=np.array([r['C_full'] for r in rows],float)
    ref=np.array([r['C_first_order_common_mixer'] for r in rows],float)
    res=np.array([r['C_residual_only'] for r in rows],float)
    common=np.array([r['C_common_only'] for r in rows],float)
    eref=np.abs(ref-full); eres=np.abs(res-full)
    corr_ref=safe_corr(ref,full); corr_res=safe_corr(res,full)
    signmatch=int(np.sum((np.sign(ref)==np.sign(full)) | (full==0)))
    wins=int(np.sum(eref<eres)); losses=int(np.sum(eref>eres))
    P1=bool(np.mean(np.abs(common))<.001)
    P2=bool(corr_ref>.95 and np.mean(eref)<.08 and signmatch>=10)
    P3=bool(np.mean(eref)<np.mean(eres) and wins>=10)
    summary=dict(
        bodies=len(rows),
        mean_absC_full=float(np.mean(np.abs(full))),
        mean_absC_common=float(np.mean(np.abs(common))),
        mean_absC_ref=float(np.mean(np.abs(ref))),
        mean_absC_residual=float(np.mean(np.abs(res))),
        corr_ref_full=float(corr_ref), corr_residual_full=float(corr_res),
        corr_advantage_ref_over_residual=float(corr_ref-corr_res),
        mean_abs_error_ref=float(np.mean(eref)), mean_abs_error_residual=float(np.mean(eres)),
        ref_better_error_bodies=wins, residual_better_error_bodies=losses,
        ref_sign_matches=signmatch,
        mean_common_to_residual_amplitude_ratio=float(np.mean([r['common_to_residual_mean_amplitude_ratio'] for r in rows])),
        registered=dict(P1=P1,P2=P2,P3=P3,all_pass=bool(P1 and P2 and P3)),
    )
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(dict(experiment='common_mode_mixer_confirm_v01',prereg='COMMON_MODE_MIXER_CONFIRM_PREREG_V01.md',seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,summary=summary,rows=rows),indent=2),encoding='utf-8')
    print('\nCOMMON-MODE MIXER CONFIRMATION')
    print(json.dumps(summary,indent=2)); print('wrote',out)

if __name__=='__main__': main()
