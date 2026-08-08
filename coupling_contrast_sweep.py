"""Hardware development sweep: reduce strong/weak coupling contrast.

The LC compiler maps k_e to inverse coupling inductance. Historical k_arbor/k_bath
is 12,500:1. This asks how quickly temporal-order contrast disappears as the
background coupling is raised toward the strong-arbor value.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad

BATHS=(0.0002,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1.0)


def weights_at(m,kb):
    wh,wv=ae.bond_weights(m,m.body);ka=float(m.cfg.k_arbor)
    wh=np.where(np.isclose(wh,ka,rtol=0,atol=1e-10),ka,float(kb))
    wv=np.where(np.isclose(wv,ka,rtol=0,atol=1e-10),ka,float(kb))
    return wh,wv


def order_metrics(m,wh,wv,lag,steps,target):
    seq=ae.source_sequence(m,target,lag,steps)
    ps,vs,E=ae.linear_forward(m,wh,wv,seq,store=True)
    sy,sx=map(int,m.soma);p=np.abs(ps[1:,sy,sx])**2
    return float(E),float(np.max(p))


def one(m,lags,steps):
    rows=[]
    base_sign={}
    for kb in BATHS:
        wh,wv=weights_at(m,kb);csE=[];csP=[]
        for lag in lags:
            ET,pT=order_metrics(m,wh,wv,lag,steps,True);ED,pD=order_metrics(m,wh,wv,lag,steps,False)
            cE=float((ET-ED)/(ET+ED+1e-30));cP=float((pT-pD)/(pT+pD+1e-30))
            csE.append(cE);csP.append(cP)
            if kb==BATHS[0]:base_sign[lag]=(np.sign(cE),np.sign(cP))
        rows.append(dict(k_bath=float(kb),contrast_ratio=float(m.cfg.k_arbor/kb),
                         mean_abs_energy_C=float(np.mean(np.abs(csE))),mean_abs_peak_C=float(np.mean(np.abs(csP))),
                         signed_energy_C=[float(x) for x in csE],signed_peak_C=[float(x) for x in csP],
                         energy_sign_retention=float(np.mean([np.sign(c)==base_sign[l][0] for c,l in zip(csE,lags)])),
                         peak_sign_retention=float(np.mean([np.sign(c)==base_sign[l][1] for c,l in zip(csP,lags)]))))
    return dict(seed=int(m.cfg.seed),rows=rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=288);ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lags',default='14,18,20,22,26');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--out',default='runs/coupling_contrast_sweep/sweep.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        assert abs(2.5/BATHS[0]-12500)<1e-9;print('selftest ok');return
    lags=[int(x) for x in a.lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    bodies=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,lags,a.steps);bodies.append(r);print('seed',seed,'done',flush=True)
    summary=[]
    for j,kb in enumerate(BATHS):
        E=np.asarray([r['rows'][j]['mean_abs_energy_C'] for r in bodies]);P=np.asarray([r['rows'][j]['mean_abs_peak_C'] for r in bodies])
        Es=np.asarray([r['rows'][j]['energy_sign_retention'] for r in bodies]);Ps=np.asarray([r['rows'][j]['peak_sign_retention'] for r in bodies])
        E0=np.asarray([r['rows'][0]['mean_abs_energy_C'] for r in bodies]);P0=np.asarray([r['rows'][0]['mean_abs_peak_C'] for r in bodies])
        q=dict(k_bath=float(kb),contrast_ratio=float(2.5/kb),mean_abs_energy_C=float(E.mean()),mean_abs_peak_C=float(P.mean()),
               energy_fraction_of_historical=float(np.mean(E/(E0+1e-30))),peak_fraction_of_historical=float(np.mean(P/(P0+1e-30))),
               mean_energy_sign_retention=float(Es.mean()),mean_peak_sign_retention=float(Ps.mean()))
        summary.append(q);print(q,flush=True)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='coupling_contrast_dev_v01',seed_start=a.seed_start,bodies=len(bodies),lags=lags,summary=summary,rows=bodies),indent=2))
if __name__=='__main__':main()
