"""Evaluate RETURN_GATE_SPECTRAL_CONFIRM_PREREG_V01.md from a JSON receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

PERIODS=(2,6,10,14,30,42,70)


def pcorr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)<2 or np.std(a)<1e-30 or np.std(b)<1e-30:return 0.0
    return float(np.corrcoef(a,b)[0,1])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--receipt',required=True)
    ap.add_argument('--out',default='runs/return_gate_spectral_confirm/verdict.json')
    a=ap.parse_args()
    p=json.loads(Path(a.receipt).read_text())
    s=p['summary']; q=s['periodic']

    damages=[1.0-q[str(P)]['mean_map_corr'] for P in PERIODS]
    k8=[q[str(P)]['K8_mean_contamination'] for P in PERIODS]
    k16=[q[str(P)]['K16_mean_contamination'] for P in PERIODS]
    periodic_r8=pcorr(k8,damages); periodic_r16=pcorr(k16,damages)
    S0=periodic_r8>.95; S1=periodic_r16>.95

    categories=[str(P) for P in PERIODS]
    mapcorr=[q[x]['mean_map_corr'] for x in categories]+[s['random']['mean_map_corr'],s['block']['mean_map_corr']]
    c8=[q[x]['K8_mean_contamination'] for x in categories]+[s['random']['K8_mean_contamination'],s['block']['K8_mean_contamination']]
    c16=[q[x]['K16_mean_contamination'] for x in categories]+[s['random']['K16_mean_contamination'],s['block']['K16_mean_contamination']]
    category_r8=pcorr(c8,1-np.asarray(mapcorr)); category_r16=pcorr(c16,1-np.asarray(mapcorr))
    S2=category_r8>.90 and category_r16>.90

    ratio8=q['42']['K8_mean_contamination']/(q['6']['K8_mean_contamination']+1e-30)
    ratio16=q['42']['K16_mean_contamination']/(q['6']['K16_mean_contamination']+1e-30)
    mapgap=q['6']['mean_map_corr']-q['42']['mean_map_corr']
    S3=ratio8>10 and ratio16>10 and mapgap>.08

    block_random8=s['block']['K8_mean_contamination']/(s['random']['K8_mean_contamination']+1e-30)
    block_random16=s['block']['K16_mean_contamination']/(s['random']['K16_mean_contamination']+1e-30)
    behavior_gap=s['random']['mean_map_corr']-s['block']['mean_map_corr']
    S4=block_random8>3 and block_random16>3 and behavior_gap>.15

    indiv8=s['K']['8']['contamination_vs_one_minus_corr']
    indiv16=s['K']['16']['contamination_vs_one_minus_corr']
    S5=indiv8>.30 and indiv16>.30

    verdict=dict(
        prereg='RETURN_GATE_SPECTRAL_CONFIRM_PREREG_V01.md',
        bodies=s['bodies'],
        criteria={
            'S0_periodic_K8':dict(pass_=bool(S0),corr=periodic_r8),
            'S1_periodic_K16':dict(pass_=bool(S1),corr=periodic_r16),
            'S2_all_categories':dict(pass_=bool(S2),K8_corr=category_r8,K16_corr=category_r16),
            'S3_fast_slow_predictor_gap':dict(pass_=bool(S3),K8_ratio=ratio8,K16_ratio=ratio16,map_corr_gap=mapgap),
            'S4_block_vs_random':dict(pass_=bool(S4),K8_ratio=block_random8,K16_ratio=block_random16,map_corr_gap=behavior_gap),
            'S5_individual_partial_prediction':dict(pass_=bool(S5),K8_corr=indiv8,K16_corr=indiv16),
        },
        summary=s,
    )
    verdict['passed']=sum(int(x['pass_']) for x in verdict['criteria'].values())
    verdict['total']=len(verdict['criteria'])
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    if verdict['passed']!=verdict['total']:
        raise SystemExit(f"registered confirmation failed: {verdict['passed']}/{verdict['total']}")

if __name__=='__main__':main()
