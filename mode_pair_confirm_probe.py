"""Held-out confirmation of soma cross-term graph-mode pair structure.

See MODE_PAIR_CONFIRM_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from graph_mode_probe import graph_laplacian_modes
from transfer_decomposition_probe import trace_single, tshift
from mode_pair_discovery_probe import concentration


def contrast_trace(pT, pD):
    pT = np.asarray(pT, float); pD = np.asarray(pD, float)
    t = float(np.max(pT)); d = float(np.max(pD))
    return float((t - d) / (t + d + 1e-30))


def body_confirm(m, lag, steps):
    hA = trace_single(m, 0, 1.0, steps)
    hB = trace_single(m, 1, 1.0, steps)
    coords, evals, evecs = graph_laplacian_modes(m.body)
    n = len(coords)
    cidx = {p: i for i, p in enumerate(coords)}
    soma = tuple(map(int, m.soma)); si = cidx[soma]
    ys = np.asarray([p[0] for p in coords], int)
    xs = np.asarray([p[1] for p in coords], int)
    ZA = np.asarray(hA[:, ys, xs], np.complex128)
    ZB = np.asarray(hB[:, ys, xs], np.complex128)
    QA = ZA @ evecs; QB = ZB @ evecs
    phi_s = np.asarray(evecs[si], float)
    uA = QA * phi_s[None, :]
    uB = QB * phi_s[None, :]

    hsA = np.asarray(hA[:, soma[0], soma[1]], np.complex128)
    hsB = np.asarray(hB[:, soma[0], soma[1]], np.complex128)
    recA = uA.sum(axis=1); recB = uB.sum(axis=1)
    src_mae = 0.5 * (np.mean(np.abs(recA-hsA)) + np.mean(np.abs(recB-hsB)))
    src_scale = 0.5 * (np.mean(np.abs(hsA)) + np.mean(np.abs(hsB))) + 1e-30
    src_rel = float(src_mae/src_scale)

    uA_l = tshift(uA, lag); uB_l = tshift(uB, lag)
    hsA_l = tshift(hsA, lag); hsB_l = tshift(hsB, lag)
    recA_l = tshift(recA, lag); recB_l = tshift(recB, lag)

    psiT = recA + recB_l
    psiD = recB + recA_l
    pT = np.abs(psiT)**2; pD = np.abs(psiD)**2
    tT = int(np.argmax(pT)); tD = int(np.argmax(pD))
    Cfull = contrast_trace(pT, pD)

    MT = 2.0*np.real(np.outer(uA[tT], np.conj(uB_l[tT])))
    MD = 2.0*np.real(np.outer(uB[tD], np.conj(uA_l[tD])))
    M = MT-MD
    direct_delta = float(2*np.real(hsA[tT]*np.conj(hsB_l[tT])) -
                         2*np.real(hsB[tD]*np.conj(hsA_l[tD])))
    pair_delta = float(M.sum())
    pair_rel = float(abs(pair_delta-direct_delta)/(abs(direct_delta)+1e-30))

    pair_rows=[]; absvals=[]; mode_inv=np.zeros(n,float); diag_abs=0.0; m0_abs=0.0
    for i in range(n):
        v=float(M[i,i]); av=abs(v)
        pair_rows.append((i,i,v,av)); absvals.append(av); mode_inv[i]+=av; diag_abs+=av
        if i==0: m0_abs+=av
        for j in range(i+1,n):
            v=float(M[i,j]+M[j,i]); av=abs(v)
            pair_rows.append((i,j,v,av)); absvals.append(av)
            mode_inv[i]+=0.5*av; mode_inv[j]+=0.5*av
            if i==0 or j==0: m0_abs+=av
    conc=concentration(absvals); total_abs=conc['total_abs']+1e-30
    mode_frac=mode_inv/(mode_inv.sum()+1e-30)
    diag_frac=float(diag_abs/total_abs)
    m0_pair_frac=float(m0_abs/total_abs)
    low18=float(mode_frac[:min(18,n)].sum())
    ti=[i for i in (18,19,20) if i<n]
    task=float(mode_frac[ti].sum()) if ti else 0.0
    expected=float(len(ti)/n) if n else 0.0
    enrich=float(task/expected) if expected>0 else float('nan')

    # Frozen exploratory mechanism ablations on the full traces.
    envT=np.abs(recA)**2 + np.abs(recB_l)**2
    envD=np.abs(recB)**2 + np.abs(recA_l)**2
    XfullT=2*np.real(recA*np.conj(recB_l)); XfullD=2*np.real(recB*np.conj(recA_l))
    X00T=2*np.real(uA[:,0]*np.conj(uB_l[:,0])); X00D=2*np.real(uB[:,0]*np.conj(uA_l[:,0]))
    nzA=uA[:,1:].sum(axis=1); nzB=uB[:,1:].sum(axis=1)
    nzA_l=tshift(nzA,lag); nzB_l=tshift(nzB,lag)
    XnzT=2*np.real(nzA*np.conj(nzB_l)); XnzD=2*np.real(nzB*np.conj(nzA_l))
    X0mixT=XfullT-X00T-XnzT; X0mixD=XfullD-X00D-XnzD
    XdiagT=2*np.real(np.sum(uA*np.conj(uB_l),axis=1))
    XdiagD=2*np.real(np.sum(uB*np.conj(uA_l),axis=1))
    XoffT=XfullT-XdiagT; XoffD=XfullD-XdiagD

    C_no0mix=contrast_trace(envT+X00T+XnzT, envD+X00D+XnzD)
    C_0mixonly=contrast_trace(envT+X0mixT, envD+X0mixD)
    C_diag=contrast_trace(envT+XdiagT, envD+XdiagD)
    C_off=contrast_trace(envT+XoffT, envD+XoffD)

    return dict(
        cells=n, soma=list(soma), soma_contrast_full=Cfull,
        source_reconstruction_rel=src_rel, pair_reconstruction_rel=pair_rel,
        concentration=conc, diagonal_abs_fraction=diag_frac,
        mode0_involvement_fraction=float(mode_frac[0]),
        mode0_pair_abs_fraction=m0_pair_frac,
        low_modes_0_17_involvement_fraction=low18,
        task_band_18_20_involvement_fraction=task,
        task_band_18_20_enrichment=enrich,
        mode_involvement_fraction=[float(x) for x in mode_frac],
        exploratory=dict(
            C_no_mode0_nonzero_crossmix=C_no0mix,
            C_mode0_nonzero_crossmix_only=C_0mixonly,
            C_diagonal_cross_only=C_diag,
            C_offdiagonal_cross_only=C_off,
            absC_loss_remove_mode0_mix=float(abs(Cfull)-abs(C_no0mix)),
            absC_loss_diagonal_only=float(abs(Cfull)-abs(C_diag)),
            absC_loss_offdiagonal_only=float(abs(Cfull)-abs(C_off)),
            signed_error_remove_mode0_mix=float(abs(Cfull-C_no0mix)),
        ),
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=84)
    ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/mode_pair_confirm/mode_pair_confirm.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    pT=np.array([1.,3.,2.]); pD=np.array([1.,1.,1.])
    assert abs(contrast_trace(pT,pD)-0.5)<1e-12
    print('selftest ok')


def main():
    a=parse_args()
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
        r=body_confirm(m,a.lag,a.steps); r['seed']=seed; rows.append(r)
        c=r['concentration']
        print(f"seed {seed}: f50={c['frac50']:.4f} f80={c['frac80']:.4f} PR={c['participation_fraction']:.4f} "
              f"diag={r['diagonal_abs_fraction']:.4f} m0inv={r['mode0_involvement_fraction']:.4f} "
              f"m0pair={r['mode0_pair_abs_fraction']:.4f} low18={r['low_modes_0_17_involvement_fraction']:.4f} "
              f"band18-20={r['task_band_18_20_enrichment']:.3f}",flush=True)
    if not rows: raise SystemExit('No valid bodies')

    def v(path):
        z=[]
        for r in rows:
            x=r
            for k in path: x=x[k]
            z.append(float(x))
        return np.asarray(z,float)

    # mean mode involvement and rank
    mat=np.asarray([r['mode_involvement_fraction'] for r in rows],float)
    mm=mat.mean(axis=0); top=int(np.argmax(mm))
    P1=bool(v(['source_reconstruction_rel']).mean()<1e-10 and v(['pair_reconstruction_rel']).mean()<1e-10)
    P2=bool(v(['concentration','frac50']).mean()<.020 and v(['concentration','frac80']).mean()<.070 and v(['concentration','participation_fraction']).mean()<.050)
    P3=bool(v(['diagonal_abs_fraction']).mean()<.060)
    P4=bool(v(['mode0_involvement_fraction']).mean()>.150 and v(['mode0_pair_abs_fraction']).mean()>.300 and top==0)
    P5=bool(v(['low_modes_0_17_involvement_fraction']).mean()>.700)
    P6=bool(v(['task_band_18_20_enrichment']).mean()<1.0)

    summary=dict(
        bodies=len(rows),
        source_reconstruction_rel_mean=float(v(['source_reconstruction_rel']).mean()),
        pair_reconstruction_rel_mean=float(v(['pair_reconstruction_rel']).mean()),
        mean_pair_fraction_50=float(v(['concentration','frac50']).mean()),
        mean_pair_fraction_80=float(v(['concentration','frac80']).mean()),
        mean_effective_pair_fraction=float(v(['concentration','participation_fraction']).mean()),
        mean_diagonal_abs_fraction=float(v(['diagonal_abs_fraction']).mean()),
        mean_mode0_involvement_fraction=float(v(['mode0_involvement_fraction']).mean()),
        mean_mode0_pair_abs_fraction=float(v(['mode0_pair_abs_fraction']).mean()),
        highest_mean_involvement_mode=top,
        mean_low_modes_0_17_involvement=float(v(['low_modes_0_17_involvement_fraction']).mean()),
        mean_task_band_18_20_enrichment=float(v(['task_band_18_20_enrichment']).mean()),
        registered=dict(P1=P1,P2=P2,P3=P3,P4=P4,P5=P5,P6=P6,all_pass=bool(P1 and P2 and P3 and P4 and P5 and P6)),
        mean_mode_involvement=[float(x) for x in mm],
        exploratory=dict(
            mean_absC_full=float(np.mean(np.abs(v(['soma_contrast_full'])))),
            mean_absC_no_mode0_mix=float(np.mean(np.abs(v(['exploratory','C_no_mode0_nonzero_crossmix'])))),
            mean_absC_mode0_mix_only=float(np.mean(np.abs(v(['exploratory','C_mode0_nonzero_crossmix_only'])))),
            mean_absC_diagonal_only=float(np.mean(np.abs(v(['exploratory','C_diagonal_cross_only'])))),
            mean_absC_offdiagonal_only=float(np.mean(np.abs(v(['exploratory','C_offdiagonal_cross_only'])))),
            mean_absC_loss_remove_mode0_mix=float(v(['exploratory','absC_loss_remove_mode0_mix']).mean()),
            mean_signed_error_remove_mode0_mix=float(v(['exploratory','signed_error_remove_mode0_mix']).mean()),
        )
    )
    payload=dict(experiment='mode_pair_confirm_v01',prereg='MODE_PAIR_CONFIRM_PREREG_V01.md',seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,summary=summary,rows=rows)
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(payload,indent=2),encoding='utf-8')

    print('\nMODE-PAIR CONFIRMATION RECEIPT')
    for k in ('P1','P2','P3','P4','P5','P6','all_pass'): print(k,summary['registered'][k])
    print(f"f50/f80/PR = {summary['mean_pair_fraction_50']:.4f}/{summary['mean_pair_fraction_80']:.4f}/{summary['mean_effective_pair_fraction']:.4f}")
    print(f"diag={summary['mean_diagonal_abs_fraction']:.4f} m0inv={summary['mean_mode0_involvement_fraction']:.4f} m0pair={summary['mean_mode0_pair_abs_fraction']:.4f}")
    print(f"low0-17={summary['mean_low_modes_0_17_involvement']:.4f} band18-20 enrich={summary['mean_task_band_18_20_enrichment']:.4f} topmode={top}")
    print('exploratory',summary['exploratory'])
    print(f'wrote {out}')

if __name__=='__main__': main()
