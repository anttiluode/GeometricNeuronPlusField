"""Map free-modal coordinates from the exact 31x31 operator back to the 70-cell body graph basis.

Descriptive mechanism audit after the matched-tuner benchmark. The weak bath creates
961 full-grid modes; this asks whether the winning free spectral coordinates are
simply the embedded body resonances already identified by the graph-mode microscope.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from collections import Counter
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from full_pole_containment import stiffness_matrix
from matched_tuner_audit import modal_setup,aggregate_modal
from graph_mode_probe import graph_laplacian_modes


def selected_full_modes(m,lags,steps,P):
    wh,wv=ae.bond_weights(m,m.body);lam,bA,bB,c=modal_setup(m,wh,wv);N=len(lam);Z=[np.zeros(N) for _ in range(4)]
    z=aggregate_modal(m.cfg,lam,bA,bB,c,lags,steps,*Z,True);G=list(z[1:5]);names='FABC'
    vals=np.concatenate([np.abs(g) for g in G]);pick=np.argsort(vals)[::-1][:P]
    return [(names[int(q//N)],int(q%N),float(G[int(q//N)][int(q%N)])) for q in pick],wh,wv


def one(m,lags,steps,P):
    selected,wh,wv=selected_full_modes(m,lags,steps,P)
    K=stiffness_matrix(wh,wv);full_lam,full_phi=np.linalg.eigh(K)
    coords,body_lam,body_phi=graph_laplacian_modes(m.body)
    w=m.body.shape[1];idx=np.asarray([y*w+x for y,x in coords],int)
    rows=[]
    for typ,j,g in selected:
        r=np.asarray(full_phi[idx,j],float);rn=float(np.linalg.norm(r));
        if rn>1e-15:r=r/rn
        ov=np.abs(body_phi.T@r);ib=int(np.argmax(ov));best=float(ov[ib])
        body_fraction=float(np.sum(full_phi[idx,j]**2))
        rows.append(dict(type=typ,full_index=j,initial_gradient=g,full_eigenvalue=float(full_lam[j]),
                         body_mode=ib,body_eigenvalue=float(body_lam[ib]),body_overlap_abs=best,
                         body_overlap_power=float(best*best),full_mode_power_on_body=body_fraction))
    return dict(seed=int(m.cfg.seed),cells=int(m.body.sum()),selected=rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=288);ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lags',default='16,20,24');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--out',default='runs/full_body_mode_map/map.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        x=np.eye(3);r=np.array([0.,1.,0.]);assert int(np.argmax(abs(x.T@r)))==1;print('selftest ok');return
    lags=[int(x) for x in a.lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,lags,a.steps,a.P);rows.append(r)
        print('seed',seed,[(x['type'],x['full_index'],x['body_mode'],round(x['body_overlap_abs'],3),round(x['full_mode_power_on_body'],3)) for x in r['selected']],flush=True)
    s=[x for r in rows for x in r['selected']];uniq={(r['seed'],x['full_index']):x for r in rows for x in r['selected']}
    u=list(uniq.values());cnt=Counter(x['body_mode'] for x in s)
    summary=dict(bodies=len(rows),selected_coordinates=len(s),unique_selected_full_modes=len(u),
                 mean_body_overlap_abs=float(np.mean([x['body_overlap_abs'] for x in u])),
                 min_body_overlap_abs=float(np.min([x['body_overlap_abs'] for x in u])),
                 mean_full_mode_power_on_body=float(np.mean([x['full_mode_power_on_body'] for x in u])),
                 min_full_mode_power_on_body=float(np.min([x['full_mode_power_on_body'] for x in u])),
                 selected_body_modes_hist={str(k):int(v) for k,v in sorted(cnt.items())},
                 fraction_selected_body_modes_0_17=float(np.mean([x['body_mode']<=17 for x in s])),
                 fraction_selected_body_modes_18_20=float(np.mean([18<=x['body_mode']<=20 for x in s])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
