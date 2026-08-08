"""How local are the oracle free-modal coordinates in the physical grid basis?

Development audit only. Reconstruct the F/A/B/C coordinates selected by the strong
matched-tuner baseline, then quantify a lower bound on the physical support needed
to realize each coordinate exactly in the 31x31 nearest-neighbour grid.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from full_pole_containment import stiffness_matrix
from matched_tuner_audit import modal_setup,aggregate_modal


def physical_metrics(phi,shape):
    p=np.asarray(phi,float)**2
    h,w=shape
    ipr=float(1.0/(np.sum(p*p)+1e-30))
    maxp=float(np.max(p))
    # For a pure modal stiffness shift DeltaK = phi phi^T, nearest-neighbour local
    # bond hardware can only populate diagonal and adjacent off-diagonal matrix entries.
    # Any Frobenius mass outside that sparsity pattern is an unavoidable realization error.
    diag=float(np.sum(p*p)); neigh=0.0
    for y in range(h):
        for x in range(w-1):neigh += 2.0*p[y*w+x]*p[y*w+x+1]
    for y in range(h-1):
        for x in range(w):neigh += 2.0*p[y*w+x]*p[(y+1)*w+x]
    local_pattern_mass=float(diag+neigh)
    matrix_resid_lb=float(np.sqrt(max(0.0,1.0-local_pattern_mass)))
    # An input/output-residue coordinate Delta b ~ phi is globally distributed in
    # the physical basis. Best single-node local approximation leaves this residual.
    vector_resid_best1=float(np.sqrt(max(0.0,1.0-maxp)))
    return dict(participation_ratio=ipr,max_node_power=maxp,
                best_single_node_vector_relative_residual=vector_resid_best1,
                local_matrix_pattern_power_fraction=local_pattern_mass,
                local_bond_matrix_relative_residual_lower_bound=matrix_resid_lb)


def one(m,lags,steps,P):
    wh,wv=ae.bond_weights(m,m.body);lam,bA,bB,c=modal_setup(m,wh,wv);N=len(lam)
    Z=[np.zeros(N) for _ in range(4)]
    z=aggregate_modal(m.cfg,lam,bA,bB,c,lags,steps,*Z,True);G=list(z[1:5]);names='FABC'
    vals=np.concatenate([np.abs(g) for g in G]);pick=np.argsort(vals)[::-1][:P]
    K=stiffness_matrix(wh,wv);lam2,Phi=np.linalg.eigh(K)
    rows=[]
    for q in pick:
        j=int(q//N);i=int(q%N);typ=names[j]
        met=physical_metrics(Phi[:,i],m.body.shape)
        rows.append(dict(type=typ,index=i,initial_gradient=float(G[j][i]),eigenvalue=float(lam2[i]),**met))
    return dict(seed=int(m.cfg.seed),N=N,selected=rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=240);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--lags',default='16,20,24');ap.add_argument('--steps',type=int,default=210);ap.add_argument('--P',type=int,default=8);ap.add_argument('--out',default='runs/spectral_coordinate_locality/locality.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:
        phi=np.ones(9)/3;m=physical_metrics(phi,(3,3));assert m['participation_ratio']>8.99;print('selftest ok');return
    lags=[int(x) for x in a.lags.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,lags,a.steps,a.P);rows.append(r)
        print('seed',seed,[(x['type'],x['index'],round(x['participation_ratio'],1),round(x['local_bond_matrix_relative_residual_lower_bound'],3)) for x in r['selected']],flush=True)
    sel=[x for r in rows for x in r['selected']]
    F=[x for x in sel if x['type']=='F'];R=[x for x in sel if x['type']!='F']
    def mean(xs,k):return float(np.mean([x[k] for x in xs])) if xs else float('nan')
    summary=dict(bodies=len(rows),selected_total=len(sel),selected_F=len(F),selected_residue=len(R),
                 mean_F_participation=mean(F,'participation_ratio'),mean_F_local_matrix_residual_lb=mean(F,'local_bond_matrix_relative_residual_lower_bound'),
                 min_F_local_matrix_residual_lb=float(min([x['local_bond_matrix_relative_residual_lower_bound'] for x in F])) if F else float('nan'),
                 mean_residue_participation=mean(R,'participation_ratio'),mean_residue_best1_residual=mean(R,'best_single_node_vector_relative_residual'),
                 min_residue_best1_residual=float(min([x['best_single_node_vector_relative_residual'] for x in R])) if R else float('nan'))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(summary=summary,rows=rows),indent=2));print(summary)
if __name__=='__main__':main()
