"""Ask whether learned quasi-active material is only a distance profile.

Development control on the material-learning discovery bodies.

Train exactly as in hcn_material_learning.py, then compare:

1. learned: full local density map;
2. radialized: replace every cell's density by the mean learned density at the
   same graph distance from the soma;
3. within-distance shuffle: preserve every distance shell's exact density
   multiset but permute values among branches/cells within that shell;
4. global shuffle: destroy both distance and branch relation.

If learned ~= radialized ~= within-distance shuffle, graph distance is the
important coordinate.  If learned beats distance-preserving controls, local
branch geometry contributes beyond path length.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites
from hcn_material_learning import body_probe as train_body, eval_profile


def reconstruct_density(m,row):
    z=np.zeros(m.body.shape,float)
    for q in row['learned_density']:
        z[tuple(q['cell'])]=float(q['density'])
    return z


def geometry_controls(m,row,omegas,tau,mu,nshuffle):
    learned=reconstruct_density(m,row)
    body=m.body.astype(bool);dist=m.graph_distance_from_soma()
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv);sites=injection_sites(m)
    rng=np.random.default_rng(int(m.cfg.seed)+818_181)

    radial=np.zeros_like(learned)
    for dd in sorted(set(int(dist[p]) for p in map(tuple,np.argwhere(body)))):
        mask=body&(dist==dd);radial[mask]=float(np.mean(learned[mask]))

    within=[]
    for _ in range(int(nshuffle)):
        z=np.zeros_like(learned)
        for dd in sorted(set(int(dist[p]) for p in map(tuple,np.argwhere(body)))):
            cells=np.argwhere(body&(dist==dd));vals=np.asarray([learned[tuple(p)] for p in cells],float)
            rng.shuffle(vals)
            for p,v in zip(cells,vals):z[tuple(p)]=v
        within.append(eval_profile(m,L,z,omegas,tau,mu,sites))

    vals=learned[body].copy();glob=[]
    for _ in range(int(nshuffle)):
        v=vals.copy();rng.shuffle(v);z=np.zeros_like(learned);z[body]=v
        glob.append(eval_profile(m,L,z,omegas,tau,mu,sites))

    full=eval_profile(m,L,learned,omegas,tau,mu,sites)
    rad=eval_profile(m,L,radial,omegas,tau,mu,sites)
    wi=dict(coherence2=float(np.mean([x['coherence2'] for x in within])),
            mean_phase_rms=float(np.mean([x['mean_phase_rms'] for x in within])))
    gl=dict(coherence2=float(np.mean([x['coherence2'] for x in glob])),
            mean_phase_rms=float(np.mean([x['mean_phase_rms'] for x in glob])))
    residual=learned[body]-radial[body]
    return dict(
        learned=full,radialized=rad,within_distance_shuffle=wi,global_shuffle=gl,
        coherence_advantage_over_radial=float(full['coherence2']-rad['coherence2']),
        coherence_advantage_over_within=float(full['coherence2']-wi['coherence2']),
        phase_rms_advantage_over_radial=float(rad['mean_phase_rms']-full['mean_phase_rms']),
        phase_rms_advantage_over_within=float(wi['mean_phase_rms']-full['mean_phase_rms']),
        rms_branch_residual=float(np.sqrt(np.mean(residual**2))),
        learned_density_std=float(np.std(learned[body])),
        branch_residual_fraction=float(np.sqrt(np.mean(residual**2))/(np.std(learned[body])+1e-30)),
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=580)
    ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2.0)
    ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005)
    ap.add_argument('--ratio',type=float,default=10.0)
    ap.add_argument('--steps',type=int,default=50)
    ap.add_argument('--step-fraction',type=float,default=.10)
    ap.add_argument('--nshuffle',type=int,default=24)
    ap.add_argument('--out',default='runs/hcn_material_geometry/dev.json')
    return ap.parse_args()


def main():
    a=parse_args();fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    omegas=[float(x) for x in a.omegas.split(',')]
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        trained=train_body(m,omegas,a.tau,a.mu,a.g0,a.ratio,a.steps,a.step_fraction,4)
        gc=geometry_controls(m,trained,omegas,a.tau,a.mu,a.nshuffle)
        rows.append(dict(seed=seed,morphology=trained['morphology'],geometry_controls=gc))
        print(f"seed {seed}: full={gc['learned']['coherence2']:.4f} radial={gc['radialized']['coherence2']:.4f} "
              f"within={gc['within_distance_shuffle']['coherence2']:.4f} global={gc['global_shuffle']['coherence2']:.4f} "
              f"branch-resid={gc['branch_residual_fraction']:.3f}",flush=True)
    def M(k):return float(np.mean([r['geometry_controls'][k] for r in rows]))
    summary=dict(
        bodies=len(rows),
        learned_coherence2=float(np.mean([r['geometry_controls']['learned']['coherence2'] for r in rows])),
        radialized_coherence2=float(np.mean([r['geometry_controls']['radialized']['coherence2'] for r in rows])),
        within_distance_shuffle_coherence2=float(np.mean([r['geometry_controls']['within_distance_shuffle']['coherence2'] for r in rows])),
        global_shuffle_coherence2=float(np.mean([r['geometry_controls']['global_shuffle']['coherence2'] for r in rows])),
        mean_advantage_over_radial=M('coherence_advantage_over_radial'),
        mean_advantage_over_within=M('coherence_advantage_over_within'),
        mean_branch_residual_fraction=M('branch_residual_fraction'),
        learned_beats_radial_bodies=int(sum(r['geometry_controls']['coherence_advantage_over_radial']>0 for r in rows)),
        learned_beats_within_bodies=int(sum(r['geometry_controls']['coherence_advantage_over_within']>0 for r in rows)),
    )
    payload=dict(experiment='learned_material_distance_preserving_controls_v01',development_only=True,
                 seed_start=a.seed_start,seeds_requested=a.seeds,omegas=omegas,summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nDISTANCE-PRESERVING CONTROLS');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
