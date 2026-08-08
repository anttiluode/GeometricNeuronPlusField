"""Hostile control for the branch-only material release.

Start from the corrected optimized radial material profile, freeze every graph-
distance shell total, and release only cell-to-cell redistribution within each
shell using the exact material gradient (via material_within_shell_release).

Then destroy the learned branch/cell placement while preserving, separately in
every distance shell:

* the shell total;
* the exact multiset/histogram of learned material values;
* the radial profile.

If the learned release still beats these within-shell shuffles, the residual is
about *where* material sits among equal-distance cells, not merely introducing
heterogeneity into a shell.

Development only; no novelty claim.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import hcn_material_learning as hml
import material_radial_vs_full_v02 as rvf2
import material_within_shell_release as wsr
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites

rvf=rvf2.rvf


def grid(shape, ids, vals):
    z=np.zeros(int(np.prod(shape)),float);z[ids]=np.asarray(vals,float)
    return z.reshape(shape)


def body_probe(m,omegas,tau,mu,g0,ratio,radial_steps,release_steps,step_fraction,nshuffle):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    body=m.body.astype(bool);h,w=body.shape
    cells=[tuple(map(int,p)) for p in np.argwhere(body)]
    ids=np.asarray([p[0]*w+p[1] for p in cells],int)
    dmap,dist,levels,groups,counts=rvf.shell_structure(m,cells,ids)
    dmax=max(float(dist.max()),1.0)
    hand=np.zeros_like(dmap,float);hand[body]=float(g0)*(1+(float(ratio)-1)*dmap[body]/dmax)
    total=float(hand[body].sum());cap=float(g0)*float(ratio)
    uniform=np.zeros_like(hand,float);uniform[body]=total/len(ids)
    sites=injection_sites(m)

    radial,_,_,_=rvf.train_shells(m,L,uniform,omegas,tau,mu,sites,ids,groups,counts,total,cap,radial_steps,step_fraction)
    released,_,_=wsr.train_release(m,L,radial,omegas,tau,mu,sites,ids,groups,cap,release_steps,step_fraction)

    er=hml.eval_profile(m,L,radial,omegas,tau,mu,sites)
    el=hml.eval_profile(m,L,released,omegas,tau,mu,sites)
    vals=np.asarray(released).ravel()[ids]
    rng=np.random.default_rng(int(m.cfg.seed)+7_331_991)
    sh=[]
    for _ in range(int(nshuffle)):
        vv=vals.copy()
        for q in groups:
            tmp=vv[q].copy();rng.shuffle(tmp);vv[q]=tmp
        sh.append(hml.eval_profile(m,L,grid(body.shape,ids,vv),omegas,tau,mu,sites))

    sR=float(np.mean([x['coherence2'] for x in sh]))
    sP=float(np.mean([x['mean_phase_rms'] for x in sh]))
    sA=float(np.mean([x['mean_median_amp'] for x in sh]))
    per={}
    for om in omegas:
        l=next(x for x in el['frequency'] if abs(x['omega']-om)<1e-12)
        sr=[]
        for x in sh:
            sr.append(next(y for y in x['frequency'] if abs(y['omega']-om)<1e-12))
        per[str(float(om))]=dict(
            learned_R2=float(l['coherence2']),
            shuffle_R2=float(np.mean([x['coherence2'] for x in sr])),
            learned_minus_shuffle_R2=float(l['coherence2']-np.mean([x['coherence2'] for x in sr])),
            learned_phase_rms=float(l['soma_phase_rms']),
            shuffle_phase_rms=float(np.mean([x['soma_phase_rms'] for x in sr])),
        )
    return dict(seed=int(m.cfg.seed),cells=len(cells),shells=len(groups),
                radial_R2=float(er['coherence2']),released_R2=float(el['coherence2']),
                shuffle_R2=sR,
                release_minus_radial=float(el['coherence2']-er['coherence2']),
                release_minus_shuffle=float(el['coherence2']-sR),
                phase_rms_gain_vs_shuffle=float(sP-el['mean_phase_rms']),
                amp_ratio_vs_shuffle=float(el['mean_median_amp']/(sA+1e-30)),
                per_frequency=per)


def summarize(rows,omegas):
    out=dict(
        bodies=len(rows),
        mean_release_minus_radial=float(np.mean([r['release_minus_radial'] for r in rows])),
        positive_vs_radial=int(sum(r['release_minus_radial']>0 for r in rows)),
        mean_release_minus_shell_shuffle=float(np.mean([r['release_minus_shuffle'] for r in rows])),
        positive_vs_shell_shuffle=int(sum(r['release_minus_shuffle']>0 for r in rows)),
        mean_phase_rms_gain_vs_shuffle=float(np.mean([r['phase_rms_gain_vs_shuffle'] for r in rows])),
        median_amp_ratio_vs_shuffle=float(np.median([r['amp_ratio_vs_shuffle'] for r in rows])),
    )
    out['frequency']={}
    for om in omegas:
        k=str(float(om));d=[r['per_frequency'][k]['learned_minus_shuffle_R2'] for r in rows]
        out['frequency'][k]=dict(mean_learned_minus_shuffle_R2=float(np.mean(d)),positive=int(sum(x>0 for x in d)))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=622);ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2);ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10)
    ap.add_argument('--radial-steps',type=int,default=160);ap.add_argument('--release-steps',type=int,default=80)
    ap.add_argument('--step-fraction',type=float,default=.10);ap.add_argument('--nshuffle',type=int,default=12)
    ap.add_argument('--out',default='runs/material_within_shell_shuffle/dev_622_627.json')
    a=ap.parse_args()
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    oms=[float(x) for x in a.omegas.split(',')];rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,oms,a.tau,a.mu,a.g0,a.ratio,a.radial_steps,a.release_steps,a.step_fraction,a.nshuffle);rows.append(r)
        print(f"seed {seed}: release-radial={r['release_minus_radial']:+.6f} release-shellshuffle={r['release_minus_shuffle']:+.6f}")
    out=dict(config=vars(a),rows=rows,summary=summarize(rows,oms))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nWITHIN-SHELL SHUFFLE CONTROL')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
