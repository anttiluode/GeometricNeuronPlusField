"""Isolate material degrees of freedom beyond readout distance.

Development-only mechanism diagnostic.

Procedure for each frozen arbor:

1. Optimize the corrected hostile radial learner: one material density per
   integer graph-distance shell, fixed total budget and cap.
2. Freeze the *total material in every shell* at that radial optimum.
3. Release only cell-to-cell redistribution inside each shell.
4. Optimize the same joint .03/.04 phase-coherence objective using the exact
   local material gradient projected to zero-sum within every shell.

The release learner cannot alter the radial profile at all.  Any improvement is
therefore a pure within-distance / branch-specific effect.

This is cleaner than comparing two independently trained solutions from uniform
because it nests the branch-only design space directly on top of the optimized
distance-only solution.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
import hcn_material_learning as hml
import material_radial_vs_full_v02 as rvf2
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites

rvf=rvf2.rvf


def within_shell_component(g,groups):
    g=np.asarray(g,float);z=g.copy()
    for q in groups:z[q]-=float(np.mean(z[q]))
    return z


def shell_totals(vals,groups):
    return np.asarray([float(np.sum(vals[q])) for q in groups],float)


def project_each_shell(v,groups,totals,cap):
    v=np.asarray(v,float);out=v.copy()
    for q,t in zip(groups,totals):
        out[q]=hml.project_capped_simplex(v[q],float(t),float(cap))
    return out


def vals_to_grid(shape,ids,vals):
    z=np.zeros(int(np.prod(shape)),float);z[ids]=np.asarray(vals,float)
    return z.reshape(shape)


def train_release(m,L,radial,omegas,tau,mu,sites,ids,groups,cap,steps,step_fraction):
    vals=np.asarray(radial,float).ravel()[ids].copy()
    totals=shell_totals(vals,groups)
    best_vals=vals.copy();best=-np.inf;hist=[]
    for it in range(int(steps)+1):
        d=vals_to_grid(m.body.shape,ids,vals)
        obj,g,det=hml.coherence_and_gradient(m,L,d,omegas,tau,mu,sites,ids,True)
        gp=within_shell_component(g,groups)
        gn=float(np.linalg.norm(gp));gm=float(np.max(np.abs(gp)))
        hist.append(dict(iteration=it,objective=float(obj),within_grad_norm=gn,
                         within_grad_max=gm,min_density=float(vals.min()),max_density=float(vals.max())))
        if obj>best:best=float(obj);best_vals=vals.copy()
        if it==steps or gm<1e-14:break
        base=float(step_fraction)*float(cap)/gm
        accepted=False
        for bt in range(16):
            trial=project_each_shell(vals+base*(.5**bt)*gp,groups,totals,cap)
            td=vals_to_grid(m.body.shape,ids,trial)
            tobj,_=hml.coherence_and_gradient(m,L,td,omegas,tau,mu,sites,ids,False)
            if tobj>=obj-1e-12:
                vals=trial;accepted=True;break
        if not accepted:break
    return vals_to_grid(m.body.shape,ids,best_vals),best,hist


def body_probe(m,omegas,tau,mu,g0,ratio,radial_steps,release_steps,step_fraction):
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

    radial,rad_best,rad_hist,shell_vals=rvf.train_shells(
        m,L,uniform,omegas,tau,mu,sites,ids,groups,counts,total,cap,radial_steps,step_fraction)
    released,rel_best,rel_hist=train_release(
        m,L,radial,omegas,tau,mu,sites,ids,groups,cap,release_steps,step_fraction)

    er=hml.eval_profile(m,L,radial,omegas,tau,mu,sites)
    ee=hml.eval_profile(m,L,released,omegas,tau,mu,sites)
    per={}
    for om in omegas:
        r=next(x for x in er['frequency'] if abs(x['omega']-om)<1e-12)
        e=next(x for x in ee['frequency'] if abs(x['omega']-om)<1e-12)
        per[str(float(om))]=dict(
            delta_R2=float(e['coherence2']-r['coherence2']),
            phase_rms_gain=float(r['soma_phase_rms']-e['soma_phase_rms']),
            amp_ratio=float(e['median_amp']/(r['median_amp']+1e-30)),
        )

    # Numerical shell-total audit.
    rv=np.asarray(radial).ravel()[ids];ev=np.asarray(released).ravel()[ids]
    shell_err=max(abs(float(np.sum(rv[q])-np.sum(ev[q]))) for q in groups)
    return dict(seed=int(m.cfg.seed),cells=len(cells),shells=len(groups),
                radial=er,released=ee,
                joint_gain=float(ee['coherence2']-er['coherence2']),
                phase_rms_gain=float(er['mean_phase_rms']-ee['mean_phase_rms']),
                amp_ratio=float(ee['mean_median_amp']/(er['mean_median_amp']+1e-30)),
                max_shell_total_error=float(shell_err),
                per_frequency=per,radial_history=rad_hist,release_history=rel_hist)


def summarize(rows,omegas):
    out=dict(
        bodies=len(rows),
        mean_joint_gain=float(np.mean([r['joint_gain'] for r in rows])),
        positive_joint=int(sum(r['joint_gain']>0 for r in rows)),
        mean_phase_rms_gain=float(np.mean([r['phase_rms_gain'] for r in rows])),
        positive_phase_rms=int(sum(r['phase_rms_gain']>0 for r in rows)),
        median_amp_ratio=float(np.median([r['amp_ratio'] for r in rows])),
        max_shell_total_error=float(max(r['max_shell_total_error'] for r in rows)),
    )
    f={}
    for om in omegas:
        k=str(float(om))
        f[k]=dict(
            mean_delta_R2=float(np.mean([r['per_frequency'][k]['delta_R2'] for r in rows])),
            positive=int(sum(r['per_frequency'][k]['delta_R2']>0 for r in rows)),
            mean_phase_rms_gain=float(np.mean([r['per_frequency'][k]['phase_rms_gain'] for r in rows])),
            median_amp_ratio=float(np.median([r['per_frequency'][k]['amp_ratio'] for r in rows])),
        )
    out['frequency']=f
    return out


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=622);ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2.0);ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10.0)
    ap.add_argument('--radial-steps',type=int,default=160);ap.add_argument('--release-steps',type=int,default=80)
    ap.add_argument('--step-fraction',type=float,default=.10)
    ap.add_argument('--out',default='runs/material_within_shell_release/dev_622_627.json')
    return ap.parse_args()


def main():
    a=parse_args();fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    oms=[float(x) for x in a.omegas.split(',')];rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,oms,a.tau,a.mu,a.g0,a.ratio,a.radial_steps,a.release_steps,a.step_fraction);rows.append(r)
        print(f"seed {seed}: branch-only gain={r['joint_gain']:+.6f}; " +
              ' '.join(f"w{om:.2f}:{r['per_frequency'][str(om)]['delta_R2']:+.6f}" for om in oms))
    out=dict(config=vars(a),rows=rows,summary=summarize(rows,oms))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nWITHIN-SHELL RELEASE')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
