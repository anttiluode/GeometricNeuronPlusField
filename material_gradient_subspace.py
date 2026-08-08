"""Decompose material credit into radial and within-shell components.

Development-only mechanism diagnostic on already-seen bodies.

The corrected radial control showed:

    omega=.03  full material does not beat optimized radial material
    omega=.04  a small full-material advantage remains

This script asks whether that frequency split is already visible in the exact
local material gradient itself.

For each cellwise gradient g, decompose in the per-cell Euclidean metric:

    global budget tangent:  g_t = g - mean(g)
    radial component:       g_r = shell_mean(g) - mean(g)
    within-shell component: g_b = g - shell_mean(g)

where a shell is one integer graph distance from the soma/readout.

`g_b` is exactly the first-order credit unavailable to a distance-only material
field.  It sums to zero inside every shell, so it also preserves total material
to first order.

We measure this decomposition at both the uniform state and the optimized
radial state, separately at omega=.03 and omega=.04.  We additionally compare
the within-shell gradient to simple morphology/topology features to see whether
it is recognizably branch structured.

No held-out confirmation or novelty claim is made here.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import adjoint_eligibility_probe as ae
import hcn_material_learning as hml
import material_radial_vs_full_v02 as rvf2  # installs corrected shell projection
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites

rvf=rvf2.rvf


def rho(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if len(a)<2 or np.std(a)<1e-14 or np.std(b)<1e-14:return float('nan')
    return float(spearmanr(a,b).statistic)


def shell_decompose(g, groups):
    g=np.asarray(g,float)
    gm=float(np.mean(g))
    shell_mean=np.zeros_like(g)
    for q in groups:
        shell_mean[q]=float(np.mean(g[q]))
    gt=g-gm
    gr=shell_mean-gm
    gb=g-shell_mean
    nt=float(np.linalg.norm(gt));nr=float(np.linalg.norm(gr));nb=float(np.linalg.norm(gb))
    return dict(
        norm_tangent=nt,norm_radial=nr,norm_within_shell=nb,
        within_fraction=float(nb/(nt+1e-30)),
        within_energy_fraction=float(nb*nb/(nt*nt+1e-30)),
        radial_energy_fraction=float(nr*nr/(nt*nt+1e-30)),
        orthogonality=float(np.dot(gr,gb)/(nr*nb+1e-30)),
    ),gb,gr


def graph_features(m,cells):
    body=m.body.astype(bool);h,w=body.shape
    index={p:k for k,p in enumerate(cells)}
    neigh=[]
    for r,c in cells:
        q=[]
        for dr,dc in ((1,0),(-1,0),(0,1),(0,-1)):
            p=(r+dr,c+dc)
            if p in index:q.append(index[p])
        neigh.append(q)
    degree=np.asarray([len(q) for q in neigh],float)

    dmap=m.graph_distance_from_soma().astype(int)
    dist=np.asarray([int(dmap[p]) for p in cells],int)

    # Deterministic shortest-path parent, sufficient as a descriptive rooted-tree
    # feature even if the occupied grid contains an occasional equal-distance cycle.
    soma_idx=index[tuple(m.soma)]
    parent=np.full(len(cells),-1,int);parent[soma_idx]=soma_idx
    for i in np.argsort(dist):
        if i==soma_idx:continue
        cand=[j for j in neigh[i] if dist[j]==dist[i]-1]
        if cand:parent[i]=min(cand)

    subtree=np.ones(len(cells),float)
    for i in np.argsort(dist)[::-1]:
        if i==soma_idx:continue
        p=parent[i]
        if p>=0 and p!=i:subtree[p]+=subtree[i]

    # Major branch identity = first shortest-path child below soma.
    branch=np.full(len(cells),-1,int);branch[soma_idx]=soma_idx
    for i in np.argsort(dist):
        if i==soma_idx:continue
        p=parent[i]
        if p<0:continue
        branch[i]=i if p==soma_idx else branch[p]

    # Distance to nearest terminal (degree <= 1, excluding soma if possible).
    terminals=[i for i,d in enumerate(degree) if d<=1 and i!=soma_idx]
    if not terminals:terminals=[soma_idx]
    dt=np.full(len(cells),np.inf);dq=deque(terminals)
    for i in terminals:dt[i]=0
    while dq:
        i=dq.popleft()
        for j in neigh[i]:
            if dt[j]>dt[i]+1:
                dt[j]=dt[i]+1;dq.append(j)

    # Shell population is another pure-distance degeneracy descriptor.
    shell_size=np.zeros(len(cells),float)
    for d in np.unique(dist):shell_size[dist==d]=np.sum(dist==d)

    return dict(distance=dist.astype(float),degree=degree,subtree=subtree,
                distance_to_terminal=dt,shell_size=shell_size,branch=branch)


def branch_eta2(values,branch):
    x=np.asarray(values,float);b=np.asarray(branch,int)
    gm=float(np.mean(x));den=float(np.sum((x-gm)**2))
    if den<1e-30:return float('nan')
    num=0.0
    for q in np.unique(b):
        z=x[b==q]
        if len(z):num+=len(z)*(float(np.mean(z))-gm)**2
    return float(num/den)


def gradient_at(m,L,density,omega,tau,mu,sites,ids):
    obj,g,det=hml.coherence_and_gradient(m,L,density,[omega],tau,mu,sites,ids,True)
    return float(obj),np.asarray(g,float),det[0]


def body_probe(m,omegas,tau,mu,g0,ratio,full_steps,radial_steps,step_fraction):
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

    full,full_best,_=hml.train_density(m,L,uniform,omegas,tau,mu,sites,ids,total,cap,full_steps,step_fraction)
    radial,rad_best,_,shell_vals=rvf.train_shells(m,L,uniform,omegas,tau,mu,sites,ids,groups,counts,total,cap,radial_steps,step_fraction)

    features=graph_features(m,cells)
    states={}
    for state_name,density in [('uniform',uniform),('radial',radial)]:
        freq={}
        for om in omegas:
            obj,g,det=gradient_at(m,L,density,om,tau,mu,sites,ids)
            dec,gb,gr=shell_decompose(g,groups)
            top={k:rho(gb,v) for k,v in features.items() if k!='branch'}
            top['branch_eta2']=branch_eta2(gb,features['branch'])
            freq[str(float(om))]=dict(objective=obj,decomposition=dec,
                                      within_topology=top,
                                      grad_min=float(g.min()),grad_max=float(g.max()))
        states[state_name]=freq

    ev_full=hml.eval_profile(m,L,full,omegas,tau,mu,sites)
    ev_rad=hml.eval_profile(m,L,radial,omegas,tau,mu,sites)
    per_freq={}
    for om in omegas:
        f=next(x for x in ev_full['frequency'] if abs(x['omega']-om)<1e-12)
        r=next(x for x in ev_rad['frequency'] if abs(x['omega']-om)<1e-12)
        per_freq[str(float(om))]=dict(full_minus_radial_R2=float(f['coherence2']-r['coherence2']),
                                      radial_minus_full_phase_rms=float(r['soma_phase_rms']-f['soma_phase_rms']))

    return dict(seed=int(m.cfg.seed),cells=len(cells),shells=len(groups),
                full_R2=ev_full['coherence2'],radial_R2=ev_rad['coherence2'],
                full_minus_radial=float(ev_full['coherence2']-ev_rad['coherence2']),
                per_frequency=per_freq,states=states)


def summarize(rows,omegas):
    out=dict(bodies=len(rows),mean_full_minus_radial=float(np.mean([r['full_minus_radial'] for r in rows])))
    freq={}
    for om in omegas:
        k=str(float(om));z={}
        z['mean_full_minus_radial_R2']=float(np.mean([r['per_frequency'][k]['full_minus_radial_R2'] for r in rows]))
        for state in ('uniform','radial'):
            d=[r['states'][state][k]['decomposition'] for r in rows]
            t=[r['states'][state][k]['within_topology'] for r in rows]
            z[state]=dict(
                mean_within_fraction=float(np.mean([q['within_fraction'] for q in d])),
                mean_within_energy_fraction=float(np.mean([q['within_energy_fraction'] for q in d])),
                mean_radial_energy_fraction=float(np.mean([q['radial_energy_fraction'] for q in d])),
                mean_branch_eta2=float(np.nanmean([q['branch_eta2'] for q in t])),
                mean_rho_degree=float(np.nanmean([q['degree'] for q in t])),
                mean_rho_subtree=float(np.nanmean([q['subtree'] for q in t])),
                mean_rho_terminal_distance=float(np.nanmean([q['distance_to_terminal'] for q in t])),
                mean_rho_shell_size=float(np.nanmean([q['shell_size'] for q in t])),
            )
        freq[k]=z
    out['frequency']=freq
    return out


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=622);ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2.0);ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10.0)
    ap.add_argument('--full-steps',type=int,default=50);ap.add_argument('--radial-steps',type=int,default=160)
    ap.add_argument('--step-fraction',type=float,default=.10)
    ap.add_argument('--out',default='runs/material_gradient_subspace/dev_622_627.json')
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
        r=body_probe(m,oms,a.tau,a.mu,a.g0,a.ratio,a.full_steps,a.radial_steps,a.step_fraction);rows.append(r)
        print(f"seed {seed}: full-radial={r['full_minus_radial']:+.5f} " +
              ' '.join(f"w{om:.2f}:{r['per_frequency'][str(om)]['full_minus_radial_R2']:+.5f}" for om in oms))
    out=dict(config=vars(a),rows=rows,summary=summarize(rows,oms))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nMATERIAL GRADIENT SUBSPACE')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
