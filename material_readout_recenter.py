"""Moved-readout control for material-adjoint learning.

The original material learner discovers density increasing with graph distance
from the soma, but the training objective is defined at the soma.  This control
moves the consequential readout to the occupied cell farthest from the soma and
repeats the *same* phase-coherence learning from a uniform material state.

No distance coordinate enters training.  After learning, compare density with
both distance from the moved readout and distance from the anatomical soma.

If the material coordinate re-centers on the moved readout, the effect is a
readout-relative property of the transfer geometry rather than a hidden soma
prior.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import splu
from scipy.stats import spearmanr

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import weighted_laplacian_sparse, circular_rms
from hcn_material_learning import bounded_simplex_project


def n4(p, shape):
    y,x=p;h,w=shape
    for q in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
        if 0<=q[0]<h and 0<=q[1]<w:yield q


def graph_distance(body,start):
    b=np.asarray(body,bool);out=np.full(b.shape,10**9,int)
    start=tuple(map(int,start));out[start]=0;q=deque([start])
    while q:
        p=q.popleft();d=out[p]+1
        for r in n4(p,b.shape):
            if b[r] and d<out[r]:out[r]=d;q.append(r)
    return out


def choose_far_readout(m):
    ds=graph_distance(m.body,m.soma);cells=np.argwhere(m.body>0)
    vals=np.asarray([ds[tuple(p)] for p in cells])
    mx=int(vals.max());cand=[tuple(map(int,p)) for p,v in zip(cells,vals) if int(v)==mx]
    cand.sort();return cand[0]


def ah_discrete(dt,tau,omega):
    lam=np.exp(1j*float(omega));return (dt/float(tau))/(lam-1.0+dt/float(tau))


def qa_operator(m,L,density,omega,tau,mu):
    c=m.cfg;dt=float(c.dt);lam=np.exp(1j*float(omega));ah=ah_discrete(dt,tau,omega)
    dyn=((lam-1.0+dt*float(c.damping))*(lam-1.0)/(dt*dt*lam)+float(c.restoring))
    local=np.asarray(density,float).ravel()*(1.0+float(mu)*ah);N=m.body.size
    return (dyn*eye(N,dtype=np.complex128,format='csr')+
            diags(local,0,shape=(N,N),dtype=np.complex128)-float(c.stiffness)*L),complex(1+float(mu)*ah)


def sites_for_readout(m,readout,min_dist=2):
    d=graph_distance(m.body,readout)
    pts=[tuple(map(int,p)) for p in np.argwhere((m.body>0)&(d>=int(min_dist)))]
    pts.sort(key=lambda p:(int(d[p]),p[0],p[1]));return pts,d


def solve_fields(m,L,density,omega,tau,mu,sites,readout):
    A,c=qa_operator(m,L,density,omega,tau,mu);lu=splu(A.tocsc())
    h,w=m.body.shape;N=h*w;ids=np.asarray([p[0]*w+p[1] for p in sites],int)
    B=np.zeros((N,len(sites)),np.complex128);B[ids,np.arange(len(sites))]=1
    X=lu.solve(B);rid=readout[0]*w+readout[1];H=X[rid,:]
    e=np.zeros(N,np.complex128);e[rid]=1;y=lu.solve(e,trans='T')
    return H,X,y,c


def objective_gradient(m,L,density,omegas,tau,mu,sites,readout,var_ids,grad=True):
    objs=[];gs=[];details=[];eps=1e-12
    for om in omegas:
        H,X,y,c=solve_fields(m,L,density,om,tau,mu,sites,readout)
        r=np.maximum(np.abs(H),eps);u=H/r;mb=np.mean(u);R2=float(abs(mb)**2);objs.append(R2)
        details.append(dict(omega=float(om),coherence2=R2,soma_phase_rms=circular_rms(np.angle(H)),median_amp=float(np.median(np.abs(H)))))
        if grad:
            dH=-c*X[var_ids,:].T*y[var_ids][None,:]
            dr=np.real(np.conj(u)[:,None]*dH);du=dH/r[:,None]-u[:,None]*dr/r[:,None]
            dm=np.mean(du,axis=0);gs.append(2*np.real(np.conj(mb)*dm))
    if grad:return float(np.mean(objs)),np.mean(gs,axis=0),details
    return float(np.mean(objs)),details


def train(m,L,initial,omegas,tau,mu,sites,readout,var_ids,total,cap,steps,step_fraction):
    d=initial.copy();best=d.copy();bestobj=-np.inf;hist=[]
    for it in range(int(steps)+1):
        obj,g,det=objective_gradient(m,L,d,omegas,tau,mu,sites,readout,var_ids,True)
        hist.append(dict(iteration=it,objective=obj,frequency=det))
        if obj>bestobj:bestobj=obj;best=d.copy()
        if it==steps:break
        gm=float(np.max(np.abs(g)))
        if gm<1e-14:break
        old=d.ravel()[var_ids].copy();base=float(step_fraction)*cap/gm;accepted=False
        for bt in range(12):
            v=bounded_simplex_project(old+base*(.5**bt)*g,total,cap)
            z=d.copy().ravel();z[var_ids]=v;z=z.reshape(d.shape)
            no,_=objective_gradient(m,L,z,omegas,tau,mu,sites,readout,var_ids,False)
            if no>=obj-1e-12:d=z;accepted=True;break
        if not accepted:break
    return best,bestobj,hist


def density_stats(m,density,readout,var_cells,var_ids):
    dr=graph_distance(m.body,readout);ds=graph_distance(m.body,m.soma)
    den=density.ravel()[var_ids];xr=np.asarray([dr[p] for p in var_cells],float);xs=np.asarray([ds[p] for p in var_cells],float)
    rr=float(spearmanr(xr,den).statistic);rs=float(spearmanr(xs,den).statistic)
    return dict(spearman_to_readout=rr,spearman_to_soma=rs,
                pearson_to_readout=float(np.corrcoef(xr,den)[0,1]),pearson_to_soma=float(np.corrcoef(xs,den)[0,1]),
                mean_density=float(np.mean(den)),zero_fraction=float(np.mean(den<1e-8)),cap_fraction=float(np.mean(den>.049999)))


def run_readout(m,L,readout,omegas,tau,mu,total,cap,steps,step_fraction):
    body=m.body.astype(bool);h,w=body.shape;cells=[tuple(map(int,p)) for p in np.argwhere(body)];ids=np.asarray([p[0]*w+p[1] for p in cells],int)
    uniform=np.zeros_like(body,float);uniform[body]=total/len(cells)
    sites,_=sites_for_readout(m,readout)
    before,_=objective_gradient(m,L,uniform,omegas,tau,mu,sites,readout,ids,False)
    learned,best,hist=train(m,L,uniform,omegas,tau,mu,sites,readout,ids,total,cap,steps,step_fraction)
    st=density_stats(m,learned,readout,cells,ids)
    return dict(readout=list(map(int,readout)),sites=len(sites),uniform_coherence2=before,learned_coherence2=best,
                gain=float(best-before),stats=st,history=hist,
                density=[dict(cell=[p[0],p[1]],density=float(learned[p]),distance_readout=int(graph_distance(m.body,readout)[p]),distance_soma=int(graph_distance(m.body,m.soma)[p])) for p in cells])


def body_probe(m,omegas,tau,mu,g0,ratio,steps,step_fraction):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv);body=m.body.astype(bool)
    ds=graph_distance(m.body,m.soma).astype(float);dmax=max(float(ds[body].max()),1)
    hand=np.zeros_like(ds,float);hand[body]=float(g0)*(1+(float(ratio)-1)*ds[body]/dmax)
    total=float(hand[body].sum());cap=float(g0)*float(ratio)
    far=choose_far_readout(m)
    soma=run_readout(m,L,tuple(m.soma),omegas,tau,mu,total,cap,steps,step_fraction)
    moved=run_readout(m,L,far,omegas,tau,mu,total,cap,steps,step_fraction)
    return dict(seed=int(m.cfg.seed),cells=int(body.sum()),far_readout=list(far),soma=soma,moved=moved)


def summarize(rows):
    def M(which,key):return float(np.mean([r[which][key] for r in rows]))
    return dict(bodies=len(rows),
                soma_mean_gain=M('soma','gain'),moved_mean_gain=M('moved','gain'),
                soma_mean_rho_readout=float(np.mean([r['soma']['stats']['spearman_to_readout'] for r in rows])),
                moved_mean_rho_readout=float(np.mean([r['moved']['stats']['spearman_to_readout'] for r in rows])),
                moved_mean_rho_anatomical_soma=float(np.mean([r['moved']['stats']['spearman_to_soma'] for r in rows])),
                moved_positive_readout_rho=int(sum(r['moved']['stats']['spearman_to_readout']>0 for r in rows)),
                moved_negative_soma_rho=int(sum(r['moved']['stats']['spearman_to_soma']<0 for r in rows)),
                soma_positive_readout_rho=int(sum(r['soma']['stats']['spearman_to_readout']>0 for r in rows)))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=580);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--omegas',default='0.03,0.04');ap.add_argument('--tau',type=float,default=2);ap.add_argument('--mu',type=float,default=.5);ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10)
    ap.add_argument('--steps',type=int,default=50);ap.add_argument('--step-fraction',type=float,default=.1);ap.add_argument('--out',default='runs/material_readout_recenter/dev.json');a=ap.parse_args()
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    omegas=[float(x) for x in a.omegas.split(',')];rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=body_probe(m,omegas,a.tau,a.mu,a.g0,a.ratio,a.steps,a.step_fraction);rows.append(r)
        print(f"seed {seed}: soma rho={r['soma']['stats']['spearman_to_readout']:+.3f} gain={r['soma']['gain']:+.3f}; "
              f"moved {r['far_readout']} rho(new)={r['moved']['stats']['spearman_to_readout']:+.3f} "
              f"rho(old soma)={r['moved']['stats']['spearman_to_soma']:+.3f} gain={r['moved']['gain']:+.3f}",flush=True)
    s=summarize(rows);payload=dict(experiment='material_readout_recenter_dev_v01',development_only=True,summary=s,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8');print('\nREADOUT RECENTER');print(json.dumps(s,indent=2))

if __name__=='__main__':main()
