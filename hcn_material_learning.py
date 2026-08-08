"""Learn quasi-active material density from a uniform start using local adjoint sensitivity.

Discovery experiment.  Nothing in the objective or update rule contains graph
distance or a desired distal gradient.

For the confirmed quasi-active material model we optimize only somatic phase
coherence across many harmonic injection locations at omega=.03/.04.  Local
material density is nonnegative, capped at the hand-gradient maximum, and the
total material budget is fixed to the amount used by the confirmed hand-drawn
profile.

After learning (never during it), we ask whether the learned density is related
to graph distance from the soma.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import splu
from scipy.stats import spearmanr

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites, circular_rms


def ah_discrete(dt,tau,omega):
    lam=np.exp(1j*float(omega))
    return (dt/float(tau))/(lam-1.0+dt/float(tau))


def qa_operator(m,L,density,omega,tau,mu):
    c=m.cfg;dt=float(c.dt);lam=np.exp(1j*float(omega))
    ah=ah_discrete(dt,tau,omega)
    dyn=((lam-1.0+dt*float(c.damping))*(lam-1.0)/(dt*dt*lam)+float(c.restoring))
    local=np.asarray(density,float).ravel()*(1.0+float(mu)*ah)
    N=m.body.size
    A=(dyn*eye(N,dtype=np.complex128,format='csr')+
       diags(local,0,shape=(N,N),dtype=np.complex128)-float(c.stiffness)*L)
    return A,complex(1.0+float(mu)*ah)


def bounded_simplex_project(v,total,cap,atol=1e-12):
    """Euclidean projection onto 0<=x<=cap, sum x=total via scalar bisection."""
    v=np.asarray(v,float)
    total=float(total);cap=float(cap)
    if total < -atol or total > cap*len(v)+atol:
        raise ValueError('infeasible bounded simplex')
    lo=float(np.min(v-cap))-1.0
    hi=float(np.max(v))+1.0
    for _ in range(100):
        mid=(lo+hi)/2
        x=np.clip(v-mid,0.0,cap)
        if x.sum()>total:lo=mid
        else:hi=mid
    x=np.clip(v-(lo+hi)/2,0.0,cap)
    # Tiny residual correction over free coordinates.
    res=total-x.sum()
    free=np.where((x>1e-14)&(x<cap-1e-14))[0]
    if abs(res)>1e-10 and len(free):
        x[free]+=res/len(free);x=np.clip(x,0,cap)
    return x


def solve_fields(m,L,density,omega,tau,mu,sites):
    A,cfac=qa_operator(m,L,density,omega,tau,mu)
    lu=splu(A.tocsc())
    h,w=m.body.shape;N=h*w
    ids=np.asarray([p[0]*w+p[1] for p in sites],int)
    B=np.zeros((N,len(sites)),np.complex128);B[ids,np.arange(len(sites))]=1.0
    X=lu.solve(B)
    soma_id=m.soma[0]*w+m.soma[1]
    H=X[soma_id,:]
    # Transpose field. A is complex symmetric in this model, but use trans='T'
    # explicitly so the code mirrors the general derivative.
    es=np.zeros(N,np.complex128);es[soma_id]=1.0
    y=lu.solve(es,trans='T')
    return H,X,y,cfac


def coherence_and_gradient(m,L,density,omegas,tau,mu,sites,var_ids,need_grad=True):
    vals=[];grads=[];details=[]
    eps=1e-12
    for om in omegas:
        H,X,y,cfac=solve_fields(m,L,density,om,tau,mu,sites)
        r=np.maximum(np.abs(H),eps);u=H/r;mbar=np.mean(u)
        R2=float(abs(mbar)**2);vals.append(R2)
        details.append(dict(omega=float(om),coherence2=R2,
                            soma_phase_rms=circular_rms(np.angle(H)),
                            median_amp=float(np.median(np.abs(H)))))
        if need_grad:
            # dH[j,i]/dd_i = -c * y_i * X[i,j]
            xv=X[var_ids,:].T
            dH=-cfac*xv*y[var_ids][None,:]
            dr=np.real(np.conj(u)[:,None]*dH)
            du=dH/r[:,None]-u[:,None]*dr/r[:,None]
            dm=np.mean(du,axis=0)
            dR=2.0*np.real(np.conj(mbar)*dm)
            grads.append(dR)
    obj=float(np.mean(vals))
    if need_grad:return obj,np.mean(np.asarray(grads),axis=0),details
    return obj,details


def finite_diff_gradient_check(m,L,density,omegas,tau,mu,sites,var_ids,rng,ncheck=5,eps=1e-6):
    obj,g,_=coherence_and_gradient(m,L,density,omegas,tau,mu,sites,var_ids,True)
    picks=rng.choice(len(var_ids),size=min(int(ncheck),len(var_ids)),replace=False)
    rows=[]
    flat=np.asarray(density,float).ravel()
    for k in picks:
        idx=int(var_ids[k]);p=flat.copy();q=flat.copy();p[idx]+=eps;q[idx]-=eps
        op,_=coherence_and_gradient(m,L,p.reshape(m.body.shape),omegas,tau,mu,sites,var_ids,False)
        om,_=coherence_and_gradient(m,L,q.reshape(m.body.shape),omegas,tau,mu,sites,var_ids,False)
        fd=(op-om)/(2*eps);an=float(g[k])
        rel=float(abs(fd-an)/(abs(fd)+abs(an)+1e-12))
        rows.append(dict(var_index=int(k),grid_index=idx,analytic=an,finite_difference=float(fd),relative_error=rel))
    return dict(objective=obj,max_relative_error=float(max(x['relative_error'] for x in rows)),rows=rows)


def train_density(m,L,initial,omegas,tau,mu,sites,var_ids,total,cap,steps,step_fraction):
    d=np.asarray(initial,float).copy();hist=[]
    best_obj=-np.inf;best=d.copy()
    for it in range(int(steps)+1):
        obj,g,det=coherence_and_gradient(m,L,d,omegas,tau,mu,sites,var_ids,True)
        hist.append(dict(iteration=it,objective=obj,frequency=det,
                         grad_max=float(np.max(np.abs(g))),density_min=float(d.ravel()[var_ids].min()),
                         density_max=float(d.ravel()[var_ids].max())))
        if obj>best_obj:best_obj=obj;best=d.copy()
        if it==steps:break
        gm=float(np.max(np.abs(g)))
        if gm<1e-14:break
        base=float(step_fraction)*float(cap)/gm
        old=d.ravel()[var_ids].copy();accepted=False
        for bt in range(12):
            trialv=bounded_simplex_project(old+base*(.5**bt)*g,total,cap)
            trial=d.copy().ravel();trial[var_ids]=trialv;trial=trial.reshape(d.shape)
            tobj,_=coherence_and_gradient(m,L,trial,omegas,tau,mu,sites,var_ids,False)
            if tobj>=obj-1e-12:
                d=trial;accepted=True;break
        if not accepted:break
    return best,best_obj,hist


def morphology_stats(m,density,var_cells,var_ids):
    distmap=m.graph_distance_from_soma();dist=np.asarray([distmap[p] for p in var_cells],float)
    den=np.asarray(density,float).ravel()[var_ids]
    rho=float(spearmanr(dist,den).statistic) if np.std(den)>1e-12 else float('nan')
    pear=float(np.corrcoef(dist,den)[0,1]) if np.std(den)>1e-12 else float('nan')
    q25,q75=np.quantile(dist,[.25,.75]);prox=dist<=q25;far=dist>=q75
    return dict(spearman_distance=rho,pearson_distance=pear,
                proximal_mean=float(np.mean(den[prox])),distal_mean=float(np.mean(den[far])),
                distal_over_proximal=float(np.mean(den[far])/(np.mean(den[prox])+1e-30)),
                min_density=float(np.min(den)),max_density=float(np.max(den)),mean_density=float(np.mean(den)))


def eval_profile(m,L,density,omegas,tau,mu,sites):
    obj,det=coherence_and_gradient(m,L,density,omegas,tau,mu,sites,np.array([],int),False)
    return dict(coherence2=obj,frequency=det,
                mean_phase_rms=float(np.mean([x['soma_phase_rms'] for x in det])),
                mean_median_amp=float(np.mean([x['median_amp'] for x in det])))


def controls(m,L,learned,uniform,hand,omegas,tau,mu,sites,var_cells,var_ids,rng,nshuffle):
    out=dict(uniform=eval_profile(m,L,uniform,omegas,tau,mu,sites),
             learned=eval_profile(m,L,learned,omegas,tau,mu,sites),
             hand=eval_profile(m,L,hand,omegas,tau,mu,sites))
    vals=learned.ravel()[var_ids].copy()
    sh=[]
    for _ in range(int(nshuffle)):
        vv=vals.copy();rng.shuffle(vv);z=np.zeros_like(learned).ravel();z[var_ids]=vv
        sh.append(eval_profile(m,L,z.reshape(learned.shape),omegas,tau,mu,sites))
    out['shuffle_learned']=dict(coherence2=float(np.mean([x['coherence2'] for x in sh])),
                                mean_phase_rms=float(np.mean([x['mean_phase_rms'] for x in sh])),
                                mean_median_amp=float(np.mean([x['mean_median_amp'] for x in sh])))
    # Same learned histogram, reversed graph-distance rank assignment.
    distmap=m.graph_distance_from_soma();order=np.argsort([distmap[p] for p in var_cells])
    sorted_vals=np.sort(vals)[::-1];rv=np.empty_like(vals);rv[order]=sorted_vals
    z=np.zeros_like(learned).ravel();z[var_ids]=rv
    out['reverse_learned']=eval_profile(m,L,z.reshape(learned.shape),omegas,tau,mu,sites)
    return out


def body_probe(m,omegas,tau,mu,g0,ratio,steps,step_fraction,nshuffle):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    body=m.body.astype(bool);h,w=body.shape
    var_cells=[tuple(map(int,p)) for p in np.argwhere(body)]
    var_ids=np.asarray([p[0]*w+p[1] for p in var_cells],int)
    dist=m.graph_distance_from_soma().astype(float);dmax=max(float(dist[body].max()),1.0)
    hand=np.zeros_like(dist,float);hand[body]=float(g0)*(1+(float(ratio)-1)*dist[body]/dmax)
    total=float(hand[body].sum());cap=float(g0)*float(ratio)
    uniform=np.zeros_like(hand);uniform[body]=total/len(var_ids)
    sites=injection_sites(m)
    rng=np.random.default_rng(int(m.cfg.seed)+606_606)
    gc=finite_diff_gradient_check(m,L,uniform,omegas,tau,mu,sites,var_ids,rng)
    learned,best,hist=train_density(m,L,uniform,omegas,tau,mu,sites,var_ids,total,cap,steps,step_fraction)
    ctrl=controls(m,L,learned,uniform,hand,omegas,tau,mu,sites,var_cells,var_ids,rng,nshuffle)
    morph=morphology_stats(m,learned,var_cells,var_ids)
    return dict(seed=int(m.cfg.seed),cells=int(body.sum()),sites=len(sites),budget=total,cap=cap,
                gradient_check=gc,best_objective=best,history=hist,controls=ctrl,morphology=morph,
                learned_density=[dict(cell=[p[0],p[1]],distance=int(dist[p]),density=float(learned[p])) for p in var_cells])


def summarize(rows):
    def M(path):
        vals=[]
        for r in rows:
            q=r
            for k in path:q=q[k]
            vals.append(float(q))
        return float(np.mean(vals))
    return dict(
        bodies=len(rows),
        max_gradient_check_relative_error=float(max(r['gradient_check']['max_relative_error'] for r in rows)),
        uniform_coherence2=M(('controls','uniform','coherence2')),
        learned_coherence2=M(('controls','learned','coherence2')),
        hand_coherence2=M(('controls','hand','coherence2')),
        shuffled_learned_coherence2=M(('controls','shuffle_learned','coherence2')),
        learned_phase_rms=M(('controls','learned','mean_phase_rms')),
        uniform_phase_rms=M(('controls','uniform','mean_phase_rms')),
        hand_phase_rms=M(('controls','hand','mean_phase_rms')),
        shuffled_learned_phase_rms=M(('controls','shuffle_learned','mean_phase_rms')),
        mean_spearman_distance=M(('morphology','spearman_distance')),
        mean_pearson_distance=M(('morphology','pearson_distance')),
        mean_distal_over_proximal=M(('morphology','distal_over_proximal')),
        positive_spearman_bodies=int(sum(r['morphology']['spearman_distance']>0 for r in rows)),
        learned_beats_uniform_bodies=int(sum(r['controls']['learned']['coherence2']>r['controls']['uniform']['coherence2'] for r in rows)),
        learned_beats_shuffled_bodies=int(sum(r['controls']['learned']['coherence2']>r['controls']['shuffle_learned']['coherence2'] for r in rows)),
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
    ap.add_argument('--nshuffle',type=int,default=12)
    ap.add_argument('--out',default='runs/hcn_material_learning/dev.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    x=bounded_simplex_project(np.array([-.2,.3,1.7]),1.0,.6)
    assert np.all(x>=-1e-12) and np.all(x<=.6+1e-12) and abs(x.sum()-1)<1e-9
    print('selftest ok')


def main():
    a=parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    omegas=[float(x) for x in a.omegas.split(',')]
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,omegas,a.tau,a.mu,a.g0,a.ratio,a.steps,a.step_fraction,a.nshuffle);rows.append(r)
        c=r['controls'];mo=r['morphology']
        print(f"seed {seed}: R2 {c['uniform']['coherence2']:.4f}->{c['learned']['coherence2']:.4f} "
              f"hand={c['hand']['coherence2']:.4f} shuf={c['shuffle_learned']['coherence2']:.4f} "
              f"rho_dist={mo['spearman_distance']:+.3f} distal/prox={mo['distal_over_proximal']:.2f} "
              f"graderr={r['gradient_check']['max_relative_error']:.2e}",flush=True)
    if not rows:raise SystemExit('No valid bodies')
    s=summarize(rows)
    payload=dict(experiment='adjoint_material_learning_discovery_v01',development_only=True,
                 objective='maximize mean squared circular coherence of soma transfer phases; no distance term',
                 constraints='nonnegative density, cap=confirmed hand-profile maximum, fixed total density budget',
                 seed_start=a.seed_start,seeds_requested=a.seeds,omegas=omegas,tau_h=a.tau,mu=a.mu,
                 g0_budget_reference=a.g0,ratio_budget_reference=a.ratio,steps=a.steps,step_fraction=a.step_fraction,
                 summary=s,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nMATERIAL LEARNING DISCOVERY')
    print(json.dumps(s,indent=2))

if __name__=='__main__':main()
