"""Development test of the reciprocity trade using a genuinely nonsymmetric operator.

A local real skew-symmetric nearest-neighbour coupling A is added to the ordinary
symmetric weighted Laplacian L:

    H_beta = L + beta A,       A^T = -A.

Therefore H_beta^T = H_-beta exactly.  The skew part is held fixed while we audit
sensitivity with respect to the ordinary symmetric bond conductances.  This cleanly
separates the credit-transport question from the structural coordinate definition.

For each body and beta we compare:
  exact      algorithmic adjoint using H_beta^T;
  same-H     time-reversed soma derivative replayed through H_beta again;
  transpose  identical replay through H_-beta (mathematical transpose positive control).

If reciprocity is the reason same-medium physical backprop works, same-H should peel
away as beta grows while transpose replay should remain exact.

This is a development probe.  The skew coupling is an abstract stable-ish local
nonreciprocity model, not a calibrated optical isolator/circulator model.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from reciprocal_adjoint_probe import retro_source_sequence, flat_pair, normalized_l2
from transfer_decomposition_probe import safe_corr

BETAS=(0.0,.02,.05,.10,.20,.30,.40,.60)


def skew_apply(u,ah,av):
    """Apply a real skew-symmetric nearest-neighbour matrix."""
    u=np.asarray(u,np.complex128);out=np.zeros_like(u)
    # Horizontal edge L,R: A_LR=+a, A_RL=-a.
    out[:,:-1]+=ah*u[:,1:]
    out[:,1:]-=ah*u[:,:-1]
    # Vertical edge U,D: A_UD=+a, A_DU=-a.
    out[:-1,:]+=av*u[1:,:]
    out[1:,:]-=av*u[:-1,:]
    return out


def spatial_apply(u,wh,wv,ah,av,beta):
    return ae.weighted_lap(u,wh,wv) + float(beta)*skew_apply(u,ah,av)


def forward(m,wh,wv,ah,av,beta,seq,store=True):
    c=m.cfg;dt=float(c.dt);stiff=float(c.stiffness);damp=float(c.damping);rest=float(c.restoring)
    psi=np.zeros(m.body.shape,np.complex128);vel=np.zeros_like(psi);energy=0.;soma=tuple(map(int,m.soma))
    if store:ps=[psi.copy()];vs=[vel.copy()]
    for src in seq:
        vel=vel+dt*(stiff*spatial_apply(psi,wh,wv,ah,av,beta)-damp*vel-rest*psi+src)
        psi=psi+dt*vel;energy+=float(abs(psi[soma])**2)
        if store:ps.append(psi.copy());vs.append(vel.copy())
    if store:return np.asarray(ps),np.asarray(vs),float(energy)
    return float(energy)


def exact_mu(m,wh,wv,ah,av,beta,ps,coeff):
    """Exact discrete adjoint mu[n] aligned to forward transitions n=0..T-1."""
    c=m.cfg;dt=float(c.dt);stiff=float(c.stiffness);damp=float(c.damping);rest=float(c.restoring);avd=1.-dt*damp
    lp=np.zeros(m.body.shape,np.complex128);lv=np.zeros_like(lp);soma=tuple(map(int,m.soma));T=len(ps)-1
    mus=np.zeros((T,)+m.body.shape,np.complex128)
    for k in range(T,0,-1):
        lp[soma]+=float(coeff)*ps[k][soma]
        mu=lv+dt*lp;mus[k-1]=mu
        lv_prev=avd*mu
        # transpose(H_beta)=H_-beta because L is symmetric and A skew-symmetric.
        lp_prev=lp+dt*(stiff*spatial_apply(mu,wh,wv,ah,av,-float(beta))-rest*mu)
        lp,lv=lp_prev,lv_prev
    return mus


def physical_mu(m,wh,wv,ah,av,beta,g):
    seq=retro_source_sequence(m,g,reverse=True)
    rp,rv,_=forward(m,wh,wv,ah,av,beta,seq,store=True);T=len(g)
    return np.asarray([rp[T-n] for n in range(T)],np.complex128)


def bond_grad_from_mu(m,ps,mu):
    """Gradient wrt the symmetric conductance coordinates; skew background fixed."""
    dt=float(m.cfg.dt);stiff=float(m.cfg.stiffness)
    gh=np.zeros((m.body.shape[0],m.body.shape[1]-1),float);gv=np.zeros((m.body.shape[0]-1,m.body.shape[1]),float)
    for n in range(len(mu)):
        prev=ps[n];q=mu[n]
        dpsi=prev[:,1:]-prev[:,:-1];dmu=q[:,:-1]-q[:,1:]
        gh+=2.*dt*stiff*np.real(np.conj(dmu)*dpsi)
        dpsi=prev[1:,:]-prev[:-1,:];dmu=q[:-1,:]-q[1:,:]
        gv+=2.*dt*stiff*np.real(np.conj(dmu)*dpsi)
    return gh,gv


def metrics(exh,exv,ah,av):
    ex=flat_pair(exh,exv);ap=flat_pair(ah,av)
    return dict(corr=float(safe_corr(ex,ap)),relative_l2=normalized_l2(ex,ap))


def body_beta(m,wh,wv,ah,av,beta,lag,steps):
    sT=ae.source_sequence(m,True,lag,steps);sD=ae.source_sequence(m,False,lag,steps)
    pT,vT,ET=forward(m,wh,wv,ah,av,beta,sT,True);pD,vD,ED=forward(m,wh,wv,ah,av,beta,sD,True)
    if not (np.isfinite(ET) and np.isfinite(ED)) or max(ET,ED)>1e12:
        return dict(stable=False,ET=float(ET),ED=float(ED))
    S=ET+ED+1e-30;aT=2.*ED/S**2;aD=-2.*ET/S**2
    gT=aT*np.asarray(pT[1:,m.soma[0],m.soma[1]],np.complex128);gD=aD*np.asarray(pD[1:,m.soma[0],m.soma[1]],np.complex128)
    xT=exact_mu(m,wh,wv,ah,av,beta,pT,aT);xD=exact_mu(m,wh,wv,ah,av,beta,pD,aD)
    sameT=physical_mu(m,wh,wv,ah,av,beta,gT);sameD=physical_mu(m,wh,wv,ah,av,beta,gD)
    trT=physical_mu(m,wh,wv,ah,av,-beta,gT);trD=physical_mu(m,wh,wv,ah,av,-beta,gD)
    exghT,exgvT=bond_grad_from_mu(m,pT,xT);exghD,exgvD=bond_grad_from_mu(m,pD,xD)
    sghT,sgvT=bond_grad_from_mu(m,pT,sameT);sghD,sgvD=bond_grad_from_mu(m,pD,sameD)
    tghT,tgvT=bond_grad_from_mu(m,pT,trT);tghD,tgvD=bond_grad_from_mu(m,pD,trD)
    exmu=np.concatenate([(xT+xD).ravel()]);smu=np.concatenate([(sameT+sameD).ravel()]);tmu=np.concatenate([(trT+trD).ravel()])
    return dict(stable=True,C=float((ET-ED)/S),
                same_field=dict(corr=float(safe_corr(exmu,smu)),relative_l2=normalized_l2(exmu,smu)),
                transpose_field=dict(corr=float(safe_corr(exmu,tmu)),relative_l2=normalized_l2(exmu,tmu)),
                same_grad=metrics(exghT+exghD,exgvT+exgvD,sghT+sghD,sgvT+sgvD),
                transpose_grad=metrics(exghT+exghD,exgvT+exgvD,tghT+tghD,tgvT+tgvD))


def one(m,lag,steps):
    wh,wv=ae.bond_weights(m,m.body);rng=np.random.default_rng(int(m.cfg.seed)+919191)
    sh=rng.choice(np.array([-1.,1.]),size=wh.shape);sv=rng.choice(np.array([-1.,1.]),size=wv.shape)
    # Scale skew couplings with the local symmetric edge scale; this sets a clean
    # dimensionless nonreciprocity strength beta.  The skew background itself is fixed.
    ah=wh*sh;av=wv*sv
    return dict(seed=int(m.cfg.seed),conditions={f'{b:g}':body_beta(m,wh,wv,ah,av,float(b),lag,steps) for b in BETAS})


def selftest():
    rng=np.random.default_rng(1);u=rng.normal(size=(4,5))+1j*rng.normal(size=(4,5));v=rng.normal(size=(4,5))+1j*rng.normal(size=(4,5))
    ah=rng.normal(size=(4,4));av=rng.normal(size=(3,5))
    Au=skew_apply(u,ah,av);Av=skew_apply(v,ah,av)
    # <v,Au> = -<Av,u> for real skew A.
    assert abs(np.vdot(v,Au)+np.vdot(Av,u))<1e-10
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=400);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--out',default='runs/nonreciprocal_adjoint/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps);rows.append(r)
        print('seed',seed,[(b, r['conditions'][f'{b:g}']['stable'], round(r['conditions'][f'{b:g}'].get('same_grad',{}).get('corr',float('nan')),4), round(r['conditions'][f'{b:g}'].get('transpose_grad',{}).get('corr',float('nan')),8)) for b in BETAS],flush=True)
    summary=dict(bodies=len(rows),betas=[float(x) for x in BETAS],conditions={})
    for b in BETAS:
        k=f'{b:g}';q=[r['conditions'][k] for r in rows if r['conditions'][k]['stable']]
        summary['conditions'][k]=dict(stable_bodies=len(q))
        if q:
            for name in ('same_field','transpose_field','same_grad','transpose_grad'):
                summary['conditions'][k][name]=dict(mean_corr=float(np.mean([x[name]['corr'] for x in q])),mean_relative_l2=float(np.mean([x[name]['relative_l2'] for x in q])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='nonreciprocal_adjoint_dev_v01',summary=summary,rows=rows),indent=2))
    print('\nNONRECIPROCAL ADJOINT DEV')
    for b in BETAS:
        q=summary['conditions'][f'{b:g}'];print('beta',b,'stable',q['stable_bodies'],'same-grad',q.get('same_grad'),'transpose-grad',q.get('transpose_grad'))
if __name__=='__main__':main()
