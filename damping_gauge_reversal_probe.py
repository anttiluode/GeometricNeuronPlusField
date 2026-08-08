"""Can uniform damping be compiled out so the internal forward field retraces physically?

Development-only hardware mechanism probe.

The mature reciprocal wave update is

    v[n+1]   = a v[n] + dt A0 psi[n] + dt source[n]
    psi[n+1] = psi[n] + dt v[n+1]

with a = 1-dt*damping and A0 = stiffness*L-restoring*I.
Eliminating velocity gives

    psi[n+1] = M psi[n] - a psi[n-1] + dt^2 source[n].

For uniform damping, let r=sqrt(a) and psi[n]=r^n z[n].  Then EXACTLY

    z[n+1] = Q z[n] - z[n-1] + u[n]
    Q = M/r
    u[n] = dt^2 r^(-(n+1)) source[n].

The transformed recurrence has unit reverse coefficient and is time-reversal
symmetric.  Damping has moved to known global temporal envelopes on source and
readout.

This probe tests four things against the original exact discrete model:

1. Gauge forward identity: r^n z[n] reproduces psi[n].
2. Gauge gradient identity: reverse-mode gradient in z coordinates reproduces
   the original bond gradient.
3. Physical retracing identity: from the final transformed state, the same Q
   recurrence plus reversed transformed forcing reconstructs z in reverse.
4. Interference identity: the exact bond gradient is reconstructed from only
   integrated branch energies of two reverse trials carrying retraced-forward
   +/- transformed-adjoint fields.  No T-sample local forward history is used
   by this reconstruction.

The simulation is not yet a hardware claim.  It assumes that a physical
implementation can (a) realize Q as a reciprocal conservative recurrence,
(b) prepare the reverse initial momentum/state exactly, (c) apply the known
source/readout exponential envelopes, and (d) repeat the same final state for
the +/- trials.  Noise/pass cost are separate questions.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
from transfer_decomposition_probe import safe_corr


def flat_pair(h,v):
    return np.concatenate([np.ravel(h),np.ravel(v)])


def rel_l2(a,b):
    a=np.asarray(a);b=np.asarray(b)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))


def edge_diffs(x):
    x=np.asarray(x)
    return x[..., :, 1:]-x[..., :, :-1], x[..., 1:, :]-x[..., :-1, :]


def params(m):
    c=m.cfg
    dt=float(c.dt);stiff=float(c.stiffness);damp=float(c.damping);rest=float(c.restoring)
    a=1.0-dt*damp
    if a<=0:
        raise ValueError(f'gauge requires 1-dt*damping>0, got {a}')
    r=math.sqrt(a)
    beta=dt*dt*stiff/r
    q0=(1.0+a-dt*dt*rest)/r
    return dt,stiff,damp,rest,a,r,beta,q0


def q_apply(m,wh,wv,x):
    dt,stiff,damp,rest,a,r,beta,q0=params(m)
    return q0*np.asarray(x)+beta*ae.weighted_lap(x,wh,wv)


def gauge_forward(m,wh,wv,seq):
    """Return z[-1..T] logically as z0..zT plus transformed forcing u[0..T-1]."""
    dt,stiff,damp,rest,a,r,beta,q0=params(m)
    zm1=np.zeros(m.body.shape,np.complex128)
    z0=np.zeros_like(zm1)
    zs=[z0.copy()]
    us=[]
    for n,src in enumerate(seq):
        u=(dt*dt)*(r**(-(n+1)))*np.asarray(src,np.complex128)
        z1=q_apply(m,wh,wv,z0)-zm1+u
        us.append(u)
        zs.append(z1.copy())
        zm1,z0=z0,z1
    return np.asarray(zs),np.asarray(us)


def reconstruct_psi(z,r):
    shape=(len(z),)+(1,)*(z.ndim-1)
    scale=(r**np.arange(len(z))).reshape(shape)
    return scale*z


def gauge_objective_sources(m,z,coeff,r):
    """Complex reverse-mode source q[k] for J=coeff*sum |psi_s[k]|^2."""
    T=len(z)-1
    q=np.zeros_like(z)
    soma=tuple(map(int,m.soma))
    for k in range(1,T+1):
        # psi_k=r^k z_k and dJ=2 Re conj(q_k) dz_k.
        q[k][soma]=float(coeff)*(r**(2*k))*z[k][soma]
    return q


def gauge_adjoint_gradient(m,wh,wv,z,coeff):
    """Exact reverse-mode bond gradient for transformed second-order recurrence.

    p[k] is total adjoint of z[k] after adding its objective derivative.
    Transition n uses z[n] -> z[n+1], so local parameter credit pairs z[n]
    with p[n+1].
    """
    dt,stiff,damp,rest,a,r,beta,q0=params(m)
    T=len(z)-1
    q=gauge_objective_sources(m,z,coeff,r)
    # p has states 0..T+1; p[T+1]=0.
    p=np.zeros((T+2,)+m.body.shape,np.complex128)
    gh=np.zeros_like(wh,float);gv=np.zeros_like(wv,float)
    for k in range(T,0,-1):
        p[k]+=q[k]
        mu=p[k]
        forward=z[k-1]
        dfh=forward[:,1:]-forward[:,:-1]
        dmh=mu[:,:-1]-mu[:,1:]
        dfv=forward[1:,:]-forward[:-1,:]
        dmv=mu[:-1,:]-mu[1:,:]
        gh += 2.0*beta*np.real(np.conj(dmh)*dfh)
        gv += 2.0*beta*np.real(np.conj(dmv)*dfv)
        p[k-1]+=q_apply(m,wh,wv,mu)
        if k-2>=0:
            p[k-2]-=mu
        # p[-1] corresponds fixed z[-1] and is irrelevant.
    return gh,gv,p,q


def reverse_forward_states(m,wh,wv,z,u):
    """Retrace transformed forward field with same Q.

    Returns w[j]=z[T-j], j=0..T.  The first two states are the terminal
    position and reversed-momentum state.  Subsequent evolution uses the
    reversed transformed forcing.
    """
    T=len(z)-1
    w=[z[T].copy(),z[T-1].copy()]
    for j in range(0,T-1):
        # w[j+2] uses the original transformed force from transition T-1-j.
        nxt=q_apply(m,wh,wv,w[j+1])-w[j]+u[T-1-j]
        w.append(nxt)
    return np.asarray(w)


def causal_adjoint_states(m,wh,wv,q):
    """Generate a[j]=p[T-j+1]?  Specifically a[j]=p[T-j+1] for j>=1.

    We return A[0..T] with A[0]=0 and A[j]=p[T-j+1] for j=1..T.
    It obeys the same Q recurrence from zero initial state, driven by q in
    reverse order.  This alignment makes A[j] pair with W[j]=z[T-j].
    """
    T=len(q)-1
    am1=np.zeros(m.body.shape,np.complex128)
    a0=np.zeros_like(am1)
    arr=[a0.copy()]
    for j in range(T):
        src=q[T-j]
        a1=q_apply(m,wh,wv,a0)-am1+src
        arr.append(a1.copy())
        am1,a0=a0,a1
    return np.asarray(arr)


def interference_gradient(m,w,a):
    """Gradient from two integrated +/- branch-energy trials.

    At reverse index j=1..T:
        w[j] = z[T-j]
        a[j] = p[T-j+1]
    which is exactly the forward/adjoint pair in the transformed gradient.

    The physical trials have branch fields w+a and w-a.  Their integrated
    square-law difference gives 4 Re(conj(Delta a) Delta w).  The sign below
    accounts for the edge orientation used by ae.adjoint_grad.
    """
    *_,beta,q0=params(m)
    wp=w[1:]+a[1:]
    wm=w[1:]-a[1:]
    ph,pv=edge_diffs(wp)
    mh,mv=edge_diffs(wm)
    cross_h=0.25*np.sum(np.abs(ph)**2-np.abs(mh)**2,axis=0)
    cross_v=0.25*np.sum(np.abs(pv)**2-np.abs(mv)**2,axis=0)
    # ae convention uses dmu = mu_left-mu_right = -Delta(mu), while
    # edge_diffs uses Delta=right-left/down-up.
    gh=-2.0*beta*cross_h
    gv=-2.0*beta*cross_v
    return gh.real,gv.real


def spectral_radius_q(m,wh,wv):
    """Small-grid dense audit of Q eigenvalue range."""
    h,w=m.body.shape;N=h*w
    eye=np.eye(N,dtype=np.complex128).reshape((N,h,w))
    cols=[]
    for i in range(N):
        cols.append(q_apply(m,wh,wv,eye[i]).ravel())
    Q=np.stack(cols,axis=1)
    ev=np.linalg.eigvalsh(Q.real)
    return float(ev.min()),float(ev.max()),float(np.max(np.abs(ev)))


def one_order(m,wh,wv,seq,coeff):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    z,u=gauge_forward(m,wh,wv,seq)
    *_,r,_,_=params(m)
    pr=reconstruct_psi(z,r)
    fwd_err=rel_l2(p,pr)

    gh0,gv0=ae.adjoint_grad(m,wh,wv,p,v,coeff)
    gh1,gv1,padj,q=gauge_adjoint_gradient(m,wh,wv,z,coeff)
    grad_err=rel_l2(flat_pair(gh0,gv0),flat_pair(gh1,gv1))
    grad_corr=float(safe_corr(flat_pair(gh0,gv0),flat_pair(gh1,gv1)))

    wrev=reverse_forward_states(m,wh,wv,z,u)
    reverse_err=rel_l2(wrev,z[::-1])
    a=causal_adjoint_states(m,wh,wv,q)
    # Compare generated causal adjoint against the backward p alignment.
    T=len(z)-1
    target=np.zeros_like(a)
    for j in range(1,T+1):
        target[j]=padj[T-j+1]
    adj_replay_err=rel_l2(target,a)

    ghi,gvi=interference_gradient(m,wrev,a)
    int_err=rel_l2(flat_pair(gh1,gv1),flat_pair(ghi,gvi))
    int_corr=float(safe_corr(flat_pair(gh1,gv1),flat_pair(ghi,gvi)))
    return dict(E=float(E),forward_gauge_relative_l2=fwd_err,
                gradient_gauge_relative_l2=grad_err,gradient_gauge_corr=grad_corr,
                reverse_retrace_relative_l2=reverse_err,
                adjoint_replay_relative_l2=adj_replay_err,
                interference_gradient_relative_l2=int_err,
                interference_gradient_corr=int_corr,
                original=(gh0,gv0),interference=(ghi,gvi))


def one(m,lag,steps):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,lag,steps)
    seqD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,seqT,store=False)
    ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30
    aT=2.0*ED/(S*S);aD=-2.0*ET/(S*S)
    T=one_order(m,wh,wv,seqT,aT)
    D=one_order(m,wh,wv,seqD,aD)
    exact_h=T['original'][0]+D['original'][0]
    exact_v=T['original'][1]+D['original'][1]
    phys_h=T['interference'][0]+D['interference'][0]
    phys_v=T['interference'][1]+D['interference'][1]
    comb_ref=flat_pair(exact_h,exact_v);comb=flat_pair(phys_h,phys_v)
    qmin,qmax,qrho=spectral_radius_q(m,wh,wv)
    dt,stiff,damp,rest,a,r,beta,q0=params(m)
    def strip(x):
        return {k:v for k,v in x.items() if k not in ('original','interference')}
    return dict(seed=int(m.cfg.seed),C=float((ET-ED)/S),
                dt=dt,damping=damp,a=a,r=r,beta=beta,q0=q0,
                Q_eigen_min=qmin,Q_eigen_max=qmax,Q_spectral_radius=qrho,
                target=strip(T),distractor=strip(D),
                combined_interference_corr=float(safe_corr(comb_ref,comb)),
                combined_interference_relative_l2=rel_l2(comb_ref,comb))


def selftest():
    # Scalar recurrence reversal independent of neuron implementation.
    rng=np.random.default_rng(5);T=20;Q=1.4
    u=rng.normal(size=T)
    z0=0.;zm1=0.;z=[z0]
    for n in range(T):
        z1=Q*z0-zm1+u[n];z.append(z1);zm1,z0=z0,z1
    z=np.asarray(z)
    w=[z[-1],z[-2]]
    for j in range(T-1):w.append(Q*w[j+1]-w[j]+u[T-1-j])
    assert np.max(np.abs(np.asarray(w)-z[::-1]))<1e-10
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472)
    ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/damping_gauge_reversal/dev.json')
    ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=one(m,a.lag,a.steps);rows.append(r)
        print(f"seed {seed}: fwd={r['target']['forward_gauge_relative_l2']:.2e} "
              f"grad={r['target']['gradient_gauge_relative_l2']:.2e} "
              f"retrace={r['target']['reverse_retrace_relative_l2']:.2e} "
              f"adj={r['target']['adjoint_replay_relative_l2']:.2e} "
              f"interf={r['combined_interference_relative_l2']:.2e} corr={r['combined_interference_corr']:.12f} "
              f"Q=[{r['Q_eigen_min']:.4f},{r['Q_eigen_max']:.4f}]",flush=True)
    if not rows:raise SystemExit('No valid bodies')
    keys=['forward_gauge_relative_l2','gradient_gauge_relative_l2','reverse_retrace_relative_l2','adjoint_replay_relative_l2','interference_gradient_relative_l2']
    summary=dict(bodies=len(rows))
    for name in ('target','distractor'):
        summary[name]={k:dict(mean=float(np.mean([r[name][k] for r in rows])),max=float(np.max([r[name][k] for r in rows]))) for k in keys}
    summary['combined']=dict(
        mean_corr=float(np.mean([r['combined_interference_corr'] for r in rows])),
        min_corr=float(np.min([r['combined_interference_corr'] for r in rows])),
        mean_relative_l2=float(np.mean([r['combined_interference_relative_l2'] for r in rows])),
        max_relative_l2=float(np.max([r['combined_interference_relative_l2'] for r in rows])),
        Q_eigen_min=float(np.min([r['Q_eigen_min'] for r in rows])),
        Q_eigen_max=float(np.max([r['Q_eigen_max'] for r in rows])),
        r=float(np.mean([r['r'] for r in rows])),
    )
    out=dict(config=vars(a),rows=rows,summary=summary)
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nDAMPING GAUGE REVERSAL')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
