"""Physically consistent residual-loss test for the damping-gauge echo compiler.

Development only on reused bodies.

The first damage probe used

    x[n+1] = Q x[n] - (1-eps) x[n-1] + u[n]

when adding residual loss.  That is not the recurrence produced by the same
velocity integrator used in the repository: changing a=1-dt*gamma also changes
the coefficient of x[n].  With stiffness/restoring held fixed the consistent
recurrence is

    x[n+1] = (Q - eps I) x[n] - (1-eps) x[n-1] + u[n].

This matters because Q has eigenvalues close to +2; the malformed v0.1 damage
could push modes outside the discrete stability region by construction.

This probe asks the more useful question: if the compiled reversible core has a
small *real* uniform residual loss, how well does an ordinary same-body echo
estimate the exact gradient of that actual lossy device?

For every eps:
  1. use the original ideal gauge only to compile the intended external source
     waveform u[n];
  2. propagate the *actual* physical core with M_eps=Q-eps I, a_eps=1-eps;
  3. evaluate the same gauged output objective on that actual trajectory;
  4. compute its exact reverse-mode bond gradient digitally (reference only);
  5. attempt a physical echo with the same lossy operator, reversed transformed
     input, and reversed objective source;
  6. read the branch overlap by integrated +/- energy;
  7. compare echo gradient with the exact gradient of the actual lossy device.

Thus this is not merely comparing a damaged device with an ideal model.  It
asks whether time-reversal echo remains a useful *gradient estimator* once the
physical core ceases to be Hamiltonian.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import damping_gauge_reversal_probe as dg
from transfer_decomposition_probe import safe_corr


def flat(h,v): return np.concatenate([np.ravel(h),np.ravel(v)])

def rel(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))

def cosine(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))


def m_apply(m,wh,wv,x,eps):
    return dg.q_apply(m,wh,wv,x) - float(eps)*np.asarray(x)


def forward_actual(m,wh,wv,u,eps):
    a=1.0-float(eps)
    xm1=np.zeros(m.body.shape,np.complex128)
    x0=np.zeros_like(xm1)
    xs=[x0.copy()]
    for src in u:
        x1=m_apply(m,wh,wv,x0,eps)-a*xm1+src
        xs.append(x1.copy())
        xm1,x0=x0,x1
    return np.asarray(xs)


def objective_sources(m,x,coeff,r_desired):
    """Derivative source of original gauged output objective on actual x."""
    T=len(x)-1
    q=np.zeros_like(x)
    soma=tuple(map(int,m.soma))
    for k in range(1,T+1):
        q[k][soma]=float(coeff)*(float(r_desired)**(2*k))*x[k][soma]
    return q


def exact_gradient_actual(m,wh,wv,x,q,eps):
    """Reverse-mode gradient of actual lossy second-order recurrence."""
    *_,beta,_=dg.params(m)
    a=1.0-float(eps)
    T=len(x)-1
    p=np.zeros((T+2,)+m.body.shape,np.complex128)
    gh=np.zeros_like(wh,float); gv=np.zeros_like(wv,float)
    for k in range(T,0,-1):
        p[k]+=q[k]
        mu=p[k]
        f=x[k-1]
        dfh=f[:,1:]-f[:,:-1]; dmh=mu[:,:-1]-mu[:,1:]
        dfv=f[1:,:]-f[:-1,:]; dmv=mu[:-1,:]-mu[1:,:]
        gh += 2.0*beta*np.real(np.conj(dmh)*dfh)
        gv += 2.0*beta*np.real(np.conj(dmv)*dfv)
        p[k-1]+=m_apply(m,wh,wv,mu,eps)
        if k-2>=0: p[k-2]-=a*mu
    return gh,gv,p


def physical_retrace(m,wh,wv,x,u,eps):
    """Naive same-loss echo from terminal canonical reversal."""
    a=1.0-float(eps); T=len(x)-1
    out=[x[T].copy(),x[T-1].copy()]
    for j in range(T-1):
        nxt=m_apply(m,wh,wv,out[j+1],eps)-a*out[j]+u[T-1-j]
        out.append(nxt)
    return np.asarray(out)


def physical_adjoint(m,wh,wv,q,eps):
    """Causal same-loss returned error field, driven by reversed q."""
    a=1.0-float(eps); T=len(q)-1
    am1=np.zeros(m.body.shape,np.complex128); a0=np.zeros_like(am1)
    out=[a0.copy()]
    for j in range(T):
        a1=m_apply(m,wh,wv,a0,eps)-a*am1+q[T-j]
        out.append(a1.copy()); am1,a0=a0,a1
    return np.asarray(out)


def interference_gradient(m,w,a):
    *_,beta,_=dg.params(m)
    wp=w[1:]+a[1:]; wm=w[1:]-a[1:]
    ph,pv=dg.edge_diffs(wp); mh,mv=dg.edge_diffs(wm)
    ch=.25*np.sum(np.abs(ph)**2-np.abs(mh)**2,axis=0)
    cv=.25*np.sum(np.abs(pv)**2-np.abs(mv)**2,axis=0)
    return (-2.0*beta*ch).real,(-2.0*beta*cv).real


def mode_stability(m,wh,wv,eps):
    h,w=m.body.shape; N=h*w
    eye=np.eye(N,dtype=np.complex128).reshape((N,h,w))
    cols=[m_apply(m,wh,wv,eye[i],eps).ravel() for i in range(N)]
    M=np.stack(cols,axis=1).real
    vals=np.linalg.eigvalsh(M)
    a=1.0-float(eps)
    # For each scalar recurrence lambda^2-q lambda+a=0.
    max_root=0.0
    for q in vals:
        roots=np.roots([1.0,-float(q),a])
        max_root=max(max_root,float(np.max(np.abs(roots))))
    return float(vals.min()),float(vals.max()),float(max_root)


def prepare_sources(m,wh,wv,seq):
    z_ideal,u=dg.gauge_forward(m,wh,wv,seq)
    *_,r,_,_=dg.params(m)
    return u,r


def one_eps(m,wh,wv,seqT,seqD,coeffT,coeffD,eps):
    uT,r=prepare_sources(m,wh,wv,seqT); uD,_=prepare_sources(m,wh,wv,seqD)
    refs=[]; phys=[]; retr=[]; adjs=[]
    for u,coeff in ((uT,coeffT),(uD,coeffD)):
        x=forward_actual(m,wh,wv,u,eps)
        q=objective_sources(m,x,coeff,r)
        gh,gv,p=exact_gradient_actual(m,wh,wv,x,q,eps)
        w=physical_retrace(m,wh,wv,x,u,eps)
        aa=physical_adjoint(m,wh,wv,q,eps)
        gi=interference_gradient(m,w,aa)
        refs.append((gh,gv)); phys.append(gi)
        retr.append(rel(x[::-1],w))
        # Ideal backward alignment target from exact reverse-mode.
        target=np.zeros_like(aa); T=len(x)-1
        for j in range(1,T+1): target[j]=p[T-j+1]
        adjs.append(rel(target,aa))
    ref=flat(refs[0][0]+refs[1][0],refs[0][1]+refs[1][1])
    got=flat(phys[0][0]+phys[1][0],phys[0][1]+phys[1][1])
    qmin,qmax,rho=mode_stability(m,wh,wv,eps)
    return dict(eps=float(eps),corr=float(safe_corr(ref,got)),cosine=cosine(ref,got),
                relative_l2=rel(ref,got),norm_ratio=float(np.linalg.norm(got)/(np.linalg.norm(ref)+1e-30)),
                retrace_relative_l2=float(np.mean(retr)),adjoint_alignment_relative_l2=float(np.mean(adjs)),
                M_eigen_min=qmin,M_eigen_max=qmax,max_characteristic_root=rho)


def one_body(m,lag,steps,epses):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,lag,steps); seqD=ae.source_sequence(m,False,lag,steps)
    # Coefficients are frozen from the original intended damped task.  This keeps
    # the outer contrast objective definition unchanged while the actual physical
    # core is damaged.
    ET=ae.linear_forward(m,wh,wv,seqT,store=False); ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30; cT=2*ED/(S*S); cD=-2*ET/(S*S)
    return [dict(seed=int(m.cfg.seed),**one_eps(m,wh,wv,seqT,seqD,cT,cD,e)) for e in epses]


def summarize(rows):
    out={}
    for e in sorted(set(r['eps'] for r in rows)):
        q=[r for r in rows if r['eps']==e]
        out[str(e)]=dict(n=len(q),mean_corr=float(np.mean([x['corr'] for x in q])),min_corr=float(np.min([x['corr'] for x in q])),
                         mean_cosine=float(np.mean([x['cosine'] for x in q])),min_cosine=float(np.min([x['cosine'] for x in q])),
                         mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
                         mean_norm_ratio=float(np.mean([x['norm_ratio'] for x in q])),
                         mean_retrace_l2=float(np.mean([x['retrace_relative_l2'] for x in q])),
                         mean_adjoint_alignment_l2=float(np.mean([x['adjoint_alignment_relative_l2'] for x in q])),
                         max_characteristic_root=float(np.max([x['max_characteristic_root'] for x in q])),
                         max_M_eigen=float(np.max([x['M_eigen_max'] for x in q])))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472); ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--lag',type=int,default=20); ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--epses',default='0,0.0001,0.0005,0.001,0.002,0.005,0.01,0.02,0.05')
    ap.add_argument('--out',default='runs/damping_gauge_residual_loss_v02/dev_472_475.json')
    a=ap.parse_args(); epses=[float(x) for x in a.epses.split(',')]
    fa=Path(a.functional_arbors).resolve(); sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); b=m.bootstrap()
        if not b.get('ok'): continue
        m.mature=True
        rr=one_body(m,a.lag,a.steps,epses); rows.extend(rr)
        print('seed',seed,[(x['eps'],round(x['corr'],6),round(x['max_characteristic_root'],8)) for x in rr],flush=True)
    out=dict(config=vars(a),rows=rows,summary=summarize(rows))
    q=Path(a.out); q.parent.mkdir(parents=True,exist_ok=True); q.write_text(json.dumps(out,indent=2))
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__': main()
