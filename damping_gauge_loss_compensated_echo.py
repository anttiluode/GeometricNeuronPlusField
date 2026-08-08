"""Exact echo compensation for a uniformly lossy compiled core.

Development-only algebraic/hardware probe.

Consider the actual physical compiled recurrence with uniform residual loss

    x[n+1] = M x[n] - a x[n-1] + u[n],   0<a<=1,

where M is symmetric.  Ordinary same-body reversal cannot reproduce x in
reverse unless a=1.

However, if r[j]=x[T-j] and we define

    y[j] = a^j r[j],

then y obeys the *same lossy physical recurrence*:

    y[j+2] = M y[j+1] - a y[j]
             + a^(j+1) u[T-1-j].

So the same lossy body can generate an exponentially attenuated exact reverse
trajectory using only:

    y[0] = x[T]
    y[1] = a x[T-1]
    reversed input with known envelope a^(j+1).

For a symmetric second-order recurrence the causal returned adjoint b[j]
already obeys the same lossy operator and equals p[T-j+1].  The exact parameter
gradient needs

    sum_j Re(conj(Delta x[T-j]) Delta b[j]).

The physical echo supplies Delta y[j] = a^j Delta x[T-j].  Therefore a shared
integration weight a^(-j) restores the exact overlap:

    sum_j a^(-j) Re(conj(Delta y[j]) Delta b[j])
      = exact gradient overlap.

This probe tests the identity against the exact reverse-mode gradient of the
actual uniformly lossy compiled device.  No local time history is stored.

The price is explicit and should not be hidden:

* a calibrated uniform residual-loss factor a;
* scaled terminal reverse state y[1]=a*x[T-1];
* a global reverse-input envelope a^(j+1);
* a global detector/integrator gain a^(-j), whose dynamic range grows with
  sequence length and loss.

This is a scaling-conjugacy/compiler identity, not claimed novel.
"""
from __future__ import annotations

import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import damping_gauge_reversal_probe as dg
import damping_gauge_residual_loss_v02 as dl
from transfer_decomposition_probe import safe_corr


def flat(h,v): return np.concatenate([np.ravel(h),np.ravel(v)])

def rel(a,b):
    a=np.asarray(a);b=np.asarray(b)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))

def cosine(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))


def compensated_retrace(m,wh,wv,x,u,eps):
    a=1.0-float(eps); T=len(x)-1
    out=[x[T].copy(),a*x[T-1].copy()]
    for j in range(T-1):
        src=(a**(j+1))*u[T-1-j]
        nxt=dl.m_apply(m,wh,wv,out[j+1],eps)-a*out[j]+src
        out.append(nxt)
    return np.asarray(out)


def weighted_interference_gradient(m,y,b,eps):
    """Integrated +/- energy difference with global a^-j gain."""
    *_,beta,_=dg.params(m)
    a=1.0-float(eps); T=len(y)-1
    yp=y[1:]+b[1:]; ym=y[1:]-b[1:]
    ph,pv=dg.edge_diffs(yp); mh,mv=dg.edge_diffs(ym)
    # one scalar global detector weight per reverse-time sample j=1..T
    wt=(a**(-np.arange(1,T+1,dtype=float)))
    whgt=wt.reshape((T,1,1))
    ch=.25*np.sum(whgt*(np.abs(ph)**2-np.abs(mh)**2),axis=0)
    cv=.25*np.sum(whgt*(np.abs(pv)**2-np.abs(mv)**2),axis=0)
    return (-2.0*beta*ch).real,(-2.0*beta*cv).real


def one_order(m,wh,wv,u,coeff,r_desired,eps):
    x=dl.forward_actual(m,wh,wv,u,eps)
    q=dl.objective_sources(m,x,coeff,r_desired)
    gh,gv,p=dl.exact_gradient_actual(m,wh,wv,x,q,eps)
    y=compensated_retrace(m,wh,wv,x,u,eps)
    b=dl.physical_adjoint(m,wh,wv,q,eps)
    gi=weighted_interference_gradient(m,y,b,eps)

    # Exact trajectory identities.
    a=1.0-float(eps); T=len(x)-1
    target=np.asarray([(a**j)*x[T-j] for j in range(T+1)])
    targetb=np.zeros_like(b)
    for j in range(1,T+1): targetb[j]=p[T-j+1]
    return dict(exact=(gh,gv),physical=gi,
                attenuated_retrace_l2=rel(target,y),
                adjoint_alignment_l2=rel(targetb,b),
                detector_gain_end=float(a**(-T)),
                reverse_source_scale_end=float(a**T))


def one_eps(m,wh,wv,seqT,seqD,cT,cD,eps):
    uT,r=dl.prepare_sources(m,wh,wv,seqT);uD,_=dl.prepare_sources(m,wh,wv,seqD)
    T=one_order(m,wh,wv,uT,cT,r,eps);D=one_order(m,wh,wv,uD,cD,r,eps)
    ref=flat(T['exact'][0]+D['exact'][0],T['exact'][1]+D['exact'][1])
    got=flat(T['physical'][0]+D['physical'][0],T['physical'][1]+D['physical'][1])
    return dict(eps=float(eps),corr=float(safe_corr(ref,got)),cosine=cosine(ref,got),
                relative_l2=rel(ref,got),norm_ratio=float(np.linalg.norm(got)/(np.linalg.norm(ref)+1e-30)),
                retrace_l2=float(np.mean([T['attenuated_retrace_l2'],D['attenuated_retrace_l2']])),
                adjoint_l2=float(np.mean([T['adjoint_alignment_l2'],D['adjoint_alignment_l2']])),
                detector_gain_end=T['detector_gain_end'],reverse_source_scale_end=T['reverse_source_scale_end'])


def one_body(m,lag,steps,epses):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,seqT,store=False);ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30;cT=2*ED/(S*S);cD=-2*ET/(S*S)
    return [dict(seed=int(m.cfg.seed),**one_eps(m,wh,wv,seqT,seqD,cT,cD,e)) for e in epses]


def summarize(rows):
    out={}
    for e in sorted(set(r['eps'] for r in rows)):
        q=[r for r in rows if r['eps']==e]
        out[str(e)]=dict(n=len(q),mean_corr=float(np.mean([x['corr'] for x in q])),min_corr=float(np.min([x['corr'] for x in q])),
                         mean_cosine=float(np.mean([x['cosine'] for x in q])),mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
                         max_retrace_l2=float(np.max([x['retrace_l2'] for x in q])),max_adjoint_l2=float(np.max([x['adjoint_l2'] for x in q])),
                         detector_gain_end=float(q[0]['detector_gain_end']),reverse_source_scale_end=float(q[0]['reverse_source_scale_end']))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--epses',default='0,0.0001,0.0005,0.001,0.002,0.005,0.01,0.02,0.05')
    ap.add_argument('--out',default='runs/damping_gauge_loss_compensated_echo/dev_472_475.json')
    a=ap.parse_args();epses=[float(x) for x in a.epses.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True
        rr=one_body(m,a.lag,a.steps,epses);rows.extend(rr)
        print('seed',seed,[(x['eps'],f"{x['corr']:.12f}",f"{x['relative_l2']:.2e}") for x in rr],flush=True)
    out=dict(config=vars(a),rows=rows,summary=summarize(rows))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
