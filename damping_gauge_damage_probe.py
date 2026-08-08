"""Damage the exact damping-gauge echo compiler on reused development bodies.

This probe does NOT touch held-out bodies.  It starts from the exact transformed
forward/adjoint construction in `damping_gauge_reversal_probe.py` and damages
only the proposed physical echo implementation.

Damage families:

1. momentum reversal gain alpha
      ideal reverse initial state has w1=z[T-1].  An imperfect momentum flip is
      represented as
          w1 = z[T] - alpha * (z[T]-z[T-1])
      with alpha=1 ideal.

2. residual reverse loss epsilon
      ideal conservative recurrence
          w[j+2] = Q w[j+1] - w[j] + u_rev
      becomes
          w[j+2] = Q w[j+1] - (1-epsilon) w[j] + u_rev.
      The same damaged recurrence is used for the returned adjoint.  This is a
      simple determinant-loss model, not a component-level circuit model.

3. reverse-operator drift sigma
      branch couplings used in the reverse trials are multiplied by independent
      log-normal perturbations with RMS scale sigma.  Forward inference remains
      at the calibrated operator.

4. local energy-readout noise eta
      after forming the exact +/- integrated branch energies, additive Gaussian
      noise is applied with RMS eta times the RMS magnitude of the true gradient
      map.  This isolates detector/integrator noise from echo error.

For every damage point we report correlation/cosine/relative-L2 of the complete
target+distractor bond-gradient map against the exact original gradient.
"""
from __future__ import annotations

import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import damping_gauge_reversal_probe as dg
from transfer_decomposition_probe import safe_corr


def flat(h,v): return np.concatenate([np.ravel(h),np.ravel(v)])

def cosine(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))

def rel(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))


def reverse_q_apply(m,wh,wv,x,loss=0.0):
    # q_apply itself is the calibrated conservative first-state coefficient.
    # Loss enters only through the coefficient on the state two steps back.
    return dg.q_apply(m,wh,wv,x)


def damaged_reverse_forward(m,wh,wv,z,u,alpha=1.0,loss=0.0):
    T=len(z)-1
    w0=z[T].copy()
    w1=z[T]-float(alpha)*(z[T]-z[T-1])
    out=[w0,w1]
    back=float(1.0-loss)
    for j in range(T-1):
        nxt=dg.q_apply(m,wh,wv,out[j+1])-back*out[j]+u[T-1-j]
        out.append(nxt)
    return np.asarray(out)


def damaged_causal_adjoint(m,wh,wv,q,loss=0.0):
    T=len(q)-1
    am1=np.zeros(m.body.shape,np.complex128);a0=np.zeros_like(am1)
    out=[a0.copy()];back=float(1.0-loss)
    for j in range(T):
        a1=dg.q_apply(m,wh,wv,a0)-back*am1+q[T-j]
        out.append(a1.copy());am1,a0=a0,a1
    return np.asarray(out)


def perturb_weights(wh,wv,sigma,rng):
    if sigma<=0:return wh,wv
    # multiplicative approximately RMS-sigma error; zero weights remain zero.
    ph=np.exp(float(sigma)*rng.standard_normal(wh.shape)-0.5*float(sigma)**2)
    pv=np.exp(float(sigma)*rng.standard_normal(wv.shape)-0.5*float(sigma)**2)
    return wh*ph,wv*pv


def prepare_order(m,wh,wv,seq,coeff):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    z,u=dg.gauge_forward(m,wh,wv,seq)
    gh,gv,padj,q=dg.gauge_adjoint_gradient(m,wh,wv,z,coeff)
    return dict(z=z,u=u,q=q,exact=(gh,gv))


def physical_order(m,base_wh,base_wv,dat,alpha,loss,sigma,rng):
    wh,wv=perturb_weights(base_wh,base_wv,sigma,rng)
    w=damaged_reverse_forward(m,wh,wv,dat['z'],dat['u'],alpha,loss)
    a=damaged_causal_adjoint(m,wh,wv,dat['q'],loss)
    gh,gv=dg.interference_gradient(m,w,a)
    return gh,gv


def one_body(m,lag,steps,alphas,losses,drifts,noises):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,seqT,store=False);ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30;aT=2*ED/(S*S);aD=-2*ET/(S*S)
    T=prepare_order(m,wh,wv,seqT,aT);D=prepare_order(m,wh,wv,seqD,aD)
    ref=flat(T['exact'][0]+D['exact'][0],T['exact'][1]+D['exact'][1])
    rr=float(np.sqrt(np.mean(ref*ref)))
    rows=[]
    seed=int(m.cfg.seed)

    def score(kind,value,alpha=1.0,loss=0.0,sigma=0.0,noise=0.0,rep=0):
        rng=np.random.default_rng(seed*1000003+rep*9176+int(round(value*1e7))+hash(kind)%100003)
        gt=physical_order(m,wh,wv,T,alpha,loss,sigma,rng)
        gd=physical_order(m,wh,wv,D,alpha,loss,sigma,rng)
        got=flat(gt[0]+gd[0],gt[1]+gd[1])
        if noise>0:
            got=got+rng.standard_normal(got.shape)*(float(noise)*rr)
        return dict(seed=seed,kind=kind,value=float(value),rep=int(rep),
                    corr=float(safe_corr(ref,got)),cosine=cosine(ref,got),relative_l2=rel(ref,got),
                    norm_ratio=float(np.linalg.norm(got)/(np.linalg.norm(ref)+1e-30)))

    for a in alphas: rows.append(score('momentum_alpha',a,alpha=a))
    for x in losses: rows.append(score('reverse_loss',x,loss=x))
    # average several independent drift/noise realizations later at summary level
    for x in drifts:
        for rep in range(4): rows.append(score('operator_drift',x,sigma=x,rep=rep))
    for x in noises:
        for rep in range(4): rows.append(score('gradient_readout_noise',x,noise=x,rep=rep))
    return rows


def summarize(rows):
    out={}
    kinds=sorted(set(r['kind'] for r in rows))
    for k in kinds:
        vals=sorted(set(r['value'] for r in rows if r['kind']==k))
        out[k]={}
        for v in vals:
            q=[r for r in rows if r['kind']==k and r['value']==v]
            out[k][str(v)]=dict(
                n=len(q),mean_corr=float(np.mean([x['corr'] for x in q])),min_corr=float(np.min([x['corr'] for x in q])),
                mean_cosine=float(np.mean([x['cosine'] for x in q])),min_cosine=float(np.min([x['cosine'] for x in q])),
                mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
                mean_norm_ratio=float(np.mean([x['norm_ratio'] for x in q])))
    return out


def parse_list(s):return [float(x) for x in s.split(',') if x.strip()]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472);ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--alphas',default='1,0.999,0.995,0.99,0.98,0.95,0.9')
    ap.add_argument('--losses',default='0,0.0001,0.0005,0.001,0.002,0.005,0.01')
    ap.add_argument('--drifts',default='0,0.001,0.005,0.01,0.02,0.05')
    ap.add_argument('--noises',default='0,0.001,0.005,0.01,0.02,0.05,0.1')
    ap.add_argument('--out',default='runs/damping_gauge_damage/dev_472_475.json')
    a=ap.parse_args();fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True
        z=one_body(m,a.lag,a.steps,parse_list(a.alphas),parse_list(a.losses),parse_list(a.drifts),parse_list(a.noises));rows.extend(z)
        print('seed',seed,'done',flush=True)
    out=dict(config=vars(a),rows=rows,summary=summarize(rows))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
