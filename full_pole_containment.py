"""B0 containment sanity check from BENCHMARK_SCOPE_V01.md.

Diagonalize the exact mature weighted grid operator used by the linear adjoint model
and show that an unconstrained modal/pole bank reconstructs the identical soma trace.
This is expected mathematics, not a performance claim.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # patches ae.weighted_lap to exact exterior-leak version


def stiffness_matrix(wh,wv):
    h=wh.shape[0]; w=wv.shape[1]; N=h*w
    K=np.zeros((N,N),float)
    def add(i,j,k):
        K[i,i]+=k;K[j,j]+=k;K[i,j]-=k;K[j,i]-=k
    for y in range(h):
        for x in range(w-1):
            add(y*w+x,y*w+x+1,float(wh[y,x]))
    for y in range(h-1):
        for x in range(w):
            add(y*w+x,(y+1)*w+x,float(wv[y,x]))
    kb=float(min(np.min(wh),np.min(wv)))
    # exact_weighted_lap leaks to zero outside the finite square
    for y in range(h):
        K[y*w+0,y*w+0]+=kb;K[y*w+w-1,y*w+w-1]+=kb
    for x in range(w):
        K[0*w+x,0*w+x]+=kb;K[(h-1)*w+x,(h-1)*w+x]+=kb
    return K


def modal_trace(m,K,seq):
    lam,Phi=np.linalg.eigh(K)
    N=len(lam);dt=float(m.cfg.dt);st=float(m.cfg.stiffness)
    damp=float(m.cfg.damping);rest=float(m.cfg.restoring)
    q=np.zeros(N,complex);v=np.zeros(N,complex)
    soma=int(m.soma[0])*m.body.shape[1]+int(m.soma[1])
    c=Phi[soma,:]
    out=[]
    for src in seq:
        f=np.asarray(src,complex).ravel()
        nz=np.flatnonzero(np.abs(f)>0)
        sm=np.zeros(N,complex)
        for i in nz:
            sm += Phi[int(i),:]*f[int(i)]
        v += dt*(-st*lam*q-damp*v-rest*q+sm)
        q += dt*v
        out.append(np.dot(c,q))
    return np.asarray(out),lam


def direct_trace(m,wh,wv,seq):
    ps,vs,E=ae.linear_forward(m,wh,wv,seq,store=True)
    soma=tuple(map(int,m.soma))
    return np.asarray(ps[1:,soma[0],soma[1]],complex),float(E)


def one(m,lag,steps):
    wh,wv=ae.bond_weights(m,m.body)
    K=stiffness_matrix(wh,wv)
    # operator construction audit
    rng=np.random.default_rng(int(m.cfg.seed)+9090)
    u=rng.normal(size=m.body.shape)+1j*rng.normal(size=m.body.shape)
    lhs=(-K@u.ravel()).reshape(u.shape)
    rhs=ae.weighted_lap(u,wh,wv)
    op_err=float(np.linalg.norm(lhs-rhs)/(np.linalg.norm(rhs)+1e-30))
    rows={}
    for name,target in [('T',True),('D',False)]:
        seq=ae.source_sequence(m,target,lag,steps)
        d,E=direct_trace(m,wh,wv,seq)
        p,lam=modal_trace(m,K,seq)
        rel=float(np.linalg.norm(p-d)/(np.linalg.norm(d)+1e-30))
        rows[name]=dict(trace_relative_l2=rel,max_abs_error=float(np.max(np.abs(p-d))),energy_direct=E,
                        energy_modal=float(np.sum(np.abs(p)**2)))
    ET=rows['T']['energy_direct'];ED=rows['D']['energy_direct']
    PT=rows['T']['energy_modal'];PD=rows['D']['energy_modal']
    Cd=float((ET-ED)/(ET+ED+1e-30));Cp=float((PT-PD)/(PT+PD+1e-30))
    return dict(seed=int(m.cfg.seed),N=int(K.shape[0]),operator_relative_error=op_err,
                target=rows['T'],distractor=rows['D'],contrast_direct=Cd,contrast_modal=Cp,
                contrast_abs_error=float(abs(Cd-Cp)),lambda_min=float(lam[0]),lambda_max=float(lam[-1]))


def selftest():
    wh=np.ones((2,1));wv=np.ones((1,2));K=stiffness_matrix(wh,wv)
    assert np.max(np.abs(K-K.T))<1e-15
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=264);ap.add_argument('--seeds',type=int,default=1)
    ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--out',default='runs/pole_containment/pole_containment.json');ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps);rows.append(r)
        print(f"seed {seed}: op={r['operator_relative_error']:.2e} T={r['target']['trace_relative_l2']:.2e} D={r['distractor']['trace_relative_l2']:.2e} C={r['contrast_abs_error']:.2e}",flush=True)
    summary=dict(bodies=len(rows),max_operator_error=float(max(r['operator_relative_error'] for r in rows)),
                 max_trace_error=float(max(max(r['target']['trace_relative_l2'],r['distractor']['trace_relative_l2']) for r in rows)),
                 max_contrast_error=float(max(r['contrast_abs_error'] for r in rows)))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(summary=summary,rows=rows),indent=2))
    print(summary)
if __name__=='__main__':main()
