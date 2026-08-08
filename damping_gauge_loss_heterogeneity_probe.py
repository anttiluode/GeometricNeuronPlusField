"""Damage the scalar-loss assumption of the compensated echo compiler.

Development only on reused seeds.

Two hostile tests:

A. Uniform loss, wrong scalar calibration.
   The actual core has eps_true, but terminal scaling, reversed-source envelope
   and detector envelope use eps_hat.  This measures calibration tolerance.

B. Spatially heterogeneous loss.
   Each cell receives eps_i around a chosen mean.  The actual recurrence is

       x[n+1] = (Q-E) x[n] - (I-E) x[n-1] + u[n]

   with diagonal E.  The exact digital adjoint of that *actual* device is the
   reference.  The physical echo is allowed only one scalar compensation factor
   equal to the known mean loss.  This asks how much non-proportional loss can
   be tolerated before the memory-free scalar compiler loses gradient direction.

No fresh confirmation bodies are touched.
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
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/(np.linalg.norm(a)+1e-30))
def cosine(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))


def map_apply(m,wh,wv,x,epsmap):
    return dg.q_apply(m,wh,wv,x)-epsmap*np.asarray(x)


def forward_map(m,wh,wv,u,epsmap):
    back=1.0-epsmap
    xm1=np.zeros(m.body.shape,np.complex128);x0=np.zeros_like(xm1);xs=[x0.copy()]
    for src in u:
        x1=map_apply(m,wh,wv,x0,epsmap)-back*xm1+src
        xs.append(x1.copy());xm1,x0=x0,x1
    return np.asarray(xs)


def objective_sources(m,x,coeff,r_desired):
    T=len(x)-1;q=np.zeros_like(x);s=tuple(map(int,m.soma))
    for k in range(1,T+1):q[k][s]=float(coeff)*(float(r_desired)**(2*k))*x[k][s]
    return q


def exact_gradient_map(m,wh,wv,x,q,epsmap):
    *_,beta,_=dg.params(m);back=1.0-epsmap;T=len(x)-1
    p=np.zeros((T+2,)+m.body.shape,np.complex128)
    gh=np.zeros_like(wh,float);gv=np.zeros_like(wv,float)
    for k in range(T,0,-1):
        p[k]+=q[k];mu=p[k];f=x[k-1]
        dfh=f[:,1:]-f[:,:-1];dmh=mu[:,:-1]-mu[:,1:]
        dfv=f[1:,:]-f[:-1,:];dmv=mu[:-1,:]-mu[1:,:]
        gh+=2*beta*np.real(np.conj(dmh)*dfh);gv+=2*beta*np.real(np.conj(dmv)*dfv)
        p[k-1]+=map_apply(m,wh,wv,mu,epsmap)
        if k-2>=0:p[k-2]-=back*mu
    return gh,gv,p


def physical_adjoint_map(m,wh,wv,q,epsmap):
    back=1.0-epsmap;T=len(q)-1
    am1=np.zeros(m.body.shape,np.complex128);a0=np.zeros_like(am1);out=[a0.copy()]
    for j in range(T):
        a1=map_apply(m,wh,wv,a0,epsmap)-back*am1+q[T-j]
        out.append(a1.copy());am1,a0=a0,a1
    return np.asarray(out)


def compensated_echo_map(m,wh,wv,x,u,epsmap,eps_hat):
    ah=1.0-float(eps_hat);back=1.0-epsmap;T=len(x)-1
    out=[x[T].copy(),ah*x[T-1].copy()]
    for j in range(T-1):
        src=(ah**(j+1))*u[T-1-j]
        nxt=map_apply(m,wh,wv,out[j+1],epsmap)-back*out[j]+src
        out.append(nxt)
    return np.asarray(out)


def weighted_interference(m,y,b,eps_hat):
    *_,beta,_=dg.params(m);ah=1.0-float(eps_hat);T=len(y)-1
    yp=y[1:]+b[1:];ym=y[1:]-b[1:]
    ph,pv=dg.edge_diffs(yp);mh,mv=dg.edge_diffs(ym)
    wt=(ah**(-np.arange(1,T+1,dtype=float))).reshape((T,1,1))
    ch=.25*np.sum(wt*(np.abs(ph)**2-np.abs(mh)**2),axis=0)
    cv=.25*np.sum(wt*(np.abs(pv)**2-np.abs(mv)**2),axis=0)
    return (-2*beta*ch).real,(-2*beta*cv).real


def one_order(m,wh,wv,u,coeff,r_desired,epsmap,eps_hat):
    x=forward_map(m,wh,wv,u,epsmap);q=objective_sources(m,x,coeff,r_desired)
    gh,gv,p=exact_gradient_map(m,wh,wv,x,q,epsmap)
    y=compensated_echo_map(m,wh,wv,x,u,epsmap,eps_hat)
    b=physical_adjoint_map(m,wh,wv,q,epsmap)
    gi=weighted_interference(m,y,b,eps_hat)
    return (gh,gv),gi


def score_case(m,wh,wv,seqT,seqD,cT,cD,epsmap,eps_hat):
    uT,r=dl.prepare_sources(m,wh,wv,seqT);uD,_=dl.prepare_sources(m,wh,wv,seqD)
    rt,pt=one_order(m,wh,wv,uT,cT,r,epsmap,eps_hat)
    rd,pd=one_order(m,wh,wv,uD,cD,r,epsmap,eps_hat)
    ref=flat(rt[0]+rd[0],rt[1]+rd[1]);got=flat(pt[0]+pd[0],pt[1]+pd[1])
    return dict(corr=float(safe_corr(ref,got)),cosine=cosine(ref,got),relative_l2=rel(ref,got),
                norm_ratio=float(np.linalg.norm(got)/(np.linalg.norm(ref)+1e-30)))


def make_epsmap(m,mean_eps,cv,rng):
    body=m.body.astype(bool);e=np.full(m.body.shape,float(mean_eps),float)
    if cv>0:
        # multiplicative lognormal field with approximately requested coefficient
        # of variation for small/moderate cv, renormalized to exact body mean.
        z=rng.standard_normal(m.body.shape)
        mult=np.exp(float(cv)*z-.5*float(cv)**2)
        mult/=float(np.mean(mult[body]))
        e[body]=float(mean_eps)*mult[body]
    e=np.clip(e,0.0,.5)
    return e


def one_body(m,lag,steps,cal_mean,cal_errors,het_means,het_cvs,reps):
    wh,wv=ae.bond_weights(m,m.body);seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,seqT,store=False);ED=ae.linear_forward(m,wh,wv,seqD,store=False);S=ET+ED+1e-30
    cT=2*ED/(S*S);cD=-2*ET/(S*S);rows=[];seed=int(m.cfg.seed)
    # A: exact uniform physical loss, wrong scalar estimate.
    for de in cal_errors:
        eh=float(cal_mean)*(1.0+float(de));emap=np.full(m.body.shape,float(cal_mean),float)
        s=score_case(m,wh,wv,seqT,seqD,cT,cD,emap,eh)
        rows.append(dict(seed=seed,kind='calibration_relative_error',mean_eps=float(cal_mean),value=float(de),rep=0,eps_hat=eh,**s))
    # B: nonuniform physical loss, compensation gets exact mean only.
    for me in het_means:
        for cv in het_cvs:
            for rep in range(int(reps)):
                rng=np.random.default_rng(seed*100003+rep*7919+int(me*1e7)+int(cv*1e6))
                emap=make_epsmap(m,me,cv,rng)
                s=score_case(m,wh,wv,seqT,seqD,cT,cD,emap,float(me))
                body=m.body.astype(bool)
                rows.append(dict(seed=seed,kind='spatial_loss_cv',mean_eps=float(me),value=float(cv),rep=rep,
                                 actual_mean=float(np.mean(emap[body])),actual_std=float(np.std(emap[body])),**s))
    return rows


def summarize(rows):
    out={}
    for kind in sorted(set(r['kind'] for r in rows)):
        out[kind]={}
        means=sorted(set(r['mean_eps'] for r in rows if r['kind']==kind))
        for me in means:
            out[kind][str(me)]={}
            vals=sorted(set(r['value'] for r in rows if r['kind']==kind and r['mean_eps']==me))
            for v in vals:
                q=[r for r in rows if r['kind']==kind and r['mean_eps']==me and r['value']==v]
                out[kind][str(me)][str(v)]=dict(n=len(q),mean_corr=float(np.mean([x['corr'] for x in q])),min_corr=float(np.min([x['corr'] for x in q])),
                                                mean_cosine=float(np.mean([x['cosine'] for x in q])),mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
                                                mean_norm_ratio=float(np.mean([x['norm_ratio'] for x in q])))
    return out


def plist(s):return [float(x) for x in s.split(',') if x.strip()]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=472);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--cal-mean',type=float,default=.005);ap.add_argument('--cal-errors',default='-0.1,-0.05,-0.02,-0.01,-0.005,0,0.005,0.01,0.02,0.05,0.1')
    ap.add_argument('--het-means',default='0.001,0.005,0.01');ap.add_argument('--het-cvs',default='0,0.01,0.05,0.1,0.2,0.5');ap.add_argument('--reps',type=int,default=3)
    ap.add_argument('--out',default='runs/damping_gauge_loss_heterogeneity/dev_472_475.json');a=ap.parse_args()
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));b=m.bootstrap()
        if not b.get('ok'):continue
        m.mature=True;rows.extend(one_body(m,a.lag,a.steps,a.cal_mean,plist(a.cal_errors),plist(a.het_means),plist(a.het_cvs),a.reps));print('seed',seed,'done',flush=True)
    out=dict(config=vars(a),rows=rows,summary=summarize(rows));q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2));print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
