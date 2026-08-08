"""Adjoint prediction of one-cell structural-event consequences.

See ADJOINT_ELIGIBILITY_DISCOVERY_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from transfer_decomposition_probe import safe_corr
from structural_interference_probe import body_probe as nonlinear_event_probe, n4


def bond_weights(m, body):
    b = np.asarray(body, bool)
    ka = float(m.cfg.k_arbor)
    kb = float(m.cfg.k_mature_bath)
    wh = np.where(b[:, :-1] & b[:, 1:], ka, kb).astype(np.float64)
    wv = np.where(b[:-1, :] & b[1:, :], ka, kb).astype(np.float64)
    return wh, wv


def weighted_lap(u, wh, wv):
    u = np.asarray(u)
    out = np.zeros_like(u, dtype=np.result_type(u, np.complex128))
    dh = u[:, 1:] - u[:, :-1]
    out[:, :-1] += wh * dh
    out[:, 1:]  -= wh * dh
    dv = u[1:, :] - u[:-1, :]
    out[:-1, :] += wv * dv
    out[1:, :]  -= wv * dv
    return out


def addsrc(a, b, shape):
    out = np.zeros(shape, np.complex128)
    if not isinstance(a, (float, int, np.floating)):
        out += np.asarray(a, np.complex128)
    if not isinstance(b, (float, int, np.floating)):
        out += np.asarray(b, np.complex128)
    return out


def source_sequence(m, target, lag, steps):
    first, second = ((0,1) if target else (1,0))
    seq=[]
    for t in range(int(steps)):
        a=m.pulse_source(first,t,False)
        b=m.pulse_source(second,t-int(lag),False)
        seq.append(addsrc(a,b,m.body.shape))
    return seq


def linear_forward(m, wh, wv, seq, store=True):
    c=m.cfg
    dt=float(c.dt); stiff=float(c.stiffness); damp=float(c.damping); rest=float(c.restoring)
    psi=np.zeros(m.body.shape,np.complex128)
    vel=np.zeros_like(psi)
    if store:
        ps=[psi.copy()]; vs=[vel.copy()]
    energy=0.0
    soma=tuple(map(int,m.soma))
    for src in seq:
        lap=weighted_lap(psi,wh,wv)
        vel=vel+dt*(stiff*lap-damp*vel-rest*psi+src)
        psi=psi+dt*vel
        energy += float(abs(psi[soma])**2)
        if store:
            ps.append(psi.copy()); vs.append(vel.copy())
    if store:
        return np.asarray(ps),np.asarray(vs),float(energy)
    return float(energy)


def contrast(ET, ED):
    return float((ET-ED)/(ET+ED+1e-30))


def linear_contrast(m, wh, wv, seqT, seqD):
    ET=linear_forward(m,wh,wv,seqT,store=False)
    ED=linear_forward(m,wh,wv,seqD,store=False)
    return contrast(ET,ED),ET,ED


def adjoint_grad(m, wh, wv, ps, vs, coeff):
    """Gradient of coeff * sum_t |psi_s(t)|^2 wrt every bond conductance.

    Complex adjoints use dJ = 2 Re <lambda, dx>.
    """
    c=m.cfg
    dt=float(c.dt); stiff=float(c.stiffness); damp=float(c.damping); rest=float(c.restoring)
    av=1.0-dt*damp
    lam_psi=np.zeros(m.body.shape,np.complex128)
    lam_vel=np.zeros_like(lam_psi)
    gh=np.zeros_like(wh,float); gv=np.zeros_like(wv,float)
    soma=tuple(map(int,m.soma))
    T=len(ps)-1
    for k in range(T,0,-1):
        # Output at state k.
        lam_psi[soma] += float(coeff)*ps[k][soma]

        # psi_k = psi_{k-1} + dt v_k, so both adjoints act through v_k.
        mu = lam_vel + dt*lam_psi
        prev = ps[k-1]

        # Local bond derivative.  For edge i-j:
        # d(L psi)_i/dw = psi_j-psi_i, and opposite at j.
        dpsi_h = prev[:,1:] - prev[:,:-1]
        dmu_h = mu[:,:-1] - mu[:,1:]
        gh += 2.0*dt*stiff*np.real(np.conj(dmu_h)*dpsi_h)

        dpsi_v = prev[1:,:] - prev[:-1,:]
        dmu_v = mu[:-1,:] - mu[1:,:]
        gv += 2.0*dt*stiff*np.real(np.conj(dmu_v)*dpsi_v)

        # Reverse the state transition.
        lam_vel_prev = av*mu
        lam_psi_prev = lam_psi + dt*(stiff*weighted_lap(mu,wh,wv)-rest*mu)
        lam_psi,lam_vel = lam_psi_prev,lam_vel_prev
    return gh,gv


def contrast_adjoint(m, body, seqT, seqD):
    wh,wv=bond_weights(m,body)
    pT,vT,ET=linear_forward(m,wh,wv,seqT,store=True)
    pD,vD,ED=linear_forward(m,wh,wv,seqD,store=True)
    S=ET+ED+1e-30
    aT=2.0*ED/(S*S)
    aD=-2.0*ET/(S*S)
    ghT,gvT=adjoint_grad(m,wh,wv,pT,vT,aT)
    ghD,gvD=adjoint_grad(m,wh,wv,pD,vD,aD)
    gh=ghT+ghD; gv=gvT+gvD
    C=contrast(ET,ED)

    # Forward-only local task-energy control, evaluated on output states.
    dhT=pT[1:,:,1:]-pT[1:,:,:-1]
    dhD=pD[1:,:,1:]-pD[1:,:,:-1]
    Fh=np.sum(np.abs(dhT)**2-np.abs(dhD)**2,axis=0).real
    dvT=pT[1:,1:,:]-pT[1:,:-1,:]
    dvD=pD[1:,1:,:]-pD[1:,:-1,:]
    Fv=np.sum(np.abs(dvT)**2-np.abs(dvD)**2,axis=0).real
    return dict(C=C,ET=ET,ED=ED,wh=wh,wv=wv,gh=gh,gv=gv,Fh=Fh,Fv=Fv,
                pT=pT,pD=pD)


def edge_lookup(h,v,p,q):
    y0,x0=map(int,p); y1,x1=map(int,q)
    if y0==y1 and abs(x0-x1)==1:
        return float(h[y0,min(x0,x1)])
    if x0==x1 and abs(y0-y1)==1:
        return float(v[min(y0,y1),x0])
    raise ValueError(f'not a 4-neighbor edge: {p} {q}')


def set_edge(h,v,p,q,value):
    y0,x0=map(int,p); y1,x1=map(int,q)
    if y0==y1 and abs(x0-x1)==1:
        h[y0,min(x0,x1)]=value; return
    if x0==x1 and abs(y0-y1)==1:
        v[min(y0,y1),x0]=value; return
    raise ValueError('not edge')


def event_edges(body,event):
    b=np.asarray(body,bool)
    p=tuple(map(int,event['cell']))
    if event['kind']=='add':
        qs=[q for q in n4(*p,b.shape) if b[q]]
        if len(qs)!=1:
            raise ValueError(f'addition does not have exactly one occupied neighbor: {p} -> {qs}')
        return [(p,qs[0],+1.0)]
    qs=[q for q in n4(*p,b.shape) if b[q]]
    return [(p,q,-1.0) for q in qs]


def event_score(arrh,arrv,body,event,delta_k):
    s=0.0
    for p,q,sgn in event_edges(body,event):
        s += sgn*float(delta_k)*edge_lookup(arrh,arrv,p,q)
    return float(s)


def mutate_body(body,event):
    b=np.asarray(body,bool).copy()
    p=tuple(map(int,event['cell']))
    b[p]=(event['kind']=='add')
    return b


def finite_diff_check(m, base, seqT, seqD, eps=1e-5):
    # Use the bond with largest absolute analytic gradient to avoid dividing by a
    # numerically null derivative.
    gh,gv=base['gh'],base['gv']
    if np.max(np.abs(gh)) >= np.max(np.abs(gv)):
        k=np.unravel_index(np.argmax(np.abs(gh)),gh.shape)
        p=(int(k[0]),int(k[1])); q=(int(k[0]),int(k[1]+1)); g=float(gh[k])
    else:
        k=np.unravel_index(np.argmax(np.abs(gv)),gv.shape)
        p=(int(k[0]),int(k[1])); q=(int(k[0]+1),int(k[1])); g=float(gv[k])
    hp=base['wh'].copy(); vp=base['wv'].copy(); hm=base['wh'].copy(); vm=base['wv'].copy()
    w0=edge_lookup(base['wh'],base['wv'],p,q)
    set_edge(hp,vp,p,q,w0+eps); set_edge(hm,vm,p,q,w0-eps)
    cp,_,_=linear_contrast(m,hp,vp,seqT,seqD)
    cm,_,_=linear_contrast(m,hm,vm,seqT,seqD)
    fd=float((cp-cm)/(2*eps))
    rel=float(abs(fd-g)/(abs(fd)+abs(g)+1e-12))
    return dict(edge=[list(p),list(q)],gradient=g,finite_difference=fd,relative_error=rel)


def sign_test_two_sided(vals):
    z=np.asarray(vals,float); z=z[np.isfinite(z)&(z!=0)]
    if len(z)==0:return float('nan'),0,0
    w=int(np.sum(z>0));l=int(np.sum(z<0));n=w+l;k=min(w,l)
    tail=sum(math.comb(n,i) for i in range(k+1))/(2**n)
    return float(min(1.0,2*tail)),w,l


def body_adjoint_probe(m, lag, steps, max_each):
    # Get the registered nonlinear counterfactual event set and observed changes.
    observed=nonlinear_event_probe(m,lag,steps,max_each)
    seqT=source_sequence(m,True,lag,steps)
    seqD=source_sequence(m,False,lag,steps)
    base=contrast_adjoint(m,m.body,seqT,seqD)
    fdcheck=finite_diff_check(m,base,seqT,seqD)

    # Laplacian construction control against FunctionalArbor itself.
    rng=np.random.default_rng(int(m.cfg.seed)+991)
    u=rng.normal(size=m.body.shape)+1j*rng.normal(size=m.body.shape)
    lap0=m._lap(u,True)
    lap1=weighted_lap(u,base['wh'],base['wv'])
    lap_rel=float(np.linalg.norm(lap1-lap0)/(np.linalg.norm(lap0)+1e-30))

    dk=float(m.cfg.k_arbor-m.cfg.k_mature_bath)
    out_events=[]
    for e in observed['events']:
        pred=event_score(base['gh'],base['gv'],m.body,e,dk)
        fwd=event_score(base['Fh'],base['Fv'],m.body,e,dk)
        b1=mutate_body(m.body,e)
        wh1,wv1=bond_weights(m,b1)
        c1,_,_=linear_contrast(m,wh1,wv1,seqT,seqD)
        dlin=float(c1-base['C'])
        out_events.append(dict(
            kind=e['kind'],cell=e['cell'],dist_soma=e['dist_soma'],
            pred_adjoint=pred,
            pred_forward_only=fwd,
            dC_lin_finite=dlin,
            dC_int=float(e['dCint']),
            dC_peak=float(e['dCpeak']),
            contrib_R=float(e['contrib_R']),
            contrib_V=float(e['contrib_V']),
        ))

    def corr(keya,keyb):return safe_corr([e[keya] for e in out_events],[e[keyb] for e in out_events])
    return dict(
        seed=int(m.cfg.seed),cells=int(m.body.sum()),n_events=len(out_events),
        base_C_lin=float(base['C']),base_C_int=float(observed['base']['Cint']),base_C_peak=float(observed['base']['Cpeak']),
        laplacian_relative_error=lap_rel,
        fd_check=fdcheck,
        r_adj_lin=float(corr('pred_adjoint','dC_lin_finite')),
        r_lin_int=float(corr('dC_lin_finite','dC_int')),
        r_adj_int=float(corr('pred_adjoint','dC_int')),
        r_forward_lin=float(corr('pred_forward_only','dC_lin_finite')),
        r_adj_peak=float(corr('pred_adjoint','dC_peak')),
        r_lin_peak=float(corr('dC_lin_finite','dC_peak')),
        improvement_adj_over_forward=float(corr('pred_adjoint','dC_lin_finite')-corr('pred_forward_only','dC_lin_finite')),
        events=out_events,
    )


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=180)
    ap.add_argument('--seeds',type=int,default=12)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--steps',type=int,default=210)
    ap.add_argument('--max-each',type=int,default=6)
    ap.add_argument('--out',default='runs/adjoint_eligibility_discovery/adjoint_eligibility.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()


def selftest():
    # Laplacian conservation/symmetry sanity check.
    wh=np.ones((3,2));wv=np.ones((2,3));u=np.arange(9,dtype=float).reshape(3,3)
    L=weighted_lap(u,wh,wv)
    assert abs(np.sum(L))<1e-12
    assert L.shape==u.shape
    print('selftest ok')


def main():
    a=parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor

    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_adjoint_probe(m,a.lag,a.steps,a.max_each);rows.append(r)
        print(f"seed {seed}: r adj/lin={r['r_adj_lin']:+.3f} lin/int={r['r_lin_int']:+.3f} "
              f"adj/int={r['r_adj_int']:+.3f} fwd/lin={r['r_forward_lin']:+.3f} "
              f"fd={r['fd_check']['relative_error']:.2e}",flush=True)

    if not rows:raise SystemExit('No valid bodies')
    def vals(k):return np.asarray([r[k] for r in rows],float)
    events=[e for r in rows for e in r['events']]
    nonint=[e for e in events if abs(e['dC_int'])>1e-5]
    nonpeak=[e for e in events if abs(e['dC_peak'])>1e-5]
    sign_int=float(np.mean([np.sign(e['pred_adjoint'])==np.sign(e['dC_int']) for e in nonint])) if nonint else float('nan')
    sign_peak=float(np.mean([np.sign(e['pred_adjoint'])==np.sign(e['dC_peak']) for e in nonpeak])) if nonpeak else float('nan')
    imp=vals('improvement_adj_over_forward')
    ip,iw,il=sign_test_two_sided(imp)

    kind_summary={}
    for kind in ('add','delete'):
        es=[e for e in events if e['kind']==kind]
        if es:
            kind_summary[kind]=dict(
                n=len(es),
                pooled_r_adj_lin=float(safe_corr([e['pred_adjoint'] for e in es],[e['dC_lin_finite'] for e in es])),
                pooled_r_adj_int=float(safe_corr([e['pred_adjoint'] for e in es],[e['dC_int'] for e in es])),
                pooled_r_adj_peak=float(safe_corr([e['pred_adjoint'] for e in es],[e['dC_peak'] for e in es])),
            )

    summary=dict(
        bodies=len(rows),total_events=len(events),
        mean_fd_relative_error=float(np.mean([r['fd_check']['relative_error'] for r in rows])),
        max_fd_relative_error=float(np.max([r['fd_check']['relative_error'] for r in rows])),
        mean_laplacian_relative_error=float(np.mean(vals('laplacian_relative_error'))),
        mean_r_adj_lin=float(np.nanmean(vals('r_adj_lin'))),
        positive_r_adj_lin_bodies=int(np.sum(vals('r_adj_lin')>0)),
        mean_r_lin_int=float(np.nanmean(vals('r_lin_int'))),
        positive_r_lin_int_bodies=int(np.sum(vals('r_lin_int')>0)),
        mean_r_adj_int=float(np.nanmean(vals('r_adj_int'))),
        positive_r_adj_int_bodies=int(np.sum(vals('r_adj_int')>0)),
        pooled_sign_adj_int=sign_int,
        mean_r_forward_lin=float(np.nanmean(vals('r_forward_lin'))),
        mean_improvement_adj_over_forward=float(np.nanmean(imp)),
        improved_over_forward_bodies=int(np.sum(imp>0)),
        worse_than_forward_bodies=int(np.sum(imp<0)),
        improvement_sign_p=float(ip),
        mean_r_adj_peak=float(np.nanmean(vals('r_adj_peak'))),
        mean_r_lin_peak=float(np.nanmean(vals('r_lin_peak'))),
        pooled_sign_adj_peak=sign_peak,
        base_C_lin_vs_Cint_corr=float(safe_corr(vals('base_C_lin'),vals('base_C_int'))),
        kind_summary=kind_summary,
    )
    summary['D0_pass']=bool(summary['mean_fd_relative_error']<1e-3)
    summary['D1_pass']=bool(summary['mean_r_adj_lin']>.70 and summary['positive_r_adj_lin_bodies']>=10)
    summary['D2_pass']=bool(summary['mean_r_lin_int']>.65 and summary['positive_r_lin_int_bodies']>=10)
    summary['D3_pass']=bool(summary['mean_r_adj_int']>.55 and summary['positive_r_adj_int_bodies']>=9 and np.isfinite(sign_int) and sign_int>.70)
    summary['D4_pass']=bool(summary['mean_improvement_adj_over_forward']>.15 and summary['improved_over_forward_bodies']>=9)

    payload=dict(experiment='adjoint_eligibility_discovery_v01',
                 prereg='ADJOINT_ELIGIBILITY_DISCOVERY_PREREG_V01.md',
                 seed_start=a.seed_start,seeds_requested=a.seeds,lag=a.lag,steps=a.steps,max_each=a.max_each,
                 summary=summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nADJOINT ELIGIBILITY DISCOVERY RECEIPT')
    for k,v in summary.items():print(f' {k}: {v}')

if __name__=='__main__':
    main()
