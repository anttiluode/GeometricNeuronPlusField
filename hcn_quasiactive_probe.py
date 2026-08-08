"""Complete minimal quasi-active HCN-gradient development probe.

The first `hcn_impedance_probe.py` intentionally tested a tiny delayed-restorative
material term.  Its morphology-indexed gradient synchronized the soma, but the
mechanism diagnostic showed the wrong distance-dependent phase sign for an HCN
interpretation.

The omission was biologically material: a linearized voltage-dependent channel
contributes BOTH

  1. resting/static membrane conductance (gamma_R), and
  2. delayed feedback (mu*m).

For an HCN-like restorative current mu > 0.  This probe therefore uses one
spatial density field d(x) to generate both terms:

  v'   = K L psi - damping*v - restoring*psi
         - d(x)*psi - mu_ratio*d(x)*z + source
  z'   = (psi-z)/tau_h

The same density histogram is compared as smooth soma->distal, shuffled,
uniform-mean, and reversed placement.  This remains a second-order wave toy with
quasi-active membrane material, NOT a conductance-based neuron.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import splu

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import (
    weighted_laplacian_sparse, build_profiles, injection_sites, circular_rms,
)


def safe_corr(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    ok=np.isfinite(a)&np.isfinite(b);a=a[ok];b=b[ok]
    if len(a)<3 or np.std(a)<1e-12 or np.std(b)<1e-12:return float('nan')
    return float(np.corrcoef(a,b)[0,1])


def circ_delta(a,b):
    return np.angle(np.exp(1j*(np.asarray(a)-np.asarray(b))))


def qa_matrix(m,L,density,omega,tau_h,mu_ratio):
    c=m.cfg;dt=float(c.dt);lam=np.exp(1j*float(omega))
    ah=(dt/float(tau_h))/(lam-1.0+dt/float(tau_h))
    dyn=((lam-1.0+dt*float(c.damping))*(lam-1.0)/(dt*dt*lam)+float(c.restoring))
    d=np.asarray(density,float).ravel()
    local=d*(1.0+float(mu_ratio)*ah)
    N=m.body.size
    return (dyn*eye(N,dtype=np.complex128,format='csr')+
            diags(local,0,shape=(N,N),dtype=np.complex128)-
            float(c.stiffness)*L)


def transfer_vectors(m,L,density,omega,tau_h,mu_ratio,sites):
    A=qa_matrix(m,L,density,omega,tau_h,mu_ratio).tocsc();lu=splu(A)
    h,w=m.body.shape;N=h*w
    idx=np.asarray([p[0]*w+p[1] for p in sites],int)
    B=np.zeros((N,len(sites)),np.complex128);B[idx,np.arange(len(sites))]=1.0
    X=lu.solve(B);si=m.soma[0]*w+m.soma[1]
    return X[si,np.arange(len(sites))],X[idx,np.arange(len(sites))]


def profile_frequency(m,L,prof,omega,tau,mu_ratio,sites):
    distmap=m.graph_distance_from_soma();dist=np.asarray([distmap[p] for p in sites],float)
    H={};Q={}
    for name in ('zero','uniform','smooth','reverse'):
        H[name],Q[name]=transfer_vectors(m,L,prof[name],omega,tau,mu_ratio,sites)
    sh=[transfer_vectors(m,L,z,omega,tau,mu_ratio,sites) for z in prof['shuffled']]
    hsh=np.asarray([x[0] for x in sh]);qsh=np.asarray([x[1] for x in sh])

    adv0=circ_delta(np.angle(H['smooth']),np.angle(H['zero']))
    advu=circ_delta(np.angle(H['smooth']),np.angle(H['uniform']))
    q25,q75=np.quantile(dist,[.25,.75]);prox=dist<=q25;far=dist>=q75
    def adv_stats(x):
        return dict(mean=float(np.mean(x)),corr_distance=safe_corr(dist,x),
                    distal_minus_prox=float(np.mean(x[far])-np.mean(x[prox])),
                    positive_fraction=float(np.mean(x>0)))

    phase={k:circular_rms(np.angle(H[k])) for k in H}
    phase['shuffle']=float(np.mean([circular_rms(np.angle(x)) for x in hsh]))
    local={k:circular_rms(np.angle(Q[k])) for k in Q}
    local['shuffle']=float(np.mean([circular_rms(np.angle(x)) for x in qsh]))
    amp={k:float(np.median(np.abs(H[k]))) for k in H}
    amp['shuffle']=float(np.mean([np.median(np.abs(x)) for x in hsh]))
    return dict(omega=float(omega),soma_phase_rms=phase,local_phase_rms=local,
                soma_amp_median=amp,advance_vs_zero=adv_stats(adv0),
                advance_vs_uniform=adv_stats(advu))


def body_probe(m,g0,ratio,tau,mu_ratio,omegas,nshuffle):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    sites=injection_sites(m);rng=np.random.default_rng(int(m.cfg.seed)+424_242)
    prof=build_profiles(m,g0,ratio,rng,nshuffle)
    rows=[profile_frequency(m,L,prof,o,tau,mu_ratio,sites) for o in omegas]
    return dict(seed=int(m.cfg.seed),cells=int(m.body.sum()),sites=len(sites),byfreq=rows)


def summarize(rows,omegas):
    freq={}
    for om in omegas:
        q=[z for r in rows for z in r['byfreq'] if abs(z['omega']-om)<1e-12]
        def M(path):
            vals=[]
            for x in q:
                y=x
                for k in path:y=y[k]
                vals.append(float(y))
            return float(np.mean(vals))
        z=dict(
            zero=M(('soma_phase_rms','zero')),
            uniform=M(('soma_phase_rms','uniform')),
            smooth=M(('soma_phase_rms','smooth')),
            shuffle=M(('soma_phase_rms','shuffle')),
            reverse=M(('soma_phase_rms','reverse')),
            smooth_local=M(('local_phase_rms','smooth')),
            zero_local=M(('local_phase_rms','zero')),
            smooth_amp=M(('soma_amp_median','smooth')),
            uniform_amp=M(('soma_amp_median','uniform')),
            shuffle_amp=M(('soma_amp_median','shuffle')),
            advance0_mean=M(('advance_vs_zero','mean')),
            advance0_corr=M(('advance_vs_zero','corr_distance')),
            advance0_distal_minus_prox=M(('advance_vs_zero','distal_minus_prox')),
            advance_uniform_corr=M(('advance_vs_uniform','corr_distance')),
            advance_uniform_distal_minus_prox=M(('advance_vs_uniform','distal_minus_prox')),
        )
        z['gain_vs_shuffle']=z['shuffle']-z['smooth']
        z['gain_vs_uniform']=z['uniform']-z['smooth']
        z['gain_vs_zero']=z['zero']-z['smooth']
        z['local_retention']=z['smooth_local']/(z['zero_local']+1e-30)
        z['amp_ratio_shuffle']=z['smooth_amp']/(z['shuffle_amp']+1e-30)
        # Development score only: force the candidate to earn BOTH synchrony and
        # the HCN-like distal lead sign.  Penalize collapse of the local field.
        lead=min(z['advance0_distal_minus_prox'],.75)
        lead_penalty=max(0.0,-z['advance0_distal_minus_prox'])
        local_penalty=max(0.0,.65-z['local_retention'])
        z['bio_score']=(z['gain_vs_shuffle']+z['gain_vs_uniform']+.25*z['gain_vs_zero']+
                        .20*max(0.0,lead)-.40*lead_penalty-1.5*local_penalty)
        freq[str(float(om))]=z
    # Best adjacent pair, not isolated best frequency.
    oms=list(map(float,omegas));pairs=[]
    for a,b in zip(oms[:-1],oms[1:]):
        za=freq[str(a)];zb=freq[str(b)]
        pairs.append(dict(omegas=[a,b],mean_score=float((za['bio_score']+zb['bio_score'])/2),
                          mean_gain_shuffle=float((za['gain_vs_shuffle']+zb['gain_vs_shuffle'])/2),
                          mean_lead=float((za['advance0_distal_minus_prox']+zb['advance0_distal_minus_prox'])/2),
                          min_local_retention=float(min(za['local_retention'],zb['local_retention']))))
    pairs.sort(key=lambda x:x['mean_score'],reverse=True)
    return dict(bodies=len(rows),frequency=freq,best_adjacent_pair=pairs[0],pairs=pairs)


def csvfloats(s):return [float(x) for x in str(s).split(',') if x.strip()]

def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=560)
    ap.add_argument('--seeds',type=int,default=4)
    ap.add_argument('--g0s',default='0.005,0.01,0.02')
    ap.add_argument('--ratios',default='5,7,10')
    ap.add_argument('--taus',default='2,4,6')
    ap.add_argument('--mus',default='0.5,1,2')
    ap.add_argument('--omegas',default='0.02,0.03,0.04,0.05,0.06,0.08')
    ap.add_argument('--nshuffle',type=int,default=3)
    ap.add_argument('--out',default='runs/hcn_quasiactive/dev.json')
    ap.add_argument('--selftest',action='store_true')
    return ap.parse_args()

def selftest():
    # Restorative branch has positive static and delayed coefficients.
    assert np.real(1+2*(.1/(.1+1j*.2)))>1
    print('selftest ok')

def main():
    a=parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    g0s=csvfloats(a.g0s);ratios=csvfloats(a.ratios);taus=csvfloats(a.taus);mus=csvfloats(a.mus);omegas=csvfloats(a.omegas)
    bodies=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if boot.get('ok'):m.mature=True;bodies.append(m)
    if not bodies:raise SystemExit('No valid bodies')
    candidates=[]
    for g0 in g0s:
      for ratio in ratios:
       for tau in taus:
        for mu in mus:
         rows=[body_probe(m,g0,ratio,tau,mu,omegas,a.nshuffle) for m in bodies]
         s=summarize(rows,omegas);rec=dict(g0=g0,ratio=ratio,tau_h=tau,mu_ratio=mu,summary=s,rows=rows)
         candidates.append(rec);bp=s['best_adjacent_pair']
         print(f"g0={g0:g} R={ratio:g} tau={tau:g} mu={mu:g}: pair={bp['omegas']} score={bp['mean_score']:+.4f} "
               f"gain={bp['mean_gain_shuffle']:+.4f} lead={bp['mean_lead']:+.4f} local={bp['min_local_retention']:.3f}",flush=True)
    candidates.sort(key=lambda x:x['summary']['best_adjacent_pair']['mean_score'],reverse=True)
    payload=dict(experiment='complete_quasiactive_hcn_gradient_dev_v01',development_only=True,
                 model_warning='second-order wave + linear quasi-active membrane proxy; not conductance-based HCN',
                 seed_start=a.seed_start,seeds_requested=a.seeds,bodies=len(bodies),omegas=omegas,
                 candidates=candidates)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nTOP COMPLETE-QA CANDIDATES')
    for x in candidates[:10]:
        b=x['summary']['best_adjacent_pair']
        print(f"g0={x['g0']:g} R={x['ratio']:g} tau={x['tau_h']:g} mu={x['mu_ratio']:g} "
              f"pair={b['omegas']} score={b['mean_score']:+.5f} gain={b['mean_gain_shuffle']:+.5f} "
              f"lead={b['mean_lead']:+.5f} local={b['min_local_retention']:.3f}")

if __name__=='__main__':main()
