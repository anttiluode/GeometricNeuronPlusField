"""Development probe for compressing the time-domain forward x adjoint correlation.

The exact causal reciprocal replay still needs the local anti-diagonal temporal
correlation

    sum_r forward[T-1-r] * conj(retro[r]).

A DFT turns that anti-diagonal correlation into a sum of independent frequency-bin
products.  This script asks whether a *small* set of lock-in/phasor accumulators per
bond can approximate the exact gradient, avoiding a full T-sample local history.

Two rankings are compared:
  oracle   rank bins by their actual total absolute gradient contribution;
  boundary rank bins using only external/source and soma-return spectra.

The all-bin reconstruction is an algebraic positive control and must reproduce the
exact reciprocal gradient to floating-point precision.  This is a development probe;
no sparse-frequency claim is preregistered in advance.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad  # patch exact mature-boundary Laplacian
from reciprocal_adjoint_probe import retro_source_sequence, flat_pair, normalized_l2
from transfer_decomposition_probe import safe_corr

KLIST=(1,2,4,8,12,16,24,32,48,64,96,128,160,210)


def edge_series(states, retro=False):
    """Return horizontal/vertical oriented edge differences over time."""
    z=np.asarray(states,np.complex128)
    if not retro:
        h=z[:,:,1:]-z[:,:,:-1]
        v=z[:,1:,:]-z[:,:-1,:]
    else:
        # Match reciprocal_adjoint overlap sign convention: left-right, up-down.
        h=z[:,:,:-1]-z[:,:,1:]
        v=z[:,:-1,:]-z[:,1:,:]
    return h,v


def anti_corr_spectral_maps(m, forward_states, retro_states):
    """Per-frequency maps whose sum equals local reversed-forward x retro overlap.

    forward_states are p[0:T]; retro_states are physical replay states rp[1:T+1].
    For length T, FFT(reverse(f))[k] = exp(+i 2pi k/T) FFT(f)[-k].
    """
    Fh,Fv=edge_series(forward_states,False)
    Rh,Rv=edge_series(retro_states,True)
    T=len(Fh); k=np.arange(T); neg=(-k)%T
    phase=np.exp(2j*np.pi*k/T)
    FH=np.fft.fft(Fh,axis=0); FV=np.fft.fft(Fv,axis=0)
    RH=np.fft.fft(Rh,axis=0); RV=np.fft.fft(Rv,axis=0)
    fac=2.0*float(m.cfg.dt)*float(m.cfg.stiffness)/T
    GH=fac*np.real(np.conj(RH)*phase[:,None,None]*FH[neg])
    GV=fac*np.real(np.conj(RV)*phase[:,None,None]*FV[neg])
    return GH,GV


def port_spectrum_score(seq, retro_seq):
    """Hardware-accessible ranking proxy from external drive and soma-return source."""
    S=np.fft.fft(np.asarray(seq,np.complex128),axis=0)
    R=np.fft.fft(np.asarray(retro_seq,np.complex128),axis=0)
    T=len(S); k=np.arange(T); neg=(-k)%T
    sn=np.sqrt(np.sum(np.abs(S[neg])**2,axis=(1,2)))
    rn=np.sqrt(np.sum(np.abs(R)**2,axis=(1,2)))
    return sn*rn


def order_data(m,wh,wv,seq,weight):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    g=weight*np.asarray(p[1:,m.soma[0],m.soma[1]],np.complex128)
    rseq=retro_source_sequence(m,g,reverse=True)
    rp,rv,_=ae.linear_forward(m,wh,wv,rseq,store=True)
    GH,GV=anti_corr_spectral_maps(m,p[:-1],rp[1:])
    eh,ev=ae.adjoint_grad(m,wh,wv,p,v,weight)
    return dict(GH=GH,GV=GV,exact=(eh,ev),port_score=port_spectrum_score(seq,rseq),E=E)


def approx_metrics(exact_h,exact_v,GH,GV,order,K):
    kk=np.asarray(order[:min(int(K),len(order))],int)
    ah=np.sum(GH[kk],axis=0); av=np.sum(GV[kk],axis=0)
    ex=flat_pair(exact_h,exact_v); ap=flat_pair(ah,av)
    mx=np.max(np.abs(ex))+1e-30
    mask=np.abs(ex)>.01*mx
    sign=float(np.mean(np.sign(ex[mask])==np.sign(ap[mask]))) if np.any(mask) else float('nan')
    return dict(K=int(len(kk)),corr=float(safe_corr(ex,ap)),relative_l2=normalized_l2(ex,ap),strong_sign_agreement=sign)


def one(m,lag,steps):
    wh,wv=ae.bond_weights(m,m.body)
    seqT=ae.source_sequence(m,True,lag,steps); seqD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,seqT,store=False); ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30; aT=2*ED/S**2; aD=-2*ET/S**2
    T=order_data(m,wh,wv,seqT,aT); D=order_data(m,wh,wv,seqD,aD)
    GH=T['GH']+D['GH']; GV=T['GV']+D['GV']
    eh=T['exact'][0]+D['exact'][0]; ev=T['exact'][1]+D['exact'][1]
    ex=flat_pair(eh,ev); allsp=flat_pair(np.sum(GH,axis=0),np.sum(GV,axis=0))
    allbin=dict(corr=float(safe_corr(ex,allsp)),relative_l2=normalized_l2(ex,allsp))

    oracle_score=np.sum(np.abs(GH),axis=(1,2))+np.sum(np.abs(GV),axis=(1,2))
    port_score=T['port_score']+D['port_score']
    oracle=np.argsort(oracle_score)[::-1]
    boundary=np.argsort(port_score)[::-1]
    ks=sorted(set(min(int(k),steps) for k in KLIST))
    curves={}
    for name,order in [('oracle',oracle),('boundary',boundary)]:
        curves[name]={str(k):approx_metrics(eh,ev,GH,GV,order,k) for k in ks}

    # How concentrated is the signed-gradient spectral mass before any truncation?
    mass=oracle_score/(np.sum(oracle_score)+1e-30); cs=np.cumsum(np.sort(mass)[::-1])
    k50=int(np.searchsorted(cs,.5)+1); k80=int(np.searchsorted(cs,.8)+1); k95=int(np.searchsorted(cs,.95)+1)
    return dict(seed=int(m.cfg.seed),C=float((ET-ED)/S),steps=int(steps),allbin=allbin,
                spectral_mass_k50=k50,spectral_mass_k80=k80,spectral_mass_k95=k95,
                oracle_top_bins=[int(x) for x in oracle[:16]],boundary_top_bins=[int(x) for x in boundary[:16]],curves=curves)


def selftest():
    rng=np.random.default_rng(4); T=17
    f=rng.normal(size=T)+1j*rng.normal(size=T); b=rng.normal(size=T)+1j*rng.normal(size=T)
    F=np.fft.fft(f); B=np.fft.fft(b); k=np.arange(T)
    x=np.sum(f[::-1]*np.conj(b))
    y=np.sum(np.exp(2j*np.pi*k/T)*F[(-k)%T]*np.conj(B))/T
    assert abs(x-y)<1e-10
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--functional-arbors',default='FunctionalArbors'); ap.add_argument('--seed-start',type=int,default=240); ap.add_argument('--seeds',type=int,default=4); ap.add_argument('--lag',type=int,default=20); ap.add_argument('--steps',type=int,default=210); ap.add_argument('--out',default='runs/spectral_correlation_compression/dev.json'); ap.add_argument('--selftest',action='store_true'); a=ap.parse_args()
    if a.selftest: selftest(); return
    fa=Path(a.functional_arbors).resolve(); sys.path.insert(0,str(fa)); from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True; r=one(m,a.lag,a.steps); rows.append(r)
        print('seed',seed,'all',round(r['allbin']['corr'],8),r['allbin']['relative_l2'],'mass',r['spectral_mass_k50'],r['spectral_mass_k80'],r['spectral_mass_k95'],
              'B8',round(r['curves']['boundary']['8']['corr'],3),'B16',round(r['curves']['boundary']['16']['corr'],3),'B32',round(r['curves']['boundary']['32']['corr'],3),flush=True)
    if not rows: raise SystemExit('No valid bodies')
    ks=sorted(set(min(int(k),a.steps) for k in KLIST))
    summary=dict(bodies=len(rows),steps=a.steps,
                 allbin_mean_corr=float(np.mean([r['allbin']['corr'] for r in rows])),
                 allbin_max_relative_l2=float(np.max([r['allbin']['relative_l2'] for r in rows])),
                 mean_k50=float(np.mean([r['spectral_mass_k50'] for r in rows])),
                 mean_k80=float(np.mean([r['spectral_mass_k80'] for r in rows])),
                 mean_k95=float(np.mean([r['spectral_mass_k95'] for r in rows])),curves={})
    for name in ('oracle','boundary'):
        summary['curves'][name]={}
        for k in ks:
            q=[r['curves'][name][str(k)] for r in rows]
            summary['curves'][name][str(k)]=dict(mean_corr=float(np.mean([x['corr'] for x in q])),
                                                  mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),
                                                  mean_strong_sign_agreement=float(np.mean([x['strong_sign_agreement'] for x in q])))
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(dict(experiment='spectral_correlation_compression_dev_v01',summary=summary,rows=rows),indent=2))
    print('\nSPECTRAL CORRELATION COMPRESSION DEV')
    print('all-bin',summary['allbin_mean_corr'],summary['allbin_max_relative_l2'],'k50/80/95',summary['mean_k50'],summary['mean_k80'],summary['mean_k95'])
    for k in ks:
        b=summary['curves']['boundary'][str(k)]; o=summary['curves']['oracle'][str(k)]
        print('K',k,'boundary r/rel',round(b['mean_corr'],4),round(b['mean_relative_l2'],4),'oracle',round(o['mean_corr'],4),round(o['mean_relative_l2'],4))

if __name__=='__main__': main()
