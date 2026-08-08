"""Device-inspired error audit for compressed reciprocal spectral gradients.

This is a development probe, not a claim that the toy's error variables are a calibrated
model of one particular silicon chip.  The values are chosen to sit near scales used
in the in-situ-photonic-backprop literature (notably phase-error studies around
0.025-0.05 rad and small tap/readout noise factors).

The ideal K-bin gradient contribution is a sum of complex local phasor products Z_k,
with the physical gradient Re(sum Z_k).  We separate three error classes:

  phase_sigma : one phase-setting error per selected frequency, shared by all bonds;
  amp_sigma   : one multiplicative amplitude-calibration error per selected frequency;
  tap_sigma   : independent local readout noise per bond/frequency, normalized to the
                RMS magnitude of that bin's local complex products.

Shared phase/amplitude errors represent controller/port calibration.  Local tap noise
represents detector/readout error.  This distinction matters: drawing a new unrelated
phase error at every bond would model a different and much harsher failure mode.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from reciprocal_adjoint_probe import retro_source_sequence, flat_pair, normalized_l2
from spectral_correlation_compression_probe import edge_series, port_spectrum_score
from transfer_decomposition_probe import safe_corr

CONDITIONS={
    'ideal':            dict(phase=0.0,   amp=0.0,   tap=0.0),
    'phase025':         dict(phase=.025,  amp=0.0,   tap=0.0),
    'phase050':         dict(phase=.050,  amp=0.0,   tap=0.0),
    'amp025':           dict(phase=0.0,   amp=.025,  tap=0.0),
    'amp050':           dict(phase=0.0,   amp=.050,  tap=0.0),
    'tap005':           dict(phase=0.0,   amp=0.0,   tap=.005),
    'tap010':           dict(phase=0.0,   amp=0.0,   tap=.010),
    'tap020':           dict(phase=0.0,   amp=0.0,   tap=.020),
    'moderate_combo':   dict(phase=.025,  amp=.025,  tap=.010),
    'high_combo':       dict(phase=.050,  amp=.050,  tap=.020),
}


def complex_spectral_maps(m,forward_states,retro_states):
    Fh,Fv=edge_series(forward_states,False)
    Rh,Rv=edge_series(retro_states,True)
    T=len(Fh); k=np.arange(T); neg=(-k)%T
    phase=np.exp(2j*np.pi*k/T)
    FH=np.fft.fft(Fh,axis=0); FV=np.fft.fft(Fv,axis=0)
    RH=np.fft.fft(Rh,axis=0); RV=np.fft.fft(Rv,axis=0)
    fac=2.0*float(m.cfg.dt)*float(m.cfg.stiffness)/T
    ZH=fac*np.conj(RH)*phase[:,None,None]*FH[neg]
    ZV=fac*np.conj(RV)*phase[:,None,None]*FV[neg]
    return ZH,ZV


def order_complex(m,wh,wv,seq,weight):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    g=weight*np.asarray(p[1:,m.soma[0],m.soma[1]],np.complex128)
    rseq=retro_source_sequence(m,g,reverse=True)
    rp,rv,_=ae.linear_forward(m,wh,wv,rseq,store=True)
    ZH,ZV=complex_spectral_maps(m,p[:-1],rp[1:])
    eh,ev=ae.adjoint_grad(m,wh,wv,p,v,weight)
    return dict(ZH=ZH,ZV=ZV,exact=(eh,ev),port_score=port_spectrum_score(seq,rseq),E=E)


def noisy_map(ZH,ZV,order,K,pars,rng):
    kk=np.asarray(order[:min(int(K),len(order))],int)
    zh=ZH[kk].copy();zv=ZV[kk].copy()
    n=len(kk)
    if pars['phase']>0:
        d=rng.normal(0.0,float(pars['phase']),size=n)
        rot=np.exp(1j*d)
        zh*=rot[:,None,None];zv*=rot[:,None,None]
    if pars['amp']>0:
        a=1.0+rng.normal(0.0,float(pars['amp']),size=n)
        zh*=a[:,None,None];zv*=a[:,None,None]
    gh=np.real(zh);gv=np.real(zv)
    if pars['tap']>0:
        for j in range(n):
            # Noise scale is local-complex-product RMS for this bin, separately H/V.
            sh=float(np.sqrt(np.mean(np.abs(zh[j])**2))+1e-30)
            sv=float(np.sqrt(np.mean(np.abs(zv[j])**2))+1e-30)
            gh[j]+=rng.normal(0.0,float(pars['tap'])*sh,size=gh[j].shape)
            gv[j]+=rng.normal(0.0,float(pars['tap'])*sv,size=gv[j].shape)
    return np.sum(gh,axis=0),np.sum(gv,axis=0),kk


def metrics(exh,exv,ah,av):
    ex=flat_pair(exh,exv);ap=flat_pair(ah,av)
    mx=float(np.max(np.abs(ex))+1e-30);mask=np.abs(ex)>.01*mx
    return dict(corr=float(safe_corr(ex,ap)),relative_l2=normalized_l2(ex,ap),
                strong_sign_agreement=float(np.mean(np.sign(ex[mask])==np.sign(ap[mask]))) if np.any(mask) else float('nan'))


def one(m,lag,steps,reps):
    wh,wv=ae.bond_weights(m,m.body)
    sT=ae.source_sequence(m,True,lag,steps);sD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,sT,store=False);ED=ae.linear_forward(m,wh,wv,sD,store=False)
    S=ET+ED+1e-30;aT=2*ED/S**2;aD=-2*ET/S**2
    T=order_complex(m,wh,wv,sT,aT);D=order_complex(m,wh,wv,sD,aD)
    ZH=T['ZH']+D['ZH'];ZV=T['ZV']+D['ZV']
    eh=T['exact'][0]+D['exact'][0];ev=T['exact'][1]+D['exact'][1]
    score=T['port_score']+D['port_score'];order=np.argsort(score)[::-1]
    out={}
    for K in (8,16):
        out[str(K)]={}
        for ci,(name,pars) in enumerate(CONDITIONS.items()):
            qs=[]
            for rep in range(int(reps)):
                rng=np.random.default_rng(int(m.cfg.seed)*1000003 + K*1009 + ci*101 + rep)
                ah,av,kk=noisy_map(ZH,ZV,order,K,pars,rng);qs.append(metrics(eh,ev,ah,av))
            out[str(K)][name]=dict(corr=float(np.mean([q['corr'] for q in qs])),relative_l2=float(np.mean([q['relative_l2'] for q in qs])),
                                    strong_sign_agreement=float(np.mean([q['strong_sign_agreement'] for q in qs])))
    return dict(seed=int(m.cfg.seed),C=float((ET-ED)/S),results=out)


def selftest():
    # Zero error must be deterministic and phase rotation must preserve shape dimensions.
    z=np.ones((3,2,2),complex);rng=np.random.default_rng(1)
    h,v,k=noisy_map(z,z,np.arange(3),2,CONDITIONS['ideal'],rng)
    assert np.allclose(h,2.0) and np.allclose(v,2.0)
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=400);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--reps',type=int,default=8);ap.add_argument('--out',default='runs/device_error/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=one(m,a.lag,a.steps,a.reps);rows.append(r)
        print('seed',seed,'K8',[(n,round(r['results']['8'][n]['corr'],4)) for n in ('ideal','phase025','phase050','tap010','tap020','moderate_combo','high_combo')],flush=True)
    summary=dict(bodies=len(rows),reps=a.reps,conditions=CONDITIONS,K={})
    for K in ('8','16'):
        summary['K'][K]={}
        for name in CONDITIONS:
            q=[r['results'][K][name] for r in rows]
            summary['K'][K][name]=dict(mean_corr=float(np.mean([x['corr'] for x in q])),mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),mean_strong_sign_agreement=float(np.mean([x['strong_sign_agreement'] for x in q])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='device_error_dev_v01',summary=summary,rows=rows),indent=2))
    print('\nDEVICE-INSPIRED ERROR DEV')
    for K in ('8','16'):
        print('K',K)
        for n,q in summary['K'][K].items():print(n,'r',round(q['mean_corr'],5),'rel',round(q['mean_relative_l2'],5),'sign',round(q['mean_strong_sign_agreement'],5))
if __name__=='__main__':main()
