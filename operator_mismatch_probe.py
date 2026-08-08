"""Development audit for forward/adjoint operator mismatch.

The exact reciprocal-adjoint identity assumes the same linear reciprocal spatial
operator is used in the forward and retro passes. Real hardware can drift between
those passes. This probe keeps each pass internally reciprocal but perturbs every
retro-pass edge conductance by a fixed fractional error relative to the forward
operator.

For one mismatch realization, the same edge error is shared by target/distractor
retro passes. The forward pass is untouched. We audit both the full time-domain
physical overlap and the already-confirmed K=8/K=16 boundary-selected spectral
compression.

This is development only; fresh-body thresholds are set only after these reused
seeds have run.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

import adjoint_eligibility_probe as ae
import adjoint_dose_probe as ad
from reciprocal_adjoint_probe import retro_source_sequence, flat_pair, normalized_l2
from device_error_probe import complex_spectral_maps
from spectral_correlation_compression_probe import port_spectrum_score
from transfer_decomposition_probe import safe_corr

# The first run found essentially no damage through 10%, so the development-only
# sweep is deliberately extended into very large mismatch to locate a transition.
SIGMAS=(0.0,.001,.0025,.005,.01,.02,.05,.10,.20,.30,.50,.75)


def one_order(m,wh,wv,seq,weight,whr,wvr):
    p,v,E=ae.linear_forward(m,wh,wv,seq,store=True)
    g=weight*np.asarray(p[1:,m.soma[0],m.soma[1]],np.complex128)
    rseq=retro_source_sequence(m,g,reverse=True)
    rp,rv,_=ae.linear_forward(m,whr,wvr,rseq,store=True)
    ZH,ZV=complex_spectral_maps(m,p[:-1],rp[1:])
    eh,ev=ae.adjoint_grad(m,wh,wv,p,v,weight)
    return dict(ZH=ZH,ZV=ZV,exact=(eh,ev),score=port_spectrum_score(seq,rseq),E=E)


def metrics(exh,exv,ah,av):
    ex=flat_pair(exh,exv);ap=flat_pair(ah,av)
    mx=float(np.max(np.abs(ex))+1e-30);mask=np.abs(ex)>.01*mx
    return dict(corr=float(safe_corr(ex,ap)),relative_l2=normalized_l2(ex,ap),
                strong_sign_agreement=float(np.mean(np.sign(ex[mask])==np.sign(ap[mask]))) if np.any(mask) else float('nan'))


def run_rep(m,lag,steps,sigma,rep):
    wh,wv=ae.bond_weights(m,m.body)
    rng=np.random.default_rng(int(m.cfg.seed)*1000003 + int(round(sigma*1e6))*101 + int(rep))
    if sigma==0:
        dh=np.zeros_like(wh);dv=np.zeros_like(wv)
    else:
        dh=rng.normal(0.,sigma,size=wh.shape);dv=rng.normal(0.,sigma,size=wv.shape)
    # Fractional pass mismatch; clip only to keep conductances positive.
    whr=np.maximum(1e-12,wh*(1.+dh));wvr=np.maximum(1e-12,wv*(1.+dv))
    seqT=ae.source_sequence(m,True,lag,steps);seqD=ae.source_sequence(m,False,lag,steps)
    ET=ae.linear_forward(m,wh,wv,seqT,store=False);ED=ae.linear_forward(m,wh,wv,seqD,store=False)
    S=ET+ED+1e-30;aT=2*ED/S**2;aD=-2*ET/S**2
    T=one_order(m,wh,wv,seqT,aT,whr,wvr);D=one_order(m,wh,wv,seqD,aD,whr,wvr)
    ZH=T['ZH']+D['ZH'];ZV=T['ZV']+D['ZV'];eh=T['exact'][0]+D['exact'][0];ev=T['exact'][1]+D['exact'][1]
    full=metrics(eh,ev,np.real(np.sum(ZH,axis=0)),np.real(np.sum(ZV,axis=0)))
    score=T['score']+D['score'];order=np.argsort(score)[::-1]
    out={'full':full}
    for K in (8,16):
        kk=np.asarray(order[:K],int);ah=np.real(np.sum(ZH[kk],axis=0));av=np.real(np.sum(ZV[kk],axis=0))
        out[str(K)]=metrics(eh,ev,ah,av)
    base=np.concatenate([wh.ravel(),wv.ravel()]);ret=np.concatenate([whr.ravel(),wvr.ravel()])
    out['operator_relative_l2']=float(np.linalg.norm(ret-base)/(np.linalg.norm(base)+1e-30))
    return out


def one_body(m,lag,steps,reps):
    conditions={}
    for s in SIGMAS:
        qs=[run_rep(m,lag,steps,float(s),r) for r in range(int(reps))]
        conditions[f'{s:g}']={name:{k:float(np.mean([q[name][k] for q in qs])) for k in ('corr','relative_l2','strong_sign_agreement')} for name in ('full','8','16')}
        conditions[f'{s:g}']['mean_operator_relative_l2']=float(np.mean([q['operator_relative_l2'] for q in qs]))
    return dict(seed=int(m.cfg.seed),conditions=conditions)


def selftest():
    assert SIGMAS[0]==0 and .5 in SIGMAS
    print('selftest ok')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=400);ap.add_argument('--seeds',type=int,default=4);ap.add_argument('--lag',type=int,default=20);ap.add_argument('--steps',type=int,default=210);ap.add_argument('--reps',type=int,default=6);ap.add_argument('--out',default='runs/operator_mismatch/dev.json');ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    if a.selftest:selftest();return
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=one_body(m,a.lag,a.steps,a.reps);rows.append(r)
        print('seed',seed,[(s,round(r['conditions'][f'{s:g}']['full']['corr'],4),round(r['conditions'][f'{s:g}']['8']['corr'],4)) for s in SIGMAS],flush=True)
    summary=dict(bodies=len(rows),reps=a.reps,sigmas=[float(x) for x in SIGMAS],conditions={})
    for s in SIGMAS:
        key=f'{s:g}';summary['conditions'][key]={}
        summary['conditions'][key]['mean_operator_relative_l2']=float(np.mean([r['conditions'][key]['mean_operator_relative_l2'] for r in rows]))
        for name in ('full','8','16'):
            q=[r['conditions'][key][name] for r in rows];summary['conditions'][key][name]=dict(mean_corr=float(np.mean([x['corr'] for x in q])),mean_relative_l2=float(np.mean([x['relative_l2'] for x in q])),mean_strong_sign_agreement=float(np.mean([x['strong_sign_agreement'] for x in q])))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(dict(experiment='operator_mismatch_dev_v02',summary=summary,rows=rows),indent=2))
    print('\nOPERATOR MISMATCH DEV')
    for s in SIGMAS:
        q=summary['conditions'][f'{s:g}'];print('sigma',s,'oprel',round(q['mean_operator_relative_l2'],4),'full',round(q['full']['mean_corr'],5),round(q['full']['mean_relative_l2'],5),'K8',round(q['8']['mean_corr'],5),round(q['8']['mean_relative_l2'],5),'K16',round(q['16']['mean_corr'],5),round(q['16']['mean_relative_l2'],5))
if __name__=='__main__':main()
