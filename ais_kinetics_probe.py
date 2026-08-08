"""Do h/n gate time constants set the active boundary's frequency-selection band?

Parameters and upstream system are frozen from AIS active v0.1/v0.2.
Only dh/dt or dn/dt kinetic speed is multiplied by 0.5 or 2.0.
See AIS_KINETICS_PREREG_V01.md.
"""
from __future__ import annotations

import argparse, json, math, sys
from pathlib import Path
import numpy as np

from ais_active_probe import AISConfig, ActiveAIS, _rates, _steady_gates, task_trace, modulated_trace, normalize_traces


class KineticAIS(ActiveAIS):
    def run_scaled(self, normalized_drive, h_scale=1.0, n_scale=1.0):
        c=self.cfg
        x=np.clip(np.asarray(normalized_drive,float),0.0,c.clip_drive)
        current=c.input_gain*x
        V=float(c.V0); m,h,n=_steady_gates(V)
        sub=max(1,int(round(c.frame_ms/c.ode_dt_ms))); dt=c.frame_ms/sub
        volts=np.zeros(len(current)); spikes=[]; tm=0.0; prev=V
        for i,I in enumerate(current):
            for _ in range(sub):
                am,bm,ah,bh,an,bn=_rates(V)
                dm=am*(1-m)-bm*m
                dh=h_scale*(ah*(1-h)-bh*h)
                dn=n_scale*(an*(1-n)-bn*n)
                INa=c.gNa*(m**3)*h*(V-c.ENa)
                IK=c.gK*(n**4)*(V-c.EK)
                IL=c.gL*(V-c.EL)
                dV=(float(I)-INa-IK-IL)/c.C
                prev=V; V += dt*dV
                m=float(np.clip(m+dt*dm,0,1)); h=float(np.clip(h+dt*dh,0,1)); n=float(np.clip(n+dt*dn,0,1))
                tm += dt
                if prev < c.spike_threshold_mv <= V: spikes.append(tm/c.frame_ms)
            volts[i]=V
        return volts,np.asarray(spikes,float)


def alloc(counts):
    z=np.asarray(counts,float); s=z.sum(); return z/s if s>0 else np.full(len(z),np.nan)

def tv(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(.5*np.abs(a-b).sum()) if np.isfinite(a).all() and np.isfinite(b).all() else float('nan')

def center(freqs,p):
    p=np.asarray(p,float)
    if not np.isfinite(p).all(): return float('nan')
    return float(2**np.sum(p*np.log2(np.asarray(freqs,float))))

def sign_two_sided(w,l):
    n=w+l
    if n==0:return float('nan')
    k=min(w,l); tail=sum(math.comb(n,i) for i in range(k+1))/2**n
    return min(1.0,2*tail)


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='../FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=0); ap.add_argument('--seeds',type=int,default=24)
    ap.add_argument('--lag',type=int,default=20); ap.add_argument('--task-steps',type=int,default=180)
    ap.add_argument('--freqs',default='0.00625,0.0125,0.025,0.05,0.0833333,0.125')
    ap.add_argument('--freq-steps',type=int,default=640); ap.add_argument('--burn',type=int,default=160)
    ap.add_argument('--out',default='runs/ais_kinetics/ais_kinetics.json')
    return ap.parse_args()


def main():
    a=parse_args(); fa=Path(a.functional_arbors).resolve()
    if not fa.exists(): raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    freqs=[float(x) for x in a.freqs.split(',') if x.strip()]
    variants={'full':(1.,1.),'h_slow':(.5,1.),'h_fast':(2.,1.),'n_slow':(1.,.5),'n_fast':(1.,2.)}
    ais=KineticAIS(AISConfig()); rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        raw={'task_T':task_trace(m,a.lag,True,a.task_steps),'task_D':task_trace(m,a.lag,False,a.task_steps)}
        for j,f in enumerate(freqs): raw[f'f{j}']=modulated_trace(m,f,a.freq_steps,source=0)
        x,scale=normalize_traces(raw)
        rec={}
        for name,(hs,ns) in variants.items():
            counts=[]
            for j,f in enumerate(freqs):
                _,sp=ais.run_scaled(x[f'f{j}'],hs,ns)
                counts.append(int(np.sum((sp>=a.burn)&(sp<a.freq_steps))))
            p=alloc(counts); band=[j for j,f in enumerate(freqs) if .05<=f<.125]
            rec[name]=dict(counts=counts,total=int(sum(counts)),allocation=p.tolist(),center_freq=center(freqs,p),
                           band_share=float(np.nansum(p[band])) if np.isfinite(p).any() else float('nan'))
        fp=rec['full']['allocation']
        for name in variants:
            if name=='full':continue
            rec[name]['tv_from_full']=tv(rec[name]['allocation'],fp)
            rec[name]['center_log2_shift']=float(np.log2(rec[name]['center_freq']/rec['full']['center_freq'])) \
                if rec[name]['center_freq']>0 and rec['full']['center_freq']>0 else float('nan')
            rec[name]['band_share_delta']=float(rec[name]['band_share']-rec['full']['band_share']) \
                if np.isfinite(rec[name]['band_share']) and np.isfinite(rec['full']['band_share']) else float('nan')
        rows.append(dict(seed=seed,cells=int(m.body.sum()),scale=scale,variants=rec))
        print(f'seed {seed:2d}: center {rec["full"]["center_freq"]:.4f} | '
              f'h .5/{2}: {rec["h_slow"]["center_freq"]:.4f}/{rec["h_fast"]["center_freq"]:.4f}  '
              f'n .5/{2}: {rec["n_slow"]["center_freq"]:.4f}/{rec["n_fast"]["center_freq"]:.4f}',flush=True)
    if not rows:raise SystemExit('No bodies')
    summary={}
    for name in variants:
        tot=np.array([r['variants'][name]['total'] for r in rows],float)
        cf=np.array([r['variants'][name]['center_freq'] for r in rows],float)
        bs=np.array([r['variants'][name]['band_share'] for r in rows],float)
        q=dict(emitting_bodies=int(np.sum(tot>0)),total_mean=float(tot.mean()),
               center_freq_mean=float(np.nanmean(cf)) if np.isfinite(cf).any() else float('nan'),
               center_freq_median=float(np.nanmedian(cf)) if np.isfinite(cf).any() else float('nan'),
               band_share_mean=float(np.nanmean(bs)) if np.isfinite(bs).any() else float('nan'))
        if name!='full':
            sh=np.array([r['variants'][name]['center_log2_shift'] for r in rows],float)
            tvv=np.array([r['variants'][name]['tv_from_full'] for r in rows],float)
            q.update(center_log2_shift_mean=float(np.nanmean(sh)) if np.isfinite(sh).any() else float('nan'),
                     shift_positive=int(np.sum(sh>0)),shift_negative=int(np.sum(sh<0)),
                     tv_mean=float(np.nanmean(tvv)) if np.isfinite(tvv).any() else float('nan'))
        summary[name]=q
    # Registered opposite-order checks for each gate among bodies where both variants emit.
    pairs={}
    for gate in ('h','n'):
        slow=gate+'_slow'; fast=gate+'_fast'; dif=[]
        for r in rows:
            cs=r['variants'][slow]['center_freq']; cf=r['variants'][fast]['center_freq']
            if np.isfinite(cs) and np.isfinite(cf) and cs>0 and cf>0:
                dif.append(float(np.log2(cf/cs)))
        arr=np.asarray(dif,float); w=int(np.sum(arr>0)); l=int(np.sum(arr<0))
        pairs[gate]=dict(valid=len(arr),mean_fast_minus_slow_log2=float(arr.mean()) if len(arr) else float('nan'),
                         fast_higher=w,fast_lower=l,sign_p=sign_two_sided(w,l))
    freq_summary=[]
    for j,f in enumerate(freqs):
        z={'freq':f}
        for name in variants:
            z[name]=dict(mean_events=float(np.mean([r['variants'][name]['counts'][j] for r in rows])),
                         mean_allocation=float(np.nanmean([r['variants'][name]['allocation'][j] for r in rows])))
        freq_summary.append(z)
    payload=dict(experiment='ais_kinetics_v01',variants=variants,freqs=freqs,
                 summary=summary,paired_order=pairs,frequency=freq_summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nAIS KINETICS RECEIPT')
    print(f' full center {summary["full"]["center_freq_mean"]:.5f}')
    for name in ('h_slow','h_fast','n_slow','n_fast'):
        s=summary[name];print(f' {name:7s}: emit {s["emitting_bodies"]}/{len(rows)} center {s["center_freq_mean"]:.5f} shift {s["center_log2_shift_mean"]:+.3f} oct TV {s["tv_mean"]:.3f}')
    for gate,q in pairs.items(): print(f' {gate}: fast-vs-slow {q["mean_fast_minus_slow_log2"]:+.3f} oct, higher/lower {q["fast_higher"]}/{q["fast_lower"]}, p={q["sign_p"]:.5g}')
    print(' wrote',out)

if __name__=='__main__':main()
