"""Which gating memory creates the active boundary's frequency-allocation shape?

Frozen from AIS active v0.1/v0.2.  No gain or conductance is retuned.
See AIS_GATE_ABLATION_PREREG_V01.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from ais_active_probe import (
    AISConfig, ActiveAIS, _rates, _steady_gates,
    task_trace, modulated_trace, normalize_traces,
)


def gate_inf(V: float):
    am, bm, ah, bh, an, bn = _rates(V)
    return am/(am+bm), ah/(ah+bh), an/(an+bn)


class GateAblatedAIS(ActiveAIS):
    def run_mode(self, normalized_drive: np.ndarray, mode: str = 'full'):
        if mode == 'full':
            return self.run(normalized_drive)
        c = self.cfg
        x = np.clip(np.asarray(normalized_drive, float), 0.0, c.clip_drive)
        current_frames = c.input_gain * x
        V = float(c.V0)
        m, h, n = _steady_gates(V)
        sub = max(1, int(round(c.frame_ms / c.ode_dt_ms)))
        dt = c.frame_ms / sub
        volts = np.zeros(len(current_frames), float)
        spikes=[]; tm=0.0; prev=V
        instant_m = mode in ('m_instant','all_instant')
        instant_h = mode in ('h_instant','all_instant')
        instant_n = mode in ('n_instant','all_instant')
        if mode not in ('m_instant','h_instant','n_instant','all_instant'):
            raise ValueError(mode)
        for i,I in enumerate(current_frames):
            for _ in range(sub):
                am,bm,ah,bh,an,bn = _rates(V)
                dm=am*(1-m)-bm*m; dh=ah*(1-h)-bh*h; dn=an*(1-n)-bn*n
                INa=c.gNa*(m**3)*h*(V-c.ENa)
                IK=c.gK*(n**4)*(V-c.EK)
                IL=c.gL*(V-c.EL)
                dV=(float(I)-INa-IK-IL)/c.C
                prev=V; V += dt*dV
                mi,hi,ni = gate_inf(V)
                m = mi if instant_m else float(np.clip(m+dt*dm,0,1))
                h = hi if instant_h else float(np.clip(h+dt*dh,0,1))
                n = ni if instant_n else float(np.clip(n+dt*dn,0,1))
                tm += dt
                if prev < c.spike_threshold_mv <= V:
                    spikes.append(tm/c.frame_ms)
            volts[i]=V
        return volts,np.asarray(spikes,float)


def allocation(counts):
    z=np.asarray(counts,float); s=z.sum()
    return z/s if s>0 else np.full(len(z),np.nan)


def tv(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if not np.isfinite(a).all() or not np.isfinite(b).all(): return float('nan')
    return float(.5*np.abs(a-b).sum())


def peak_freq(freqs,p):
    p=np.asarray(p,float)
    if not np.isfinite(p).any(): return float('nan')
    return float(freqs[int(np.nanargmax(p))])


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='../FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=0)
    ap.add_argument('--seeds',type=int,default=24)
    ap.add_argument('--lag',type=int,default=20)
    ap.add_argument('--task-steps',type=int,default=180)
    ap.add_argument('--freqs',default='0.00625,0.0125,0.025,0.05,0.0833333,0.125')
    ap.add_argument('--freq-steps',type=int,default=640)
    ap.add_argument('--burn',type=int,default=160)
    ap.add_argument('--out',default='runs/ais_gate_ablation/ais_gate_ablation.json')
    return ap.parse_args()


def main():
    a=parse_args(); fa=Path(a.functional_arbors).resolve()
    if not fa.exists(): raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    freqs=[float(x) for x in a.freqs.split(',') if x.strip()]
    modes=['full','m_instant','h_instant','n_instant','all_instant']
    ais=GateAblatedAIS(AISConfig())
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed)); boot=m.bootstrap()
        if not boot.get('ok'): continue
        m.mature=True
        # Keep v0.1/v0.2 normalization battery exactly: task pair + all frequency traces.
        raw={'task_T':task_trace(m,a.lag,True,a.task_steps),
             'task_D':task_trace(m,a.lag,False,a.task_steps)}
        for j,f in enumerate(freqs): raw[f'f{j}']=modulated_trace(m,f,a.freq_steps,source=0)
        x,scale=normalize_traces(raw)
        rec={}
        for mode in modes:
            counts=[]
            for j,f in enumerate(freqs):
                _,sp=ais.run_mode(x[f'f{j}'],mode)
                counts.append(int(np.sum((sp>=a.burn)&(sp<a.freq_steps))))
            alloc=allocation(counts)
            band=[j for j,f in enumerate(freqs) if 0.05 <= f < 0.125]
            rec[mode]=dict(counts=counts,total=int(sum(counts)),allocation=alloc.tolist(),
                           band_share=float(np.nansum(alloc[band])) if np.isfinite(alloc).any() else float('nan'),
                           peak_freq=peak_freq(freqs,alloc))
        fp=np.asarray(rec['full']['allocation'],float)
        for mode in modes[1:]:
            rec[mode]['tv_from_full']=tv(rec[mode]['allocation'],fp)
            rec[mode]['band_share_delta']=float(rec[mode]['band_share']-rec['full']['band_share']) \
                if np.isfinite(rec[mode]['band_share']) and np.isfinite(rec['full']['band_share']) else float('nan')
            rec[mode]['peak_agrees']=bool(rec[mode]['peak_freq']==rec['full']['peak_freq']) \
                if np.isfinite(rec[mode]['peak_freq']) and np.isfinite(rec['full']['peak_freq']) else False
        rows.append(dict(seed=seed,cells=int(m.body.sum()),scale=scale,modes=rec))
        print(f'seed {seed:2d}: full total {rec["full"]["total"]:3d} peak {rec["full"]["peak_freq"]:.5f} | '
              + ' '.join(f'{q}:TV={rec[q]["tv_from_full"]:.2f}' if np.isfinite(rec[q]['tv_from_full']) else f'{q}:void'
                         for q in modes[1:]),flush=True)
    if not rows: raise SystemExit('No valid bodies')
    summary={}
    for mode in modes:
        total=np.asarray([r['modes'][mode]['total'] for r in rows],float)
        bshare=np.asarray([r['modes'][mode]['band_share'] for r in rows],float)
        item=dict(total_events_mean=float(total.mean()),total_events_median=float(np.median(total)),
                  emitting_bodies=int(np.sum(total>0)),
                  band_share_mean=float(np.nanmean(bshare)) if np.isfinite(bshare).any() else float('nan'))
        if mode!='full':
            tt=np.asarray([r['modes'][mode]['tv_from_full'] for r in rows],float)
            dd=np.asarray([r['modes'][mode]['band_share_delta'] for r in rows],float)
            item.update(tv_mean=float(np.nanmean(tt)) if np.isfinite(tt).any() else float('nan'),
                        tv_median=float(np.nanmedian(tt)) if np.isfinite(tt).any() else float('nan'),
                        tv_valid=int(np.isfinite(tt).sum()),
                        band_share_delta_mean=float(np.nanmean(dd)) if np.isfinite(dd).any() else float('nan'),
                        peak_agreement_count=int(sum(r['modes'][mode]['peak_agrees'] for r in rows)))
        summary[mode]=item
    # Mean natural event count and allocation by frequency for readability.
    freq_summary=[]
    for j,f in enumerate(freqs):
        q={'freq':f}
        for mode in modes:
            q[mode]=dict(mean_events=float(np.mean([r['modes'][mode]['counts'][j] for r in rows])),
                         mean_allocation=float(np.nanmean([r['modes'][mode]['allocation'][j] for r in rows])))
        freq_summary.append(q)
    payload=dict(experiment='ais_gate_ablation_v01',modes=modes,freqs=freqs,
                 interpretation='gating-memory ablation only; active parameters frozen from v0.1/v0.2',
                 summary=summary,frequency=freq_summary,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print('\nAIS GATING-MEMORY ABLATION')
    print(f'  full total events mean {summary["full"]["total_events_mean"]:.2f}, band share {summary["full"]["band_share_mean"]:.3f}')
    for mode in modes[1:]:
        s=summary[mode]
        print(f'  {mode:11s}: emit {s["emitting_bodies"]}/{len(rows)}  TV {s["tv_mean"]:.3f}  '
              f'band delta {s["band_share_delta_mean"]:+.3f}  peak agree {s["peak_agreement_count"]}/{len(rows)}')
    print(f'  wrote {out}')

if __name__=='__main__': main()
