"""Two-tone cross-frequency envelope/beat probe.

Purpose
-------
Return to the gamma/theta intuition without the cheap finite-burst confound.
A source waveform contains exactly two continuous complex tones w1,w2.  In the
linear spectrum there is no component at Delta=|w2-w1|.  A quadratic/power
readout generates the difference-frequency phasor

    B_j = H_j(w1) * conj(H_j(w2))

for source location j.  Its phase is the carrier phase difference and therefore
measures envelope/group-delay alignment across source locations.

The strong signature is NOT simply high beat coherence.  It is:

  beat phase coherent across locations
  while individual carrier phase remains substantially less coherent.

We compare passive/zero material, uniform material, the confirmed hand gradient,
and the self-organized material learned ONLY on the low direct frequencies
(.03,.04).  Thus high-carrier beat performance is out-of-objective transfer.

This is a development probe.  Frequencies are simulation angular frequencies;
calling the upper pairs "gamma-like" is only a scale-separation analogy, not a
Hz mapping.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import adjoint_eligibility_probe as ae
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites
from hcn_material_learning import train_density, solve_fields


def circ_coherence(z):
    z=np.asarray(z,np.complex128)
    u=z/np.maximum(np.abs(z),1e-30)
    return float(abs(np.mean(u))**2)


def profile_transfers(m,L,density,omegas,tau,mu,sites):
    out={}
    for om in sorted(set(float(x) for x in omegas)):
        H,_,_,_=solve_fields(m,L,density,om,tau,mu,sites)
        out[om]=H
    return out


def pair_metrics(T,w1,w2):
    d=abs(float(w2)-float(w1));h1=T[float(w1)];h2=T[float(w2)];hd=T[d]
    beat=h1*np.conj(h2)
    c1=circ_coherence(h1);c2=circ_coherence(h2);cb=circ_coherence(beat);cd=circ_coherence(hd)
    return dict(
        w1=float(w1),w2=float(w2),delta=float(d),
        carrier1_R2=c1,carrier2_R2=c2,mean_carrier_R2=float((c1+c2)/2),
        beat_R2=cb,direct_delta_R2=cd,
        beat_minus_carrier=float(cb-(c1+c2)/2),
        beat_minus_direct_delta=float(cb-cd),
        median_beat_amplitude=float(np.median(np.abs(beat))),
        median_carrier_amplitude=float(np.median((np.abs(h1)+np.abs(h2))/2)),
    )


def build_materials(m,L,train_omegas,tau,mu,g0,ratio,steps,step_fraction,nshuffle):
    body=m.body.astype(bool);h,w=body.shape
    cells=[tuple(map(int,p)) for p in np.argwhere(body)]
    ids=np.asarray([p[0]*w+p[1] for p in cells],int)
    dist=m.graph_distance_from_soma().astype(float);dmax=max(float(dist[body].max()),1.0)
    hand=np.zeros_like(dist,float);hand[body]=float(g0)*(1+(float(ratio)-1)*dist[body]/dmax)
    total=float(hand[body].sum());cap=float(g0)*float(ratio)
    uniform=np.zeros_like(hand);uniform[body]=total/len(ids)
    zero=np.zeros_like(hand)
    sites=injection_sites(m)
    learned,_,_=train_density(m,L,uniform,train_omegas,tau,mu,sites,ids,total,cap,steps,step_fraction)
    rng=np.random.default_rng(int(m.cfg.seed)+771_771)
    sh=[];vals=learned[body].copy()
    for _ in range(int(nshuffle)):
        vv=vals.copy();rng.shuffle(vv);z=np.zeros_like(learned);z[body]=vv;sh.append(z)
    return dict(zero=zero,uniform=uniform,hand=hand,learned=learned,shuffle_learned=sh),sites


def body_probe(m,pairs,train_omegas,tau,mu,g0,ratio,steps,step_fraction,nshuffle):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    mats,sites=build_materials(m,L,train_omegas,tau,mu,g0,ratio,steps,step_fraction,nshuffle)
    needed=[]
    for a,b in pairs:needed.extend([a,b,abs(b-a)])
    prof={}
    for name in ('zero','uniform','hand','learned'):
        T=profile_transfers(m,L,mats[name],needed,tau,mu,sites)
        prof[name]=[pair_metrics(T,a,b) for a,b in pairs]
    sh=[]
    for z in mats['shuffle_learned']:
        T=profile_transfers(m,L,z,needed,tau,mu,sites)
        sh.append([pair_metrics(T,a,b) for a,b in pairs])
    shmean=[]
    for k,(a,b) in enumerate(pairs):
        keys=['carrier1_R2','carrier2_R2','mean_carrier_R2','beat_R2','direct_delta_R2','beat_minus_carrier','beat_minus_direct_delta','median_beat_amplitude','median_carrier_amplitude']
        q=dict(w1=float(a),w2=float(b),delta=float(abs(b-a)))
        for key in keys:q[key]=float(np.mean([x[k][key] for x in sh]))
        shmean.append(q)
    prof['shuffle_learned']=shmean
    return dict(seed=int(m.cfg.seed),cells=int(m.body.sum()),sites=len(sites),profiles=prof)


def summarize(rows,pairs):
    out={}
    for pi,(a,b) in enumerate(pairs):
        rec={}
        for name in ('zero','uniform','hand','learned','shuffle_learned'):
            q=[r['profiles'][name][pi] for r in rows]
            rec[name]={k:float(np.mean([x[k] for x in q])) for k in ('mean_carrier_R2','beat_R2','direct_delta_R2','beat_minus_carrier','beat_minus_direct_delta','median_beat_amplitude')}
        rec['learned_advantage_vs_uniform']=float(rec['learned']['beat_R2']-rec['uniform']['beat_R2'])
        rec['learned_advantage_vs_shuffle']=float(rec['learned']['beat_R2']-rec['shuffle_learned']['beat_R2'])
        rec['hand_advantage_vs_uniform']=float(rec['hand']['beat_R2']-rec['uniform']['beat_R2'])
        out[f'{a:.6g},{b:.6g}']=rec
    # Broad signatures across all carrier pairs.
    lr=[];ur=[];sr=[];sep=[];dirgap=[]
    for rec in out.values():
        lr.append(rec['learned']['beat_R2']);ur.append(rec['uniform']['beat_R2']);sr.append(rec['shuffle_learned']['beat_R2'])
        sep.append(rec['learned']['beat_minus_carrier']);dirgap.append(rec['learned']['beat_minus_direct_delta'])
    return dict(
        bodies=len(rows),pairs=len(pairs),pair=out,
        mean_learned_beat_R2=float(np.mean(lr)),mean_uniform_beat_R2=float(np.mean(ur)),
        mean_shuffle_beat_R2=float(np.mean(sr)),
        mean_learned_advantage_uniform=float(np.mean(np.asarray(lr)-np.asarray(ur))),
        mean_learned_advantage_shuffle=float(np.mean(np.asarray(lr)-np.asarray(sr))),
        mean_learned_beat_minus_carrier=float(np.mean(sep)),
        mean_learned_beat_minus_direct_delta=float(np.mean(dirgap)),
    )


def parse_pairs(s):
    out=[]
    for tok in str(s).split(';'):
        if not tok.strip():continue
        a,b=tok.split(',');out.append((float(a),float(b)))
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--functional-arbors',default='FunctionalArbors');ap.add_argument('--seed-start',type=int,default=616);ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--pairs',default='0.10,0.13;0.16,0.19;0.24,0.27;0.36,0.39;0.12,0.16;0.20,0.24;0.32,0.36;0.48,0.52')
    ap.add_argument('--train-omegas',default='0.03,0.04');ap.add_argument('--tau',type=float,default=2);ap.add_argument('--mu',type=float,default=.5);ap.add_argument('--g0',type=float,default=.005);ap.add_argument('--ratio',type=float,default=10)
    ap.add_argument('--steps',type=int,default=50);ap.add_argument('--step-fraction',type=float,default=.10);ap.add_argument('--nshuffle',type=int,default=8);ap.add_argument('--out',default='runs/cross_frequency_beat/dev.json');a=ap.parse_args()
    pairs=parse_pairs(a.pairs);train=[float(x) for x in a.train_omegas.split(',')]
    fa=Path(a.functional_arbors).resolve();sys.path.insert(0,str(fa));from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True;r=body_probe(m,pairs,train,a.tau,a.mu,a.g0,a.ratio,a.steps,a.step_fraction,a.nshuffle);rows.append(r)
        vals=[q['beat_R2'] for q in r['profiles']['learned']];cars=[q['mean_carrier_R2'] for q in r['profiles']['learned']]
        print(f"seed {seed}: learned beat mean={np.mean(vals):.3f} carriers={np.mean(cars):.3f} delta={np.mean(np.asarray(vals)-np.asarray(cars)):+.3f}",flush=True)
    s=summarize(rows,pairs);payload=dict(experiment='two_tone_cross_frequency_beat_dev_v01',development_only=True,
        linear_input_statement='two continuous tones only; no Fourier line at difference frequency before quadratic/power observation',
        beat_definition='B_j=H_j(w1)*conj(H_j(w2)); phase measures carrier phase difference / envelope timing',
        train_omegas=train,pairs=pairs,summary=s,rows=rows)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2),encoding='utf-8');print('\nCROSS-FREQUENCY BEAT DEV');print(json.dumps(s,indent=2))
if __name__=='__main__':main()
