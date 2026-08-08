"""Decompose what information the material-adjoint learner is actually using.

Exploratory mechanism diagnostic.

The confirmed learner updates quasi-active material with the exact local gradient,
which contains both the local forward field and the readout-launched transpose
field.  Older neuroscience provides two dangerous near-neighbours:

* a returning bAP/Ca signal can implicitly encode soma distance and establish
  synaptic democracy with local homeostatic plasticity;
* activity-dependent local Ca regulation can spontaneously create nonuniform
  intrinsic conductance maps.

This probe therefore asks whether the final material placement can be explained
mostly by a coordinate available from the transpose/return field alone.

For each body, learn the exact map from uniform material exactly as in
hcn_material_learning.py.  Then keep the *exact learned density histogram* and
reassign those same values by several scalar coordinates measured at the
uniform state:

    distance_rank       explicit graph-distance oracle
    return_delay_rank   two-frequency group delay of transpose field
    return_amp_rank     attenuation of transpose-field amplitude
    forward_amp_rank    mean local forward-field amplitude
    random_shuffle      histogram-matched null

Because the histogram is identical, only placement differs.  We also report
Spearman correlations among each coordinate, graph distance, exact learned
material, and the initial exact gradient.

No novelty claim is attached to this exploratory file.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import adjoint_eligibility_probe as ae
import hcn_material_learning as hml
from hcn_impedance_probe import weighted_laplacian_sparse, injection_sites


def rho(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if len(a)<2 or np.std(a)<1e-14 or np.std(b)<1e-14:
        return float('nan')
    return float(spearmanr(a,b).statistic)


def rank_assign(template_values, coordinate, larger_coordinate_gets_larger=True):
    """Assign exact template histogram monotonically by a coordinate."""
    vals=np.sort(np.asarray(template_values,float))
    c=np.asarray(coordinate,float)
    order=np.argsort(c)
    if not larger_coordinate_gets_larger:
        order=order[::-1]
    out=np.empty_like(vals)
    out[order]=vals
    return out


def uniform_coordinates(m,L,uniform,omegas,tau,mu,sites,var_ids):
    ys=[];Xs=[]
    for om in omegas:
        H,X,y,cfac=hml.solve_fields(m,L,uniform,om,tau,mu,sites)
        ys.append(y[var_ids])
        Xs.append(X[var_ids,:])

    yarr=np.asarray(ys)                       # [freq,node]
    # Return attenuation coordinate: distal locations generally see weaker
    # transpose amplitude, so minus log amplitude should increase outward.
    return_amp=-np.log(np.mean(np.abs(yarr),axis=0)+1e-30)

    # Return group-delay coordinate from the two frozen frequencies.  For a
    # local response y(omega) ~ exp(-i omega tau),
    # angle(y(w0)*conj(y(w1))) / (w1-w0) ~= tau.
    if len(omegas)>=2:
        dw=float(omegas[1]-omegas[0])
        z=yarr[0]*np.conj(yarr[1])
        return_delay=np.angle(z)/dw
    else:
        return_delay=-np.angle(yarr[0])/float(omegas[0])

    # Purely forward, no knowledge of readout: average local field magnitude
    # over all injection locations and frozen frequencies.  Negating log makes
    # low-access/attenuated regions large, matching the orientation used for
    # the return-amplitude coordinate.
    fmag=np.mean(np.stack([np.mean(np.abs(X),axis=1) for X in Xs],axis=0),axis=0)
    forward_amp=-np.log(fmag+1e-30)

    return dict(return_amp=return_amp,return_delay=return_delay,forward_amp=forward_amp)


def make_grid(shape,var_ids,vals):
    z=np.zeros(int(np.prod(shape)),float);z[var_ids]=np.asarray(vals,float)
    return z.reshape(shape)


def body_probe(m,omegas,tau,mu,g0,ratio,steps,step_fraction,nshuffle):
    wh,wv=ae.bond_weights(m,m.body);L=weighted_laplacian_sparse(wh,wv)
    body=m.body.astype(bool);h,w=body.shape
    cells=[tuple(map(int,p)) for p in np.argwhere(body)]
    ids=np.asarray([p[0]*w+p[1] for p in cells],int)
    distmap=m.graph_distance_from_soma().astype(float)
    dist=np.asarray([distmap[p] for p in cells],float)
    dmax=max(float(dist.max()),1.0)

    hand=np.zeros_like(distmap,float)
    hand[body]=float(g0)*(1+(float(ratio)-1)*distmap[body]/dmax)
    total=float(hand[body].sum());cap=float(g0)*float(ratio)
    uniform=np.zeros_like(hand);uniform[body]=total/len(ids)
    sites=injection_sites(m)

    # Exact initial objective/gradient and exact confirmed learner.
    obj0,g0_exact,_=hml.coherence_and_gradient(m,L,uniform,omegas,tau,mu,sites,ids,True)
    learned,best,hist=hml.train_density(m,L,uniform,omegas,tau,mu,sites,ids,total,cap,steps,step_fraction)
    learned_vals=learned.ravel()[ids]

    coords=uniform_coordinates(m,L,uniform,omegas,tau,mu,sites,ids)
    coords['distance']=dist

    # Histogram-matched monotone placements.  For return_delay and distance,
    # larger coordinate gets larger material.  return_amp/forward_amp were
    # oriented above so larger also means more attenuated/far-like.
    placements={}
    for name,c in coords.items():
        vv=rank_assign(learned_vals,c,True)
        placements[name+'_rank']=hml.eval_profile(m,L,make_grid(body.shape,ids,vv),omegas,tau,mu,sites)

    # Existing exact learned map and uniform baseline.
    placements['exact']=hml.eval_profile(m,L,learned,omegas,tau,mu,sites)
    placements['uniform']=hml.eval_profile(m,L,uniform,omegas,tau,mu,sites)

    # Histogram-matched global shuffles.
    rng=np.random.default_rng(int(m.cfg.seed)+919_191)
    sh=[]
    for _ in range(int(nshuffle)):
        vv=learned_vals.copy();rng.shuffle(vv)
        sh.append(hml.eval_profile(m,L,make_grid(body.shape,ids,vv),omegas,tau,mu,sites))
    placements['shuffle']=dict(
        coherence2=float(np.mean([q['coherence2'] for q in sh])),
        mean_phase_rms=float(np.mean([q['mean_phase_rms'] for q in sh])),
        mean_median_amp=float(np.mean([q['mean_median_amp'] for q in sh])),
    )

    correlations={
        'learned_vs_distance':rho(learned_vals,dist),
        'learned_vs_return_delay':rho(learned_vals,coords['return_delay']),
        'learned_vs_return_amp':rho(learned_vals,coords['return_amp']),
        'learned_vs_forward_amp':rho(learned_vals,coords['forward_amp']),
        'gradient0_vs_distance':rho(g0_exact,dist),
        'gradient0_vs_return_delay':rho(g0_exact,coords['return_delay']),
        'gradient0_vs_return_amp':rho(g0_exact,coords['return_amp']),
        'gradient0_vs_forward_amp':rho(g0_exact,coords['forward_amp']),
        'return_delay_vs_distance':rho(coords['return_delay'],dist),
        'return_amp_vs_distance':rho(coords['return_amp'],dist),
        'forward_amp_vs_distance':rho(coords['forward_amp'],dist),
    }

    return dict(seed=int(m.cfg.seed),cells=int(body.sum()),sites=len(sites),
                initial_objective=float(obj0),best_objective=float(best),
                correlations=correlations,placements=placements)


def summarize(rows):
    names=['exact','uniform','distance_rank','return_delay_rank','return_amp_rank','forward_amp_rank','shuffle']
    p={}
    for n in names:
        p[n]=dict(
            coherence2=float(np.mean([r['placements'][n]['coherence2'] for r in rows])),
            phase_rms=float(np.mean([r['placements'][n]['mean_phase_rms'] for r in rows])),
            median_amp=float(np.mean([r['placements'][n]['mean_median_amp'] for r in rows])),
            exact_minus=float(np.mean([r['placements']['exact']['coherence2']-r['placements'][n]['coherence2'] for r in rows])),
        )
    keys=list(rows[0]['correlations']) if rows else []
    corr={k:float(np.nanmean([r['correlations'][k] for r in rows])) for k in keys}
    wins={n:int(sum(r['placements']['exact']['coherence2']>r['placements'][n]['coherence2'] for r in rows)) for n in names if n!='exact'}
    return dict(bodies=len(rows),placements=p,mean_correlations=corr,exact_wins_bodies=wins)


def parse_args():
    ap=argparse.ArgumentParser()
    ap.add_argument('--functional-arbors',default='FunctionalArbors')
    ap.add_argument('--seed-start',type=int,default=622)
    ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--omegas',default='0.03,0.04')
    ap.add_argument('--tau',type=float,default=2.0)
    ap.add_argument('--mu',type=float,default=.5)
    ap.add_argument('--g0',type=float,default=.005)
    ap.add_argument('--ratio',type=float,default=10.0)
    ap.add_argument('--steps',type=int,default=50)
    ap.add_argument('--step-fraction',type=float,default=.10)
    ap.add_argument('--nshuffle',type=int,default=12)
    ap.add_argument('--out',default='runs/material_credit_decomposition/dev_622_627.json')
    return ap.parse_args()


def main():
    a=parse_args();fa=Path(a.functional_arbors).resolve()
    if not fa.exists():raise SystemExit(f'FunctionalArbors not found at {fa}')
    sys.path.insert(0,str(fa))
    from v09_causal_eligibility.eligibility_arbor import V09Config,CausalEligibilityArbor
    oms=[float(x) for x in a.omegas.split(',')]
    rows=[]
    for seed in range(a.seed_start,a.seed_start+a.seeds):
        m=CausalEligibilityArbor(V09Config(seed=seed));boot=m.bootstrap()
        if not boot.get('ok'):continue
        m.mature=True
        r=body_probe(m,oms,a.tau,a.mu,a.g0,a.ratio,a.steps,a.step_fraction,a.nshuffle)
        rows.append(r)
        c=r['correlations'];p=r['placements']
        print(f"seed {seed}: exact={p['exact']['coherence2']:.4f} dist={p['distance_rank']['coherence2']:.4f} "
              f"retDelay={p['return_delay_rank']['coherence2']:.4f} retAmp={p['return_amp_rank']['coherence2']:.4f} "
              f"rho learned/dist={c['learned_vs_distance']:+.3f} learned/retDelay={c['learned_vs_return_delay']:+.3f}")
    out=dict(config=vars(a),rows=rows,summary=summarize(rows))
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(json.dumps(out,indent=2))
    print('\nMATERIAL CREDIT DECOMPOSITION')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':main()
