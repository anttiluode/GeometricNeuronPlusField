"""Frozen evaluator for MATERIAL_WITHIN_SHELL_CONFIRM_PREREG_V01.md."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--receipt',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    d=json.loads(Path(a.receipt).read_text());rows=d.get('rows',[]);n=len(rows)
    frac=lambda k: (sum(bool(k(r)) for r in rows)/n if n else 0.0)
    mean=lambda f: (float(np.mean([float(f(r)) for r in rows])) if rows else float('nan'))

    radial_gain=mean(lambda r:r['radial_minus_uniform'])
    release_gain=mean(lambda r:r['release_minus_radial'])
    shuffle_gain=mean(lambda r:r['release_minus_shuffle'])
    phase_gain=mean(lambda r:r['phase_rms_gain_vs_radial'])
    amp=[float(r['amp_ratio_vs_radial']) for r in rows]
    med_amp=float(np.median(amp)) if amp else float('nan')

    freq={}
    for om in ('0.03','0.04'):
        dr=[float(r['per_frequency'][om]['learned_minus_radial_R2']) for r in rows]
        ds=[float(r['per_frequency'][om]['learned_minus_shuffle_R2']) for r in rows]
        freq[om]=dict(mean_vs_radial=float(np.mean(dr)) if dr else float('nan'),
                      mean_vs_shuffle=float(np.mean(ds)) if ds else float('nan'),
                      positive_vs_radial=int(sum(x>0 for x in dr)),
                      positive_vs_shuffle=int(sum(x>0 for x in ds)))

    criteria={
      'W0_usable_population':dict(pass_=n>=10,bodies=n),
      'W1_radial_useful':dict(pass_=(n>=10 and radial_gain>.06 and frac(lambda r:r['radial_minus_uniform']>0)>=.75),
                             mean_gain=radial_gain,positive_fraction=frac(lambda r:r['radial_minus_uniform']>0)),
      'W2_release_beats_radial':dict(pass_=(n>=10 and release_gain>.002 and frac(lambda r:r['release_minus_radial']>0)>=.75),
                                     mean_gain=release_gain,positive_fraction=frac(lambda r:r['release_minus_radial']>0)),
      'W3_exact_placement_matters':dict(pass_=(n>=10 and shuffle_gain>.0025 and frac(lambda r:r['release_minus_shuffle']>0)>=.75),
                                        mean_gain=shuffle_gain,positive_fraction=frac(lambda r:r['release_minus_shuffle']>0)),
      'W4_phase_rms_agrees':dict(pass_=(n>=10 and phase_gain>.005 and frac(lambda r:r['phase_rms_gain_vs_radial']>0)>=2/3),
                                 mean_gain=phase_gain,positive_fraction=frac(lambda r:r['phase_rms_gain_vs_radial']>0)),
      'W5_amplitude_sane':dict(pass_=(n>=10 and .75<=med_amp<=1.25 and frac(lambda r:.60<=r['amp_ratio_vs_radial']<=1.40)>=.75),
                               median_ratio=med_amp,in_range_fraction=frac(lambda r:.60<=r['amp_ratio_vs_radial']<=1.40)),
      'W6_both_frequencies':dict(pass_=(n>=10 and all(freq[o]['mean_vs_radial']>0 and freq[o]['mean_vs_shuffle']>0 for o in ('0.03','0.04'))),frequency=freq),
      'W7_distance_dominates':dict(pass_=(n>=10 and radial_gain>0 and release_gain/(radial_gain+1e-30)<.15),
                                   radial_gain=radial_gain,release_gain=release_gain,release_over_radial=release_gain/(radial_gain+1e-30)),
    }
    passed=sum(bool(x['pass_']) for x in criteria.values());total=len(criteria)
    verdict=dict(prereg='MATERIAL_WITHIN_SHELL_CONFIRM_PREREG_V01.md',bodies=n,criteria=criteria,
                 passed=passed,total=total,all_pass=passed==total)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True);Path(a.out).write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    if not verdict['all_pass']:
        raise SystemExit('Held-out within-shell confirmation failed one or more frozen criteria')

if __name__=='__main__':main()
