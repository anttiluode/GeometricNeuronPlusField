"""Evaluate MATERIAL_READOUT_RECENTER_CONFIRM_PREREG_V01."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np


def evaluate(p):
    rows=p['rows'];n=len(rows);frac=.75
    mg=np.asarray([float(r['moved']['gain']) for r in rows])
    sg=np.asarray([float(r['soma']['gain']) for r in rows])
    nr=np.asarray([float(r['moved']['stats']['spearman_to_readout']) for r in rows])
    old=np.asarray([float(r['moved']['stats']['spearman_to_soma']) for r in rows])
    sr=np.asarray([float(r['soma']['stats']['spearman_to_readout']) for r in rows])
    delta=nr-old
    crit={}
    crit['R0_moved_learning']={'pass_':bool(mg.mean()>.05 and np.mean(mg>0)>=frac),'mean_gain':float(mg.mean()),'positive':int(np.sum(mg>0)),'fraction_positive':float(np.mean(mg>0))}
    crit['R1_new_readout_coordinate']={'pass_':bool(nr.mean()>.40 and np.mean(nr>0)>=frac),'mean_rho_new':float(nr.mean()),'positive':int(np.sum(nr>0)),'fraction_positive':float(np.mean(nr>0))}
    crit['R2_new_beats_old_coordinate']={'pass_':bool(delta.mean()>.40 and np.mean(delta>0)>=frac),'mean_delta_rho':float(delta.mean()),'positive_delta':int(np.sum(delta>0)),'fraction_positive':float(np.mean(delta>0))}
    crit['R3_not_soma_centered']={'pass_':bool(old.mean()<.15),'mean_rho_old_soma':float(old.mean())}
    crit['R4_soma_reproduces']={'pass_':bool(sg.mean()>.05 and sr.mean()>.40),'mean_soma_gain':float(sg.mean()),'mean_soma_rho':float(sr.mean())}
    ratio=float(mg.mean()/(sg.mean()+1e-30))
    crit['R5_gain_same_order']={'pass_':bool(ratio>.45),'moved_over_soma_gain':ratio}
    passed=sum(int(x['pass_']) for x in crit.values())
    return dict(prereg='MATERIAL_READOUT_RECENTER_CONFIRM_PREREG_V01.md',bodies=n,criteria=crit,passed=passed,total=len(crit),all_pass=bool(passed==len(crit)),pooled=dict(mean_moved_gain=float(mg.mean()),mean_soma_gain=float(sg.mean()),mean_rho_new=float(nr.mean()),mean_rho_old_soma=float(old.mean()),mean_delta_rho=float(delta.mean()),mean_soma_rho=float(sr.mean())))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--receipt',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();p=json.loads(Path(a.receipt).read_text(encoding='utf-8'));s=evaluate(p);Path(a.out).parent.mkdir(parents=True,exist_ok=True);Path(a.out).write_text(json.dumps(s,indent=2),encoding='utf-8');print(json.dumps(s,indent=2));
    if not s['all_pass']:raise SystemExit('Moved-readout confirmation failed one or more frozen criteria')
if __name__=='__main__':main()
