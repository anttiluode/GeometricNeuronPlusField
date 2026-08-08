"""WHY is the soma the best place? Test the trivial explanation first.

GeometricNeuronPlusField reports that per-cell temporal-order selectivity is much
higher at the soma (|C| ~ 0.228) than at an average body cell (~0.054), top
quartile in 24/24 bodies, and above the median local aperture at every radius.
It caveats correctly that the soma is a designated convergence root.

But there is a sharper, cheaper explanation that should be ruled out before the
result is read as being about the soma at all:

    C = (max_t y_target - max_t y_distractor) / (sum)

is a RATIO. Wherever one source dominates the response, the ratio saturates and
becomes insensitive to order. Selectivity should therefore be largest where the
two sources arrive with COMPARABLE AMPLITUDE -- the balance point. That is a
coincidence-detector statement, and it is exactly the Jeffress arrangement the
task was designed around.

So: measure at every body cell
    b(x)  = min(pA, pB) / max(pA, pB)      amplitude balance, single sources
    dl(x) = |dist(x,A) - dist(x,B)|        path-length difference
    |C|(x)                                 measured order selectivity

and ask which predicts |C|. If balance does, then "the soma is informative" is
"the soma is the balance point" -- true, useful, and a consequence of how the
body was grown, not an independent finding about somata.
"""
from __future__ import annotations
import sys, json, math, argparse
from collections import deque
import numpy as np

sys.path.insert(0, '.')
from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor
from v06_ephaptic_growth.ephaptic_arbor import n4


def gdist(body, start):
    b = np.asarray(body, bool); d = np.full(b.shape, -1, int)
    if not b[start]: return d
    d[start] = 0; q = deque([tuple(start)])
    while q:
        p = q.popleft()
        for r in n4(*p, *b.shape):
            if b[r] and d[r] < 0:
                d[r] = d[p] + 1; q.append(r)
    return d


def peak_map(m, lag=None, target=None, which=None, steps=None):
    """Per-cell max over time of |psi|^2 for a given drive."""
    steps = int(steps or m.cfg.probe_steps); m.reset_fast(True)
    acc = np.zeros(m.body.shape)
    if which is not None:
        for t in range(steps):
            m.advance(m.pulse_source(which, t, False), False, True, 'none')
            acc = np.maximum(acc, np.abs(m.psi) ** 2)
    else:
        first, second = (0, 1) if target else (1, 0)
        for t in range(steps):
            a = m.pulse_source(first, t, False); b = m.pulse_source(second, t - lag, False)
            src = b if isinstance(a, float) else (a if isinstance(b, float) else a + b)
            m.advance(src, False, True, 'none')
            acc = np.maximum(acc, np.abs(m.psi) ** 2)
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=12); ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--out', type=str, default='soma_why.json')
    a = ap.parse_args()

    rows = []
    for seed in range(a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        boot = m.bootstrap()
        if not boot.get('ok'): continue
        m.mature = True
        body = m.body.astype(bool)
        pA = peak_map(m, which=0); pB = peak_map(m, which=1)
        T = peak_map(m, lag=a.lag, target=True); D = peak_map(m, lag=a.lag, target=False)
        C = (T - D) / (T + D + 1e-30)
        srcA, srcB = m.sources[0], m.sources[1]
        dA = gdist(m.body, srcA); dB = gdist(m.body, srcB)

        cells = np.argwhere(body)
        idx = tuple(cells.T)
        absC = np.abs(C[idx])
        bal = np.minimum(pA[idx], pB[idx]) / (np.maximum(pA[idx], pB[idx]) + 1e-30)
        dl = np.abs(dA[idx] - dB[idx]).astype(float)
        energy = (pA[idx] + pB[idx])
        dsom = gdist(m.body, m.soma)[idx].astype(float)
        si = int(np.flatnonzero((cells[:, 0] == m.soma[0]) & (cells[:, 1] == m.soma[1]))[0])

        pct = float((absC < absC[si]).mean())
        rows.append(dict(seed=seed, n=len(cells), soma_absC=float(absC[si]),
                         mean_absC=float(absC.mean()), soma_pct=pct,
                         soma_balance=float(bal[si]), mean_balance=float(bal.mean()),
                         balance_pct=float((bal < bal[si]).mean()),
                         r_bal=float(np.corrcoef(bal, absC)[0, 1]),
                         r_dl=float(np.corrcoef(dl, absC)[0, 1]),
                         r_energy=float(np.corrcoef(np.log(energy + 1e-30), absC)[0, 1]),
                         r_dsoma=float(np.corrcoef(dsom, absC)[0, 1])))
        r = rows[-1]
        print(f"seed {seed:2d}  soma |C| {r['soma_absC']:.3f} vs mean {r['mean_absC']:.3f} "
              f"(pctile {pct:.3f})  |  soma balance {r['soma_balance']:.3f} vs mean "
              f"{r['mean_balance']:.3f} (pctile {r['balance_pct']:.3f})  |  "
              f"corr(|C|, balance) {r['r_bal']:+.3f}", flush=True)

    n = len(rows)
    print('\n' + '=' * 88)
    print(f'WHY IS THE SOMA THE INFORMATIVE PLACE? — {n} frozen bodies, every body cell measured')
    print('=' * 88)
    sc = np.array([r['soma_absC'] for r in rows]); mc = np.array([r['mean_absC'] for r in rows])
    print(f'reproduces the reported effect:  soma |C| {sc.mean():.4f}  vs mean body cell {mc.mean():.4f}')
    print(f'   soma selectivity percentile   {np.mean([r["soma_pct"] for r in rows]):.3f} '
          f'(top quartile in {sum(r["soma_pct"]>=.75 for r in rows)}/{n} bodies)')
    print(f'\nAND the soma is also the AMPLITUDE BALANCE POINT:')
    print(f'   soma balance {np.mean([r["soma_balance"] for r in rows]):.3f} vs mean cell '
          f'{np.mean([r["mean_balance"] for r in rows]):.3f}, percentile '
          f'{np.mean([r["balance_pct"] for r in rows]):.3f}')
    print(f'\nwhich predictor explains per-cell selectivity?  (mean corr with |C| across bodies)')
    for k, lab in [('r_bal', 'amplitude balance min/max'), ('r_dl', '|path(A) - path(B)|'),
                   ('r_energy', 'log total energy'), ('r_dsoma', 'graph distance from soma')]:
        v = np.array([r[k] for r in rows])
        print(f'   {lab:28s} {v.mean():+.3f}   (range {v.min():+.3f} .. {v.max():+.3f})')
    json.dump(dict(seeds=n, rows=rows), open(a.out, 'w'))
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
