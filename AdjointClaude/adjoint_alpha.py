"""Is the adjoint wrong, or is one cell simply outside the linear regime?

The gradient check says dJ/dk is exact to ~1e-7. Yet the first-order prediction of
a one-cell ADDITION does not track the exact counterfactual. The obvious suspect:
adding a cell takes a bond from k_bath = 2e-4 to k_arbor = 2.5, a factor of 12500.
That is not a small perturbation.

So partially add the cell. Set the new bond to k = KB + alpha*(KA-KB) and sweep
alpha. If the adjoint prediction alpha*dJ/dk*(KA-KB) tracks the exact dJ at small
alpha and decays as alpha grows, the adjoint is fine and the discrete edit is
simply far outside its radius of validity -- which names the fix (density
relaxation / SIMP), rather than leaving a bare null.
"""
import sys, numpy as np
sys.path.insert(0,'.')
import adjoint_eligibility as A
from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor
from v06_ephaptic_growth.ephaptic_arbor import n4

STEPS, LAG, SEEDS, EV = 150, 20, 8, 10
ALPHAS = [0.001, 0.01, 0.05, 0.2, 0.5, 1.0]
res = {al: [] for al in ALPHAS}

for seed in range(SEEDS):
    m = CausalEligibilityArbor(V09Config(seed=seed))
    if not m.bootstrap().get('ok'): continue
    m.mature = True
    body = m.body.astype(bool); n = body.shape[0]; N = n*n
    soma = m.soma[0]*n + m.soma[1]
    srcs = [int(np.ravel_multi_index(m.source_terminal(w),(n,n))) for w in (0,1)]
    ea, eb, k = A.bonds_of(body)
    run, J0, _ = A.run_pair(k, ea, eb, N, srcs, soma, LAG, STEPS)
    grad = A.adjoint_grad(k, ea, eb, N, run, soma, STEPS, 'energy')
    bidx = {(min(i,j),max(i,j)):t for t,(i,j) in enumerate(zip(ea,eb))}

    occ = set(map(tuple, np.argwhere(body)))
    cand = {}
    for cy,cx in occ:
        for r in n4(cy,cx,n,n):
            if r not in occ and 0<r[0]<n-1 and 0<r[1]<n-1:
                nb=[q for q in n4(*r,n,n) if q in occ]
                if len(nb)==1: cand[r]=nb[0]
    rng = np.random.default_rng(seed+77)
    items = list(cand.items())
    if len(items)>EV: items=[items[i] for i in rng.choice(len(items),EV,replace=False)]

    for al in ALPHAS:
        pr, ex = [], []
        for cell, nb in items:
            i=cell[0]*n+cell[1]; j=nb[0]*n+nb[1]
            t=bidx.get((min(i,j),max(i,j)))
            if t is None: continue
            k2=k.copy(); k2[t]=A.KB+al*(A.KA-A.KB)
            _,J1,_=A.run_pair(k2,ea,eb,N,srcs,soma,LAG,STEPS)
            pr.append(grad[t]*al*(A.KA-A.KB)); ex.append(J1-J0)
        pr=np.array(pr); ex=np.array(ex)
        if len(pr)>=4 and pr.std()>0 and ex.std()>0:
            res[al].append((float(np.corrcoef(pr,ex)[0,1]),
                            float(np.median(np.abs(pr-ex)/(np.abs(ex)+1e-30)))))
    print('seed',seed,'done',flush=True)

print('\n'+'='*78)
print('DOES THE ADJOINT PREDICTION HOLD AS THE EDIT GETS SMALLER?')
print('='*78)
print(f'{"alpha":>8} {"k of new bond":>14} {"mean corr":>11} {"positive":>10} {"median rel err":>16}')
print('-'*78)
for al in ALPHAS:
    v=np.array([x[0] for x in res[al]]); e=np.array([x[1] for x in res[al]])
    kk=A.KB+al*(A.KA-A.KB)
    print(f'{al:8.3f} {kk:14.4f} {v.mean():+11.3f} {int((v>0).sum()):>6}/{len(v):<3} {np.median(e):16.3e}')
print('-'*78)
print('alpha=1.0 is the real one-cell addition. If corr is high at small alpha and')
print('collapses toward alpha=1, the adjoint is exact and the discrete edit is simply')
print('outside its radius of validity.')
