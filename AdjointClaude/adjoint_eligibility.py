"""ADJOINT ELIGIBILITY — can one backward pass give the local credit signal?

STRUCTURAL_INTERFERENCE_V01 ends by naming the right next question: can the effect
of a local structural event on the soma objective be predicted from a LOCAL
sensitivity signal, without building the event and recomputing the whole arbor?

For a linear reciprocal wave operator that is the classical adjoint-state method,
and it is exactly what the whole v0.7-v0.9 credit arc was missing. This script
implements it against the real FunctionalArbor physics and checks it three ways.

THE SETUP.  Their update, minus the negligible saturation term (already certified
negligible by the r = 0.99999984 superposition test):

    v^{n+1} = v^n + dt*( c*L(k) psi^n - gamma v^n - rho psi^n + s^n )
    psi^{n+1} = psi^n + dt*v^{n+1}

with L(k) psi = sum_e k_e (psi_b - psi_a) on each bond, i.e. L = -sum_e k_e (d_a-d_b)(d_a-d_b)^T.
c, gamma, rho are real, so the real and imaginary parts of psi evolve independently
under the same real operator — everything below is real arithmetic on two fields.

THE OBJECTIVE.  A smooth energy contrast at the soma:

    J = sum_n |psi_T(soma,n)|^2  -  sum_n |psi_D(soma,n)|^2

THE ADJOINT.  With z = (psi, v), z^{n+1} = M z^n + f^n and

    M = [[ I + dt*A , dt*B ],
         [ A        , B    ]],    A = dt*(c L - rho I),  B = (1-gamma*dt) I

run p^n = M^T p^{n+1} + dJ/dz^n backward. Then the sensitivity of J to the
conductance of ANY bond is a purely LOCAL product across that bond:

    dJ/dk_e = -dt^2 c * sum_n (psi_a^n - psi_b^n) * (P_a^{n+1} - P_b^{n+1})
    P = dt*p_psi + p_v

One forward pass and one backward pass give this for EVERY bond at once. It is
local, it is a coincidence of a forward field with a retrograde field, and it is
the eligibility trace the tag experiments were trying to guess.

THREE CHECKS, in increasing strength:
  1. GRADIENT CHECK      adjoint dJ/dk_e vs central finite differences
  2. PREDICTION          adjoint prediction of dJ for adding one cell vs the exact
                         before/after counterfactual of the same objective
  3. TASK RELEVANCE      does the adjoint prediction track the change in the
                         historical PEAK contrast the lineage actually optimises
"""
from __future__ import annotations
import sys, json, math, argparse
from collections import deque
import numpy as np

sys.path.insert(0, '.')
from v09_causal_eligibility.eligibility_arbor import V09Config, CausalEligibilityArbor
from v06_ephaptic_growth.ephaptic_arbor import n4

CFG = V09Config()
DT, C, GAM, RHO = CFG.dt, CFG.stiffness, CFG.damping, CFG.restoring
KA, KB = CFG.k_arbor, CFG.k_mature_bath
OMEGA, PF, AMP = CFG.carrier_omega, CFG.pulse_frames, CFG.source_amp


# ───────────────────────── the reduced, exactly-differentiable model ──────────
def bonds_of(body):
    """All 4-neighbour bonds on the grid, with conductance k_arbor inside the body
    and k_mature_bath elsewhere — exactly the mature bond_fields rule."""
    n = body.shape[0]
    ea, eb, k = [], [], []
    for y in range(n):
        for x in range(n):
            i = y * n + x
            for dy, dx in ((0, 1), (1, 0)):
                yy, xx = y + dy, x + dx
                if yy >= n or xx >= n: continue
                j = yy * n + xx
                ea.append(i); eb.append(j)
                k.append(KA if (body[y, x] and body[yy, xx]) else KB)
    return np.array(ea), np.array(eb), np.array(k, float)


def pulse(which, srcs, N, steps):
    """Complex source train, split into the two independent real fields."""
    sr = np.zeros((steps, N)); si = np.zeros((steps, N))
    p = srcs[which]
    for q in range(PF):
        env = math.sin(math.pi * (q + 1) / (PF + 1)) ** 2
        sr[q, p] += AMP * env * math.cos(OMEGA * q)
        si[q, p] += AMP * env * math.sin(OMEGA * q)
    return sr, si


def lap(k, ea, eb, u, N):
    f = k * (u[eb] - u[ea])
    out = np.zeros(N)
    np.add.at(out, ea, f); np.add.at(out, eb, -f)
    return out


def forward(k, ea, eb, N, src, steps, store=True):
    psi = np.zeros(N); vel = np.zeros(N)
    P = np.zeros((steps + 1, N)) if store else None
    if store: P[0] = psi
    for n in range(steps):
        vel = vel + DT * (C * lap(k, ea, eb, psi, N) - GAM * vel - RHO * psi + src[n])
        psi = psi + DT * vel
        if store: P[n + 1] = psi
    return P if store else psi


def run_pair(k, ea, eb, N, srcs, soma, lag, steps):
    """Both stimulus orders. Returns stored psi_r/psi_i histories and the objectives."""
    out = {}
    for name, order in (('T', (0, 1)), ('D', (1, 0))):
        ar, ai = pulse(order[0], srcs, N, steps)
        br, bi = pulse(order[1], srcs, N, steps)
        sr = ar.copy(); si = ai.copy()
        if lag < steps:
            sr[lag:] += br[:steps - lag]; si[lag:] += bi[:steps - lag]
        Pr = forward(k, ea, eb, N, sr, steps)
        Pi = forward(k, ea, eb, N, si, steps)
        E = float((Pr[1:, soma] ** 2 + Pi[1:, soma] ** 2).sum())
        pk = float((Pr[1:, soma] ** 2 + Pi[1:, soma] ** 2).max())
        out[name] = dict(Pr=Pr, Pi=Pi, sr=sr, si=si, E=E, pk=pk)
    J = out['T']['E'] - out['D']['E']
    Cpk = (out['T']['pk'] - out['D']['pk']) / (out['T']['pk'] + out['D']['pk'] + 1e-30)
    return out, J, Cpk


def adjoint_grad(k, ea, eb, N, run, soma, steps, objective='energy'):
    """One backward pass per field per order -> dJ/dk_e for every bond at once.

    objective='energy'  J = sum_n |psi_soma|^2 (T) - same (D)      smooth, used for
                        the gradient check.
    objective='peak'    J = (pkT-pkD)/(pkT+pkD), the historical lineage objective.
                        dJ/dpsi is then supported only at the argmax step.
    """
    grad = np.zeros(len(ea))
    if objective == 'peak':
        pT, pD = run['T']['pk'], run['D']['pk']
        den = (pT + pD) ** 2 + 1e-30
        wT, wD = 2 * pD / den, -2 * pT / den
        nT = int(np.argmax(run['T']['Pr'][1:, soma] ** 2 + run['T']['Pi'][1:, soma] ** 2)) + 1
        nD = int(np.argmax(run['D']['Pr'][1:, soma] ** 2 + run['D']['Pi'][1:, soma] ** 2)) + 1
        peak = {'T': (wT, nT), 'D': (wD, nD)}
    for name, sgn in (('T', +1.0), ('D', -1.0)):
        for fld in ('Pr', 'Pi'):
            P = run[name][fld]
            pp = np.zeros(N); pv = np.zeros(N)      # p^{steps}
            if objective == 'energy':
                pp[soma] += 2.0 * P[steps][soma]
            else:
                w, nstar = peak[name]
                if nstar == steps: pp[soma] += w * 2.0 * P[steps][soma]
            acc = np.zeros(len(ea))
            for n in range(steps - 1, -1, -1):
                Pn = DT * pp + pv                   # the p^{n+1} combination
                # d(L psi)/dk_e = (psi_b-psi_a) at a and (psi_a-psi_b) at b, so the
                # contraction with Pn is -(psi_a-psi_b)(Pn_a-Pn_b)
                acc -= (P[n][ea] - P[n][eb]) * (Pn[ea] - Pn[eb])
                Lp = lap(k, ea, eb, Pn, N)          # p^n = M^T p^{n+1} + dJ/dz^n
                new_pp = pp + (DT * C) * Lp - (DT * RHO) * Pn
                new_pv = (1 - GAM * DT) * Pn      # p_v^n = B^T (dt p_psi + p_v)^{n+1}
                if objective == 'energy':
                    if n >= 1: new_pp[soma] += 2.0 * P[n][soma]
                else:
                    w, nstar = peak[name]
                    if n == nstar: new_pp[soma] += w * 2.0 * P[n][soma]
                pp, pv = new_pp, new_pv
            grad += (sgn if objective == 'energy' else 1.0) * (DT * C) * acc
    return grad


# ───────────────────────── graph helpers ──────────────────────────────────────
def connected(body):
    n = body.shape[0]; occ = np.argwhere(body)
    if not len(occ): return False
    s = tuple(occ[0]); seen = {s}; q = deque([s])
    while q:
        p = q.popleft()
        for r in n4(*p, n, n):
            if body[r] and r not in seen: seen.add(r); q.append(r)
    return len(seen) == len(occ)


def sign_p(d):
    d = np.asarray(d, float); d = d[np.abs(d) > 1e-14]; m = len(d)
    if m == 0: return 1.0, 0, 0
    pos = int((d > 0).sum()); kk = min(pos, m - pos)
    return min(1.0, 2 * sum(math.comb(m, i) for i in range(kk + 1)) / 2 ** m), pos, m


# ───────────────────────── main ───────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=12)
    ap.add_argument('--steps', type=int, default=150)
    ap.add_argument('--lag', type=int, default=20)
    ap.add_argument('--events', type=int, default=10)
    ap.add_argument('--out', type=str, default='adjoint_eligibility.json')
    a = ap.parse_args()

    rows = []; gradcheck = []
    for seed in range(a.seeds):
        m = CausalEligibilityArbor(V09Config(seed=seed))
        if not m.bootstrap().get('ok'): continue
        m.mature = True
        body = m.body.astype(bool); n = body.shape[0]; N = n * n
        soma = m.soma[0] * n + m.soma[1]
        srcs = [int(np.ravel_multi_index(m.source_terminal(w), (n, n))) for w in (0, 1)]
        ea, eb, k = bonds_of(body)
        run, J0, Cpk0 = run_pair(k, ea, eb, N, srcs, soma, a.lag, a.steps)
        grad = adjoint_grad(k, ea, eb, N, run, soma, a.steps, 'energy')
        gradC = adjoint_grad(k, ea, eb, N, run, soma, a.steps, 'peak')

        # ---- CHECK 1: gradient check against central finite differences
        rng = np.random.default_rng(seed + 5150)
        big = np.argsort(-np.abs(grad))[:60]
        pick = rng.choice(big, size=6, replace=False)
        for e in pick:
            h = max(1e-4, abs(k[e]) * 1e-3)
            kp = k.copy(); kp[e] += h
            km = k.copy(); km[e] -= h
            _, Jp, _ = run_pair(kp, ea, eb, N, srcs, soma, a.lag, a.steps)
            _, Jm, _ = run_pair(km, ea, eb, N, srcs, soma, a.lag, a.steps)
            fd = (Jp - Jm) / (2 * h)
            gradcheck.append(dict(seed=seed, e=int(e), adj=float(grad[e]), fd=float(fd)))

        # ---- CHECKS 2 & 3: one-cell additions, adjoint prediction vs exact
        occ = set(map(tuple, np.argwhere(body)))
        cand = []
        for cy, cx in occ:
            for r in n4(cy, cx, n, n):
                if r not in occ and 0 < r[0] < n - 1 and 0 < r[1] < n - 1:
                    nb = [q for q in n4(*r, n, n) if q in occ]
                    if len(nb) == 1: cand.append((r, nb[0]))
        cand = list({c[0]: c for c in cand}.values())
        if len(cand) > a.events:
            cand = [cand[i] for i in rng.choice(len(cand), a.events, replace=False)]

        bidx = {}
        for t, (i, j) in enumerate(zip(ea, eb)): bidx[(min(i, j), max(i, j))] = t

        evs = []
        for (cell, nb) in cand:
            b2 = body.copy(); b2[cell] = True
            if not connected(b2): continue
            i = cell[0] * n + cell[1]; j = nb[0] * n + nb[1]
            t = bidx.get((min(i, j), max(i, j)))
            if t is None: continue
            # first-order adjoint prediction: this bond goes from bath to arbor
            pred = float(grad[t] * (KA - KB))
            predC = float(gradC[t] * (KA - KB))
            k2 = k.copy()
            for r in n4(*cell, n, n):
                if b2[r]:
                    ii = cell[0] * n + cell[1]; jj = r[0] * n + r[1]
                    tt = bidx.get((min(ii, jj), max(ii, jj)))
                    if tt is not None: k2[tt] = KA
            _, J1, Cpk1 = run_pair(k2, ea, eb, N, srcs, soma, a.lag, a.steps)
            evs.append(dict(pred=pred, predC=predC, dJ=float(J1 - J0), dC=float(Cpk1 - Cpk0)))
        if len(evs) >= 4:
            P = np.array([e['pred'] for e in evs]); DJ = np.array([e['dJ'] for e in evs])
            DC = np.array([e['dC'] for e in evs])
            r_pred = float(np.corrcoef(P, DJ)[0, 1]) if P.std() > 0 and DJ.std() > 0 else np.nan
            PC = np.array([e['predC'] for e in evs])
            r_task = float(np.corrcoef(PC, DC)[0, 1]) if PC.std() > 0 and DC.std() > 0 else np.nan
            sg = float((np.sign(P) == np.sign(DJ)).mean())
            sgc = float((np.sign(PC) == np.sign(DC)).mean())
            rows.append(dict(seed=seed, n_events=len(evs), r_pred=r_pred, r_task=r_task,
                             sign_dJ=sg, sign_dC=sgc, events=evs))
            print(f'seed {seed:2d}  {len(evs):2d} events  corr(adjoint, exact dJ) {r_pred:+.3f}  '
                  f'sign {sg:.2f}  |  corr(adjoint, d peak-contrast) {r_task:+.3f}  sign {sgc:.2f}',
                  flush=True)

    # ── report ──
    ad = np.array([g['adj'] for g in gradcheck]); fd = np.array([g['fd'] for g in gradcheck])
    rel = np.abs(ad - fd) / (np.abs(fd) + 1e-12)
    print('\n' + '=' * 86)
    print(f'ADJOINT ELIGIBILITY — {len(rows)} bodies, {a.steps} steps, lag {a.lag}')
    print('=' * 86)
    print('CHECK 1  gradient check: adjoint dJ/dk vs central finite differences')
    print(f'   {len(gradcheck)} bonds   corr {np.corrcoef(ad,fd)[0,1]:.8f}   '
          f'median relative error {np.median(rel):.2e}   max {rel.max():.2e}')
    print(f'   {"IMPLEMENTATION VERIFIED" if np.median(rel)<1e-3 else "*** ADJOINT IS WRONG ***"}')

    rp = np.array([r['r_pred'] for r in rows]); rt = np.array([r['r_task'] for r in rows])
    sj = np.array([r['sign_dJ'] for r in rows]); sc = np.array([r['sign_dC'] for r in rows])
    p1, n1, m1 = sign_p(rp); p2, n2, m2 = sign_p(rt)
    print('\nCHECK 2  one backward pass predicts the exact one-cell counterfactual')
    print(f'   mean corr(adjoint prediction, exact dJ)   {rp.mean():+.3f}   '
          f'median {np.median(rp):+.3f}   positive {n1}/{m1}   sign p {p1:.5f}')
    print(f'   mean sign agreement                        {sj.mean():.3f}')
    print("\nCHECK 3  a peak-objective adjoint predicts the lineage's own contrast")
    print(f'   mean corr(adjoint prediction, d peak C)   {rt.mean():+.3f}   '
          f'median {np.median(rt):+.3f}   positive {n2}/{m2}   sign p {p2:.5f}')
    print(f'   mean sign agreement                        {sc.mean():.3f}')

    json.dump(dict(rows=rows, gradcheck=gradcheck), open(a.out, 'w'))
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
