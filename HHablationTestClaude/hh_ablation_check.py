"""Is "h instantaneous -> 0/24 bodies spike" a mechanism result, or a dead arm?

The AIS probe made each gate instantaneous (gate = its voltage-dependent
equilibrium value at every step) while preserving the equilibrium curve, and
found:

    m instantaneous   24/24 bodies still fire
    h instantaneous    0/24
    n instantaneous    0/24

and concluded that h(t) and n(t) history is "dynamically necessary".

That conclusion may not be available from this ablation, for the same reason the
v0.8 trophic arm and the v0.9 timing arm were not valid nulls: the arm may be
dead by construction rather than dead because the memory mattered. Concretely,
in Hodgkin-Huxley the spike exists only because tau_m << tau_h, tau_n. Set
h = h_inf(V) and sodium inactivation cancels the upstroke as fast as it is
generated. Set n = n_inf(V) and potassium shunts it. Neither is a statement
about frequency selection; both are the textbook reason the equations spike.

The test that distinguishes the two readings: sweep the drive amplitude over
orders of magnitude and ask whether ANY operating point restores spiking. If the
ablated membrane is silent everywhere, the ablation removed excitability, not a
memory -- and the 0/24 result carries no information about what h and n do to
the frequency profile.

Standard HH (Hodgkin & Huxley 1952) at 6.3 C. No parameter is tuned to the answer.
"""
from __future__ import annotations
import numpy as np

C_M, G_NA, G_K, G_L = 1.0, 120.0, 36.0, 0.3
E_NA, E_K, E_L = 50.0, -77.0, -54.387
V0 = -65.0


def rates(V):
    am = 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10)) if abs(V + 40) > 1e-7 else 1.0
    bm = 4.0 * np.exp(-(V + 65) / 18)
    ah = 0.07 * np.exp(-(V + 65) / 20)
    bh = 1.0 / (1 + np.exp(-(V + 35) / 10))
    an = 0.01 * (V + 55) / (1 - np.exp(-(V + 55) / 10)) if abs(V + 55) > 1e-7 else 0.1
    bn = 0.125 * np.exp(-(V + 65) / 80)
    return am, bm, ah, bh, an, bn


def inf(a, b):
    return a / (a + b)


def run(I, dt=0.005, T=200.0, instant=()):
    n_steps = int(T / dt)
    V = V0
    am, bm, ah, bh, an, bn = rates(V)
    m, h, n = inf(am, bm), inf(ah, bh), inf(an, bn)
    spikes, above = 0, False
    for _ in range(n_steps):
        am, bm, ah, bh, an, bn = rates(V)
        m = inf(am, bm) if 'm' in instant else m + dt * (am * (1 - m) - bm * m)
        h = inf(ah, bh) if 'h' in instant else h + dt * (ah * (1 - h) - bh * h)
        n = inf(an, bn) if 'n' in instant else n + dt * (an * (1 - n) - bn * n)
        I_ion = G_NA * m ** 3 * h * (V - E_NA) + G_K * n ** 4 * (V - E_K) + G_L * (V - E_L)
        V = V + dt * (I - I_ion) / C_M
        if not np.isfinite(V):
            return -1
        if V > 0 and not above:
            spikes += 1; above = True
        elif V < -20:
            above = False
    return spikes


if __name__ == '__main__':
    arms = [((), 'intact'), (('m',), 'm instantaneous'), (('h',), 'h instantaneous'),
            (('n',), 'n instantaneous'), (('h', 'n'), 'h and n instantaneous'),
            (('m', 'h', 'n'), 'all instantaneous')]
    drives = [0.5, 1, 2, 4, 6, 8, 10, 15, 20, 30, 50, 80, 120, 200, 400, 800]
    print('Spike count in 200 ms of constant drive, standard HH, drive swept 0.5 -> 800 uA/cm^2')
    print('=' * 96)
    print(f'{"arm":24s} ' + ' '.join(f'{d:>4g}' for d in drives) + '   ANY?')
    print('-' * 96)
    for instant, label in arms:
        counts = [run(I, instant=instant) for I in drives]
        ever = any(c > 0 for c in counts)
        print(f'{label:24s} ' + ' '.join(f'{c:>4d}' for c in counts) +
              f'   {"fires" if ever else "SILENT EVERYWHERE"}')
    print('-' * 96)
    print('If an ablated arm is silent at every drive across three orders of magnitude,')
    print('it lost excitability by construction. Its 0/24 is a dead arm, not a mechanism.')
    print('The interpretable version re-tunes conductance to a MATCHED event rate first,')
    print('then compares frequency profiles -- the same rate-matching rule the v0.2')
    print('encoder controls already use, applied to the ablations too.')
