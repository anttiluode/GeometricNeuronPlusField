# Exact damping-gauge reversal — development v0.1

Date: 2026-08-08

## Result in one sentence

> **For the repository's exact uniformly damped discrete wave recurrence, a scalar exponential gauge converts the entire finite-time dynamics into an exactly reversible conservative second-order recurrence; on reused development arbors, the gauged body retraces its internal forward field and reconstructs the complete original transient bond gradient from integrated +/- branch energies to machine precision, without storing the forward field history.**

This is a development mechanism result, not a novelty claim and not yet a hardware result.

---

## The old wall

`TIME_DOMAIN_IN_SITU_WALL.md` correctly established that the direct physical adjoint protocol for the *damped* body has a local-history problem.

The exact bond gradient can be written as a delayed local correlation between forward and adjoint edge differences. During a causal backward/adjoint replay, the needed adjoint sample is present locally, but the matching quantity is the **time-reversed local forward field**.

Simply replaying the original boundary drive backwards through the same damped body does not regenerate that internal history: ordinary time reversal changes the sign of viscous damping.

The spectral branch then compressed the history to K~8/16 coherent bins, and `LOCKIN_OVERLAP_CONFIRM_V01.md` showed that the retained-bin complex overlap itself can be read out with two global interference states. But the origin/alignment of the correctly reversed local forward field remained the hardware wall.

---

# Exact discrete gauge

The repository's linear wave step is

```text
v[n+1]   = a v[n] + dt A0 psi[n] + dt source[n]
psi[n+1] = psi[n] + dt v[n+1]

a = 1 - dt*damping
A0 = stiffness*L - restoring*I
```

Eliminating velocity gives the exact second-order recurrence

```text
psi[n+1] = M psi[n] - a psi[n-1] + dt^2 source[n]

M = (1+a-dt^2*restoring) I + dt^2*stiffness*L.
```

For uniform scalar damping with `a>0`, define

```text
r = sqrt(a)
psi[n] = r^n z[n].
```

Then, with no approximation,

```text
z[n+1] = Q z[n] - z[n-1] + u[n]

Q = M/r
u[n] = dt^2 r^(-(n+1)) source[n].
```

The coefficient of the reversed state is now exactly `-1`.

The interior recurrence is therefore a reversible / determinant-one second-order map. Dissipation has not been undone dynamically; it has been moved into known temporal scale factors on the boundary/source and readout coordinates.

This is the discrete analogue of standard integrating-factor / conformally symplectic factorizations for linearly damped Hamiltonian systems. The transformation itself is not claimed novel.

---

## Objective transformation

For the original energy-style objective

```text
J = coeff * sum_k |psi_s[k]|^2
```

we have

```text
J = coeff * sum_k r^(2k) |z_s[k]|^2.
```

Thus the transformed adjoint/error source at the readout is simply weighted by the known global envelope

```text
q_s[k] = coeff * r^(2k) z_s[k].
```

No spatially distributed correction is required for uniform damping.

---

# Physical retracing identity

Let the transformed forward trajectory be

```text
z[0], z[1], ..., z[T].
```

At the end of the forward pass, the reversible body already contains its terminal state.

If a time-reversal operation reverses the momentum-like quadrature while retaining the position-like quadrature, the same `Q` dynamics and the transformed forcing replayed in reverse generate

```text
w[j] = z[T-j]
```

throughout the body.

In the discrete recurrence this means starting the return with

```text
w[0] = z[T]
w[1] = z[T-1]
```

and then applying the same recurrence with reversed `u`.

This is the internal forward history reconstructed by dynamics rather than stored in a local tape.

---

# Adjoint alignment

The transformed adjoint obeys the same `Q` recurrence.

During the reverse-time experiment one can generate the causal returned field

```text
a[j] = p[T-j+1]
```

so that at the same physical reverse index `j`:

```text
retraced forward branch field  = Delta z[T-j]
returned adjoint branch field  = Delta p[T-j+1].
```

These are exactly the two factors required by the transformed bond gradient.

---

# Full-band integrated interference

For a trainable branch, run two reverse trials with local branch fields

```text
Delta w + Delta a
Delta w - Delta a.
```

A square-law energy accumulator produces

```text
E+ = sum_j |Delta w[j] + Delta a[j]|^2
E- = sum_j |Delta w[j] - Delta a[j]|^2.
```

Then

```text
(E+ - E-) / 4
    = sum_j Re[ conj(Delta w[j]) Delta a[j] ].
```

Up to the known coupling prefactor and edge orientation, that is the **complete time-domain gradient**.

Crucially this sums the full transient bandwidth automatically.

There is no local FFT, no K complex products, and no T-sample field history in this idealized protocol.

---

# Development audit — reused seeds 472–475

`damping_gauge_reversal_probe.py` compares every stage against the existing exact damped forward/adjoint implementation.

Across all four bodies and both target/distractor trajectories:

```text
forward damped <-> gauged reconstruction
    relative L2 ~ 2e-13

original bond gradient <-> gauged reverse-mode gradient
    relative L2 ~ 5e-13
    correlation 1.000000000000

same-Q physical reverse retrace <-> z[::-1]
    relative L2 ~ 1e-14

causal transformed adjoint replay <-> reverse-mode adjoint
    relative L2 ~ 1e-15 or below

+/- integrated branch-energy gradient <-> exact full transient gradient
    relative L2 ~ 4e-13

combined target+distractor gradient
    correlation 1.000000000000
    relative L2 ~ 4e-13
```

So all identities are numerically exact to ordinary floating-point accuracy.

---

## Conservative stability audit

For the transformed recurrence

```text
z[n+1] = Q z[n] - z[n-1],
```

a scalar eigenmode is oscillatory/stable when its `Q` eigenvalue lies in `[-2,2]`.

On development bodies 472–475:

```text
minimum Q eigenvalue   ~1.8096
maximum Q eigenvalue   ~1.99965
```

So the compiled operator lies inside the conservative recurrence's stability band for the present frozen model.

This is important: the exact algebra has not produced an exponentially unstable anti-damped interior.

---

## Boundary dynamic range in the present task

For the frozen model

```text
r = .996694536957...
T = 210
```

The transformed source envelope grows as `r^(-n)`.

At the end of the 210-step window:

```text
r^(-210)  ~ 2.0043
r^(210)   ~ .49893
r^(420)   ~ .24893
```

Thus the present compiler trades distributed damping for only about a 2x boundary-drive dynamic range over the full task window.

For much longer sequences or stronger damping this cost grows exponentially and may become the dominant limitation.

---

# What this changes

The old statement

```text
damping prevents physical retracing,
therefore a local forward-history representation is unavoidable
```

is too strong **when the intended damping is spatially uniform and known**.

A more precise statement is:

```text
if the damped dynamics can be compiled by a scalar similarity transform
into a stable reversible physical recurrence,
then the internal history may be regenerated by an echo instead of stored.
```

The forward-history wall has therefore moved from a mathematical-memory requirement to an **engineering time-reversal requirement**.

---

# Prior-art boundary

This development result immediately intersects several mature ideas.

## Conformal / integrating-factor transformations

Factoring uniform linear damping out of a Hamiltonian/conformally symplectic flow is standard mathematical territory. The scalar gauge itself is not a novelty claim.

## Hamiltonian Echo Backpropagation

López-Pastor & Marquardt (Phys. Rev. X 2023) introduced Hamiltonian Echo Backpropagation for time-reversible Hamiltonian physical systems. It uses a physical time-reversal operation plus an output error perturbation so the echo dynamics carry gradient information.

Pourcel & Ernoult (NeurIPS 2025) introduced Recurrent Hamiltonian Echo Learning, a discrete-time Hamiltonian method equivalent to BPTT for their Hamiltonian recurrent units and requiring three forward-style passes independent of sequence length.

Therefore this repository does **not** own:

```text
memory-free temporal credit via a Hamiltonian echo
constant-pass reversible physical BPTT
```

## Scattering Backpropagation

Dal Cin, Marquardt & Wanjura (2025 preprint) extend physics-based training to driven-dissipative nonlinear optical systems using two scattering experiments and approximate reciprocity, but their training protocol is formulated around a stable **steady-state** scattering response.

Therefore this repository also does not own the broad statement

```text
physics-based backpropagation in a dissipative wave system.
```

---

# Narrow question opened by this result

The actual candidate bridge is:

> **Can a finite-time uniformly damped reciprocal scattering task be compiled exactly into a stable reversible physical core via a known boundary-time gauge, so that Hamiltonian-style echo/adjoint training acquires the original dissipative transient gradient without local history storage?**

That is a much narrower statement than physical backpropagation or Hamiltonian echo learning.

I have not established that this bridge is novel.

---

# Hardware assumptions still unpaid

The machine-precision development identity assumes all of the following:

1. damping in the intended model is uniform enough to admit one scalar gauge;
2. the physical compiled core can realize `Q` with very low residual dissipation;
3. the terminal momentum-like quadrature can be reversed sufficiently accurately;
4. the transformed input can be replayed backward with correct timing;
5. the transformed output/error envelope `r^(2k)` can be applied at the readout;
6. the `+` and `-` reverse trials begin from reproducibly equivalent terminal states;
7. local branch-energy integrators have adequate dynamic range and noise floor;
8. pass-to-pass drift does not destroy the interference difference.

These are now the correct kill tests.

---

# Candidate pass budget

For one forward trajectory, an implementation with no terminal-state snapshot would plausibly require:

```text
1  transformed forward inference pass
1  transformed forward recreation pass for each additional reverse phase state
2  reverse interference trials (+ and -)
```

which is roughly 4 physical traversals if the terminal state must be regenerated for the second interference trial.

If the terminal state can be captured/restored or if a one-run lock-in modulation can be adapted to the echo protocol, the count may fall toward 3.

This pass budget is independent of the number of trainable bonds, the transient length T, and the retained spectral-bin count K.

The cost has moved into time-reversal fidelity, boundary modulation, and local scalar integration.

---

## Next wall

Damage the ideal compiler before any fresh confirmation:

```text
terminal momentum-reversal error
residual physical damping
pass-to-pass operator drift
source-envelope gain/timing error
finite detector/integrator noise
longer T / stronger intended damping
nonuniform intended damping
```

If gradient direction collapses under modest imperfections, the exact identity is only an elegant algebraic compiler.

If it remains useful under realistic errors, then it deserves a held-out confirmation and a direct comparison with Hamiltonian Echo / Scattering Backpropagation assumptions.

## Wall sentence

> **Uniform damping need not force local history storage: in this exact discrete wave model it can be moved to known boundary-time envelopes, leaving a stable reversible core whose physical echo reconstructs the entire internal forward history and whose +/- branch-energy interference returns the full broadband gradient. The remaining question is no longer mathematical memory; it is whether a real low-loss time-reversal core can execute that compiler robustly enough to matter.**
