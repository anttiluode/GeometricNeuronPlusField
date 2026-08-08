# Graded bond response v0.1 — useful structural couplings are often analog, not binary

`ADJOINT_DOSE_V01.md` established that the exact soma-conditioned adjoint is locally correct but cannot be extrapolated across a full bath-to-arbor conductance jump.

The next question was whether that merely means "take many small steps until the bond is fully mature," or whether **intermediate bond strength is itself computationally preferred**.

The answer is: **the interior optimum is real and held out.** A stronger claim that a majority of locally favored additions become actively harmful at the full binary endpoint narrowly failed its held-out threshold, so the correct statement is "often," not "usually."

## Setup

For each frozen body, legal tip-like additions and safe deletions were generated exactly as in the structural-interference / adjoint-dose line.

A base-state adjoint gave each event a directional derivative toward its corresponding binary structural change. Primary analysis kept only **gradient-favored additions**:

```text
g_e = grad C_lin · DeltaK_e > 0.
```

Each selected event was then evaluated exactly on the fixed conductance path

```text
alpha = 0,
.005, .01, .02, .03, .05, .075,
.10, .15, .20, .30, .40, .50, .60,
.75, 1.0
```

with no gradient recomputation along the path.

## Discovery — seeds 216-227

Gradient-favored additions:

```text
n                                      31
interior optimum                       26 / 31 = 0.8387
harmful at binary endpoint             17 / 31 = 0.5484
mean binary regret                     0.01795
regret > 1e-5                          26 / 31 = 0.8387
median interior alpha_best             0.04
slope-sign reversal                    28 / 31 = 0.9032
mean best partial gain                 +0.00957
mean forced-binary gain                -0.00838
```

All discovery D1-D4 criteria passed.

## Held-out confirmation — seeds 228-239

Gradient-favored additions:

```text
n                                      41
interior optimum                       33 / 41 = 0.8049
harmful at binary endpoint             18 / 41 = 0.4390
mean binary regret                     0.01833
regret > 1e-5                          33 / 41 = 0.8049
median interior alpha_best             0.10
slope-sign reversal                    35 / 41 = 0.8537
mean best partial gain                 +0.00999
mean forced-binary gain                -0.00834
```

### C1 — interior optima

Registered threshold `>75%`.

```text
33 / 41 = 80.5%
```

**C1 PASS.**

### C2 — binary endpoint becomes harmful

Registered threshold `>45%`.

```text
18 / 41 = 43.9%
```

**C2 FAIL, narrowly.**

One additional event would have crossed the threshold. That is not a reason to move the goalpost. The held-out result supports "binary forcing is frequently harmful," but not the preregistered stronger frequency claim.

### C3 — binary regret

```text
mean regret                           0.01833   > 0.010
positive regret                       80.5%    >= 75%
```

**C3 PASS.**

### C4 — interior scale is genuinely graded

Registered median interval `[0.02, 0.20]`.

```text
median alpha_best among interior      0.10
```

**C4 PASS.**

### C5 — nonmonotonic response

Registered threshold `>80%`.

```text
35 / 41 = 85.4%
```

**C5 PASS.**

So the confirmation is **4/5 registered criteria**, with the one failure being the exact frequency of endpoint sign reversal, not the existence of interior optima or binary regret.

## Deletion control

Gradient-favored deletions were less consistently interior:

```text
held-out n                            28
interior optimum                      15 / 28 = 53.6%
median interior alpha_best            0.60
harmful at binary endpoint            14 / 28 = 50.0%
slope reversal                        21 / 28 = 75.0%
```

This asymmetry is useful. A newly strengthened tip-like coupling often behaves as a tunable impedance element with a low/intermediate optimum, whereas weakening an established multi-bond cell more often prefers a larger move.

## Why this matters

The naive engineering prescription after the adjoint result was:

```text
relax binary topology -> optimize continuous density -> threshold back to binary.
```

The held-out response curves say the final thresholding step is **not free**. In many cases the intermediate coupling is where the computation lives.

This model therefore supports a different structural picture:

```text
field is already weakly coupled everywhere
        ↓
geometry is a spatial modulation of coupling strength
        ↓
local conductances tune transfer histories / impedance
        ↓
soma square-law readout exposes the resulting interference
```

The binary arbor is one limiting representation of that modulation, not obviously the unique natural state space.

## Interferometer language, carefully

The soma's quadratic readout really is interferometric in the signal-processing sense. But `mode 0` is not an externally supplied frequency-offset local oscillator, so **heterodyne** is stronger than the data justify.

A more precise phrase is:

> **reference-assisted square-law detection / self-homodyne-like mode mixing.**

Mode 0 is nearly order-blind alone but acts as a large common reference component; the actual computation remains redundant across the residual mode system.

Likewise, the connection to Moiré work is real at the level of **geometry making weak relational differences visible through interference**, but mode-0-plus-residual mixing is not literally a Moiré pattern unless one additionally establishes the near-spatial-frequency structure that defines Moiré interference.

## Current architecture

The line now reduces to:

```text
analog structural coupling field rho_e
        -> propagation operator A(rho)
        -> geometry-shaped h_A, h_B
        -> amplitude visibility + directional complex compatibility
        -> local quadratic/interferometric soma readout
        -> scalar task consequence
        -> soma-seeded adjoint field
        -> local forward × adjoint bond sensitivity
        -> small analog structural update
        -> changed propagation operator
```

That is a closed computational loop in differential form.

## Wall sentence

> **The soma is not merely where the field is read; its objective defines the backward sensitivity field. Forward propagation and backward credit meet locally on each bond, while the bond itself is best treated as a graded coupling whose useful value may lie far from either binary endpoint.**

## Next clean test

The next experiment should stop auditing one bond at a time and actually learn several frontier conductances jointly by **relinearized projected gradient ascent**:

1. compute the soma adjoint;
2. update frontier bond strengths by a small step;
3. recompute the field and adjoint;
4. repeat;
5. compare against a frozen-gradient and shuffled-gradient control.

If relinearized analog learning reliably improves the held-out objective while frozen/shuffled credit does not, then the v0.7-v0.9 eligibility wall has been solved in the model's continuous structural state space.
