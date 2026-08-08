# Reciprocal credit-wave identity preregistration v0.1

## Question

`ANALOG_FRONTIER_LEARNING_V01.md` leaves one physical-plausibility wall: the exact adjoint is written as an algorithmic reverse pass.

For the mature linear FunctionalArbor wave, however, the spatial operator is reciprocal/symmetric. Eliminating velocity from the semi-implicit update reveals a stronger fact:

> **The exact bond-credit field should be generatable by the same damped wave medium, driven from the soma by the time-reversed objective-derivative waveform.**

This is a mathematical identity claim for the linear model, not a biological claim.

## Derivation fixed before the run

Write

```text
a = 1 - gamma dt
H = c L(k) - rho I
```

for the reciprocal spatial operator. The forward update

```text
v[n+1]   = a v[n] + dt H psi[n] + dt s[n]
psi[n+1] = psi[n] + dt v[n+1]
```

eliminates to

```text
psi[n+1] = (1+a) psi[n] - a psi[n-1]
           + dt^2 H psi[n] + dt^2 s[n].
```

For an objective `J`, let `p_psi,p_v` be the exact discrete adjoint and define the bond-credit field used by the conductance gradient

```text
mu[n] = dt p_psi[n] + p_v[n].
```

The reverse recurrence is

```text
mu[n] = (1+a) mu[n+1] - a mu[n+2]
        + dt^2 H mu[n+1] + dt g[n]
```

where

```text
g[n] = dJ/dpsi[n]
```

is nonzero only at the soma for the present objective.

Reverse the time index. The recurrence is then **the same forward damped wave equation** provided the retrograde source is

```text
s_retro[r] = g[T-r] / dt.
```

Thus a zero-state forward run of the same reciprocal medium with the reversed soma derivative waveform should reconstruct `mu` exactly (up to index convention), without changing damping sign or using a separately engineered transport graph.

The bond gradient is then the local overlap of the ordinary forward field difference and this retrograde field difference.

## Objective

Use the same smooth normalized integrated soma-energy contrast as the analog-learning experiments:

```text
C = (E_T-E_D)/(E_T+E_D).
```

Target and distractor trajectories each provide their own soma derivative waveform and retrograde wave; the bond gradients add.

## Fresh bodies

```text
seeds 264-275
```

Same mature v0.9 bootstrap, exact mature-boundary linear operator, lag 20, 210 steps.

## Registered tests

### R1 — reciprocal-wave bond map equals the algorithmic adjoint

For every horizontal and vertical grid bond compare the exact adjoint conductance gradient with the forward×retrograde overlap reconstructed from the physical replay.

PASS if across bodies:

```text
mean map correlation > 0.999999
mean relative L2 error < 1e-8.
```

### R2 — frontier event scores are identical

Generate the same legal tip-like frontier candidates used by analog learning and compare directional event derivatives.

PASS if pooled frontier-score correlation `>0.999999` and maximum normalized absolute error `<1e-7`.

### R3 — one learning step is numerically identical

Starting from the same candidate `rho`, normalize both exact-adjoint and reciprocal-wave frontier gradients by their own max absolute value and take `eta=.01` projected steps.

PASS if the maximum absolute difference in the updated `rho` vector is `<1e-8`.

## Ablation: time order of the retrograde source

As a descriptive control, run the same soma derivative waveform **without reversing it in time** and report its gradient-map correlation with the exact adjoint.

No pass/fail threshold is preregistered. This asks whether the time reversal is merely cosmetic or carries essential temporal credit structure.

## Interpretation fixed in advance

If R1-R3 pass, then in this reciprocal linear wave model the algorithmic adjoint has an exact wave implementation:

```text
forward task field
        ×
same-medium, time-reversed soma credit wave
        -> local bond sensitivity.
```

The old retrograde-carrier idea would therefore return in a more specific form. What must travel backward is **not a scalar reward** but a task-conditioned waveform containing the soma's temporal/complex derivative information.

This still requires storage/reversal (or some physical mechanism that produces the conjugate/reversed waveform), and therefore does not by itself establish biological plausibility.
