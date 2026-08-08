# Material-adjoint learning confirmation preregistration v0.2

Date: 2026-08-08

This is a new held-out confirmation after the formal `7/8 FAIL` in `MATERIAL_LEARNING_CONFIRM_V01.md`.

The v0.1 failure is preserved and not reinterpreted. Its sole failed criterion was the derivative-audit guard: a mathematically near-zero analytic derivative was compared with a `~1.7e-10` finite-difference residue using a relative-error denominator that was itself near zero, producing a meaningless relative error near one.

All seven substantive v0.1 learning criteria passed strongly on seeds 584–595.

v0.2 changes **only the numerical derivative-audit statistic**. The learning algorithm, objective, material physics, constraints, controls and all substantive thresholds remain frozen unchanged.

---

## Frozen physics and learner

Exactly as v0.1:

```text
omega               .03 and .04
tau_h                2
mu                   .5
budget reference g0  .005
budget ratio         10
local density cap    .05
uniform initialization
50 projected-gradient steps
step_fraction        .10
fixed total material budget
```

Objective:

```text
for each source location j:
  u_j = H_j / |H_j|

R^2 = |mean_j u_j|^2

maximize mean R^2 over omega=.03,.04
```

No graph-distance, branch, source-terminal, hand-gradient or distal-enrichment term enters training.

Controls remain:

```text
uniform
learned
hand linear distance profile
shuffled learned values
```

---

## Fresh held-out bodies

```text
seeds 604–615
12 requested bodies
```

These bodies were not used in the development set, v0.1 held-out set, or moved-readout held-out set.

No threshold below will be altered after inspecting them.

---

# Corrected L0 — zero-safe derivative audit

For every sampled analytic/finite-difference derivative pair define

```text
scale = |analytic| + |finite_difference|
abs_error = |analytic - finite_difference|
```

### Nonzero-scale samples

If

```text
scale >= 1e-6
```

require

```text
abs_error / (scale + 1e-30) < 1e-5.
```

### Near-zero samples

If

```text
scale < 1e-6
```

relative error is not used.

Require instead

```text
abs_error < 1e-8.
```

### L0 pass rule

Across all sampled derivatives on all fresh bodies:

```text
max nonzero-scale relative error < 1e-5
and
max near-zero absolute error < 1e-8.
```

If no samples fall in one regime, that regime is vacuously satisfied and the count is reported.

This rule was frozen after diagnosing v0.1 but before inspecting seeds 604–615.

For context only, the v0.1 held-out worst absolute derivative discrepancy was `~5.44e-10`, so the new near-zero threshold is about 18× looser than that observed residue.

---

# Substantive criteria — unchanged from v0.1

## L1 — learning improves the frozen phase objective

```text
mean(learned R^2 - uniform R^2) > .05
```

and learned must beat uniform on at least 75% of bodies.

## L2 — local placement matters

```text
mean(learned R^2 - shuffled-learned R^2) > .08
```

and learned must beat the globally shuffled learned map on at least 75% of bodies.

## L3 — learning beats the hand-drawn gradient

```text
mean(learned R^2 - hand R^2) > .025
```

and learned must beat the hand profile on at least 75% of bodies.

## L4 — readout-distance organization emerges without being trained

After optimization only:

```text
mean Spearman rho(density, graph distance from soma/readout) > .40
```

and positive rho on at least 75% of bodies.

## L5 — amplitude remains nondegenerate

For each body compute

```text
A_ratio = learned mean-median soma amplitude
          / uniform mean-median soma amplitude.
```

Require pooled median

```text
.40 < median(A_ratio) < 1.50
```

and at least 75% of bodies within

```text
.30 < A_ratio < 2.0.
```

## L6 — both frozen frequencies improve independently

For each of `omega=.03` and `.04`:

```text
mean(learned R^2 - uniform R^2) > .02.
```

## L7 — circular phase spread agrees with the optimized metric

```text
mean(uniform phase RMS - learned phase RMS) > .08 rad.
```

---

## Confirmation rule

v0.2 is held-out confirmed only if **all L0–L7 pass**.

A failure remains a failure.

The purpose of v0.2 is not to obtain an all-pass label by relaxing the learning claim. It is to replace a numerically undefined relative-error test at zero with a scale-aware derivative audit while leaving all substantive predictions untouched.

---

## Meaning of a pass

A pass would support the following model-level result with a numerically valid implementation audit:

> **Starting from uniform electrical material and without being given a morphology coordinate, the reciprocal material-adjoint learner repeatedly discovers a readout-distance-organized quasi-active material distribution that improves somatic phase coordination beyond uniform, shuffled and hand-drawn profiles under the same fixed material budget.**

The separate moved-readout confirmation already establishes that the emergent coordinate re-centers around consequence rather than being anatomically soma-hard-coded.
