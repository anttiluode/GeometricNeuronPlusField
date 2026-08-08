# Material branch-specific residual confirmation — preregistration v0.1

Date frozen: 2026-08-08

## Claim under test

Prior-art review and development controls indicate that most of the confirmed material-learning effect is explainable as readout-distance compensation.

The remaining candidate result is narrower:

> **A full cell-by-cell material learner may exploit branch/topological information beyond anything available to a material profile constrained to be only a function of graph distance from the readout.**

This preregistration tests that residual on untouched FunctionalArbor bodies.

No fresh-body result has been inspected before this file is committed.

---

## Frozen bodies

```text
seeds 628 through 639 inclusive
12 requested bodies
```

Any body that fails FunctionalArbor bootstrap is reported and omitted. The test requires at least 10 usable bodies.

---

## Frozen physical/material model

Same as the confirmed material learner:

```text
omega               .03 and .04
tau_h                2
mu                   .5
budget reference g0  .005
budget ratio          10
local per-cell cap    .05
uniform initialization
fixed total material budget
step_fraction         .10
```

Objective at the anatomical soma/readout:

```text
u_j = H_j / |H_j|
R^2 = |mean_j u_j|^2
```

maximize mean `R^2` over `.03,.04` and the same injection-site set.

---

## Learners compared

### Full learner

One independent material density for every occupied cell.

```text
50 projected-gradient steps
exact analytic material gradient
```

### Optimized radial learner

One shared material density for every integer graph-distance shell from the soma/readout.

All cells at the same graph distance must have the same density.

```text
120 projected-gradient steps
exact analytic shell gradient
same total material budget
same per-cell cap
same uniform start
same objective
```

The radial learner is deliberately given more iterations because it is the hostile control.

This is not the hand-drawn linear HCN gradient. It is a freely optimized distance-only solution.

---

# Frozen criteria

## B0 — distance-only material remains useful

The radial learner must improve over uniform:

```text
mean(radial R^2 - uniform R^2) > .06
positive on at least 9 / 12 usable bodies
```

Purpose: confirm that the dominant distance-compensation effect is present on the fresh population.

---

## B1 — full local material has a reproducible branch-specific residual

```text
mean(full R^2 - radial R^2) > .005
full > radial on at least 9 / 12 usable bodies
```

This is the central criterion.

Development margin was approximately `+.0124 R^2`; the frozen threshold deliberately asks for less than half that mean effect.

---

## B2 — direct circular phase spread agrees

Define

```text
phase_margin = radial mean soma phase RMS - full mean soma phase RMS
```

Require:

```text
mean phase_margin > .005 rad
positive on at least 8 / 12 usable bodies
```

Purpose: the residual must not exist only in the complex-coherence statistic.

---

## B3 — residual is not an attenuation trick

For each body define

```text
amp_ratio = full mean median soma amplitude / radial mean median soma amplitude
```

Require:

```text
pooled median amp_ratio between .50 and 1.50
at least 9 / 12 bodies between .40 and 1.80
```

---

## B4 — both frozen frequencies contribute

For each frequency separately compute the mean full-minus-radial coherence difference across bodies.

Require:

```text
omega=.03   mean(full R^2 - radial R^2) > 0
omega=.04   mean(full R^2 - radial R^2) > 0
```

No minimum effect size beyond positive sign is imposed per frequency.

---

## B5 — distance remains the dominant first-order coordinate

Let

```text
G_radial = mean(radial - uniform)
G_full   = mean(full - uniform)
```

Require:

```text
G_radial / G_full > .60
```

provided `G_full > 0`.

This criterion prevents us from rewriting the story as “distance was irrelevant” if the full learner wins. The current hypothesis is hierarchical:

```text
large first-order effect     distance from consequence
smaller residual             branch/cell-specific geometry
```

---

# Verdict rule

All six criteria must pass.

```text
6 / 6 -> held-out support for a branch-specific residual beyond optimized distance-only material
otherwise -> FAIL; keep the dominant result as distance compensation and report which residual criteria failed
```

No threshold changes after fresh seeds are examined.

---

## Interpretation boundary if it passes

A pass would **not** establish novelty.

It would establish only this model-level fact:

> A freely optimized distance-shell material field captures most of the computational gain, but a full local forward×transpose learner reproducibly extracts additional useful structure tied to the irregular arbor beyond graph distance alone.

Prior art would still own the broad ideas of dendritic democracy, local intrinsic plasticity, adjoint channel gradients, and gradient-optimized channel distributions.

## Interpretation boundary if it fails

If the full learner does not reproducibly beat optimized radial material, the strongest honest description becomes:

> The present local material-adjoint mechanism mostly rediscovers an established readout-distance compensation principle in a reciprocal wave-arbor implementation.
