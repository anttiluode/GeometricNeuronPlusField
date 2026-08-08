# Within-shell material placement confirmation — preregistration v0.1

Date frozen: 2026-08-08

## Claim under test

The corrected development hierarchy is:

```text
large effect:
    optimized graph-distance material profile

small residual:
    exact redistribution among cells at the same graph distance
```

The residual is tested with a nested design: first optimize the strongest distance-only material field, then freeze every distance shell's total material exactly and release only within-shell placement.

This preregistration tests whether that small residual generalizes to untouched FunctionalArbor bodies.

No seed in the held-out set below has been examined for this test or for the superseded branch-residual preregistration.

---

## Frozen bodies

```text
seeds 628 through 639 inclusive
12 requested bodies
```

Require at least 10 usable bootstrapped bodies.

---

## Frozen physics and objective

```text
omega               .03 and .04
tau_h                2
mu                   .5
budget reference g0  .005
budget ratio          10
local density cap     .05
uniform initialization
fixed total material budget
step_fraction         .10
```

Objective:

```text
u_j = H_j / |H_j|
R^2 = |mean_j u_j|^2
```

maximize mean `R^2` over `.03,.04` and the frozen injection-site set.

---

## Stage 1 — optimized radial learner

One shared density per integer graph-distance shell from the soma/readout.

```text
160 projected-gradient steps
exact analytic shell gradient
ordinary Euclidean shell-parameter projection
same cap and total material budget
```

This is the hostile distance-only control.

---

## Stage 2 — within-shell release

Starting from the radial optimum:

1. record total material in every graph-distance shell;
2. freeze each shell total;
3. release independent cell densities inside each shell;
4. take 80 projected-gradient steps on the same joint objective;
5. each update is projected separately inside every shell, so no radial material can move between shells.

The release therefore tests only cell/branch placement beyond graph distance.

---

## Stage 3 — same-shell placement shuffle

For each learned release map, perform 12 random controls per body.

Within every graph-distance shell separately preserve:

```text
shell total
exact multiset / histogram of learned densities
```

but randomly permute those values among equal-distance cells.

This destroys exact branch/cell placement while preserving the radial profile and the amount of heterogeneity inside each shell.

---

# Frozen criteria

## W0 — usable population

```text
at least 10 / 12 requested bodies bootstrap successfully
```

---

## W1 — distance-only material remains the dominant useful component

Require:

```text
mean(radial R2 - uniform R2) > .06
positive on at least 75% of usable bodies
```

This must pass before interpreting a branch residual.

---

## W2 — branch-only release improves the optimized radial solution

Require:

```text
mean(released R2 - radial R2) > .002
positive on at least 75% of usable bodies
```

Development value was `+.00456`, positive `6/6`.

---

## W3 — exact within-shell placement matters

Require:

```text
mean(released R2 - within-shell-shuffle R2) > .0025
positive on at least 75% of usable bodies
```

Development value was `+.00570`, positive `6/6`.

This is the central placement criterion.

---

## W4 — direct circular phase spread agrees

Using released versus optimized radial:

```text
mean(radial phase RMS - released phase RMS) > .005 rad
positive on at least 2/3 of usable bodies
```

Development value was `+.01217 rad`, positive `6/6`.

---

## W5 — residual is not an attenuation trick

For each body:

```text
amp_ratio = released mean median soma amplitude / radial mean median soma amplitude
```

Require:

```text
pooled median amp_ratio between .75 and 1.25
at least 75% of bodies between .60 and 1.40
```

Development median was `.996`.

---

## W6 — both frozen frequencies retain a non-radial placement signal

For each frequency separately require both mean differences to be positive:

```text
mean(released - radial) R2 > 0
mean(released - within-shell shuffle) R2 > 0
```

at both `.03` and `.04`.

No minimum per-frequency effect size beyond positive sign is imposed.

---

## W7 — distance remains quantitatively dominant

Define:

```text
G_radial  = mean(radial - uniform)
G_release = mean(released - radial)
```

Require:

```text
G_radial > 0
G_release / G_radial < .15
```

Development ratio was about `.044`.

This criterion prevents a passing residual from being retold as though graph distance were not the main coordinate.

---

# Verdict rule

All eight criteria `W0` through `W7` must pass.

```text
8 / 8 PASS
    -> held-out support for a small within-distance placement residual

otherwise FAIL
    -> report the failed criteria and keep distance compensation as the dominant result
```

No thresholds change after seeds 628–639 are examined.

---

## Interpretation boundary if it passes

A pass establishes only this model-level statement:

> **After the optimal readout-distance material profile is frozen, exact local material credit can reproducibly find a small additional arrangement among equal-distance cells that improves phase coordination, and the exact within-shell placement performs better than shell-histogram-matched shuffles.**

It does not establish novelty.

Prior art already owns broad versions of branch-specific conductance heterogeneity, morphology-dependent channel maps, dendritic democracy, intrinsic plasticity, adjoint channel gradients, and differentiable branch-conductance optimization.

## Interpretation boundary if it fails

Then the honest material conclusion becomes even simpler:

> The present learner mainly rediscovers a readout-distance compensation profile; branch/cell corrections observed in development are not robust enough to promote.
