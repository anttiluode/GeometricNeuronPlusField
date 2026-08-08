# Material-adjoint learning confirmation preregistration v0.1

Date: 2026-08-08

This file freezes the learning algorithm, material constraints, untouched bodies and acceptance criteria **before** seeds 584–595 are examined.

## Claim being tested

The development result in `MATERIAL_LEARNING_DISCOVERY_V01.md` suggests that a local material-adjoint learner, given only a somatic phase-coherence objective and no morphology coordinate, can redistribute a fixed budget of quasi-active material so that:

1. somatic transfer phases become more coherent;
2. the learned placement matters, not just the density histogram;
3. the learned solution beats the previously confirmed hand-drawn linear distance gradient;
4. a positive density-vs-distance-from-readout organization emerges without being specified;
5. this is not obtained by arbitrarily extinguishing transfer amplitude.

The moved-readout re-centering result is a separate development mechanism control and is **not** part of this primary confirmation.

---

## Frozen physics

Use the complete minimal quasi-active material model from `hcn_quasiactive_probe.py` / `hcn_material_learning.py`:

```text
v' = K L psi - damping*v - restoring*psi
     - d(x)*psi
     - mu*d(x)*z
     + source

z' = (psi-z)/tau_h
```

Frozen parameters:

```text
omega      .03 and .04
tau_h       2
mu          .5
```

This remains a second-order reciprocal wave model plus one-state quasi-active material, not a conductance-based biological HCN model.

---

## Frozen material budget and optimizer

The fixed total material budget is the sum of the previously confirmed hand profile

```text
g0                  .005
soma->distal ratio  10x
```

but the learner starts from **uniform density** and is never shown that profile.

Constraints:

```text
d_i >= 0
sum_i d_i = fixed total material budget
d_i <= .05
```

Optimizer:

```text
50 projected-gradient steps
step_fraction = .10
backtracking acceptance exactly as in hcn_material_learning.py
```

No graph-distance, branch, source-terminal, HCN-gradient or distal-enrichment term enters the objective or update.

---

## Frozen objective

At each frequency independently drive many occupied sites and measure soma transfer `H_j`.

Normalize each location to its unit phase phasor

```text
u_j = H_j / |H_j|
```

and maximize

```text
R^2 = |mean_j u_j|^2.
```

Training objective:

```text
mean R^2 over omega=.03,.04.
```

Absolute transfer amplitude does not weight the phase objective.

---

## Controls

After learning compare the identical frozen body with:

1. `uniform` — initial uniform density;
2. `learned` — optimized local density;
3. `hand` — confirmed hand-drawn linear soma-to-distal profile under the same total budget;
4. `shuffle_learned` — learned density values globally shuffled over occupied cells.

Graph distance is computed only **after** learning to characterize the learned solution.

---

## Development receipt motivating thresholds

Seeds 580–583:

```text
mean R^2
uniform          .48438
learned          .61750
hand             .55241
shuffled learned .45339

learned - uniform per body
+.1229 +.1121 +.1105 +.1870

learned - hand per body
+.0777 +.0569 +.0581 +.0678

learned - shuffled per body
+.1689 +.1243 +.1416 +.2216

Spearman learned density vs soma graph distance
+.907 +.676 +.901 +.609
mean +.773

learned/uniform median-soma-amplitude ratio by body
.780 .595 .842 .602
```

Per-frequency learned-minus-uniform coherence gain:

```text
omega=.03 mean +.0843
omega=.04 mean +.1820
```

The cap-robustness control subsequently showed positive distance correlation on all four development bodies at caps `.035`, `.05`, `.075`, `.10`; the held-out test nevertheless freezes the original `.05` cap.

---

## Held-out bodies

Fresh FunctionalArbor bodies:

```text
seeds 584–595
12 requested bodies
```

No threshold will be changed after these results are inspected.

---

# Frozen criteria

## L0 — analytic material gradient remains numerically correct

Each body includes the existing finite-difference audit.

```text
maximum relative derivative error < 1e-5
```

This is an implementation integrity criterion rather than a biological claim.

## L1 — learning improves the frozen phase objective

```text
mean(learned R^2 - uniform R^2) > .05
```

and learned must beat uniform on at least `9/12` bodies. If fewer than 12 bodies bootstrap, require at least 75% positive bodies.

## L2 — local placement matters

```text
mean(learned R^2 - shuffled-learned R^2) > .08
```

and learned must beat its globally shuffled density map on at least 75% of bodies.

## L3 — learning beats the hand-drawn gradient

```text
mean(learned R^2 - hand R^2) > .025
```

and learned must beat the hand profile on at least 75% of bodies.

This criterion is intentionally strong: a pass means the discovered local material arrangement is doing more than merely approximating the original biological prior.

## L4 — readout-distance organization emerges without being trained

After optimization only, compute Spearman correlation between local learned density and graph distance from the soma/readout.

Registered:

```text
mean Spearman rho > .40
positive rho on at least 75% of bodies
```

The threshold is well below the development mean `.773`.

## L5 — amplitude remains in a nondegenerate range

For each body, use the mean over the two frequencies of the median soma-transfer amplitude over injection sites.

Compute

```text
A_ratio = learned / uniform.
```

Registered pooled median:

```text
.40 < median(A_ratio) < 1.50
```

and at least 75% of bodies must have

```text
.30 < A_ratio < 2.0.
```

This does not force amplitude preservation. It only rejects a phase solution that effectively annihilates or explosively amplifies the transfer channel.

## L6 — both frozen frequencies improve independently

For **each** of `omega=.03` and `.04`:

```text
mean(learned R^2 - uniform R^2) > .02.
```

No per-body frequency threshold is registered.

## L7 — phase spread itself moves in the expected direction

Although `R^2` is the optimized quantity, the directly interpretable circular phase-spread measure should agree.

```text
mean(uniform phase RMS - learned phase RMS) > .08 rad.
```

---

## Confirmation rule

The primary material-learning claim is called **held-out confirmed only if all L0–L7 pass**.

A failure remains a failure. In particular:

- positive distance correlation without objective gain is not enough;
- objective gain without defeating shuffled placement is not enough;
- objective gain that cannot beat the hand profile is not enough for the stronger self-organization claim;
- a profile that appears distance-organized only because transfer amplitude collapses is not enough.

---

## What a pass would mean

A pass supports this model-level statement:

> **Starting from uniform electrical material and without being given a morphology coordinate, the reciprocal material-adjoint learner repeatedly discovers a readout-distance-organized quasi-active material distribution that improves somatic phase coordination beyond uniform, shuffled and hand-drawn profiles under the same material budget.**

It would not show that real HCN channel trafficking implements this gradient rule.

---

## What comes after a pass

The next confirmation should target the development-only moved-readout result in `MATERIAL_READOUT_RECENTER_V01.md`:

```text
move the consequential readout
-> transpose field changes
-> learned material coordinate should re-center around the new readout.
```

Only after that should the project return to the gamma/theta cross-frequency branch.
