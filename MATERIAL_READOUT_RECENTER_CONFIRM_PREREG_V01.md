# Material readout re-centering confirmation preregistration v0.1

Date: 2026-08-08

This preregistration freezes the moved-readout mechanism test before fresh seeds 596–603 are examined.

## Development observation

On seeds 580–583, the material learner was run twice on the same frozen body:

1. objective measured at the anatomical soma;
2. objective moved to the occupied cell farthest in graph distance from the soma.

The learner never received graph distance as an input.

Development result:

```text
soma-readout
  mean coherence gain                         +.1331
  mean rho(density, distance-to-readout)      +.7733

moved-readout
  mean coherence gain                         +.1137
  mean rho(density, distance-to-new-readout)  +.7399
  mean rho(density, distance-to-old-soma)     -.4055

new-readout rho positive                       4 / 4
new-readout rho > old-soma rho                 4 / 4
```

The result suggests that the emergent material coordinate follows the site where consequence is measured rather than remaining tied to the anatomical soma.

## Frozen intervention

Reuse `material_readout_recenter.py` unchanged with:

```text
omega          .03,.04
tau_h           2
mu              .5
budget ref g0   .005
budget ratio    10
local cap       .05
steps           50
step_fraction   .10
```

For each body:

- run the original soma-centered learner;
- choose the occupied cell at maximum graph distance from the soma as the moved readout;
- restart from uniform material;
- train the identical objective around the moved readout;
- compute density correlations only after learning.

## Held-out bodies

```text
seeds 596–603
8 requested bodies
```

## Frozen criteria

### R0 — moved-readout learning remains useful

```text
mean moved-readout coherence gain > .05
positive gain on >= 75% of bodies
```

### R1 — density organizes by distance from the moved readout

```text
mean Spearman rho(density, distance-to-new-readout) > .40
positive rho on >= 75% of bodies
```

### R2 — the new readout explains density better than the old soma

For every body define

```text
delta_rho = rho_to_new_readout - rho_to_anatomical_soma.
```

Registered:

```text
mean delta_rho > .40
delta_rho > 0 on >= 75% of bodies
```

### R3 — the material map is not still soma-centered in the population

```text
mean rho(density, distance-to-anatomical-soma) < .15
```

This does not require every individual body to be negatively soma-correlated.

### R4 — the original soma-centered result reproduces in the same fresh bodies

```text
mean soma-centered coherence gain > .05
mean soma-centered rho(density, distance-to-soma) > .40
```

### R5 — re-centering does not require a weaker learning effect

The moved objective may differ in difficulty, but the gain should remain of the same order:

```text
mean moved gain / mean soma gain > .45
```

No upper ratio is registered.

## Confirmation rule

The readout-recentering claim is held-out confirmed only if all R0–R5 pass.

## Meaning of a pass

A pass would support the model-level statement:

> **The emergent material-distance coordinate is consequence-relative rather than soma-hard-coded: moving the readout changes the transpose field and re-centers the learned electrical-material organization around the new consequential site.**

It would not imply that biological neurons literally relocate their soma or that HCN trafficking follows this exact optimizer. The moved readout is a mechanism control for the geometry of the learning rule.
