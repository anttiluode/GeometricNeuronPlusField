# Material readout re-centering control v0.1

Date: 2026-08-08

## Question

The material learner in `hcn_material_learning.py` starts from uniform density and discovers strong positive density correlation with graph distance from the soma, even though graph distance is absent from the objective.

But the soma is also the location where the objective is measured.

So there is an obvious hostile control:

> **Move the consequential readout to the far side of the same frozen arbor. Does the learned electrical-material coordinate remain tied to the anatomical soma, or re-center around the new readout?**

This is a development mechanism control on the same seeds 580–583 used for material-learning discovery. It is not a held-out confirmation.

---

## Intervention

For each frozen body:

1. choose the occupied cell with maximum graph distance from the anatomical soma as a second readout;
2. initialize the same uniform quasi-active material density under the same fixed total budget and `.05` local cap;
3. optimize the same two-frequency soma-style phase-coherence objective, now evaluated at the moved readout;
4. never expose graph distance, soma identity or a desired spatial profile to the optimizer;
5. after learning, compute density correlation with:
   - graph distance from the **moved readout**;
   - graph distance from the **anatomical soma**.

For comparison, rerun the original soma-centered learner under the same code path.

---

## Development result

### Original soma readout

```text
mean learning gain in coherence R^2      +.1331
mean Spearman density vs readout distance +.7733
positive readout-distance rho             4 / 4
```

This reproduces the previous material-learning discovery.

### Moved readout

Per body:

```text
seed 580
  moved readout [6,3]
  rho(density, distance from moved readout)   +.776
  rho(density, distance from old soma)        -.088
  learning gain                               +.124

seed 581
  moved readout [5,28]
  rho(new readout distance)                   +.582
  rho(old soma distance)                      -.861
  gain                                        +.104

seed 582
  moved readout [6,28]
  rho(new readout distance)                   +.631
  rho(old soma distance)                      +.106
  gain                                        +.081

seed 583
  moved readout [28,27]
  rho(new readout distance)                   +.971
  rho(old soma distance)                      -.779
  gain                                        +.146
```

Pooled:

```text
mean moved-readout learning gain              +.1137
mean rho to moved-readout distance             +.7399
positive rho to moved readout                   4 / 4

mean rho to anatomical-soma distance           -.4055
negative rho to anatomical soma                 3 / 4
```

---

## Interpretation

The distance organization is not anchored to a privileged hard-coded soma coordinate.

When the consequential observation point is moved, the learned material map reorganizes around the new readout.

The cleanest earned statement is therefore not

```text
material becomes distal from the soma.
```

It is

```text
material tends to increase with transfer distance from the consequential readout.
```

The biological CA1 soma happens to be one important consequential readout, so distal HCN-like enrichment can be understood in this model as one instance of a more general readout-relative compensation geometry.

This also sharpens the old soma result from the Geometric Neuron line. The soma need not be a mathematically privileged spatial coordinate. Its role is defined by where distributed transfer relations are required to become consequential.

---

## Why this matters for the adjoint story

The local material derivative is

```text
dH/dd_i = -c(omega) y_i x_i,
```

where the transpose field `y` is launched from the current readout.

Move the readout and the transpose field changes automatically.

Therefore the local sensitivity landscape changes without giving each material tuner a new global coordinate system.

In the current reciprocal model, the readout itself physically defines the coordinate by launching the return/transpose field.

This gives a compact interpretation:

```text
forward field says: what is happening here?
return field says: how does here matter to consequence over there?
local overlap says: how should this electrical material change?
```

The emerging distance relation is not explicitly measured by the local tuner. It is implicit in the geometry of those two fields.

---

## Current wall sentence

> **The learned electrical-material coordinate follows the site of consequence: moving the readout across the same frozen arbor re-centers the emergent density-vs-distance organization around the new readout rather than the anatomical soma.**

This is still a development result. A future held-out material-learning confirmation should include a moved-readout criterion if re-centering is part of the final claim.
