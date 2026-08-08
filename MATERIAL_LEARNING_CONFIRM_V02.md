# Material-adjoint learning confirmation v0.2 — held-out 8/8 pass

Date: 2026-08-08

## Result in one sentence

> **Starting from uniform quasi-active electrical material and without receiving any morphology coordinate, a fixed-budget local material-adjoint learner repeatedly discovers a readout-distance-organized material distribution that improves somatic phase coordination beyond uniform, globally shuffled, and hand-drawn distance-gradient controls; all 8/8 preregistered criteria passed on 12 fresh bodies.**

This is a computational result in the present reciprocal wave-arbor model. It does not establish a biological HCN trafficking rule.

---

## Why v0.2 was necessary

The first held-out confirmation, recorded in `MATERIAL_LEARNING_CONFIRM_V01.md`, formally failed `7/8` because its derivative-audit statistic used relative error even at a true zero derivative.

One analytic derivative on seed 584 was approximately `9.5e-17`; central finite difference left a harmless `~1.7e-10` residue. Dividing that difference by a denominator that was itself effectively zero produced a reported relative error near one.

All seven substantive v0.1 learning predictions passed strongly.

v0.2 therefore changed exactly one thing:

```text
nonzero derivative scale -> relative-error audit
near-zero derivative scale -> absolute-error audit
```

All material physics, learning code, constraints, controls and substantive thresholds were left unchanged.

The correction and fresh seed set were frozen in `MATERIAL_LEARNING_CONFIRM_PREREG_V02.md` before seeds 604–615 were examined.

---

# Frozen learner

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
u_j = H_j / |H_j|
R^2 = |mean_j u_j|^2
```

maximize mean `R^2` over the two frozen frequencies.

No graph distance, branch identity, hand gradient or distal-enrichment target enters the objective or update.

Controls:

```text
uniform
learned
hand linear distance profile
globally shuffled learned values
```

---

# Held-out confirmation — seeds 604–615

All 12 requested bodies were available.

Mean coherence `R^2`:

```text
uniform             .4649858
learned             .5977747
hand                .5372109
shuffled learned    .4229199
```

Mean circular soma phase RMS:

```text
uniform             .9196602
learned             .7069986
hand                .7753345
shuffled learned    .9704812
```

---

## L0 — zero-safe analytic gradient audit

Registered:

```text
if derivative scale >= 1e-6:
    relative error < 1e-5

if derivative scale < 1e-6:
    absolute error < 1e-8
```

Observed:

```text
nonzero-scale samples             60
near-zero samples                  0
max nonzero relative error      1.211e-7
max absolute error overall      7.408e-10
```

Worst nonzero sample:

```text
seed 605
analytic          .0002479535511
finite difference .0002479534911
absolute error    6.006e-11
relative error    1.211e-7
```

**L0 PASS.**

The corrected audit confirms that the analytic material-adjoint derivative is numerically consistent with finite differences on the fresh set.

---

## L1 — learning improves the frozen objective

Registered:

```text
mean learned-uniform > .05
positive on >=75% bodies
```

Observed:

```text
mean gain    +.1327890
positive      12 / 12
```

**L1 PASS.**

---

## L2 — learned placement matters

Registered:

```text
mean learned-shuffled > .08
positive on >=75% bodies
```

Observed:

```text
mean gain    +.1748548
positive      12 / 12
```

**L2 PASS.**

The density histogram alone is insufficient. Where the material is placed on the frozen arbor matters strongly.

---

## L3 — learning beats the hand-drawn distance gradient

Registered:

```text
mean learned-hand > .025
positive on >=75% bodies
```

Observed:

```text
mean gain    +.0605639
positive      12 / 12
```

**L3 PASS.**

The learned solution therefore does more than approximate the simple linear biological prior.

---

## L4 — a readout-distance organization emerges without being trained

Registered after learning only:

```text
mean Spearman rho(density, graph distance from soma/readout) > .40
positive on >=75% bodies
```

Observed:

```text
mean rho       +.7647426
minimum rho    +.4897622
positive        12 / 12
```

**L4 PASS.**

Every fresh body independently learned a positive density-versus-readout-distance organization.

This coordinate is not present in the loss or update rule.

---

## L5 — transfer amplitude remains nondegenerate

Registered pooled median learned/uniform amplitude ratio between `.40` and `1.50`, with at least 75% of bodies between `.30` and `2.0`.

Observed:

```text
median ratio     .6431974
mean ratio       .6536222
in range          12 / 12
```

**L5 PASS.**

The learned solution attenuates median transfer more than the confirmed hand-gradient material, but it does not obtain phase coherence by turning the transfer channel off.

---

## L6 — both frozen frequencies improve independently

Registered mean learned-uniform `R^2` gain > `.02` at each frequency.

Observed:

```text
omega=.03    +.1025495
omega=.04    +.1630284
```

**L6 PASS.**

---

## L7 — direct circular phase spread agrees

Registered:

```text
mean(uniform phase RMS - learned phase RMS) > .08 rad
```

Observed:

```text
+.2126616 rad
```

All 12 bodies improved.

**L7 PASS.**

---

# Confirmation verdict

```text
L0 PASS   zero-safe material-gradient audit
L1 PASS   learned objective improves
L2 PASS   learned placement beats shuffled values
L3 PASS   learned map beats hand distance gradient
L4 PASS   readout-distance organization emerges
L5 PASS   transfer amplitude remains nondegenerate
L6 PASS   both frozen frequencies improve
L7 PASS   circular phase spread improves

8 / 8 PASS
```

---

## Relationship to the v0.1 formal failure

The history remains:

```text
v0.1 formal confirmation   7 / 8 FAIL
    sole failure = ill-conditioned relative-error audit at zero

v0.2 fresh confirmation   8 / 8 PASS
    corrected scale-aware audit
    substantive criteria unchanged
```

Nothing in v0.1 is deleted or retroactively relabeled.

---

## Combined with the moved-readout confirmation

`MATERIAL_READOUT_RECENTER_CONFIRM_V01.md` independently showed on fresh bodies that moving consequence away from the soma causes the learned material coordinate to re-center around the new readout:

```text
moved-readout density vs new-readout distance
mean Spearman rho     +.668
positive               8 / 8

same learned density vs old soma distance
mean rho              -.648

new-vs-old coordinate delta
mean                  +1.316
positive               8 / 8
```

So the confirmed abstraction is not merely

```text
HCN-like material likes distal dendrites.
```

It is

```text
local material learning discovers a coordinate
relative to where distributed dynamics become consequential.
```

---

## Current geometric interpretation

The project began with anatomy as geometry.

The confirmed material branch now requires a broader object:

```text
morphology / connectivity G
        +
spatial electrical material M(x,omega)
        +
moving field psi(x,t)
        +
consequential readout R
```

The readout launches the transpose field that supplies local sensitivity.

That makes the effective geometry **relational**:

```text
not just where a cell is in the arbor,

but how that location sits between
local forward activity
and a site of consequence.
```

The strong distance-from-readout trend is the first-order coordinate of that relation. The earlier distance-shell controls showed smaller but reproducible branch-specific corrections on top.

---

## Biological boundary

The biological parallels remain constraints, not identifications.

Known CA1 biology provides three relevant observations:

1. HCN channels are strongly enriched in distal dendrites;
2. their inductive gradient can compensate location-dependent timing at the soma;
3. HCN1 dendritic localization is activity dependent and reversible.

The model now independently reproduces a computational reason for a readout-relative material gradient and shows that local adjoint sensitivity can generate one.

It has **not** shown that biological HCN trafficking computes this adjoint or optimizes the present circular-coherence objective.

---

## Next wall

The morphology/material/readout branch is now sufficiently constrained to stop tuning it.

The next interesting return to theta/gamma should be hostile from the start:

- distinguish slow/theta spectral content already present in a gamma burst envelope from genuinely generated cross-frequency components;
- compare linear field, local nonlinear/synaptic conversion, and soma quadratic mixing;
- ask whether the confirmed quasi-active material preferentially transfers the slow consequential component without pretending that simple beat-frequency algebra is a biological discovery.

## Wall sentence

> **Uniform electrical material can self-organize under a local forward×transpose rule into a consequence-centered spatial impedance field: the dominant learned coordinate is distance from the readout, it re-centers when consequence moves, and it outperforms both shuffled material and a hand-drawn distal gradient on held-out arbors.**
