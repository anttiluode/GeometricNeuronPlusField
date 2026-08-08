# Quasi-active HCN-gradient confirmation preregistration v0.1

Date: 2026-08-08

This file freezes the candidate and acceptance criteria **before** examining seeds 568–579.

## Why this exists

The first morphology-indexed delayed-restorative proxy in `hcn_impedance_probe.py` reduced somatic phase spread but failed the biological mechanism diagnostic: the distance-dependent phase shift had the wrong sign for an HCN-like inductive compensation story.

The correction was to use the complete minimal quasi-active structure motivated by linearized voltage-dependent current theory: one channel-density field contributes both a static/resting conductance and a delayed restorative component.

The model in `hcn_quasiactive_probe.py` is

```text
v' = K L psi - damping*v - restoring*psi
     - d(x)*psi - mu*d(x)*z + source

z' = (psi-z)/tau_h
```

with `mu > 0` restorative.

This is still a **second-order wave toy plus quasi-active membrane material**, not a conductance-based HCN neuron and not evidence that CA1 neurons implement this exact equation.

## Frozen candidate

Development bodies used for candidate selection: seeds 560–563.

The frozen candidate is **not** the highest aggregate development score. It was selected for body-by-body robustness.

```text
g0                  0.005
soma->distal ratio  10
channel tau         2
mu/static ratio     0.5
frequencies         omega = 0.03, 0.04
```

The smooth density field is linear in graph distance from the soma and reaches `10 * g0` at the maximum occupied graph distance.

Controls use the same frozen arbor geometry:

1. `zero`: no added quasi-active material;
2. `uniform`: same mean density everywhere;
3. `smooth`: soma-to-distal density gradient;
4. `shuffle`: identical smooth-profile density values shuffled across occupied cells;
5. `reverse`: identical density histogram assigned in reverse graph-distance order.

The shuffled control therefore asks whether **morphology-indexed placement**, rather than merely the amount or histogram of quasi-active material, matters.

## Development receipt motivating the freeze

Across four development bodies × two frequencies = eight observations:

```text
smooth beats shuffled soma phase spread       8 / 8
smooth beats uniform soma phase spread        8 / 8
distal-minus-proximal phase advance > 0       8 / 8

pooled mean gain vs shuffled                  +0.1757 rad
pooled mean gain vs uniform                   +0.1559 rad
pooled mean distal-minus-prox phase advance   +1.5399 rad
pooled mean local phase retention              1.170
pooled mean smooth/shuffle amplitude ratio     0.961
```

Frequency-level development values:

```text
omega=.03
  gain vs shuffle       +0.2211 rad
  gain vs uniform       +0.1623 rad
  distal phase lead     +2.3557 rad
  local phase retention  1.542

omega=.04
  gain vs shuffle       +0.1303 rad
  gain vs uniform       +0.1494 rad
  distal phase lead     +0.7242 rad
  local phase retention  0.769
```

These values are reported to make the selection history explicit. The held-out thresholds below are deliberately much weaker.

## Held-out set

Fresh frozen FunctionalArbor bodies:

```text
seeds 568–579
12 requested bodies
24 body-frequency observations if all bootstrap successfully
```

No threshold will be changed after those results are observed.

## Primary quantities

For each body and frequency, inject the same unit harmonic drive independently at many occupied arbor locations and measure the complex transfer at the soma.

### Somatic phase spread

```text
sigma_soma = circular RMS of soma transfer phase over injection locations
```

Smaller means location-dependent inputs arrive in a more coordinated somatic phase relation.

Define

```text
gain_shuffle = sigma_shuffle - sigma_smooth
gain_uniform = sigma_uniform - sigma_smooth
gain_reverse = sigma_reverse - sigma_smooth
```

Positive is favorable to the morphology-indexed gradient.

### Distal phase-advance sign

For every injection location,

```text
advance(x) = wrapped_phase(H_smooth(x)) - wrapped_phase(H_zero(x))
```

Compare the most distal distance quartile with the most proximal quartile:

```text
lead = mean_distal advance - mean_proximal advance
```

The HCN-like mechanism prediction is `lead > 0`: distal inputs receive more phase advance than proximal inputs.

### Local phase retention

```text
retention = local_phase_RMS_smooth / local_phase_RMS_zero
```

This guards against the cheap solution in which the intervention simply flattens the local dendritic field everywhere.

### Amplitude control

Use the median soma transfer amplitude over injection sites. The primary amplitude ratios are

```text
A_smooth / A_shuffle
A_smooth / A_uniform
```

The synchronization effect should not depend on extreme amplification or attenuation.

## Frozen criteria

### Q0 — morphology-indexed placement beats shuffled and uniform material

Across all body-frequency observations:

```text
mean gain_shuffle > 0.05 rad
mean gain_uniform > 0.05 rad
```

### Q1 — placement advantage is body-level robust

Out of 24 expected observations:

```text
gain_shuffle > 0 in at least 18
and
gain_uniform > 0 in at least 18
```

If fewer than 24 observations are available because a body fails bootstrap, use a required fraction of at least `0.75` for each control.

### Q2 — the correction has the HCN-like distance sign

```text
mean distal-minus-proximal lead > 0.30 rad
lead > 0 in at least 18 / 24 observations
```

With fewer observations, require positive lead in at least `0.75` of them.

### Q3 — the local field is not globally phase-flattened

```text
mean local phase retention > 0.70
```

This criterion does not require local phase spread to increase; it only rejects severe collapse as the mechanism of apparent somatic synchrony.

### Q4 — amplitude is not the dominant explanation

Across all observations, both pooled median ratios must satisfy

```text
0.5 < median(A_smooth / A_shuffle) < 2.0
0.5 < median(A_smooth / A_uniform) < 2.0
```

### Q5 — both members of the frozen low-frequency pair contribute

At **each** of `omega=.03` and `.04`:

```text
mean gain_shuffle > 0.02 rad
mean gain_uniform > 0.02 rad
mean distal-minus-proximal lead > 0.10 rad
```

This rejects a result driven by one isolated frequency.

### Q6 — gradient direction matters

Pooled over both frequencies:

```text
mean(sigma_reverse - sigma_smooth) > 0.05 rad
```

The reverse profile contains exactly the same density values, so this is a directional morphology/material control.

## Confirmation rule

The model-level quasi-active HCN bridge is called **confirmed on this held-out test only if all Q0–Q6 pass**.

A failure is retained as a failure. In particular:

- strong somatic synchronization with the wrong distal phase sign is not an HCN-like positive;
- strong synchronization caused by severe local phase collapse is not an HCN-like positive;
- an effect that does not beat shuffled same-material placement is not a geometry-indexed material positive.

## What a pass would mean

A pass would support the following narrow computational statement:

> In this frozen wave-arbor model, a morphology-indexed distribution of a complete minimal quasi-active restorative material can compensate location-dependent transfer phase at the soma better than uniform, shuffled, or reversed placement while preserving substantial local phase structure.

It would **not** show that biological HCN channels optimize this model objective, that the numerical parameter ratios correspond directly to CA1 channel densities, or that neurons use physical adjoint learning.

## Next experiment if confirmed

Do not hand-tune another spatial gradient.

The next experiment should make `d(x)` trainable under a fixed total material budget and ask whether local forward × returned sensitivity can discover a distal-enriched profile from an initially uniform material distribution.

That experiment is motivated by the known activity dependence of CA1 HCN1 spatial distribution, but the learning rule must earn its own result rather than being assumed from biology.
