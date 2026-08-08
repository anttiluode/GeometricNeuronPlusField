# Coarse biological return-code confirmation preregistration v0.1

## Status

Frozen before running fresh bodies **498–509**.

Development bodies 492–497 asked how much the exact reciprocal-adjoint soma waveform could be degraded toward a biologically coarse return event while preserving the local bond-gradient direction.

All transformed return codes are L2-dose-matched separately for target and distractor. The test therefore concerns **information in the temporal return code**, not total return energy.

This is a model/mechanism experiment, not a claim that neurons implement an adjoint.

## Frozen setup

- FunctionalArbor bootstrap bodies, no learning;
- fresh seeds `498..509`;
- lag `20`, `T=210`;
- same passive reciprocal return operator as `RECIPROCAL_ADJOINT_V01`;
- only the soma derivative waveform is transformed;
- metrics: bond-map correlation, relative L2, and strong-sign agreement.

## Return codes

### Exact

Full complex derivative waveform. Positive control.

### Real wave

```text
Re[g(t)]
```

renormalized to the same L2 dose. This removes the explicit complex quadrature while preserving the real temporal voltage-like waveform.

### Phase only

All samples are given constant amplitude while retaining their complex phase.

### Envelope + task sign

```text
sign(task coefficient) * |g(t)|
```

This preserves the temporal amplitude envelope and one global target/distractor sign but discards carrier phase.

### Sparse phase-bearing events

Keep the N strongest separated temporal peaks, make their amplitudes equal, and retain each selected sample's complex phase. Primary sparse condition: `N=32`.

### Sparse signed events

Same N peak times, equal real amplitudes, but retain only the target/distractor sign.

### Sparse positive events

Same peak times and equal positive amplitudes; even the target/distractor sign is removed. This is a destructive control.

### Periodic gate

Retain the exact return waveform only during 50%-duty periodic windows, renormalize the surviving waveform to the original L2 dose, and scan every gate phase offset.

Primary periods:

```text
P=14 frames   fast gate
P=42 frames   slow gate
```

No biological frequency calibration is implied by those frame counts.

## Registered criteria

### R0 — exact positive control

Require

```text
mean exact corr > .999999
mean exact relative L2 < 1e-10
```

### R1 — real-valued return preserves structural direction

Require

```text
mean real-wave corr > .995
mean strong-sign agreement > .99
```

This tests whether explicit complex quadrature is necessary for the useful return-map direction in this model.

### R2 — amplitude envelope + consequence sign remains informative

Require

```text
mean envelope-signed corr > .80
median envelope-signed corr > .80
>= 9/12 bodies corr > .75
```

### R3 — fast 50%-duty gating is nearly transparent

For each body scan all 14 gate offsets. Require

```text
mean across bodies of median-offset corr > .99
mean across bodies of worst-offset corr  > .98
mean median-offset strong-sign agreement > .96
```

This is deliberately phase-hostile: the result must not depend on choosing a favorable gate phase.

### R4 — slower gating exposes phase sensitivity

Require the fast gate to outperform the slow P=42 gate by

```text
mean(median_corr_P14 - median_corr_P42) > .03
mean(worst_corr_P14 - worst_corr_P42)   > .15
```

The claim is scale separation, not monotonic frequency preference.

### R5 — 32 sparse phase-bearing events retain substantial direction

Require

```text
mean sparse_phase_32 corr > .85
>= 8/12 bodies corr > .80
```

This is intentionally weaker than the real-wave criterion.

### R6 — consequence sign is load-bearing in the sparse limit

Require

```text
mean sparse_signed_32 corr - mean sparse_positive_32 corr > .50
abs(mean sparse_positive_32 corr) < .25
```

This asks whether event timing alone can substitute for the objective sign.

### R7 — envelope information beats phase-only coding

Require

```text
mean(envelope_signed corr - phase_only corr) > .08
```

This was a development clue that the return code is more tolerant to losing carrier phase than to flattening its temporal envelope.

## Kill conditions

The biological bridge is narrowed if:

- real-only return fails to preserve the gradient direction;
- fast periodic gating works only at specially selected phases;
- sparse event timing works even after task sign is removed, implying the positive control was not discriminating;
- or the envelope/phase ordering does not replicate.

## Claim boundary

A pass would support only:

> **Within the reciprocal passive-arbor model, the exact analog adjoint waveform contains substantial redundancy. Useful structural direction can survive removal of explicit quadrature, 50%-duty temporal gating, and some event sparsification, while the sign of consequence remains load-bearing.**

It would not show that biological back-propagating action potentials encode an adjoint or that theta/gamma rhythms implement these gates.