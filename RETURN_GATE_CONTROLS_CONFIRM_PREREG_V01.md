# Return-gating hostile controls — held-out preregistration v0.1

## Status

Frozen before running fresh bodies **516–527**.

Development bodies 510–515 showed that the strong 50%-duty return-gating effect is not explained by generic temporal subsampling. Regular fast interleaving preserved the exact gradient direction far better than random masks or one contiguous half-window at the same duty cycle.

This confirmation tests that distinction.

## Frozen setup

Use [`return_gate_controls.py`](return_gate_controls.py) unchanged:

- fresh FunctionalArbor bootstrap bodies: seeds `516..527`;
- lag `20`, `T=210`;
- exact target/distractor soma derivative waveforms;
- every gated waveform separately L2-dose-matched to its exact parent;
- all masks have 50% duty;
- periodic square-wave masks at `P=2,6,10,14,30,42,70`;
- all phase offsets tested for each periodic mask;
- 24 independent random exactly-half masks per body;
- 24 circular contiguous half-window offsets per body;
- no learning.

Simulation frame periods are **not** mapped to biological Hz.

## Registered criteria

### G0 — fast regular interleaving is nearly transparent

For `P=6`, across bodies require

```text
mean of phase-mean correlations > .999
mean of phase-min correlations  > .999
worst body phase-min correlation > .995
```

### G1 — P14 remains robust but is measurably below the fastest gates

Require

```text
P14 mean of phase-mean corr > .995
P14 mean of phase-min corr  > .990

P6 mean corr - P14 mean corr > .0005
```

The final inequality prevents the test from collapsing to "all distributed masks work identically."

### G2 — periodic P14 beats random 50% subsampling

Require

```text
P14 mean-of-mean corr - random mean-of-mean corr > .008
P14 mean-of-min  corr - random mean-of-min  corr > .05
```

The second comparison is the stronger hostile control: random masks may be good on average but should show much worse unlucky realizations.

### G3 — random distributed subsampling beats a contiguous half-window

Require

```text
random mean-of-mean corr - block mean-of-mean corr > .15
random mean-of-min  corr - block mean-of-min  corr > .30
```

This tests whether temporal coverage itself matters.

### G4 — regular gate period orders the degradation regime

Require the group mean correlations to satisfy

```text
P6 > P14 > P30 > P42 > P70
```

and additionally

```text
P14 - P42 > .10
P42 - P70 > .05
```

This is a scale-separation criterion, not a claim of universal monotonicity for every possible period.

### G5 — the alternating comb is not merely equivalent to random half sampling

`P=2` is the extreme regular comb. Require

```text
comb mean-of-mean corr - random mean-of-mean corr > .01
comb mean-of-min  corr - random mean-of-min  corr > .05
```

## Interpretation boundary

A pass supports:

> **For this finite-time reciprocal gradient map, equal-duty temporal masks are not equivalent. Fast regularly interleaved return support preserves structural direction far better than random or contiguous support, and performance degrades as the periodic gate approaches the task-sensitive temporal scale.**

It does not by itself prove a Nyquist theorem, identify the exact relevant bandwidth, or establish a theta/gamma mechanism.

If confirmed, the next experiment must work in frequency space and predict gate performance from mask sidebands and the return-to-gradient sensitivity spectrum.