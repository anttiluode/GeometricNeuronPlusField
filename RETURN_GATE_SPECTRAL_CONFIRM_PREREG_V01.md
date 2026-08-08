# Return-gate spectral mechanism — held-out preregistration v0.1

## Status

Frozen before running fresh bodies **540-551**.

Development bodies 534-539 found that sideband contamination landing in the previously established boundary-selected K=8/K=16 frequency bins strongly predicts the **masking regime** produced by fast periodic, slow periodic, random, and contiguous 50%-duty return masks.

The development result was much weaker for individual mask realizations. This confirmation therefore preregisters a regime-level claim and only a modest individual-level criterion.

## Frozen setup

Use [`return_gate_spectral_mechanism.py`](return_gate_spectral_mechanism.py) unchanged:

- FunctionalArbor bootstrap bodies;
- seeds `540..551`;
- lag `20`, `T=210`;
- K=8 and K=16 boundary-selected bins chosen before each gated internal map is evaluated;
- periodic masks `P=2,6,10,14,30,42,70`, all offsets;
- 12 random exactly-half masks per body;
- 12 contiguous circular half-window masks per body;
- every gated target/distractor return separately L2-dose matched to its exact parent.

For each mask

```text
m(t) = 1/2 + r(t)
```

and the registered predictor is the weighted ratio of

```text
FFT[r(t) g(t)]
```

sideband energy landing in the selected port bins to the desired

```text
1/2 G
```

energy in those bins.

## Registered criteria

### S0 — periodic regime damage is predicted by K8 contamination

Across the seven periodic period means `P=2,6,10,14,30,42,70`, require

```text
corr(K8 mean contamination, 1 - mean map corr) > .95
```

### S1 — the result is not unique to K8

For K=16 periodic regime means require

```text
corr(K16 mean contamination, 1 - mean map corr) > .95
```

### S2 — the broad mask classes are also ordered at the boundary

Across the nine category means

```text
P2,P6,P10,P14,P30,P42,P70,random,block
```

require for both K=8 and K=16

```text
corr(contamination, 1 - map corr) > .90
```

### S3 — fast/slow separation is large in the predictor itself

Require

```text
K8 contamination P42 / P6 > 10
K16 contamination P42 / P6 > 10

map corr P6 - map corr P42 > .08
```

### S4 — a contiguous half-window is spectrally dirtier than random half sampling

Require for both K values

```text
block mean contamination > 3 * random mean contamination
```

and behaviorally

```text
random mean map corr - block mean map corr > .15.
```

### S5 — individual-mask prediction remains only moderate but nonzero

Pool every individual periodic offset, random draw, and block draw across all bodies.

Require

```text
K8 corr(contamination, 1-map_corr) > .30
K16 corr(contamination, 1-map_corr) > .30
```

This deliberately does **not** demand the `.50` development value. A pass supports only partial realization-level prediction.

## Kill / narrowing conditions

The boundary-spectrum mechanism is narrowed if:

- periodic period ordering is not captured on fresh bodies;
- K8 works but K16 does not, suggesting an accidental ranking choice;
- block/random differences are not visible in the port contamination score;
- or the individual-level correlation collapses to near zero.

Even a full pass does not mean the K-bin score reconstructs the exact gated gradient. It means it predicts the temporal **regime** in which masking starts to corrupt the structural direction.

## Claim boundary

A pass supports:

> **The temporal-rate boundary of return-wave multiplexing can be estimated from soma-port spectra: fast masks preserve the gradient when their sidebands avoid the compact frequency set already associated with useful transient-gradient information, while slow/block masks contaminate that set.**

No claim is made that this exact spectral rule is implemented by theta/gamma biology.