# Return-gate spectral mechanism v0.1 — held-out boundary prediction

## Result

The temporal-rate boundary of return-wave multiplexing is strongly visible in **boundary spectra**.

Development bodies 534–539 motivated [`RETURN_GATE_SPECTRAL_CONFIRM_PREREG_V01.md`](RETURN_GATE_SPECTRAL_CONFIRM_PREREG_V01.md). The frozen criteria were tested on fresh bodies **540–551**.

```text
6 / 6 registered criteria PASS
```

The useful claim is deliberately narrower than an exact per-mask model:

> **sideband contamination landing in the same compact soma-port frequency set that already predicted the transient gradient strongly predicts which masking regime will damage the structural direction.**

## Boundary decomposition

For exact return `g(t)` and 50%-duty mask `m(t)`:

```text
m(t) = 1/2 + r(t)
```

therefore

```text
FFT[m g]
    = 1/2 G
      + FFT[r g].
```

The first term is the desired return spectrum scaled by one half.

The second term is the modulation/sideband contamination.

Before evaluating each gated internal gradient, the probe measured how much of that contamination landed in the previously established boundary-selected K=8 or K=16 frequency bins.

## Fresh periodic ladder

Held-out group means:

```text
period   map corr    K8 contamination   K16 contamination
P2       .999959       .00377              .00410
P6       .999843       .00639              .00702
P10      .999548       .01068              .01586
P14      .998406       .02450              .03782
P30      .933305       .14617              .24593
P42      .847204       .26827              .36689
P70      .710651       .52658              .56873
```

The break in structural fidelity and the rise in boundary contamination track each other closely.

## S0 — K8 predicts periodic-regime damage

Across P2/P6/P10/P14/P30/P42/P70:

```text
corr(K8 contamination, 1-map_corr)
    = 0.99856346
```

Registered requirement: `> .95`.

**PASS.**

## S1 — K16 independently predicts the periodic regime

```text
corr(K16 contamination, 1-map_corr)
    = 0.98181385
```

**PASS.**

The result is therefore not peculiar to one exact K=8 ranking.

## S2 — periodic, random and block classes are jointly predicted

Add the two hostile control classes:

```text
random
map corr            .978209
K8 contamination    .137944
K16                 .154869

contiguous block
map corr            .743307
K8 contamination    .810746
K16                 .823912
```

Across all nine categories:

```text
K8  corr(contamination, 1-map_corr) = .933449
K16                                      .949639
```

Registered requirement for both: `> .90`.

**PASS.**

## S3 — fast/slow separation is already large at the ports

```text
K8  P42 / P6 contamination ratio    41.98 x
K16 P42 / P6                        52.28 x

map corr P6 - P42                   +0.15264
```

**PASS.**

## S4 — block versus random is also visible spectrally

```text
block/random contamination ratio
K8     5.88 x
K16    5.32 x

random map corr - block map corr
       +0.23490
```

**PASS.**

So the fact that a contiguous half-window is much worse than random half-sampling is not hidden in the interior of the mesh. It is already reflected in what the mask does to the consequential boundary spectrum.

## S5 — individual-mask prediction remains partial

Pooling every periodic phase, random draw, block draw, body and class:

```text
K8  corr(contamination, 1-map_corr)   .464386
K16                                    .470639
```

Registered minimum: `.30`.

**PASS.**

But this is the important limitation.

The K-bin boundary contamination score does **not** explain every individual mask well. It explains the rate/class regime much better than realization-level variance.

That missing variance is real and should not be hidden.

## What is established now

The earlier statement

```text
fast regular return gating works
```

can now be sharpened to

```text
fast regular return gating
    -> modulation sidebands stay mostly outside
       the compact port spectrum associated with useful gradient information
    -> structural direction survives

slow / block gating
    -> modulation sidebands invade that port spectrum
    -> structural direction degrades
```

Random masks occupy an intermediate position because their spectral contamination is broad but not concentrated into one long low-frequency interruption.

## Why this matters architecturally

This makes the temporal-multiplexing problem locally measurable.

A device does not necessarily need to inspect every internal coupling to decide how fast a return/error channel must be chopped or interleaved.

The already available boundary signals provide an estimate of the consequential bandwidth.

That suggests a hardware protocol:

```text
1. estimate compact task-relevant return spectrum at ports
2. choose a modulation/interleaving rhythm whose sidebands avoid it
3. send the return/adjoint channel through that temporal schedule
4. accumulate local structural overlap
```

This is ordinary modulation theory, but it is useful because the distributed mesh itself remains the processor and gradient router.

## Relation to the biological rhythm question

This result does **not** establish that gamma is an error channel or that theta/gamma frequencies are selected by this rule in neurons.

It changes the form of the biological hypothesis:

```text
not:
    fast rhythm carries a gradient

but possibly:
    fast rhythm temporally multiplexes/gates processes
    while keeping their modulation products away from
    the frequencies in which a dendritic/somatic computation is sensitive
```

That is much closer to an ordinary reason for nested rhythms to exist.

The biology has independent evidence for theta/gamma temporal organization and for frequency-dependent dendritic impedance, but the mapping remains open.

## Wall sentence

> **The rhythm is not the message. Its useful job can be spectral housekeeping: interleave a return/consequence channel fast enough that the modulation sidebands miss the compact boundary spectrum carrying the structural information.**
