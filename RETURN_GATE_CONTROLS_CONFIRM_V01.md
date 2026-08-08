# Return-gating hostile controls v0.1 — regular temporal interleaving is load-bearing

## Result

The strong fast-gating effect is **not** explained by generic 50% temporal subsampling.

Development bodies 510–515 motivated [`RETURN_GATE_CONTROLS_CONFIRM_PREREG_V01.md`](RETURN_GATE_CONTROLS_CONFIRM_PREREG_V01.md). The frozen criteria were then tested on fresh bodies **516–527**.

```text
6 / 6 registered criteria PASS
```

The result is unusually orderly:

```text
periodic P=2       0.999969 mean corr
periodic P=6       0.999886
periodic P=10      0.999708
periodic P=14      0.998405
periodic P=30      0.972457
periodic P=42      0.851114
periodic P=70      0.720260

random 50%         0.980956
contiguous 50%     0.652009
```

All masks retain half of the soma return samples and are L2-dose matched to the exact return waveform.

So **where the retained samples occur in time is load-bearing**.

## G0 — fast regular interleaving is nearly transparent

For P=6:

```text
mean phase-mean corr       0.99988615
mean phase-min corr        0.99971675
worst-body phase-min       0.99928599
mean strong-sign           0.99957711
```

**PASS.**

## G1 — P14 remains robust but is measurably below P6

```text
P14 mean corr              0.99840530
P14 mean phase-min         0.99681378
worst-body phase-min       0.98616854
P6 - P14 mean              +0.00148085
```

**PASS.**

## G2 — periodic P14 beats random 50% subsampling

Random exactly-half masks:

```text
mean-of-mean corr          0.98095614
mean-of-min corr           0.90927539
worst-body random min      0.83052421
```

Registered P14 advantages:

```text
mean advantage             +0.01744916
min advantage              +0.08753839
```

**PASS.**

Random masks are often good, confirming substantial redundancy in the return waveform, but unlucky random masks are much less reliable than regular interleaving.

## G3 — distributed random samples beat one contiguous half-window

Contiguous circular half-window:

```text
mean-of-mean corr           0.65200919
mean-of-min corr           -0.02213074
worst-body min             -0.66922635
```

Random minus block:

```text
mean advantage              +0.32894696
min advantage               +0.93140612
```

**PASS.**

So simple sample count is nowhere near enough. Coverage across the temporal episode matters enormously.

## G4 — periodic gate period orders the degradation regime

Held-out group means:

```text
P6      0.999886
P14     0.998405
P30     0.972457
P42     0.851114
P70     0.720260
```

The preregistered strict ordering passed.

Additional differences:

```text
P14 - P42    +0.14729141
P42 - P70    +0.13085419
```

**PASS.**

The simulation has therefore exposed a real temporal-scale boundary. No biological frequency calibration is implied.

## G5 — the alternating comb is much better than random half sampling

The extreme regular comb, retaining every other return sample, gives

```text
mean corr                  0.99996924
mean phase-min             0.99996901
worst-body min             0.99992291
```

Comb minus random:

```text
mean advantage             +0.01901309
min advantage              +0.09069362
```

**PASS.**

## What this rules out

The earlier result cannot be reduced to

```text
the return waveform is redundant, so any half is enough.
```

That prediction is false.

Nor is the effect merely choosing a lucky phase of a periodic gate: every phase was scanned, and the fast-period worst cases remain extremely strong.

## The natural signal-processing explanation

Let `A` denote the linear map from a soma return code `g(t)` to the bond-gradient map.

For a 50%-duty periodic gate

```text
m(t) = 1/2 + r(t)
```

where `r` has zero mean. Ignoring the harmless dose-normalization scalar for the moment,

```text
A[m g]
  = 1/2 A[g] + A[r g].
```

The first term is a scaled copy of the correct structural direction.

The second is the error term.

In frequency space, multiplication by periodic `r(t)` shifts copies of the return spectrum by the gate harmonics. A fast regular gate has widely separated harmonic lines. A random mask has broadband spectral content. A long contiguous block has strong low-frequency mask content.

That gives a natural candidate mechanism for the observed ordering:

```text
fast regular gate
    -> scaled desired/baseband map
       + shifted sidebands mostly outside the consequential band

random gate
    -> spectral contamination broadly distributed

slow/block gate
    -> distortion overlaps the low/task-relevant temporal structure directly
```

This is standard modulation/sampling mathematics. What remains to be established is whether the measured return-to-gradient sensitivity spectrum quantitatively predicts the observed break between P14, P30, P42 and P70.

## Why this is interesting for the rhythm question

The result does **not** identify gamma, theta, chandelier cells, or any biological oscillator.

It establishes a narrower functional property:

> **A fast regular temporal gate can multiplex a consequence-return waveform while preserving almost all of its structural direction, whereas equally sparse unstructured or slow support cannot.**

That is exactly the kind of job for which nested biological rhythms are worth considering — temporal separation and reliable interleaving — without pretending the rhythm carries the gradient itself.

## Next wall

Measure the return-to-gradient temporal sensitivity spectrum and predict each gate from its sidebands.

A successful mechanism should explain, before seeing gated-map performance, why

```text
P2/P6/P10/P14   are nearly transparent
P30             begins to hurt
P42/P70         strongly interfere
```

and why random and contiguous masks fall where they do.

Only after that should the spatially varying dendritic phase experiment be interpreted.

## Wall sentence

> **The return code is not merely redundant; it is regularly subsampleable. Fast temporal interleaving preserves its structural direction almost exactly, while random and slow masks contaminate the consequential temporal band.**
