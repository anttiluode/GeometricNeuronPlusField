# Return-gate spectral mechanism — development v0.1

## Question

Why do equal-duty return masks behave so differently?

Held-out controls already established:

```text
fast regular interleaving    nearly exact gradient direction
random 50% masks             moderately worse
slow periodic masks          much worse
one contiguous half-window   often disastrous
```

This development probe asks whether that regime ordering can be predicted **at the soma ports**, without inspecting the gated internal gradient first.

Bodies: seeds **534-539**.

## Boundary decomposition

For an exact soma return `g(t)` and a 50% mask `m(t)`, write

```text
m(t) = 1/2 + r(t).
```

Therefore

```text
FFT[m g]
    = 1/2 G
      + FFT[r g].
```

The first term is a scaled copy of the desired return spectrum.

The second term is the sideband contamination created by the non-DC part of the mask.

Previous held-out work had already identified a common small set of **boundary-selected frequency bins** whose source/return spectral products preserve most of the full internal gradient direction.

So before evaluating each gated gradient, this probe measured

```text
contamination ratio
    = sideband energy landing in those selected bins
      / desired half-amplitude energy in those bins.
```

K=8 and K=16 versions were tested.

## Development result — the regime ordering is visible at the boundary

Group means for K=8:

```text
mask          map corr      selected-bin contamination
P2            0.999987             0.00334
P6            0.999950             0.00560
P10           0.999869             0.00893
P14           0.999513             0.01980
P30           0.965404             0.09349
P42           0.853237             0.24577
P70           0.734473             0.54657
random        0.974343             0.12405
block         0.689126             0.81942
```

K=16 contamination shows the same qualitative ordering:

```text
P2      0.00356
P6      0.00601
P10     0.01278
P14     0.03358
P30     0.19821
P42     0.31849
P70     0.56772
random  0.13615
block   0.82557
```

Across the nine mask regimes, K=8 contamination correlates with `1-map_corr` at approximately

```text
r = 0.979
```

and across only the seven periodic gate periods at approximately

```text
r = 0.993.
```

So the previously qualitative sideband story has a strong quantitative development result:

> **the mask begins to damage the structural direction when its non-DC modulation products are driven back into the boundary frequency bins that already carried the useful transient-gradient information.**

## But the predictor is not a complete per-mask model

Pooling every individual mask realization rather than regime means gives only moderate prediction:

```text
K=8
contamination vs 1-corr        r = 0.501
contamination vs relative L2   r = 0.495

K=16
contamination vs 1-corr        r = 0.505
contamination vs relative L2   r = 0.500
```

That distinction matters.

The simple K-bin boundary score is good at predicting

```text
which temporal masking regime is dangerous
```

but not

```text
exactly how much one particular random mask will damage one particular body.
```

That missing variance may live in:

- important bins outside the K=8/K=16 compression;
- relative phase among contaminated bins;
- target/distractor-specific cancellation;
- the spatial pattern of each frequency contribution;
- or nonlinear effects of separate L2 dose matching of the two task returns.

Do not claim more than the result supports.

## Interpretation

The confirmed fast-gating result now has a concrete signal-processing picture.

```text
FAST REGULAR MASK
m = 1/2 + narrow harmonic comb
          |
          +---- shifted copies mostly miss important return bins
          v
scaled desired boundary spectrum survives
          v
same structural direction

SLOW / BLOCK MASK
non-DC spectrum moves toward low / task-relevant offsets
          |
          +---- shifted copies contaminate important return bins
          v
gradient direction changes

RANDOM MASK
broadband mask spectrum
          |
          +---- moderate contamination spread across many bins
          v
usually good, but less reliable
```

This is ordinary modulation/sampling mathematics applied to the transient adjoint return. The potentially useful architectural result is that the relevant bandwidth can be estimated from **boundary signals** rather than from omniscient inspection of every internal bond.

## Next test

Freeze the K-bin port ranking and contamination score, then use fresh bodies to predict the group ordering of periodic, random, and block masks.

The held-out claim should remain a **regime predictor**, not a per-realization predictor.

## Wall sentence

> **Fast return gating works because its modulation sidebands mostly miss the same compact boundary spectrum that carries useful gradient information. The boundary spectrum predicts the masking regime very well, but not every individual mask.**
