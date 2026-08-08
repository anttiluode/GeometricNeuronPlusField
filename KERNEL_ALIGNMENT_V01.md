# Point-kernel alignment v0.1 — soma privilege is mostly balance, not a uniquely aligned kernel

This report follows directly from `MODE_PAIR_V01.md`.

The question was:

> **Why is the soma's quadratic point-readout kernel better aligned with the source interaction than the kernels of most other cells?**

The answer is more nuanced than the question assumed: **the soma is not unusually special in normalized task-kernel alignment by itself.** Its location privilege is mainly the amplitude-balance/convergence effect already identified by `SomaWhyClaude`. A separate interaction-alignment factor modifies the body-wide selectivity map and can improve a fixed balance-only predictor, but that incremental improvement did not meet the strict held-out sign-test criterion.

## Geometry-only point kernel

For every occupied cell `x`, the complete graph eigenbasis gives a row vector

```text
v_x = [phi_0(x), ..., phi_{N-1}(x)]
```

and point power is the quadratic form

```text
y_x = q^* K_x q
K_x = v_x v_x^T.
```

The two source histories define at every time an order-sensitive cross-source operator

```text
D(t) = H_T(t) - H_D(t)
```

where `H_T` and `H_D` are the Hermitian target/distractor source-interaction matrices.

The local cross-source order difference is exactly

```text
Delta_cross_x(t) = Re[v_x^T D(t) v_x].
```

To remove the global instantaneous strength of `D(t)`, the preregistered normalized alignment score was

```text
A_x = RMS_t( Delta_cross_x(t) / ||D(t)||_F ).
```

Amplitude balance retained the previous definition

```text
B_x = min(P_A,P_B) / max(P_A,P_B)
```

and the fixed, coefficient-free combined score was

```text
Q_x = B_x * A_x.
```

## Discovery — seeds 108-119

Registered predictions D1 and D2 failed; D3 passed.

```text
mean corr(A, |C|)                 0.40885
positive r_A bodies               12 / 12
D1 threshold mean > .50           FAIL

median soma A percentile          0.66429
soma A above median                7 / 12
D2                                 FAIL

mean corr(B, |C|)                 0.76830
mean corr(B*A, |C|)               0.84741
mean improvement                  0.07911
improved / worse                  11 / 1
sign p                            0.00635
D3                                 PASS
```

The soma itself was extreme in balance but only ordinary-to-moderately-high in normalized alignment:

```text
mean soma percentile
A alignment                        0.64286
B balance                          0.96429
B*A combined                       0.92619
|C| task selectivity               0.89167
```

`corr(A,B)` averaged only `0.3755`, so the two maps are related but far from identical.

## Held-out confirmation — seeds 120-131

Before these bodies were run, the refined interpretation and four confirmation criteria were frozen in `KERNEL_ALIGNMENT_CONFIRM_PREREG_V01.md`.

```text
mean corr(A, |C|)                 0.40258
mean corr(B, |C|)                 0.72813
mean corr(B*A, |C|)               0.78898
mean corr(A,B)                    0.38898
```

### C1 — incremental product improvement

```text
mean (r_Q - r_B)                  +0.06085
positive / negative bodies         9 / 3
two-sided sign p                  0.145996
```

The mean-effect and count thresholds were met, but the preregistered sign-test criterion `p < .05` was not.

**C1 FAIL.**

This matters. The discovery result cannot be promoted to a universal complementary-factor law from this 12-body confirmation.

### C2 — fixed product remains a strong predictor

```text
mean corr(B*A, |C|)               0.78898
registered threshold              > 0.78
```

**C2 PASS.**

The coefficient-free product remains a strong map-level predictor, even though its incremental advantage over balance is not body-universal.

### C3 — combined score keeps soma high

```text
median soma B*A percentile        0.90714
somata above 0.75                 12 / 12
```

**C3 PASS.**

### C4 — alignment alone is not the soma privilege

```text
median soma A percentile          0.69286
somata above 0.75                  4 / 12
```

**C4 PASS.**

This is the key negative-structure confirmation. The soma is not sitting at an exceptionally high normalized `K_x`/`D(t)` alignment location.

## A useful descriptive clue

The *unnormalized* cross-difference RMS correlated with `|C|` much more strongly:

```text
discovery mean r                  0.74544
confirmation mean r               0.72407
```

That is expected because it contains both interaction orientation and interaction magnitude. The normalized `A` intentionally removes the global operator magnitude and therefore asks the stricter question about kernel orientation alone.

The modal/coordinate equality check was effectively exact (`~6e-15` relative), so these findings are not a projection-reconstruction artifact.

## Current interpretation

The previous question contained a false premise. It asked why the soma's quadratic kernel was unusually aligned.

The data say:

```text
soma/root convergence
       -> unusually strong A/B amplitude balance          YES

normalized task-kernel orientation A
       -> body-wide modifier of selectivity               YES, moderate
       -> uniquely large at soma                          NO

balance * alignment
       -> strong body-wide predictor                      YES
       -> universally better than balance alone           NOT CONFIRMED
```

So the most economical mechanism is still:

```text
geometry
  -> source-specific transfer histories
  -> convergence makes both histories locally available at comparable amplitude
  -> local quadratic readout mixes their geometry-shaped modes
  -> the exact temporal/complex compatibility of those histories determines the final order effect
```

Amplitude balance is the robust location-level explanation for why the root is privileged. Modal interaction structure explains what is computed once both inputs have leverage there.

## Why the simple product is only approximate

`B_x*A_x` was intentionally coefficient-free, but it is not the exact algebra of square-law interference. At a point,

```text
|a+b|^2 = |a|^2 + |b|^2 + 2 |a||b| cos(delta_phi).
```

The natural amplitude factor is therefore an interference-visibility term such as

```text
V = 2 |a||b| / (|a|^2 + |b|^2)
```

rather than `min/max`, and the temporal factor is the lag-dependent relative complex relation, not the RMS kernel score alone.

That suggests the next experiment should stop multiplying two descriptive maps and test the actual local interference factorization directly.

## Next clean question

> **Can temporal-order selectivity be predicted by an explicit interference-visibility × lagged-complex-compatibility decomposition of the two single-source transfer histories, without any fitted coefficients?**

If yes, the amplitude-balance result and the mode-pair result collapse into the standard local algebra of coherent superposition, with geometry's role being to shape the two transfer histories before they meet.
