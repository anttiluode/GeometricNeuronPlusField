# Kernel-alignment held-out confirmation preregistration v0.1

Discovery on fresh seeds 108-119 gave a mixed but informative result:

```text
mean corr(alignment A, |C|)          0.4089     discovery D1 FAIL
median soma A percentile             0.6643     discovery D2 FAIL
mean corr(balance B, |C|)            0.7683
mean corr(B*A, |C|)                  0.8474
mean improvement (B*A)-B             0.0791
improved bodies                      11 / 12
discovery D3 PASS, sign p            0.00635
```

So the refined hypothesis is **not** that the soma is privileged because its normalized quadratic kernel is exceptionally aligned by itself.

Instead:

> **Amplitude balance supplies the soma's location privilege; task-kernel alignment supplies a partially independent temporal-compatibility factor. Their fixed product predicts the selectivity landscape better than balance alone.**

No discovery body will be reused in the confirmation.

## Held-out bodies

```text
seeds 120-131
```

All simulation parameters, lag, steps, graph basis, definitions of `A`, `B`, `Q=B*A`, and historical `C` remain exactly those in `kernel_alignment_probe.py` and `KERNEL_ALIGNMENT_DISCOVERY_PREREG_V01.md`.

No fitted coefficients are allowed.

## Registered confirmation criteria

### C1 — complementary-factor improvement replicates

For each body:

```text
Delta_r = corr(B*A, |C|) - corr(B, |C|).
```

PASS if all are true:

```text
mean Delta_r > 0.04
positive bodies >= 9/12
two-sided sign-test p < 0.05
```

### C2 — the fixed product is a strong body-wide predictor

PASS if:

```text
mean corr(B*A, |C|) > 0.78
```

No coefficient fitting or per-body calibration is permitted.

### C3 — the combined score retains the soma/root as a high-percentile point

PASS if:

```text
median soma percentile in B*A > 0.90
and at least 10/12 somata are above the 0.75 percentile.
```

### C4 — normalized alignment alone is not the soma's main location privilege

This is a registered negative-structure prediction from discovery.

PASS if:

```text
median soma alignment percentile < 0.75
and fewer than 10/12 somata exceed the 0.75 alignment percentile.
```

This prevents a later reinterpretation that `A` alone was the special soma property.

## Descriptive quantities

Also report without additional pass/fail thresholds:

- mean `corr(A,|C|)`, `corr(B,|C|)`, `corr(A,B)`;
- soma percentiles for `A`, `B`, `B*A`, and `|C|`;
- unnormalized cross-difference RMS correlation with `|C|`;
- modal/coordinate reconstruction error.

## Interpretation table fixed in advance

- `C1/C2 pass`, `C4 pass`: **balance + temporal compatibility** is the supported decomposition; soma is balance-special, not alignment-special.
- `C1/C2 fail`: the discovery product was unstable; retain amplitude balance as the simpler explanation.
- `C4 fail`: revisit whether the soma actually has an independently special task-kernel orientation.
