# Interference factorization held-out confirmation preregistration v0.1

Discovery on fresh seeds 132-143 supported the explicit local decomposition:

```text
mean corr(C_int, C_peak)                 0.7826
mean corr(|C_int|, |C_peak|)             0.7841
mean corr(V, |C_peak|)                   0.6297
mean corr(V*|Delta_rho|, |C_peak|)       0.8388
mean improvement over V                  0.2091
improved / worse                         11 / 1
sign p                                   0.00635
median soma visibility percentile        0.9286
median soma |Delta_rho| percentile       0.6929
```

All discovery predictions passed. This document freezes the held-out criteria before any new body is run.

## Held-out bodies

```text
seeds 144-155
```

Definitions, lag `20`, trace length `210`, zero-padding convention, and code in `interference_factorization_probe.py` remain unchanged. No fitted weights, lag search, or per-body calibration.

## Registered confirmation criteria

### C1 — integrated interference tracks signed peak-order computation

PASS if:

```text
mean corr(C_int, C_peak) > 0.72
positive bodies >= 11/12
```

### C2 — magnitude relation replicates

PASS if:

```text
mean corr(|C_int|, |C_peak|) > 0.72
```

### C3 — directional compatibility adds substantial information beyond visibility alone

Let

```text
N = V * |Delta_rho|.
Delta_r = corr(N, |C_peak|) - corr(V, |C_peak|).
```

PASS if all are true:

```text
mean Delta_r > 0.12
positive bodies >= 9/12
two-sided sign-test p < 0.05
```

### C4 — soma privilege is primarily amplitude visibility, not unusually extreme compatibility

PASS if both are true:

```text
median soma V percentile > 0.88
median soma |Delta_rho| percentile < 0.80
```

This is a structural prediction: convergence should make the root unusually balanced/visible, while the directional complex compatibility can be only moderate and body-specific.

## Exact-identity control

The zero-padded integrated target/distractor energy contrast must match the closed-form visibility/coherence expression to relative error `<1e-10` averaged across bodies.

## Interpretation fixed in advance

If C1-C4 pass, the strongest current local mechanism becomes:

```text
geometry -> h_A,h_B
         -> high interference visibility where amplitudes are balanced
         -> lag-directional complex compatibility
         -> coherent cross-source order signal
         -> point square-law readout
```

This does not make `C_int` identical to the historical peak-power objective. It says a standard, coefficient-free integrated interference statistic extracted from the same two transfer histories captures the same spatial order-selectivity structure.
