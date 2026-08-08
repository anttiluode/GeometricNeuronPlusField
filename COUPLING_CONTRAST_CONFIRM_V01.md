# Coupling-contrast confirmation v0.1 — the temporal code does not require 12,500:1 coupling contrast

`HARDWARE_COMPILER_V01.md` exposed a practical issue hidden by the abstract model: the historical parameters use

```text
k_arbor / k_bath = 12,500 : 1.
```

A development sweep suggested that much smaller contrast could preserve the temporal-order computation. `COUPLING_CONTRAST_CONFIRM_PREREG_V01.md` therefore froze **100:1** (`k_bath=0.025`) as the fresh primary condition.

Fresh bodies: seeds **300-311**. No learning or retuning. Five lags per body: `14,18,20,22,26`.

## C1 — peak-direction retention

Registered:

```text
median body sign retention >= 0.80
at least 9/12 bodies >= 0.80.
```

Observed:

```text
median retention                1.00
bodies >= 0.80                  9 / 12
mean over bodies                0.90
```

**C1 PASS.**

Nine bodies preserve the historical peak-contrast direction at all five lags. Three preserve 3/5.

## C2 — smooth energy-direction retention

Observed:

```text
median retention                1.00
bodies >= 0.80                  9 / 12
mean over bodies                0.7667
```

**C2 PASS**, exactly at the registered body-count threshold.

The lower mean is caused by three bodies with substantial direction reorganization; the median/body-count criterion was preregistered precisely to avoid treating 60 within-body lag observations as independent replicates.

## C3 — peak-information magnitude

For each body:

```text
R_peak = mean_lag |C_peak(100:1)| / mean_lag |C_peak(12500:1)|.
```

Registered median `>=0.80` and at least 9/12 bodies `>=0.80`.

Observed:

```text
median R_peak                  1.1918
bodies >= 0.80                 9 / 12
```

**C3 PASS.**

The pooled mean absolute peak contrast actually rises slightly:

```text
12,500:1     0.2318
100:1        0.2599
```

This is not evidence that 100:1 is universally better; body-level ratios vary substantially. It does show that the information magnitude does not depend on the historical extreme contrast.

## C4 — smooth energy-information magnitude

Observed:

```text
median R_energy                1.6157
bodies >= 0.80                11 / 12

mean |C_energy|
12,500:1                       0.1374
100:1                          0.2192
```

**C4 PASS.**

## Formal result

```text
C1 peak direction        PASS
C2 energy direction      PASS
C3 peak magnitude        PASS
C4 energy magnitude      PASS
```

So all four registered 100:1 tests pass on fresh bodies.

> **The observed temporal-order computation is not dependent on the historical 12,500:1 strong/weak coupling ratio. Reducing the ratio to 100:1 without retraining preserves the original direction on most bodies/lags and preserves or increases order-selective magnitude on most bodies.**

## Descriptive lower-ratio ladder

The full preregistered secondary ladder shows a useful transition rather than a sudden collapse.

Mean sign retention over body-level lag fractions:

```text
ratio        energy     peak
12500        1.000      1.000
2500         1.000      0.983
500          1.000      0.950
250          0.950      0.950
100          0.767      0.900   <-- confirmatory point
50           0.767      0.817
25           0.700      0.800
10           0.617      0.783
5            0.683      0.533
2.5          0.683      0.533
```

Absolute contrast remains substantial even when direction relative to the original mapping starts to reorganize. Increasing bath coupling is therefore not simply “adding noise until computation disappears.” It changes the operator and can create a different temporal computation.

## Hardware interpretation

This is encouraging for the LC/transmission-line compiler. A finite **100:1** reactive-coupling ratio is still demanding but qualitatively different from 12,500:1.

The next device-level question is no longer whether the task requires an almost disconnected bath. It is:

> **Can a candidate reciprocal coupler technology realize the required ~100:1 controllable range with acceptable loss, calibration error, and local monitoring cost?**

That question should be answered in the native parameters of the chosen platform, not by importing an MZI phase-error number into `k_e`.
