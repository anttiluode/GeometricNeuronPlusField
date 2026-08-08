# Coupling-contrast confirmation preregistration v0.1

`HARDWARE_COMPILER_V01.md` maps model coupling to inverse reactive coupling in an LC/transmission-line implementation. The historical strong/weak ratio is extreme:

```text
k_arbor / k_bath = 2.5 / 0.0002 = 12,500 : 1.
```

A development sweep on already-opened bodies 288-299 found that order-selective magnitude did not collapse as the bath coupling was increased. At **100:1** (`k_bath=0.025`) most body/lag directions were still preserved and contrast was often larger.

This fresh test asks whether that specific relaxed ratio replicates.

## Frozen protocol

Fresh bodies: seeds **300-311**.

No learning, morphology change, or retuning.

Lags:

```text
14, 18, 20, 22, 26
```

Compare the same frozen body at:

```text
historical: k_bath = 0.0002   ratio 12,500:1
primary:    k_bath = 0.025    ratio    100:1
```

All body-body arbor bonds stay at `k_arbor=2.5`. Every non-arbor bond, including exterior leakage in the exact linear surrogate, uses the tested `k_bath`.

For each body and lag measure both:

```text
energy contrast
peak soma-power contrast.
```

Direction is defined relative to that body's historical-ratio sign for the same lag.

## Registered tests

### C1 — peak-direction retention

For each body compute the fraction of the five lags whose peak-contrast sign matches the historical system.

Pass if:

```text
median body retention >= 0.80
and at least 9/12 bodies have retention >= 0.80.
```

### C2 — smooth energy-direction retention

Same criterion for energy contrast:

```text
median body retention >= 0.80
and at least 9/12 bodies have retention >= 0.80.
```

### C3 — peak-information magnitude does not require 12,500:1

For each body compute

```text
R_peak = mean_lag |C_peak(100:1)| / mean_lag |C_peak(12500:1)|.
```

Pass if:

```text
median R_peak >= 0.80
and at least 9/12 bodies have R_peak >= 0.80.
```

### C4 — energy-information magnitude

Same criterion for the smooth energy contrast ratio `R_energy`.

## Secondary descriptive sweep

The workflow may also calculate the already-fixed full contrast ladder

```text
12500, 2500, 500, 250, 100, 50, 25, 10, 5, 2.5 : 1
```

on the fresh bodies. Only the 100:1 criteria above are confirmatory. No lower-ratio threshold will be promoted after inspection.

## Interpretation

If C1-C4 pass, it is fair to say that the observed order-selective wave computation is **not dependent on the historical 12,500:1 strong/weak contrast**; a 100:1 contrast preserves most of the original directional code and magnitude without retraining.

It would not establish that 100:1 is easy for a specific fabricated platform. That depends on the chosen coupler technology and loss/Q constraints.
