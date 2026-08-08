# Device-inspired compressed-gradient learning — held-out preregistration v0.1

## Question

The broadband reciprocal learner now uses only K=8 or K=16 boundary-selected spectral correlation channels. The next hardware-facing question is whether the learner survives measurement/calibration errors at scales motivated by in-situ photonic-backpropagation studies.

This is **not** a calibrated device model of Pai et al.'s silicon chip. The noise variables below are explicit abstractions chosen near published phase-error/tap-noise study scales.

## Frozen error model

Script: `device_error_learning_probe.py`

Fresh bodies:

```text
seeds              436-447
lag                 20
steps               210
candidate bonds     up to 8
eta                 .01
iterations          40
```

Arms:

```text
exact       full discrete adjoint
K8_ideal    boundary-selected K=8, no added measurement error
K8_mod      sigma_phi=.025 rad, sigma_amp=.025, local tap sigma=.010
K8_high     sigma_phi=.050 rad, sigma_amp=.050, local tap sigma=.020
K16_high    same high errors, K=16
```

Error structure:

- phase error is one Gaussian phase-setting offset per absolute spectral bin, shared across all bonds and held fixed throughout a body's training run;
- amplitude error is one Gaussian multiplicative calibration factor per absolute spectral bin, shared across bonds and fixed through training;
- tap/readout noise is independent per bond/bin and redrawn on every gradient measurement, with standard deviation equal to the stated fraction of that bin's RMS local complex-product magnitude.

Thus systematic port/controller error is kept distinct from local measurement noise.

## Registered criteria

### C0 — exact learner gate

```text
mean exact DeltaC > .015
at least 10 / 12 bodies improve
```

If C0 fails, the hardware-error comparison is not interpreted as preservation of a useful learner.

### C1 — high-error K8 still learns

```text
mean K8_high DeltaC > .015
at least 10 / 12 bodies improve
```

Both required.

### C2 — high-error K8 preserves most exact gain

Using group means:

```text
mean(K8_high DeltaC) / mean(exact DeltaC) >= .75
mean(K8_high DeltaC - exact DeltaC) > -.010
```

Both required.

### C3 — high-error K16 preserves most exact gain

```text
mean(K16_high DeltaC) / mean(exact DeltaC) >= .80
mean(K16_high DeltaC - exact DeltaC) > -.008
```

Both required.

### C4 — gradient direction remains recognizable under errors throughout training

Average within each body over all relinearization iterations, then average bodies:

```text
K8_high  mean map correlation > .975
K16_high mean map correlation > .990
```

Both required.

### C5 — added device-inspired error does not dominate the existing K8 compression error

Compare the K8_high and K8_ideal arms directly:

```text
mean(K8_high DeltaC - K8_ideal DeltaC) > -.007
```

This tests the incremental cost of the modeled device error after the K=8 spectral truncation has already been paid.

## Descriptive only

- K8_mod arm;
- monotone-step fraction;
- final material sum;
- number of noisy arms beating exact;
- individual-body failures/outliers.

They cannot rescue a failed registered criterion.

## Interpretation boundary

If C0-C5 pass, the earned statement is:

> **Under this explicit device-inspired error model, the already-compressed broadband in-situ gradient remains useful for closed-loop training at 0.05-rad shared phase-setting error, 5% shared spectral-amplitude error, and 2%-of-local-bin-RMS tap/readout noise.**

Do not translate that directly into a specific fabrication tolerance or chip-level accuracy without a hardware-specific mapping of these abstract error variables.
