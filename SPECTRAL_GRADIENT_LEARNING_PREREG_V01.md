# Spectral gradient learning — held-out preregistration v0.1

## Question

`SPECTRAL_CORRELATION_COMPRESSION_V01.md` confirmed that a boundary-selected 8- or 16-bin phasor representation closely approximates the exact internal bond-gradient map at the base state.

The stronger question is whether that approximate physical gradient remains useful **after the structure starts changing**.

This test closes the loop. The compressed gradient is recomputed from the current arm state at every iteration and is used to update the same graded frontier conductances as the exact adjoint learner.

## Frozen implementation

Script: `spectral_gradient_learning_probe.py`

Fresh bodies:

```text
seeds              412-423
lag                 20
steps               210
candidate bonds     up to 8, same deterministic set in every arm
eta                 .01
iterations          40
```

Arms:

```text
exact    full discrete adjoint
K4       4 boundary-selected spectral bins   (descriptive stress arm)
K8       8 boundary-selected spectral bins   (primary compressed arm)
K16      16 boundary-selected spectral bins  (secondary compressed arm)
```

Every arm starts from `rho=0` on the identical candidate set. Each arm then follows its own state trajectory. At every iteration, K8/K16 recompute the port-only spectral ranking from that arm's current external-source and soma-return spectra. No internal oracle ranking is used.

The update normalization and clipping are exactly the same in all arms.

## Registered criteria

### C0 — exact learner positive control

The exact arm must still be a useful local optimizer on these fresh bodies:

```text
mean exact DeltaC > .015
at least 10 / 12 bodies improve
```

If this fails, compressed-vs-exact comparison is not interpretable as a preserved learner.

### C1 — K8 learns

```text
mean K8 DeltaC > .015
at least 10 / 12 bodies improve
```

Both required.

### C2 — K8 preserves most exact learning gain

Using group means, not the unstable mean of per-body ratios:

```text
mean(K8 DeltaC) / mean(exact DeltaC) >= .85
mean(K8 DeltaC - exact DeltaC) > -.007
```

Both required.

### C3 — K16 preserves most exact learning gain

```text
mean(K16 DeltaC) / mean(exact DeltaC) >= .85
mean(K16 DeltaC - exact DeltaC) > -.007
```

Both required.

### C4 — gradient-map fidelity survives the changing trajectory

Averaging the per-iteration map correlations within each body and then across bodies:

```text
K8  mean map correlation > .980
K16 mean map correlation > .990
```

Both required.

## Descriptive quantities, not pass/fail criteria

- monotone-step fraction;
- total final frontier material `sum rho`;
- K4 learning;
- number of bodies where a compressed arm happens to beat exact;
- per-body compressed/exact gain ratios.

These will be reported but cannot rescue failed registered criteria.

## Interpretation

If C0-C4 pass, the earned statement is:

> **A small boundary-selected phasor bank preserves most of the useful closed-loop local learning obtained from the exact time-domain adjoint, even when the gradient is relinearized after every structural update.**

This would move spectral compression from a bond-map reconstruction result to a training result.

It still would not prove an intensity-only chip implementation. Each selected bin is presently an abstract complex local phasor correlation. The next hardware step would be to realize that complex product through phase-stepped/coherent intensity measurements and count the required physical trials/sensors.
