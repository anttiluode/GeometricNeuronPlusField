# Reciprocal pass-mismatch learning — held-out preregistration v0.1

## Question

The exact same-medium adjoint requires the forward and retro passes to use the same reciprocal operator. Development showed that small random fractional edge mismatch barely changes the gradient map, while a transition appears only at very large mismatch: around 30% the K=8 learner begins to lose gain, and around 50% map/learning quality degrades strongly.

This test asks whether useful broadband compressed training survives a deliberately large but fixed **forward/retro operator mismatch** on fresh bodies.

This is not a calibrated fabrication-tolerance statement. `sigma=.20` means the retro-pass conductance of every edge is multiplied by `1 + delta_e`, with `delta_e ~ N(0,.20)`, clipped only to remain positive. The same random fractional mismatch is fixed for the whole training run and is shared by target/distractor retro passes. Each pass is internally reciprocal; the two passes simply do not have identical operators.

## Frozen implementation

Script: `operator_mismatch_learning_probe.py`

Fresh bodies:

```text
seeds              448-459
lag                 20
steps               210
candidate bonds     up to 8, same deterministic set in every arm
eta                 .01
iterations          40
```

Primary arms:

```text
exact       exact discrete adjoint
K8_ideal    K=8 compressed physical gradient, no operator mismatch
K8_m20      K=8, sigma=.20 fixed retro-pass edge mismatch
K8_m30      K=8, sigma=.30 fixed retro-pass edge mismatch
K16_m30     K=16, sigma=.30 fixed retro-pass edge mismatch
K8_m50      K=8, sigma=.50 stress arm
```

Boundary spectral bins are reselected from the current arm's port spectra at every iteration. No internal oracle is used.

## Registered criteria

### C0 — exact learner gate

```text
mean exact DeltaC > .015
at least 10 / 12 bodies improve
```

### C1 — K8 survives 20% pass mismatch

```text
mean K8_m20 DeltaC > .015
at least 10 / 12 bodies improve
```

Both required.

### C2 — K8 at 20% preserves most exact learning gain

Using group means:

```text
mean(K8_m20 DeltaC) / mean(exact DeltaC) >= .75
mean(K8_m20 DeltaC - exact DeltaC) > -.010
```

Both required.

### C3 — K16 remains useful at 30% pass mismatch

```text
mean K16_m30 DeltaC > .015
mean(K16_m30 DeltaC) / mean(exact DeltaC) >= .70
```

Both required.

### C4 — compressed gradient direction remains recognizable along the mismatched trajectory

Average the per-iteration map correlation within each body, then across bodies:

```text
K8_m20  mean map correlation > .970
K16_m30 mean map correlation > .975
```

Both required.

## Descriptive stress quantities

`K8_m30` and `K8_m50` gain, map correlation, monotonicity and body failures are reported, but no transition threshold is registered. The development sample was too small to preregister a sharp universal failure wall responsibly.

## Interpretation boundary

If C0-C4 pass, the earned statement is:

> **In this toy transient mesh, broadband compressed reciprocal training tolerates substantial fixed pass-to-pass operator mismatch: K=8 remains useful under 20% RMS edge mismatch, and K=16 remains useful under 30% RMS edge mismatch.**

Do not call `20%` or `30%` a silicon-photonic fabrication tolerance. The perturbation is a model-space conductance mismatch, not a direct mapping from heater phase error, optical loss, S-parameter drift, or nonreciprocal hardware.

This test also does not model true directed/nonreciprocal couplings within one pass. That requires a nonsymmetric forward operator and a separate test of whether the physical reverse device implements the transpose/adjoint operator.
