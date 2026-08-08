# Coordinate-curvature confirmation preregistration v0.1

The scale-invariant trust benchmark confirmed small-step parity and a shift toward spectral coordinates as trust radius grows, but its strong large-step dominance criterion narrowly failed. A reused-seed mechanism probe (240-243) then measured the missing quantity directly: **finite-step departure from each coordinate's own first-order task tangent**.

Development showed substantially larger curvature for local bond coordinates than for direct spectral coordinates.

This fresh run tests that mechanism on new bodies.

## Frozen protocol

Fresh bodies: seeds **324-335**.

Select `P=8` local and `P=8` free spectral coordinates exactly as in the scale-invariant trust benchmark, using training lags

```text
16,20,24
```

and requiring enough one-sided parameter range for the largest one-coordinate trust step.

Evaluate finite-step response over the seven-lag vector

```text
C = [C(14), C(16), C(18), C(20), C(22), C(24), C(26)].
```

For each selected coordinate separately and each training-task trust radius

```text
delta = .001, .0025, .005, .010
```

choose a parameter displacement whose **initial three-training-lag tangent** predicts RMS change `delta`. Apply only that coordinate from the common base state, then compare the actual seven-lag change `Delta C_actual` with the first-order prediction `Delta C_pred`.

Primary curvature metric:

```text
E = ||Delta C_actual - Delta C_pred||_2 / ||Delta C_pred||_2.
```

No coordinate is jointly optimized in this experiment.

## Registered tests

### C1 — local coordinates are already more curved at the smallest step

At `delta=.001`, pooled across selected coordinates:

```text
median E_graph / median E_free > 5
and median E_graph > .015.
```

### C2 — large-step curvature gap

At `delta=.010`:

```text
median E_graph / median E_free > 5
and median E_graph > .15.
```

### C3 — graph curvature grows substantially with trust radius

```text
median E_graph(.010) - median E_graph(.001) > .10.
```

### C4 — free spectral coordinates remain comparatively straight

```text
median E_free(.010) < .10.
```

### C5 — no bound-clipping explanation

Across the primary selected coordinates and tested deltas, fewer than 5% of one-coordinate probes may be clipped by parameter bounds in either arm. If clipping exceeds that threshold, the curvature comparison is considered confounded and C1-C4 are not promoted.

## Interpretation boundary

If the curvature gap confirms, the result is not that local physical coordinates point in bad directions. The task-space compiler and small-trust benchmark already say otherwise.

The intended statement is:

> **Local physical bond coordinates and direct spectral coordinates can have similar first-order task directions, but the local coordinate chart bends much faster away from the base state. Spectral coordinates are straighter optimization coordinates.**

That would explain why small matched trust steps produce parity while larger trust steps increasingly favor direct spectral tuning.
