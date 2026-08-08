# Adjoint conductance-step held-out confirmation preregistration v0.1

Discovery on fresh seeds 192-203 cleanly separated an infinitesimal/local sensitivity regime from the full binary structural jump.

Discovery dose curve:

```text
alpha        mean event correlation: base adjoint vs exact finite linear change
1e-4         0.99997
3e-4         0.99976
1e-3         0.99745
3e-3         0.98117
1e-2         0.87246
3e-2         0.30670
1e-1        -0.34697
3e-1        -0.53963
1           -0.32787
```

Additional discovery facts:

```text
median largest alpha with r >= .70       0.01
corrected Laplacian relative error        ~1.4e-11
adjoint small-epsilon FD error            ~4.8e-9
full exact linear event vs nonlinear C_int mean r  0.9810
```

All discovery D0-D4 criteria passed. This document freezes confirmation before any new body is run.

## Held-out bodies

```text
seeds 204-215
```

Same event generator, same dose grid, same corrected exact mature-boundary Laplacian, same lag 20 and 210-step trajectories. The adjoint is computed once at `alpha=0` and never updated along a dose path.

## Registered confirmation criteria

### C0 — implementation remains exact

PASS if:

```text
mean custom-vs-FunctionalArbor Laplacian relative error < 1e-9
mean small-epsilon adjoint finite-difference relative error < 1e-3.
```

### C1 — small structural steps are accurately predicted

At `alpha=1e-3`:

```text
mean within-body r > 0.95
all 12/12 bodies positive.
```

### C2 — a useful non-infinitesimal neighborhood remains

At `alpha=1e-2`:

```text
mean within-body r > 0.75
at least 10/12 bodies positive.
```

This is the practical part of the hypothesis: the useful range must extend beyond an epsilon-scale derivative check.

### C3 — the same base derivative fails as the event approaches a binary jump

PASS if both:

```text
mean r(alpha=1e-3) - mean r(alpha=1e-1) > 0.80
mean r(alpha=1) < 0.30.
```

The sign of the full-step mean is not preregistered; only substantial loss of predictive validity is.

### C4 — the radius-of-validity scale replicates

For each body let `alpha_max` be the largest tested dose with within-body `r >= .70`.

PASS if:

```text
median alpha_max >= 0.003
median alpha_max <= 0.10.
```

This rules out both "only infinitesimal" and "the derivative stays useful essentially all the way to the binary jump."

### C5 — the linear finite-event model still matches the nonlinear interference counterfactual

At the full binary event (`alpha=1`):

```text
mean corr(dC_lin_finite, dC_int) > 0.90
all 12/12 bodies positive.
```

This isolates the failure to first-order step size rather than to the linear wave surrogate itself.

## Descriptive outputs

Also report:

- complete mean/median dose curve and sign agreement;
- additions versus deletions;
- full-step correlation with historical `dC_peak`;
- fraction of events whose exact finite response reverses sign relative to the base adjoint as alpha grows.

## Interpretation fixed in advance

If C0-C5 pass, the first adjoint failure is classified as a **finite-step failure**. In this model, task-conditioned local sensitivity is accurate for gradual conductance maturation but cannot be extrapolated across an abrupt bath-to-arbor jump. The next justified experiment would therefore update sensitivity repeatedly while bonds mature in small increments.

If C1 fails, the adjoint mechanism itself is not stable across bodies. If C1 passes but C2 fails, the useful neighborhood is likely too narrow for a practical gradual-growth rule under this parameterization.
