# True nonreciprocity and physical adjoint transport — held-out preregistration v0.1

## Question

The reciprocal FunctionalArbor satisfies `H^T = H`, which is why a time-reversed objective-derivative waveform can be replayed through the same physical medium to generate the exact adjoint.

The pass-mismatch experiment showed that approximate operator equality can be surprisingly tolerant. That does not answer the more basic causality question:

> **Does same-medium physical backpropagation fail specifically when the forward operator is genuinely nonsymmetric, while propagation through the transpose operator restores the exact adjoint?**

## Frozen nonreciprocal operator

Script: `nonreciprocal_adjoint_probe.py`

Start from the usual symmetric weighted Laplacian `L` and add a local real skew-symmetric nearest-neighbour coupling `A`:

```text
H_beta = L + beta A
A^T = -A
```

Therefore

```text
H_beta^T = H_-beta.
```

The random local skew pattern is fixed per body. The ordinary symmetric bond conductances remain the structural coordinates whose gradient is audited; the skew background is held fixed.

This is an abstract local nonreciprocity model, not a calibrated isolator/circulator model.

Fresh bodies:

```text
seeds   460-471
lag     20
steps   210
betas   0, .02, .05, .10, .20, .30, .40, .60
```

For each beta:

```text
exact
  algorithmic discrete adjoint uses H_beta^T = H_-beta

same-H replay
  reverse the soma derivative waveform in time, but send it through H_beta again

transpose replay
  send the same reversed waveform through H_-beta
```

The primary observable is the gradient map with respect to the symmetric bond conductances.

## Registered criteria

### C0 — reciprocal baseline

At beta=0:

```text
mean same-H gradient correlation > .999999
mean same-H relative L2 < 1e-10
```

### C1 — transpose replay stays exact after nonreciprocity is added

At both beta=.10 and beta=.20:

```text
mean transpose-gradient correlation > .999999
maximum body transpose-gradient relative L2 < 1e-10
```

All four conditions are required.

### C2 — same-H replay is already measurably wrong at beta=.10

```text
mean same-H gradient correlation < .95
at least 10 / 12 bodies have same-H correlation < .98
```

Both required.

### C3 — same-H failure strengthens at beta=.20

```text
mean same-H gradient correlation < .88
mean same-H relative L2 > .50
```

Both required.

### C4 — transpose specifically rescues the gradient

At beta=.20:

```text
mean[ corr(transpose) - corr(same-H) ] > .10
```

At beta=.10:

```text
mean[ corr(transpose) - corr(same-H) ] > .03
```

Both required.

## Descriptive only

The beta=.02, .05, .30, .40 and .60 conditions describe the shape of the failure curve but cannot rescue a failed registered criterion.

## Interpretation boundary

If C0-C4 pass, the earned statement is:

> **Same-medium replay works because of the transpose symmetry, not merely because a wave was sent backward in time. A genuinely nonsymmetric forward operator breaks the same-H physical adjoint, while replay through the actual transpose operator restores the exact gradient.**

This directly supports the engineering statement:

> **Reciprocity buys a zero-extra-operator physical adjoint.**

It does **not** mean nonreciprocal systems cannot be trained by adjoints. They can; they simply require access to the transpose/adjoint operator by some other physical or algorithmic mechanism.
