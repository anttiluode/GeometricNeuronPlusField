# Matched-tuner v0.1 interpretation correction — scalar count was matched, coordinate scale was not

`MATCHED_TUNER_CONFIRM_V01.md` remains a valid numerical report of the optimizer that was actually run. Under that optimizer, the selected free-modal coordinates outperform local bond coordinates on all 12 held-out bodies.

However, `SPECTRAL_TO_LOCAL_COMPILER_V01.md` exposed a normalization flaw in the interpretation.

Both arms used

```text
eta = 0.01
```

but `0.01` means different things in the two parameterizations:

```text
local coordinate       rho in [0,1]
free pole coordinate   log-stiffness-like variable
free residue coordinate relative modal gain variable.
```

The same numerical step is not a matched physical or task-space step.

Across the 96 winning free coordinates, a single feasible local-bond tangent reproduces the free tangent's **25-lag functional shape** with median relative residual `0.0276`, but the median amplitude coefficient is `14.9`. Thus the two optimizers were often taking very different-sized moves through transfer-function space.

## What remains valid

- The exact same base transfer was used.
- The free optimizer achieved much better held-out objective under the declared update rules.
- The winning free coordinates map to low body modes 1-16 and align with the soma interference mechanism.
- Pole-only and residue-only free optimizers each beat the local optimizer under those same declared update rules.

## What is no longer earned

Do **not** cite v0.1 as proving:

> “free spectral coordinates are intrinsically better per trainable scalar.”

Scalar count alone does not define a fair optimizer comparison when coordinate scales differ.

## Canonical status

Treat `MATCHED_TUNER_CONFIRM_V01.md` as:

> **an optimizer-coordinate result under one explicit normalization, not a definitive parameter-efficiency benchmark.**

The canonical next benchmark must normalize updates by their induced multi-task output change (task-space trust region) or otherwise use a coordinate-scale-invariant optimization protocol.

Until that is run, the architecture question remains open.
