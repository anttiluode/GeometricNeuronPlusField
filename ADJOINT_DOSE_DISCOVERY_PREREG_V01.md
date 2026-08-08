# Adjoint conductance-step discovery preregistration v0.1

## Why this experiment exists

The first adjoint discovery (`ADJOINT_ELIGIBILITY_DISCOVERY_PREREG_V01.md`, seeds 180-191) produced a sharp split:

```text
small-epsilon adjoint finite-difference check     ~4.5e-9 relative error  PASS
exact finite linear event -> nonlinear C_int      mean r ~0.986          PASS
base adjoint -> full binary event                 mean r ~-0.417         FAIL
```

So the adjoint mathematics was locally correct and the linear wave surrogate was excellent, yet a full bath->arbor conductance jump was far outside the useful first-order regime.

A minor implementation discrepancy was also found after that run: the surrogate omitted the tiny mature-bath coupling from the four outer grid boundaries to the zero-valued exterior used by FunctionalArbor's shift operator. The observed Laplacian mismatch was only ~`4e-5` relative, too small to plausibly explain the sign reversal, but this experiment corrects it before proceeding.

> **Does the adjoint correctly predict gradual structural maturation, with prediction quality breaking down only when the local conductance change becomes too large?**

## Bodies

Fresh bodies only:

```text
seeds 192-203
```

Same legal event generator as the structural-interference experiments: up to six additions and six deletions per body, fixed source terminals, lag 20, 210 steps.

## Exact linear surrogate

Use the same mature wave coefficients as FunctionalArbor, no saturation, and now include the exact weak bath-to-zero exterior boundary terms. Verify the custom Laplacian against `m._lap(..., mature=True)`.

## Conductance path

For each registered structural event define `alpha` as the fraction of the full conductance change applied to every affected bond:

```text
alpha = 0       original body
alpha = 1       full binary event
```

Fixed dose grid:

```text
1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1
```

The base adjoint is computed **once at alpha=0**. No gradient is recomputed after moving along the path.

For each alpha compare

```text
pred_adjoint(alpha) = alpha * grad_C · DeltaK_event
```

with the exact finite linear change

```text
Delta C_lin(alpha) = C_lin(alpha) - C_lin(0).
```

## Registered discovery predictions

### D0 — corrected implementation identity

PASS if both:

```text
mean custom-vs-FunctionalArbor Laplacian relative error < 1e-10
mean small-epsilon adjoint finite-difference relative error < 1e-3.
```

### D1 — very small maturation steps are in the adjoint regime

At `alpha=1e-4`, within-body event correlations

```text
r(alpha) = corr(pred_adjoint(alpha), Delta C_lin(alpha))
```

must have mean `>0.95`, with all `12/12` bodies positive.

### D2 — prediction quality collapses over the full binary jump

PASS if

```text
mean r(1e-4) - mean r(1) > 0.60.
```

No sign is preregistered for `r(1)` itself.

### D3 — there is a nontrivial useful step-size range

For each body define the largest tested `alpha` whose correlation remains `>=0.70`. PASS if the median such alpha is at least `1e-3`.

This asks whether the usable region is larger than an infinitesimal numerical perturbation.

### D4 — full finite linear events still reproduce nonlinear interference changes

At `alpha=1`, compare exact finite `Delta C_lin(1)` with the nonlinear counterfactual `dC_int` from the same event.

PASS if mean within-body correlation `>0.90` and all `12/12` bodies are positive.

## Descriptive outputs

- mean/median r(alpha) over the entire dose grid;
- pooled sign agreement by alpha;
- additions versus deletions;
- full-step relation to `dC_peak`;
- fraction of events whose finite response changes sign relative to the base derivative as alpha grows.

## Interpretation fixed in advance

If D0-D4 pass, the failed binary adjoint experiment is a **step-size failure**, not a credit-coordinate failure. The next model change worth testing would be gradual bond/branch maturation with repeated local sensitivity updates.

If D1 fails after the boundary correction, the adjoint implementation or objective derivation remains suspect.

If D1 passes but D3 fails, the useful neighborhood is so tiny that it is unlikely to help a structural-growth model without a substantially different parameterization.
