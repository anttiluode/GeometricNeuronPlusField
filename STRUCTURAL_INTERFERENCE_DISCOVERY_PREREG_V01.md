# Structural-event interference discovery preregistration v0.1

## Question

The confirmed local mechanism is

```text
geometry -> h_A,h_B -> visibility V + directional compatibility -> order contrast.
```

The old FunctionalArbor learning wall is structural credit assignment. We can now ask a more causal question:

> **When one local structural event changes the arbor, does the resulting change in temporal-order task value follow the change in the visibility/compatibility interference statistic?**

This is not yet a biologically local learning rule. It is a counterfactual audit intended to identify the *right target quantity* for eligibility.

## Bodies

Fresh bodies only:

```text
seeds 156-167
```

Same mature v0.9 wave physics, lag `20`, 210-step probes.

## Structural events

For each frozen body construct two legal one-cell event classes while keeping soma and source terminals fixed:

1. **tip-like addition** — add one empty window cell with exactly one occupied 4-neighbor;
2. **safe deletion** — remove one occupied non-source, non-soma cell only if the body remains connected and both source terminals remain exactly the same.

Cells inside either source target mask are excluded from both classes.

Sample up to six events of each class per body from a seed-fixed RNG, for at most 12 events/body.

No growth/credit dynamics occur; each event is evaluated as a frozen before/after counterfactual.

## Soma metrics

From independently measured single-source soma histories compute:

- historical peak contrast `C_peak`;
- energy visibility `V`;
- `rho_plus`, `rho_minus`;
- exact zero-padded integrated interference contrast

```text
C_int = V*(Re rho_plus - Re rho_minus) /
        [2 + V*(Re rho_plus + Re rho_minus)].
```

For event `e` define

```text
dC_peak(e) = C_peak_after - C_peak_before
dC_int(e)  = C_int_after  - C_int_before.
```

## Exact two-factor change decomposition

Write

```text
f(V,R) = C_int
R = (Re rho_plus, Re rho_minus).
```

Use the two-factor Shapley decomposition:

```text
contrib_V = 1/2 [f(V1,R0)-f(V0,R0) + f(V1,R1)-f(V0,R1)]
contrib_R = 1/2 [f(V0,R1)-f(V0,R0) + f(V1,R1)-f(V1,R0)]
```

so exactly

```text
contrib_V + contrib_R = dC_int.
```

`contrib_R` is the structural effect mediated through directional complex compatibility; `contrib_V` is the effect mediated through amplitude opportunity/visibility.

## Registered discovery predictions

### D1 — interference change predicts peak-task change

Within each body correlate event-wise `dC_int` with `dC_peak`.

Prediction:

```text
mean body correlation > 0.50
positive-correlation bodies >= 10/12.
```

### D2 — event direction is often predicted correctly

Across all non-negligible events (`|dC_peak| > 1e-5`), compare signs of `dC_int` and `dC_peak`.

Prediction: pooled sign agreement `> 0.65`.

### D3 — full interference change beats visibility change alone

Within each body compare

```text
corr(dC_int, dC_peak)
vs
corr(dV, dC_peak).
```

Prediction: mean improvement `> 0.10` and at least `9/12` bodies improve.

## Descriptive mechanism questions

No discovery pass/fail threshold is set for these:

- fraction of `|dC_int|` attributable to `|contrib_V|` versus `|contrib_R|`;
- whether additions and deletions differ;
- relation to graph distance from soma;
- magnitude of `dC_peak` and `dC_int`;
- Shapley reconstruction error.

If one factor clearly dominates, a fresh held-out confirmation will be frozen before new bodies are run.
