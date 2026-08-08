# Analog frontier learning discovery preregistration v0.1

## Question

`BOND_RESPONSE_V01.md` showed that gradient-favored frontier additions commonly have nonmonotonic conductance-response curves and reproducible interior optima. The exact adjoint is accurate only locally, so any usable learning rule must update the sensitivity as the geometry changes.

> **Can repeated small, relinearized adjoint steps actually learn several frontier conductances jointly?**

This is the first experiment in this line that uses the adjoint as an optimizer rather than only auditing it.

## Bodies

Fresh FunctionalArbor bodies:

```text
seeds 240-251
```

Same exact linear mature-wave surrogate, fixed source terminals, lag 20, and 210-step target/distractor trajectories.

## Design variables

For each body sample up to eight legal tip-like frontier additions using the same geometry-only candidate rule as the structural-interference experiments.

Each candidate is one weak bath bond from the existing body to an empty tip cell. Give it a continuous maturation variable

```text
rho_e in [0,1]

k_e = k_bath + rho_e (k_arbor-k_bath).
```

All `rho_e` begin at zero. The candidate set is frozen for the run; this experiment does not yet expand the frontier recursively.

## Objective

Use the smooth integrated soma-energy contrast

```text
C_lin = (E_target - E_distractor) / (E_target + E_distractor).
```

The adjoint gives `dC_lin/dk_e` for every grid bond in one backward pass.

## Three arms

All arms start from the identical body and `rho=0`.

### A — relinearized adjoint

At every iteration recompute the exact adjoint at the current conductances. For candidate bonds compute

```text
g_e = (k_arbor-k_bath) dC/dk_e.
```

Normalize by the maximum absolute candidate gradient and take a projected step

```text
rho <- clip(rho + eta * g/max|g|, 0, 1).
```

### B — frozen adjoint

Compute `g` only once at the initial body and reuse it for every iteration. Same normalization, step size, and projection.

### C — shuffled relinearized control

Recompute the current adjoint every iteration, but randomly permute the candidate gradients across candidate bonds before applying the same projected update. This preserves the instantaneous gradient-value distribution while destroying spatial credit assignment.

The shuffle RNG is seed-fixed.

## Fixed optimization schedule

```text
eta        = 0.01
iterations = 40
max candidate additions = 8
```

No line search, momentum, Adam, material penalty, or post-hoc tuning.

## Registered discovery predictions

### D1 — relinearized adjoint improves the objective

PASS if:

```text
mean final Delta C_lin > 0.005
and at least 10/12 bodies improve.
```

### D2 — updating the sensitivity matters

Compare final improvement body-by-body:

```text
Delta C_relinearized - Delta C_frozen.
```

PASS if mean paired difference `>0.003` and at least `9/12` bodies favor relinearization.

### D3 — spatial credit assignment matters

Compare

```text
Delta C_relinearized - Delta C_shuffled.
```

PASS if mean paired difference `>0.003` and at least `9/12` bodies favor the true adjoint assignment.

### D4 — the local rule usually climbs rather than wanders

For each body compute the fraction of relinearized iterations with

```text
C_{t+1} >= C_t - 1e-8.
```

PASS if the median body fraction is `>0.80`.

## Descriptive outputs

Also report:

- starting and final `C_lin` for all arms;
- final candidate `rho` values and total frontier material `sum rho`;
- objective trajectories;
- fraction of candidate bonds active (`rho>1e-4`) and fraction strongly matured (`rho>0.25`);
- correlations between the initial gradient and final `rho` in the relinearized arm;
- how often an initially positive candidate later acquires a negative gradient and is backed away from.

## Interpretation fixed in advance

If D1-D4 pass, the old credit-assignment wall is solved **for this model's continuous frontier state space**: a soma-conditioned backward sensitivity field can guide multiple local conductance variables through repeated small updates, while stale or spatially shuffled credit cannot.

This would still not establish a biological mechanism. The exact adjoint is an algorithmic backward pass. The next biological/physical question would be whether reciprocity or a retrograde wave can approximate it.

If D1 passes but D2 fails, a single initial sensitivity map is already sufficient over this short range. If D1 passes but D3 fails, the candidate frontier may be too homogeneous for spatial assignment to matter. If D1 fails, the one-bond local success does not extend to joint optimization under this simple update rule.
