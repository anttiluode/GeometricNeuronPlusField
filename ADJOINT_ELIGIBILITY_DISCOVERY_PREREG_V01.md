# Adjoint structural-eligibility discovery preregistration v0.1

## Motivation

`STRUCTURAL_INTERFERENCE_V01.md` identified a causal target but not a local credit signal. One-cell edits change the confirmed integrated interference statistic mainly through the lag-directional compatibility of the two source transfer histories. Exact before/after recomputation is omniscient.

For a reciprocal linear wave operator, however, sensitivity of a scalar soma objective to a local bond change can be computed by an adjoint field. The local derivative is an overlap between the forward field and a task-conditioned field propagated backward from the readout.

> **Can an adjoint sensitivity computed on the unedited arbor predict which one-cell structural events will improve or worsen the soma interference objective?**

This is a computational mechanism test, not a claim that real neurons perform mathematical backpropagation.

## Bodies and events

Fresh FunctionalArbor bodies only:

```text
seeds 180-191
```

Use the exact event generator from `structural_interference_probe.py`:

- up to six legal one-cell tip-like additions;
- up to six safe one-cell deletions;
- soma and exact source terminals fixed;
- lag 20;
- 210 steps.

The observed nonlinear event effects `dC_int` and `dC_peak` are measured exactly as in the structural-interference experiment.

## Linearized mature wave

The adjoint model uses the same mature bond conductances, `dt`, stiffness, damping, restoring term, source pulses and source locations as FunctionalArbor, but freezes out the weak saturation nonlinearity:

```text
v_{t+1}   = v_t + dt [ stiffness L_K psi_t
                       - damping v_t
                       - restoring psi_t
                       + source_t ]
psi_{t+1} = psi_t + dt v_{t+1}.
```

`L_K` is built from the exact mature arbor/bath bond weights.

The scalar adjoint objective is the zero-padded integrated target-vs-distractor soma-energy contrast

```text
C_lin = (E_T - E_D) / (E_T + E_D).
```

For the linear system this is the same class of coherent interference statistic as `C_int`.

## Adjoint derivative

A reverse pass through the linear recurrence returns the gradient of `C_lin` with respect to every local bond conductance.

For a bond `(i,j)` the instantaneous local contribution has the form

```text
2 dt * stiffness * Re[ conj(mu_i - mu_j) * (psi_j - psi_i) ]
```

summed over time and over the target/distractor trajectories, where `mu` is the task-conditioned adjoint state emitted from the soma objective.

For a one-cell addition, the predicted event effect is the gradient of its single new arbor bond times `(k_arbor-k_bath)`. For a deletion it is the sum over all downgraded incident arbor bonds times `(k_bath-k_arbor)`.

No after-edit simulation is used to compute the adjoint prediction.

## Validation layers

For each event also compute the exact finite change of `C_lin` by actually changing the relevant bonds in the **linear** simulator. This separates:

1. adjoint first-order error;
2. linearization/finite-event error;
3. mismatch to the historical peak task.

A forward-only control is fixed in advance:

```text
F_event = Delta k * sum_edges sum_t
          ( |psi_T,i-psi_T,j|^2 - |psi_D,i-psi_D,j|^2 ).
```

It sees local task-conditioned forward activity but has no backward/readout sensitivity.

## Registered discovery predictions

### D0 — implementation check

A small-`epsilon` finite-difference check of one base bond per body must match the adjoint gradient with mean relative error `<1e-3`.

### D1 — adjoint predicts exact finite changes in the linear system

Within each body:

```text
r_adj_lin = corr(pred_adjoint, dC_lin_finite).
```

PASS if mean `r_adj_lin > 0.70` and at least `10/12` bodies are positive.

### D2 — linear finite-event changes track the nonlinear interference counterfactual

```text
r_lin_int = corr(dC_lin_finite, dC_int).
```

PASS if mean `r_lin_int > 0.65` and at least `10/12` bodies are positive.

### D3 — unedited-arbor adjoint predicts nonlinear structural interference changes

```text
r_adj_int = corr(pred_adjoint, dC_int).
```

PASS if mean `r_adj_int > 0.55`, at least `9/12` bodies are positive, and pooled sign agreement for non-negligible `dC_int` events is `>0.70`.

### D4 — the backward/readout factor matters

Compare within-body prediction of `dC_lin_finite`:

```text
corr(pred_adjoint, dC_lin_finite)
  - corr(F_event, dC_lin_finite).
```

PASS if mean improvement `>0.15` and at least `9/12` bodies improve.

## Historical task relevance — descriptive

Report, without a discovery pass/fail threshold:

- correlations/sign agreement of adjoint prediction with `dC_peak`;
- additions versus deletions;
- effect versus graph distance;
- first-order performance separately for additions and deletions.

## Interpretation

If D0-D4 pass, the old eligibility wall has a principled solution **in the model**: the correct local structural sensitivity is not recent activity alone but a forward-field × soma-adjoint overlap. The remaining biological/computational question would then be whether the explicit adjoint can be approximated by a physically local retrograde signal.

If D1 passes but D2/D3 fail, the adjoint mathematics is correct but the linear surrogate is not adequate for the actual structural events.

If D4 fails, a simpler forward-only task trace may already contain most of the needed information and the adjoint interpretation should not be privileged.
