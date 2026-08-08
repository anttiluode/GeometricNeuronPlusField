# Kernel-alignment discovery preregistration v0.1

## Question

Every occupied cell `x` has a geometry-only modal sampling vector

```text
v_x = [phi_0(x), ..., phi_{N-1}(x)]
```

and therefore a rank-1 point-power kernel

```text
K_x = v_x v_x^T.
```

The source-specific transfer histories define an order-sensitive cross-source interaction operator in modal coordinates. This experiment asks:

> **Is the soma's point-power kernel unusually well aligned with that task interaction, and does this alignment explain selectivity beyond amplitude balance alone?**

This is a frozen-body diagnostic. No growth, credit, readout tuning, source relocation, or active boundary is added.

## Bodies

Discovery uses previously unused FunctionalArbor bodies:

```text
seeds 108-119
```

with the same v0.9 bootstrap, target lag `20`, and `210` probe steps used by the recent transfer/mode-pair experiments.

## Task interaction operator

Record the two single-source complex fields and project them into the full graph-Laplacian basis:

```text
q_A(t), q_B(t).
```

For target A->B and distractor B->A define at each time

```text
H_T(t) = q_A(t) q_B(t-tau)^* + q_B(t-tau) q_A(t)^*
H_D(t) = q_B(t) q_A(t-tau)^* + q_A(t-tau) q_B(t)^*
D(t)   = H_T(t) - H_D(t).
```

For a real graph eigenbasis, the order-sensitive cross-source power difference seen by cell `x` is

```text
Delta_cross_x(t) = Re[v_x^T D(t) v_x].
```

To separate kernel orientation from the global instantaneous interaction magnitude, normalize each valid time slice by the Frobenius norm of `D(t)` and define

```text
A_x = RMS_t( Delta_cross_x(t) / ||D(t)||_F ).
```

`A_x` is the preregistered **kernel-alignment score**.

## Other maps

Historical coherent temporal-order contrast is reconstructed directly from the two single-source histories:

```text
psi_T(x,t) = h_A(x,t) + h_B(x,t-tau)
psi_D(x,t) = h_B(x,t) + h_A(x,t-tau)
C_x = (max |psi_T|^2 - max |psi_D|^2) /
      (max |psi_T|^2 + max |psi_D|^2).
```

Amplitude balance is the existing SomaWhyClaude variable

```text
B_x = min(P_A(x),P_B(x)) / max(P_A(x),P_B(x))
P_A = max_t |h_A|^2,  P_B = max_t |h_B|^2.
```

A simple combined opportunity/alignment score is fixed before looking:

```text
Q_x = B_x * A_x.
```

No fitted weights are allowed in this discovery probe.

## Registered discovery predictions

### D1 — kernel alignment tracks the task map

Within each body compute Pearson correlation across occupied cells:

```text
r_A = corr(A_x, |C_x|).
```

Prediction: mean `r_A > 0.50` and at least `9/12` bodies have `r_A > 0`.

### D2 — soma is high in the alignment landscape

Rank the soma's `A_x` among all occupied cells.

Prediction: median soma alignment percentile `> 0.75`, with at least `9/12` bodies above the median (`>0.5`).

### D3 — alignment supplies information missing from amplitude balance

For each body compute

```text
r_B = corr(B_x, |C_x|)
r_Q = corr(Q_x, |C_x|).
```

Prediction: mean `(r_Q - r_B) > 0` and at least `8/12` bodies improve.

This is deliberately a discovery threshold, not a final confirmation claim. If the qualitative structure survives, exact held-out confirmation criteria will be frozen before new bodies are run.

## Descriptive controls

Also report:

- soma percentiles for `|C|`, `B`, and `Q`;
- `corr(A,B)` to quantify redundancy;
- correlations of the unnormalized cross-difference RMS with `|C|`;
- reconstruction equality between modal `Re[v^T D v]` and the direct coordinate-space cross difference;
- field/body size and valid-time counts.

## Failure interpretation

- If `A` does not track `|C|`, the mode-pair result does not explain why one local point is privileged.
- If `A` tracks `|C|` but soma is ordinary in `A`, the soma advantage remains primarily a convergence/balance construction effect.
- If `Q` does not improve over `B`, amplitude balance and kernel alignment may be redundant descriptions of the same local geometry rather than complementary factors.
