# Graded-bond response discovery preregistration v0.1

## Question

`ADJOINT_DOSE_V01.md` confirmed that the soma-conditioned adjoint is accurate for small conductance changes but frequently points the wrong way by the time a bond is forced all the way from bath conductance to arbor conductance.

The held-out dose data contained a stronger descriptive clue: among gradient-favored additions, 26/29 had their best sampled conductance strictly inside `(0,1)` and 18/29 became harmful by the binary endpoint.

> **Is a partially matured bond itself a reproducible computational optimum, rather than merely a temporary step on the way to a binary bond?**

This experiment is frozen before fresh bodies are run.

## Bodies

Fresh FunctionalArbor bodies:

```text
seeds 216-227
```

Same mature wave coefficients, fixed source terminals, lag 20, 210-step target/distractor trajectories, and legal event generator as the adjoint-dose experiments.

## Event selection

For each body generate up to six legal tip-like additions and six legal safe deletions.

Compute the exact base-state adjoint. An event is called **gradient-favored** if its directional base derivative toward the corresponding binary event is positive:

```text
g_e = grad C_lin · DeltaK_e > 0.
```

The registered primary analysis is on gradient-favored additions. Gradient-favored deletions are a secondary control because earlier descriptive data suggested additions show the stronger interior-optimum effect.

## Fixed conductance grid

For every event evaluate the exact linear objective at

```text
alpha =
0,
0.005, 0.01, 0.02, 0.03, 0.05, 0.075,
0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60,
0.75, 1.00
```

with

```text
k(alpha) = k_bath + alpha (k_arbor-k_bath)
```

for additions, and the corresponding reverse interpolation for deletions.

No gradient is recomputed along the path. This is a response-curve audit of one fixed structural degree of freedom.

## Registered discovery predictions

### D1 — gradient-favored additions commonly prefer an interior conductance

For each gradient-favored addition, let `alpha_best` maximize exact `C_lin(alpha)` on the fixed grid.

PASS if more than `70%` have

```text
0 < alpha_best < 1.
```

### D2 — forcing the bond to binary frequently destroys the locally predicted gain

PASS if more than `40%` of gradient-favored additions have

```text
Delta C_lin(alpha=1) < 0.
```

although their base derivative is positive by selection.

### D3 — stopping at the best graded conductance materially outperforms forced binarization

For each gradient-favored addition define

```text
regret_binary = max_alpha Delta C_lin(alpha) - Delta C_lin(1).
```

PASS if:

```text
mean regret_binary > 0.005
and at least 75% of gradient-favored additions have regret_binary > 1e-5.
```

### D4 — the optimum is not concentrated at numerical epsilon

Among gradient-favored additions whose optimum is interior, PASS if median `alpha_best >= 0.03`.

This distinguishes a real graded-coupling optimum from a trivial infinitesimal improvement.

## Secondary/control outputs

Report descriptively:

- the same four quantities for gradient-favored deletions;
- distribution of `alpha_best`;
- fraction of response curves with at least one finite-difference slope sign reversal;
- full-endpoint relation to nonlinear `dC_int` and historical `dC_peak`;
- response-curve examples identified only by seed/cell, without post-hoc cherry-picking into pass/fail claims.

## Interpretation fixed in advance

If D1-D4 pass, the natural design variable in this model is not merely binary occupancy. A bond can function as an **analog coupling / impedance parameter** whose intermediate value is computationally meaningful. In that case, SIMP-style "relax then binarize" is not automatically justified; binarization would require an explicit continuation/penalty and a demonstration that it preserves the task.

If D1 passes but D2/D3 fail, intermediate conductance is useful but mainly as a continuation path toward binary anatomy.

If D1 fails, the previous interior-optimum pattern was unstable and ordinary gradual maturation toward binary bonds remains the simpler interpretation.
