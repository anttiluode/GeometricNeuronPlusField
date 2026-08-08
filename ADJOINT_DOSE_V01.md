# Adjoint dose v0.1 — the local sensitivity is right; the binary edit is too large

This result independently converges with `AdjointClaude/adjoint_alpha_result.txt`.

Claude's sweep on eight bodies found the adjoint nearly exact for small partial bond maturation and useless or sign-reversed by the full bath-to-arbor jump:

```text
alpha      mean corr
0.001      +1.000
0.010      +0.974
0.050      +0.544
0.200      -0.309
1.000      -0.428
```

Our independent implementation first reproduced the local-adjoint identity, then extended the dose grid and froze a held-out confirmation before new bodies were run.

## Important wording correction

The result does **not** prove that "no first-order local rule can ever predict a discrete structural edit." What it establishes in this model is narrower and cleaner:

> **The exact first derivative at the current conductance state cannot be extrapolated across the full `k_bath -> k_arbor` jump.**

The local gradient is correct. The finite binary event lies far outside its radius of validity.

## Corrected exact linear surrogate

The second sweep also corrected a tiny omission in the first independent surrogate: FunctionalArbor's shift-based Laplacian weakly couples the four outer boundaries to zero-valued exterior bath. Including that term reduced the custom-vs-FunctionalArbor Laplacian discrepancy to about `1e-11` relative.

The adjoint's central finite-difference derivative remained exact to about `1e-8` relative on held-out bodies.

## Discovery — seeds 192-203

Fixed dose grid:

```text
1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1
```

Mean within-body event correlation between the **single base-state adjoint prediction** and the exact finite linear event at each dose:

```text
alpha        mean r
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

Median largest tested `alpha` with body-level `r >= .70` was `0.01`.

The exact full-step **linear** structural changes still matched the nonlinear integrated-interference counterfactual very strongly (`mean r = 0.9810`). That isolated the failure to first-order extrapolation, not to the linear wave surrogate.

All registered discovery criteria passed.

## Held-out confirmation — seeds 204-215

The confirmation criteria were frozen before these bodies were run.

### C0 — implementation identity

```text
mean Laplacian relative error      1.27e-11
mean adjoint FD relative error     1.14e-08
```

**PASS.**

### C1 — alpha = 0.001 remains a clean gradient regime

```text
mean r                              0.99684
median r                            0.99953
positive bodies                    12 / 12
```

**PASS.**

### C2 — alpha = 0.01 is still practically useful

```text
mean r                              0.85405
median r                            0.94456
positive bodies                    12 / 12
mean sign agreement                0.95139
```

**PASS.**

### C3 — the base derivative collapses before the binary endpoint

```text
alpha 0.1 mean r                  -0.42012
alpha 1.0 mean r                  -0.49157
r(0.001) - r(0.1)                 1.41696
```

**PASS.**

### C4 — radius scale replicates

```text
median largest alpha with r >= .70    0.01
```

The registered acceptable interval was `[0.003, 0.10]`.

**PASS.**

### C5 — exact finite linear events still match nonlinear interference events

```text
mean corr(dC_lin(alpha=1), dC_int)     0.96324
positive bodies                        12 / 12
```

**PASS.**

So the full confirmation is C0-C5 PASS.

## The dose curve itself is the result

Held-out mean correlations:

```text
alpha        mean r      mean sign agreement
1e-4         0.99996     1.000
3e-4         0.99966     1.000
1e-3         0.99684     0.986
3e-3         0.98235     0.986
1e-2         0.85405     0.951
3e-2         0.05990     0.806
1e-1        -0.42012     0.750
3e-1        -0.65214     0.694
1           -0.49157     0.424
```

The useful local neighborhood is not merely numerical epsilon. It extends to roughly one percent of the full bath-to-arbor conductance change, sometimes farther, but it is nowhere near the binary endpoint.

## One more clue hidden in the held-out dose data

A naive next sentence would be: "therefore use continuous density relaxation and binarize at the end."

The data warn that the last clause is not yet earned.

Across the 144 held-out events, the sign of the exact finite response relative to the base derivative had reversed in:

```text
alpha 0.01       4.9%
alpha 0.03      19.4%
alpha 0.10      25.0%
alpha 0.30      30.6%
alpha 1.00      57.6%
```

Among the 66 events whose **base derivative points in the improving direction**, 37/66 ended with a *negative* full binary change. Their mean best sampled partial-conductance gain was positive (`+0.00629`), while their mean binary-endpoint change was negative (`-0.00741`).

For gradient-favored **additions** specifically:

```text
29 events total
26 / 29 had their best sampled conductance strictly inside (0,1)
18 / 29 became harmful by alpha=1
median best sampled alpha = 0.10
```

This suggests something stronger than "take smaller steps":

> **The computation may actually prefer graded bond strengths. Forcing a useful partially matured bond all the way to a binary arbor bond can undo the gain.**

That is now the next question.

## Wall sentence

> **The soma-conditioned adjoint gives the correct local structural derivative, but anatomy is strongly nonlinear over a full binary bond change. The system has a real differential regime around the current geometry; beyond it, the objective bends, often turns, and frequently reverses sign.**

## Next clean experiment

Before adopting SIMP-style "optimize continuously, binarize later," map the conductance-response curve densely on fresh bodies and ask whether gradient-favored additions genuinely have reproducible interior optima.

If they do, the model's natural structural variable may be **graded conductance / material density**, not merely binary occupancy.
