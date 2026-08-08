# Graded-bond response held-out confirmation preregistration v0.1

Discovery on fresh seeds 216-227 strongly supported the hypothesis that a locally improving new bond often has a **non-binary optimum**.

Discovery, gradient-favored additions only:

```text
n                                      31
interior best alpha                    26 / 31 = 0.8387
harmful at alpha=1                     17 / 31 = 0.5484
mean binary regret                     0.01795
positive binary regret                 26 / 31 = 0.8387
median alpha_best among interior       0.04
response curves with slope reversal    28 / 31 = 0.9032
mean best partial gain                 +0.00957
mean forced-binary gain                -0.00838
```

All discovery D1-D4 criteria passed. This confirmation is frozen before any new body is run.

## Held-out bodies

```text
seeds 228-239
```

Use exactly the same event generator, exact linear wave, base-state adjoint, lag 20, 210 steps, and fixed conductance grid from `bond_response_probe.py`. No gradient recomputation along the one-bond path.

Primary analysis remains **gradient-favored additions** (`g_e > 0` at alpha=0).

## Registered confirmation criteria

### C1 — interior optima replicate

PASS if more than `75%` of gradient-favored additions have

```text
0 < alpha_best < 1.
```

### C2 — binary forcing commonly reverses the locally improving direction

PASS if more than `45%` of gradient-favored additions have

```text
Delta C_lin(alpha=1) < 0.
```

### C3 — binary regret is substantial

PASS if both:

```text
mean [max_alpha Delta C_lin(alpha) - Delta C_lin(1)] > 0.010
at least 75% of favored additions have regret > 1e-5.
```

### C4 — the favored interior scale is genuinely graded

Among favored additions whose optimum is interior, PASS if

```text
0.02 <= median alpha_best <= 0.20.
```

This brackets the discovery median (`0.04`) without requiring exact replication of one grid point.

### C5 — nonmonotonicity is the rule rather than a rare exception

PASS if more than `80%` of gradient-favored addition response curves contain at least one finite-difference slope sign reversal on the fixed alpha grid.

## Secondary control

Report the same quantities for gradient-favored deletions, but do not use them to rescue or reject the primary addition hypothesis.

## Interpretation fixed in advance

If C1-C5 pass, the next model should treat structural strength as an analog state variable and **must not assume that continuous optimization can be safely thresholded to binary anatomy afterward**. Binarization would become a separate engineering/biophysical constraint requiring its own continuation or penalty.

If the interior optimum replicates, a useful conceptual reframe is:

```text
weak field coupling -> graded structural maturation -> tuned impedance/coupling -> field computation
```

rather than

```text
absent edge -> instantaneous full edge.
```
