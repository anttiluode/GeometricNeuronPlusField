# Scale-invariant matched-tuner confirmation v0.1 — small-step parity confirms; the crossover confirms, but large-step dominance is weaker than development

`MATCHED_TUNER_TRUST_CONFIRM_PREREG_V01.md` was frozen after the task-space compiler exposed a normalization flaw in the first matched-tuner benchmark. Instead of assigning both coordinate systems the same numerical `eta`, this benchmark gives both arms the same **predicted RMS change of the three training-lag outputs** at every iteration.

Fresh bodies: seeds **312-323**.

```text
P = 8 coordinates
40 relinearized iterations
train lags = 16,20,24
test lags  = 14,18,22,26
```

Trust-radius ladder:

```text
0.001, 0.0025, 0.005, 0.010
```

## C0 — base identity

Registered maximum graph/modal base discrepancy `<1e-8`.

Observed:

```text
2.35e-11
```

**C0 PASS.**

## The held-out curve

```text
delta      graph test    free test     free-graph    free wins
.001       0.06497       0.05019       -0.01478      4/12
.0025      0.09716       0.09527       -0.00189      7/12
.005       0.12830       0.13496       +0.00666      7/12
.010       0.13306       0.17621       +0.04315      8/12
```

Mean improvement over the common base (`C=0.03346`):

```text
delta      graph Delta   free Delta
.001       +0.03151      +0.01673
.0025      +0.06370      +0.06181
.005       +0.09485      +0.10150
.010       +0.09960      +0.14276
```

## C1 — small-trust parity

Registered:

```text
|D(.001)| < .025
```

where `D = free - graph`.

Observed:

```text
D(.001) = -0.01478
```

**C1 PASS.**

At the smallest matched functional step, local bonds are not worse. They are slightly ahead on average and win 8/12 bodies.

## C2 — intermediate-trust parity

Registered:

```text
|D(.0025)| < .035
```

Observed:

```text
D(.0025) = -0.00189
```

**C2 PASS.**

The two parameterizations are essentially tied.

## C3 — strong large-trust spectral advantage

Registered at `delta=.010`:

```text
D > .050
and free wins >= 9/12 bodies.
```

Observed:

```text
D(.010) = +0.04315
free wins = 8/12
```

**C3 FAIL.**

The direction and qualitative development trend replicate, but the preregistered strength criterion does not. Do not promote “large-step spectral dominance” as a held-out result from this experiment.

## C4 — crossover magnitude

Registered:

```text
D(.010) - D(.001) > .050.
```

Observed:

```text
+0.05793
```

**C4 PASS.**

Eight of twelve bodies move in the predicted direction; the median per-body crossover is `+0.0531`.

So coordinate choice becomes materially more favorable to the spectral parameterization as the allowed trust radius grows, even though the endpoint advantage itself missed the stronger C3 threshold.

## C5 — both parameterizations remain useful

At `delta=.001` both arms had positive mean held-out improvement over base.

**C5 PASS.**

## Formal confirmation

```text
C0  PASS   same base transfer
C1  PASS   small-step parity
C2  PASS   intermediate-step parity
C3  FAIL   large-step spectral advantage weaker than threshold
C4  PASS   crossover magnitude
C5  PASS   both learn
```

Strictly: **5/6 including the identity/sanity criterion, or 4/5 substantive registered predictions.**

## What this settles

The original `eta=.01` benchmark cannot be used as a theorem that free spectral coordinates are intrinsically better per scalar. Once the optimizer is invariant to arbitrary coordinate scaling, that advantage disappears in the local regime.

The strongest held-out statement is now:

> **At matched small task-space trust radii, eight local physical bond coordinates and eight direct spectral coordinates are similarly effective on this temporal task. As the permitted functional step grows, performance shifts toward the spectral coordinates, confirming a crossover, but the preregistered large-step dominance threshold was not reached.**

This fits the independent structural evidence:

```text
exact adjoint is a local derivative
random bond response is strongly nonmonotone
interior coupling optima are common
local task tangents can shadow oracle spectral tangents
```

The natural next mechanism question is therefore **curvature**, not another leaderboard:

> Do local-bond coordinates depart from their own first-order task-space tangent faster than direct modal coordinates as the step grows?

That can be measured directly by a finite-step tangent-fidelity experiment.
