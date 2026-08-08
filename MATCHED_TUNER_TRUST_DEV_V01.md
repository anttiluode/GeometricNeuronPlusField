# Scale-invariant matched-tuner development v0.1 — local and spectral coordinates tie at small trust radius, then diverge

`SPECTRAL_TO_LOCAL_COMPILER_V01.md` exposed the normalization flaw in the original `eta=0.01` benchmark: equal numerical coordinate steps were not equal functional steps.

`matched_tuner_trust.py` fixes that by normalizing each selected coordinate by the RMS norm of its current three-lag task tangent and scaling the joint update so both arms receive the same **predicted task-space RMS trust radius** per iteration.

This is development on reused bodies 240-243 only.

## Fixed development sweep

```text
P = 8 coordinates
40 iterations
train lags = 16,20,24
test lags  = 14,18,22,26

delta = .001, .0025, .005, .010
```

At every step, both arms achieved the requested pre-clipping predicted RMS task move to floating-point precision.

## Result

```text
delta      graph test    free test     graph-free    graph wins
.001       0.18558       0.18002       +0.00557      3/4
.0025      0.21753       0.21846       -0.00093      2/4
.005       0.24431       0.27331       -0.02900      1/4
.010       0.22197       0.33246       -0.11049      0/4
```

Mean held-out improvement over the same base:

```text
delta      graph Delta   free Delta
.001       +0.03708      +0.03151
.0025      +0.06903      +0.06996
.005       +0.09580      +0.12480
.010       +0.07347      +0.18396
```

## The important correction

At small trust radii the supposed free-coordinate superiority **vanishes**.

At `delta=.001`, the graph is slightly ahead on average and wins 3/4 bodies. At `.0025`, the arms are essentially tied.

The large free advantage appears only as the step size grows.

This matches the independent compiler result: on the 25-lag task manifold, the tangent shape of 90.6% of winning free spectral directions can be shadowed to <10% relative error by one feasible local bond direction.

So the local physical coordinates are not pointing in fundamentally worse task directions near the current state.

## Why the curves separate

The structural response studies already told us that local bond coordinates are highly curved:

- exact adjoint prediction degrades quickly as a bond moves away from the base state;
- random frontier response curves have slope reversals in 87.5% of events;
- intermediate coupling optima are common.

Direct pole/residue coordinates are much closer to the natural coordinates of the linear dynamics. They can move farther before the intended task-space direction is scrambled by coordinated changes in many modal quantities.

Thus the emerging interpretation is:

> **Local bond coordinates and free spectral coordinates have similar useful tangent directions, but very different curvature away from the tangent point.**

The free coordinates are straighter optimization coordinates. The local coordinates are natural physical actuators whose effects bend rapidly through spectral space.

## Connection to the original alpha sweep

This makes the alpha sweep much more than a warning that binary edits are “too large.” It measured the curvature scale of the physical coordinate chart.

The benchmark now sees the same phenomenon at the system level:

```text
small task-space move
    local ~= spectral

larger task-space move
    local reparameterization bends / hits nonmonotonic response
    direct spectral coordinates remain better aligned
```

That is a coherent geometric story.

## Fresh confirmation

`MATCHED_TUNER_TRUST_CONFIRM_PREREG_V01.md` freezes the entire four-radius crossover on fresh bodies 312-323. It does not cherry-pick the radius at which either arm wins.
