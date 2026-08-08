# Scale-invariant matched-tuner confirmation preregistration v0.1

`SPECTRAL_TO_LOCAL_COMPILER_V01.md` showed that the original matched-tuner benchmark used the same numerical `eta` in differently scaled coordinate systems. The free spectral arm therefore often took much larger functional steps. `matched_tuner_trust.py` removes that ambiguity by matching the **predicted RMS change of the multi-lag training-output vector** at every iteration.

A reused-seed development sweep (240-243) showed a specific crossover: local and spectral coordinates were nearly tied at small task-space trust radii, while spectral coordinates became much better at larger radii.

This fresh run tests that shape rather than choosing one favorable radius after the fact.

## Frozen protocol

Fresh bodies: seeds **312-323**.

Both arms:

```text
P = 8 trainable real coordinates
40 relinearized iterations
training lags = 16,20,24
test lags = 14,18,22,26
```

Trust-radius ladder, fixed before fresh bodies:

```text
delta = 0.001, 0.0025, 0.005, 0.010
```

At each iteration and for each arm:

1. compute the Jacobian of the three training-lag contrasts with respect to the selected coordinates;
2. normalize each coordinate by the RMS norm of its current task-space tangent;
3. take the ascent direction in that normalized coordinate system;
4. scale the whole update so the **predicted RMS change of the three training outputs equals `delta`** before box clipping;
5. apply declared parameter bounds;
6. recompute the field/Jacobian.

This update is invariant to independent positive rescaling of coordinate units.

### G8 local arm

Candidate pool: all 1860 nearest-neighbour bonds in the 31x31 medium. At the binary base state weak bonds may strengthen and strong bonds may weaken. Select 8 once by normalized improving alignment with the mean training objective, subject to enough available parameter range to support one trust-radius step.

### F8 spectral arm

Candidate pool: direct modal

```text
F pole/stiffness
A source-A residue
B source-B residue
C soma/output residue
```

coordinates of the exact same base operator. Select 8 once by the same normalized improving alignment criterion and one-step feasibility rule.

## Registered outcomes

For each delta:

```text
D(delta) = mean_bodies [ C_free(test) - C_graph(test) ].
```

### C0 — base identity

Maximum graph/modal base training-objective discrepancy `<1e-8`.

### C1 — small-trust parity

At `delta=0.001`:

```text
|D(0.001)| < 0.025.
```

This tests the compiler result that the useful local and free tangent directions are often nearly collinear once step scale is matched.

### C2 — intermediate-trust parity remains approximate

At `delta=0.0025`:

```text
|D(0.0025)| < 0.035.
```

### C3 — large-trust spectral advantage

At `delta=0.010`:

```text
D(0.010) > 0.050
and free beats graph in at least 9/12 completed bodies.
```

### C4 — crossover magnitude

```text
D(0.010) - D(0.001) > 0.050.
```

This is the main curvature prediction: coordinate choice matters increasingly as the optimizer leaves the local differential neighborhood.

### C5 — both local and spectral coordinates remain useful

At `delta=0.001` both arms must have positive mean held-out improvement over base. This is a sanity condition, not a superiority claim.

## Interpretation if the crossover confirms

The result will not be summarized as “graph wins” or “spectral wins.” The intended statement is:

> **At matched infinitesimal task-space step size, local physical bond coordinates and direct spectral coordinates are similarly effective on this task. Spectral coordinates gain a strong advantage only as the permitted trust radius grows, consistent with the much stronger curvature/nonmonotonicity of local structural response.**

That would connect three previously separate observations:

```text
adjoint exact only locally for structural change
random bond response curves are strongly nonmonotone
local task tangents shadow spectral task tangents but at different scales
```

into one geometric statement about optimization coordinates.

No thresholds, deltas, tuner counts, or lag sets will be changed on seeds 312-323 after inspection.
