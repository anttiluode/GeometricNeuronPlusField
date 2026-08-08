# Matched-tuner confirmation preregistration v0.1

## Why this is a negative-capable benchmark

`POLE_CONTAINMENT_V01.md` established that the complete frozen reciprocal wave is itself exactly a pole/residue bank. Therefore a fully free spectral model is a superset in abstract function space. The only plausible advantage for local geometry is **structured parameterization per physical tuner**, not greater linear expressivity.

Two development comparisons on reused seeds 240-243 both favored direct spectral coordinates. The second symmetry audit gave the graph access to every local bond in the 31x31 medium and strengthened the free arm with soma/output-residue coordinates; the free arm still won 4/4 development bodies.

The held-out experiment is therefore explicitly allowed to confirm a **negative architectural result**.

## Frozen protocol

Fresh bodies: seeds **288-299** (12 requested).

Task:

```text
training lags      16, 20, 24
test lags          14, 18, 22, 26
probe steps        210
```

Both arms start from the exact same mature reciprocal transfer function.

Both receive:

```text
P = 8 trainable real tuners
eta = 0.01 normalized projected-gradient step
40 iterations
same train objective
same test lags
```

### G8 — local physical coupling tuners

Candidate pool is every horizontal/vertical bond in the exact 31x31 weighted medium (`1860` candidates).

At the base state:
- weak bath bonds can strengthen;
- strong arbor bonds can weaken.

Rank by projected feasible first-order improvement and select the top 8 once. Their normalized conductances are then jointly optimized for 40 relinearized adjoint steps with clipping to `[0,1]`.

### F8 — direct free-modal tuners

Diagonalize the exact same base operator. Candidate coordinates are, for every mode:

```text
F   log modal stiffness / pole-frequency coordinate
A   source-A modal residue coordinate
B   source-B modal residue coordinate
C   soma/output modal residue coordinate
```

This gives `4N` candidate scalar coordinates (`3844` for the 961-state system). Select the 8 largest initial absolute train-objective gradients once, then jointly optimize those coordinates for the same 40 normalized steps.

This baseline is intentionally strong. It is not required to correspond to eight spatially local hardware knobs; it represents the unconstrained spectral alternative against which any claimed geometry-per-tuner advantage must survive.

## Registered outcomes

Primary quantity:

```text
Delta_test_arm = mean C_arm(test lags) - mean C_base(test lags)
```

and paired difference

```text
D = Delta_test_free - Delta_test_graph.
```

### C0 — simulator identity

Maximum base graph-vs-modal train-objective discrepancy `< 1e-8`.

### C1 — direct spectral coordinates beat local geometry on held-out lags

Confirm if:

```text
mean D > 0.05
and F8 beats G8 in at least 9/12 completed bodies.
```

### C2 — both parameterizations are trainable

Descriptive sanity gate, not needed for C1:

```text
G8 improves held-out test objective in at least 9/12 bodies
F8 improves held-out test objective in at least 9/12 bodies.
```

### C3 — no graph-performance claim if C1 passes

If C1 passes, the repo will explicitly reject the statement that the current local reciprocal geometry provides superior temporal-task performance per added trainable scalar against this free spectral baseline.

If C1 fails, we report the fixed result; no alternate baseline, tuner count, learning rate, or lag set will be substituted on these seeds.

## Interpretation boundary

Even a strong F8 win does **not** invalidate the hardware case for the spatial medium. F8 coordinates are globally defined spectral knobs. A real implementation may make them much more expensive or impossible to tune independently.

A free-modal win would instead move the hardware question to:

```text
local fabrication and control
in-situ gradient measurement
robustness to defects/drift
energy/latency
wiring/monitor cost
```

Those require explicit hardware-cost models rather than abstract parameter counting.
