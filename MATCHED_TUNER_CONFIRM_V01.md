# Matched-tuner confirmation v0.1 — free spectral coordinates decisively beat local bond tuners on this temporal task

`MATCHED_TUNER_CONFIRM_PREREG_V01.md` froze a negative-capable benchmark after two development audits had favored direct spectral coordinates.

Fresh bodies: seeds **288-299**.

Both arms started from the exact same mature reciprocal transfer function and received:

```text
8 trainable real coordinates
eta = 0.01 normalized projected-gradient step
40 iterations
training lags = 16,20,24
test lags = 14,18,22,26
```

Graph arm `G8` searched every horizontal/vertical bond in the 31x31 medium (1860 physical candidates), allowing weak bonds to strengthen and strong bonds to weaken.

Free arm `F8` searched all direct modal pole-frequency, source-A residue, source-B residue, and soma/output-residue coordinates (`4N = 3844` mathematical candidates for `N=961`) and selected the eight largest initial train-objective gradients.

## C0 — same base simulator

Registered maximum base graph-vs-modal train-objective discrepancy `<1e-8`.

Observed:

```text
max discrepancy = 3.18e-11
```

**C0 PASS.**

The comparison is not caused by different base simulators.

## Main held-out result

Across 12 fresh bodies:

```text
mean base test C                         -0.020657
mean G8 test C                           +0.028492
mean F8 test C                           +0.168527

mean G8 Delta test                       +0.049149
mean F8 Delta test                       +0.189184

mean (F8 Delta - G8 Delta)               +0.140035
F8 beats G8                               12 / 12
```

The preregistered C1 criterion was:

```text
mean (F8 Delta - G8 Delta) > 0.05
and F8 beats G8 in at least 9/12 bodies.
```

**C1 PASS, decisively.**

The paired difference is positive in all 12 bodies (descriptive exact two-sided sign p = 0.000488).

Per-body `G8 - F8` test objective differences:

```text
288  -0.1495
289  -0.1097
290  -0.0913
291  -0.1705
292  -0.0895
293  -0.3078
294  -0.1723
295  -0.0573
296  -0.0881
297  -0.1167
298  -0.2738
299  -0.0539
```

There is no near-tie body in this held-out set.

## C2 — both parameterizations are trainable

Registered sanity gate: each arm improves the held-out test objective in at least 9/12 bodies.

Observed:

```text
G8 improves test objective   12 / 12
F8 improves test objective   12 / 12
```

**C2 PASS.**

So the result is not “local geometry cannot be trained.” It can. The result is that the selected direct spectral coordinates are substantially more effective per **mathematical scalar coordinate** on this fixed-source multi-lag task.

## The claim we now reject

The preregistration said that if C1 passed, the repo would explicitly reject the statement that the current local reciprocal geometry provides superior temporal-task performance per added trainable scalar against this free spectral baseline.

It passed.

> **The current local-bond parameterization does not beat direct free-modal tuning per added mathematical scalar on this task. The free spectral parameterization wins strongly and consistently.**

Do not rescue a graph win by weakening the spectral baseline, changing tuner count, changing learning rate, or moving the lag set after this result.

## What this does and does not settle

It settles the abstract `8 scalars vs 8 scalars` performance question for this benchmark.

It does **not** establish equal physical hardware cost. A direct modal pole or residue coordinate is a globally defined transformation in the spatial basis, whereas a graph bond tuner is one local physical element. `spectral_coordinate_locality.py` separately audits that distinction.

Therefore the hardware case, if there is one, has moved to a different axis:

```text
locality of control
number of physical actuators needed to realize a spectral coordinate
in-situ gradient sensing cost
fabrication complexity
robustness to drift / defects
energy and latency
multi-port / spatial generalization
```

Those are not loopholes in this result. They are simply different benchmarks.

## Architectural conclusion

The strongest clean picture is now:

1. the reciprocal spatial wave is exactly representable as a pole/residue bank;
2. local bond tuning is a constrained structured parameterization of that bank;
3. the constraint is a disadvantage on the present fixed-source temporal objective when cost is counted only as number of abstract trainable scalars;
4. the reason to build the spatial medium must therefore come from **physical implementation economics or useful inductive bias under a task where locality matters**, not from superior unconstrained spectral expressivity.

That is a useful negative result because it tells us exactly where the next serious hardware question begins.
