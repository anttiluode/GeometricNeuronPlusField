# Spectral gradient learning v0.1 — 8–16 boundary-selected phasors preserve closed-loop training

`SPECTRAL_CORRELATION_COMPRESSION_V01.md` established that the exact local time-domain forward×adjoint correlation is spectrally sparse on this task. The held-out map test was strong, but map correlation alone does not guarantee useful optimization after the structure moves.

This experiment therefore put the compressed gradient inside the learner.

## Frozen setup

Fresh bodies: seeds **412-423**.

Every body used the same deterministic frontier candidate set in all arms:

```text
up to 8 graded frontier bonds
eta = .01
40 relinearized iterations
lag = 20
210-frame task window
```

Arms:

```text
exact    full discrete adjoint
K4       4 boundary-selected spectral bins
K8       8 boundary-selected spectral bins
K16      16 boundary-selected spectral bins
```

The K arms never use an internal oracle to choose bins. At every iteration they recompute the common frequency ranking from the current arm's **external/source spectrum × soma-return-source spectrum**, build a local gradient using only those bins, update, and relinearize again.

`SPECTRAL_GRADIENT_LEARNING_PREREG_V01.md` froze all criteria before these bodies ran.

## C0 — exact learner positive control

Registered:

```text
mean exact DeltaC > .015
>= 10 / 12 bodies improve
```

Observed:

```text
mean exact DeltaC              +.028787
improved                         11 / 12
```

Seed 420 was effectively flat rather than positive.

**C0 PASS.**

## C1 — K8 actually learns

Registered:

```text
mean K8 DeltaC > .015
>= 10 / 12 bodies improve
```

Observed:

```text
mean K8 DeltaC                 +.024688
improved                         12 / 12
```

**C1 PASS.**

## C2 — K8 preserves most of the exact learner's gain

Registered on group means:

```text
mean(K8 DeltaC) / mean(exact DeltaC) >= .85
mean(K8 - exact DeltaC) > -.007
```

Observed:

```text
group gain ratio                 .8576
mean K8 - exact                -.004099
```

**C2 PASS.**

This is close to the registered 85% boundary rather than a giant-margin result. The threshold should not be rewritten after seeing it.

## C3 — K16 preserves most exact learning

Observed:

```text
mean K16 DeltaC                +.030252
mean exact DeltaC              +.028787
group gain ratio                1.0509
mean K16 - exact              +.001464
improved                         11 / 12
```

Registered requirements were gain ratio `>= .85` and mean difference `> -.007`.

**C3 PASS.**

K16 happens to beat the exact arm in mean final gain because the approximate trajectory can take a different normalized path; seed 423 contributes strongly (`+.0785` K16 versus `+.0648` exact). This is **not** evidence that an approximate gradient is intrinsically better than the exact gradient. The optimizer uses max-normalized finite steps and box constraints, so slightly different directions can land on different finite-step trajectories.

## C4 — map fidelity survives structural motion

Registered:

```text
K8  mean map correlation > .980
K16 mean map correlation > .990
```

Observed across all relinearization steps:

```text
K8  mean map correlation          .99254
     mean relative L2             .11191

K16 mean map correlation          .99720
     mean relative L2             .06681
```

**C4 PASS.**

So the sparse representation is not only a base-state coincidence. It continues to track the moving local gradient while conductances change.

## Formal verdict

```text
C0 PASS   exact learner remains positive
C1 PASS   K8 learns on all 12 fresh bodies
C2 PASS   K8 retains 85.8% of exact mean gain
C3 PASS   K16 retains/exceeds exact mean gain in this run
C4 PASS   compressed map fidelity remains high during learning
```

**5 / 5 registered criteria pass.**

## Descriptive stress arm: K4

Four phasors are already useful but clearly less robust:

```text
mean K4 DeltaC                  +.017644
improved                         11 / 12
mean map correlation              .96886
mean relative L2                  .25037
```

One body, seed 423, is especially damaging to K4: exact gain `+.0648`, K4 only `+.0101`. That is a useful warning against turning the spectacular spectral sparsity into a universal four-number story.

## The cost picture now

The original exact time-domain implementation effectively asked every tunable bond to retain the whole 210-frame forward history until the reciprocal credit replay arrived.

The confirmed alternative is:

```text
FORWARD TASK
  each tuner accumulates only K selected complex lock-in / DFT phasors
                     ↓
SOMA OBJECTIVE
  creates time-reversed reciprocal return source
                     ↓
RETRO REPLAY
  each tuner accumulates the matching K return phasors
                     ↓
LOCAL PHASOR PRODUCTS
  sum K contributions -> approximate local gradient
                     ↓
GRADED UPDATE
```

For this registered task:

```text
K=8   12/12 learn, 85.8% of exact group gain
K=16  11/12 improve, mean gain slightly above exact in this finite-step run
```

That is now a **training** result, not just a gradient-map reconstruction result.

## What remains before calling it an in-situ hardware protocol

The current script computes the complex phasor products numerically after accumulating them. A physical implementation still needs to decide how those products are measured or multiplied locally.

The coherent-intensity identity is straightforward if the required two local phasors can be superposed:

```text
|U + exp(i theta)V|^2 - |U - exp(i theta)V|^2
    = 4 Re[exp(-i theta) conj(U) V].
```

That is the same square-law/interference trick that underlies photonic in-situ gradient measurement. But our broadband finite-time setting adds a real systems question: the forward phasor was accumulated during the earlier task pass, whereas the matching return phasor arrives later. Either a small analog phasor memory/multiplier is retained at each tuner, or the selected spectral components must be physically replayed/co-propagated in an equivalent phase-stepped protocol.

So the old 210-sample memory wall is gone, but it has become a much smaller and more concrete **K-phasor local state / coherent-readout problem**.

## Relation to existing physical backpropagation

This architecture now sits directly beside the in-situ-backpropagation literature rather than merely resembling it. Hughes et al. derived exact photonic gradients through forward/adjoint interference, and Pai et al. experimentally trained a silicon photonic mesh using forward/backward optical interference.

Our distinct wrinkle is the **broadband transient task**: the exact reciprocal adjoint is time-domain, and this experiment shows that its delayed local correlation can be compressed to a small port-selected spectral bank without destroying learning.

That is the piece to compare against real device noise, coherent phase error, attenuation, detector noise, and drift next.

## Wall sentence

> **The reciprocal mesh can be trained without storing a full local time trace: on held-out bodies, eight boundary-selected spectral correlation channels retain most of the exact adjoint learner's gain, and sixteen essentially preserve it.**
