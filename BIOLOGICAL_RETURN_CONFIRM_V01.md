# Coarse biological return-code v0.1 — held-out result

## Result

The exact reciprocal-adjoint soma waveform contains substantial redundancy, but not every coarse coding hypothesis survives held-out testing.

Development bodies 492–497 motivated [`BIOLOGICAL_RETURN_CONFIRM_PREREG_V01.md`](BIOLOGICAL_RETURN_CONFIRM_PREREG_V01.md). The frozen criteria were then evaluated on fresh bodies **498–509**.

Formal result:

```text
7 / 8 registered criteria PASS
1 / 8 FAIL
```

The failed criterion is retained as a failure. No threshold is moved post hoc.

The strongest surviving results are:

1. the **real-valued return waveform** preserves the exact bond-gradient direction almost perfectly;
2. a **fast 50%-duty periodic gate** can discard half of the return samples with almost no loss, even under the worst tested gate phase;
3. a slower gate becomes strongly phase-sensitive;
4. sparse event codes retain useful but body-dependent structural direction;
5. the **sign of consequence is load-bearing** in the sparse limit.

The attractive envelope-only result from development narrowly missed its registered held-out mean threshold.

## Frozen setup

Only the soma return code was altered. The forward field, reciprocal passive arbor, task, and structural coordinates were unchanged.

Every transformed target and distractor return waveform was separately L2-dose-matched to the corresponding exact waveform, so differences cannot be explained by total returned signal energy alone.

Metrics were bond-gradient map correlation, relative L2 error, and strong-coordinate sign agreement.

## R0 — exact positive control

```text
mean map correlation      1.0000000000
mean relative L2          2.65e-15
```

**PASS.**

## R1 — explicit complex quadrature is not needed for map direction

Replace the exact complex return waveform by

```text
Re[g(t)]
```

and renormalize its L2 dose.

Held-out result:

```text
mean map correlation      0.99999857
median correlation        0.99999861
mean strong-sign agreement 1.00000000
mean relative L2          0.07885091
```

Registered requirements were `corr > .995` and strong-sign agreement `> .99`.

**PASS.**

This does **not** mean the real waveform reproduces the exact gradient magnitude. The roughly 7.9% relative L2 difference shows a scale/amplitude change. It means the **structural direction** is almost unchanged.

That is important for the biological bridge because a membrane-voltage-like return need not carry an explicit analytic-signal quadrature to preserve the useful direction in this model.

## R2 — envelope + task sign narrowly fails

The envelope-only code was

```text
sign(task coefficient) * |g(t)|
```

with carrier phase discarded.

Held-out:

```text
mean correlation          0.79524471
median correlation        0.85361367
bodies corr > .75         9 / 12
```

The preregistered mean threshold was

```text
mean correlation > .80
```

so this is a **FAIL**.

It is a narrow numerical miss, but it remains a miss. The correct interpretation is that the temporal amplitude envelope plus one consequence sign is often informative but is **not robust enough across bodies for the registered claim**.

The body dependence is visible directly: some held-out bodies were around `.95`, while seeds 498 and 508 fell to about `.53` and `.36`.

## R3 — fast 50%-duty gating is almost transparent

The exact return waveform was retained only during alternating on/off windows, then dose-renormalized. For period `P=14`, every possible gate phase offset was evaluated.

Held-out means across bodies:

```text
best-offset correlation        0.99961432
median-offset correlation      0.99878518
worst-offset correlation       0.99762266
median-offset strong-sign      0.98604712
```

Registered requirements:

```text
median corr > .99
worst corr  > .98
median sign > .96
```

**PASS.**

This is stronger than a favorable-phase effect. In this model, the fast gate can remove half of the return time support and still preserve almost the same structural direction for **every tested phase offset**.

Simulation frame units are not biologically calibrated. `P=14` is not being called a gamma frequency.

## R4 — slower gating exposes phase sensitivity

At `P=42`:

```text
best-offset correlation        0.96891020
median-offset correlation      0.83826683
worst-offset correlation       0.52554174
median-offset strong-sign      0.84073543
```

Fast-versus-slow registered differences:

```text
median advantage P14-P42       +0.16051835
worst advantage P14-P42        +0.47208091
```

Requirements were `> .03` and `> .15`.

**PASS.**

So the robust variable is not merely duty cycle. Temporal scale relative to the consequential waveform matters strongly.

## R5 — 32 sparse phase-bearing events retain substantial direction

Keep only 32 separated high-amplitude return times, make their amplitudes equal, and retain each selected sample's complex phase.

```text
mean correlation               0.86420230
median                          0.86884054
bodies corr > .80              8 / 12
mean strong-sign agreement     0.82998614
```

**PASS.**

This is useful but nowhere near exact. Sparse eventization is a substantial information loss and is strongly body-dependent.

## R6 — consequence sign is load-bearing

The sparse-event comparison is especially clean.

At `N=32`:

```text
fixed-amplitude events retaining target/distractor sign
    mean corr                   +0.85007683

same event times, all events forced positive
    mean corr                   -0.08718019

difference                     +0.93725703
```

**PASS.**

The same pattern appears all the way down the sparse curve. Event timing by itself does not carry the structural objective. Some signed consequence / valence variable is essential in this construction.

Interestingly, for 1–16 events the phase-bearing and sign-only sparse codes are nearly identical:

```text
N      phase-bearing      sign-only
1         .7653             .7648
2         .7703             .7699
4         .7700             .7697
8         .7969             .7967
16        .8181             .8178
32        .8642             .8501
```

So at sparse peaks, detailed carrier phase contributes surprisingly little until the code becomes denser. The global consequence sign contributes enormously.

## R7 — envelope beats phase-only coding

Held-out:

```text
phase-only mean corr            0.58276272
envelope + sign mean corr       0.79524471
difference                      +0.21248199
```

Registered difference: `> .08`.

**PASS.**

Even though the envelope arm itself missed R2, the relative ordering replicates strongly: flattening the temporal amplitude structure is more damaging here than discarding the carrier phase while preserving the envelope and task sign.

## Descriptive delay result

The primary sparse code was not preregistered for a delay claim, but the curve is informative:

```text
sparse phase N=8, no added delay    0.7969
+1 frame                            0.7942
+2                                  0.7907
+4                                  0.7813
+8                                  0.7501
+16                                 0.5957
```

Small timing errors are tolerated gradually rather than causing an immediate collapse. Larger delays eventually destroy substantial structural direction.

This should be preregistered separately if pursued.

## What survives the biological thought experiment

The held-out test does **not** support the strong statement

```text
one stereotyped back-propagating spike is an adjoint signal
```

A very sparse return is only a moderate approximation, and performance varies greatly across bodies.

What survives is subtler:

```text
exact analog return waveform
        |
        | explicit complex quadrature can be removed
        v
real temporal return waveform       ~ same direction
        |
        | fast periodic time support can be removed
        v
50%-duty fast-gated waveform         ~ same direction
        |
        | severe event sparsification
        v
moderate, body-dependent direction
        |
        | remove consequence sign
        v
collapse
```

The current biological candidate is therefore **not phase inversion** and not a literal backpropagation algorithm.

It is a more general architecture:

```text
forward dendritic history
        |
        v
soma/AIS consequence
        |
        v
returning dendritic event / waveform
        |
        +---- fast temporal gating / multiplexing
        |
        +---- signed consequence / valence
        v
local overlap with recent forward state
        v
structural change
```

Passive cable reciprocity makes this worth testing inside one dendritic tree, but real active dendrites remain a major missing piece.

## The strongest new question

The `P=14` result is now too robust to leave as an observation.

Why can a 50%-duty gate discard half the exact return waveform at **any phase** while preserving the gradient map almost perfectly, whereas slower gating becomes strongly phase-sensitive?

The next experiment should derive the gated-gradient error in frequency space and ask whether the gate shifts the discarded/modulated component into temporal frequencies to which the forward/return structural sensitivity has little overlap.

That would connect this result directly to the earlier lock-in spectral-collision result without invoking biological rhythm names.

## Formal verdict

```text
R0 exact positive control                PASS
R1 real-valued return                    PASS
R2 envelope + task sign                  FAIL
R3 fast 50%-duty gate                    PASS
R4 fast/slow scale separation            PASS
R5 sparse phase N=32                     PASS
R6 consequence sign load-bearing         PASS
R7 envelope > phase-only                 PASS

TOTAL                                    7 / 8
```

## Wall sentence

> **The exact analog return contains more information than the structural update needs. In this reciprocal-arbor model the update direction survives loss of complex quadrature and even half of the return time support when that support is removed rapidly, but it does not survive loss of consequence sign.**
