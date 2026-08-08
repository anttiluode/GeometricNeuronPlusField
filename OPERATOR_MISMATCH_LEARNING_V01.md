# Reciprocal pass-mismatch learning v0.1 — the compressed learner tolerates large fixed operator drift before failing

The exact same-medium adjoint identity relies on a strong symmetry: the forward and retro passes use the same reciprocal linear operator.

`DEVICE_ERROR_LEARNING_V01.md` showed that modest phase/amplitude/tap readout errors are not especially damaging in this toy. The more fundamental hardware question is whether learning survives when the **medium itself is different on the return pass**.

`OPERATOR_MISMATCH_LEARNING_PREREG_V01.md` froze the fresh-body test before seeds 448–459 ran.

## Mismatch model

The forward pass uses the current nominal edge conductances `k_e`.

For each body/arm, the retro pass gets one fixed random fractional edge-error pattern:

```text
k_e,retro = k_e,forward * (1 + delta_e)
delta_e ~ Normal(0, sigma)
```

with conductances clipped only to remain positive.

The mismatch pattern stays fixed for all 40 learning iterations and is shared by target/distractor retro passes. As the learner changes a nominal coupling, both passes track that change, but the return pass retains its multiplicative edge error.

Important: **each pass is still internally reciprocal**. This is a pass-to-pass operator mismatch test, not a true nonsymmetric/nonreciprocal-network test.

Arms:

```text
exact       exact discrete adjoint
K8_ideal    K=8 boundary-selected compressed physical gradient, sigma=0
K8_m20      K=8, sigma=.20
K8_m30      K=8, sigma=.30
K16_m30     K=16, sigma=.30
K8_m50      K=8, sigma=.50 stress arm
```

All use the same deterministic frontier candidates, `eta=.01`, 40 relinearized iterations, lag20 and the 210-frame task.

## Static development wall

Before the held-out learning run, reused bodies 400–403 were used only to locate the mismatch scale.

The full physical gradient map was almost unchanged through surprisingly large random edge mismatch:

```text
sigma       full-map corr     full rel-L2     K8 corr      K16 corr
0            1.00000           0.00000        .99597       .99709
.05           .99992           .01453         .99594       .99702
.10           .99972           .02500         .99565       .99677
.20           .99695           .07195         .99282       .99394
.30           .98722           .15222         .98272       .98391
.50           .82829           .51951         .82518       .82573
.75           .56324          1.16783          .55388       .55512
```

So there is no knife-edge “perfect reciprocity or nothing” behavior in this particular task. Exact equality gives the mathematical identity, but useful gradient direction survives substantial approximate equality.

## Held-out confirmation — fresh seeds 448–459

### C0 — exact learner gate

Registered:

```text
mean exact DeltaC > .015
>= 10 / 12 bodies improve
```

Observed:

```text
mean exact DeltaC               +.038369
improved                          12 / 12
```

**C0 PASS.**

### C1 — K8 survives 20% RMS edge mismatch

Registered:

```text
mean K8_m20 DeltaC > .015
>= 10 / 12 bodies improve
```

Observed:

```text
mean K8_m20 DeltaC              +.034901
improved                          12 / 12
```

**C1 PASS.**

### C2 — K8 at 20% preserves most exact learning gain

Registered:

```text
mean(K8_m20 DeltaC) / mean(exact DeltaC) >= .75
mean(K8_m20 - exact) > -.010
```

Observed:

```text
group gain ratio                  .90964
mean K8_m20 - exact             -.003467
```

**C2 PASS.**

So under this model-space perturbation the K=8 learner retains about **91.0%** of the exact group-mean gain.

### C3 — K16 remains useful at 30% mismatch

Registered:

```text
mean K16_m30 DeltaC > .015
group gain ratio >= .70
```

Observed:

```text
mean K16_m30 DeltaC             +.033609
improved                          12 / 12
group gain ratio                  .87596
```

**C3 PASS.**

### C4 — map direction remains recognizable along the changing trajectory

Registered:

```text
K8_m20  mean map corr > .970
K16_m30 mean map corr > .975
```

Observed:

```text
K8_m20  mean map corr             .983845
         mean relative L2          .157936

K16_m30 mean map corr             .987073
         mean relative L2          .167057
```

**C4 PASS.**

## Formal verdict

```text
C0 PASS   exact learner is useful on fresh bodies
C1 PASS   K8 learns under sigma=.20 pass mismatch
C2 PASS   K8_m20 retains 91.0% of exact mean gain
C3 PASS   K16_m30 retains 87.6% of exact mean gain
C4 PASS   gradient-map direction remains high through training
```

**5 / 5 registered criteria pass.**

## Descriptive stress arms

The unregistered K8_m30 arm also remained useful:

```text
mean DeltaC                       +.032180
improved                           12 / 12
group gain ratio                    .83872
mean map corr                       .97588
mean relative L2                    .21110
```

So 30% RMS edge mismatch did not kill K=8 learning either, although monotone-step fraction had fallen to about `.533`.

At 50% RMS edge mismatch the picture changes sharply:

```text
K8_m50 mean DeltaC                +.007946
improved                            9 / 12
group gain ratio                    .20711
mean map corr                       .76774
mean relative L2                    .73130
monotone-step fraction              .48958
```

Three bodies had negative final gain in this stress arm.

This gives a useful qualitative wall in the toy: deterioration begins well before final learning disappears, and around 50% random pass mismatch the compressed learner is no longer a reliable approximation to the exact one.

## What the result means — and does not mean

The registered result is surprisingly robust, but `20%` and `30%` here are **not fabrication tolerances**.

They mean RMS fractional error on every abstract edge conductance in this model. A real platform's heater phase error, optical attenuation, S-parameter drift, mechanical detuning, transmission-line reactance error, or nonreciprocal element cannot be relabeled as this percentage without an explicit device model.

The likely reasons for the broad tolerance are ordinary rather than magical:

- the task gradient is redundant across many edges;
- the K-bin representation already discards some detail while preserving dominant direction;
- random edge perturbations partly average through distributed propagation;
- max-normalized finite gradient steps care more about direction/order than exact magnitude;
- the learner relinearizes after every step instead of integrating one stale wrong gradient.

Those are hypotheses, not yet separated experimentally.

## The reciprocity statement becomes more precise

The mathematical identity still requires the transpose/adjoint operator.

For a reciprocal medium,

```text
H^T = H
```

so the same physical geometry supplies it automatically.

This experiment shows that **approximate** same-operator replay can remain useful even when that equality is imperfect.

The next decisive test is therefore not more symmetric drift. It is true nonreciprocity:

```text
H != H^T
```

Construct a stable directed wave operator and compare:

```text
1. exact algorithmic adjoint using H^T
2. physical return through the same directed H
3. physical return through an explicitly transposed/reversed device H^T
```

If (2) fails while (3) restores the exact gradient, the slogan

> reciprocity buys cheap physical backpropagation

will have been turned into a direct causal experiment rather than an inference from the reciprocal case.

## Current hardware wall

The current line is now:

```text
broadband transient task
      -> reciprocal physical adjoint
      -> K=8/16 port-selected spectral correlation channels
      -> coherent-readout error tolerance
      -> substantial pass-to-pass reciprocal-operator drift tolerance
      -> local graded update
      -> successful closed-loop training
```

The next two engineering questions are sharply separated:

1. **operator directionality:** what happens when the forward medium is genuinely nonreciprocal?
2. **local measurement:** how are the K complex phasor products realized as actual intensity/phase-stepped measurements, and what is the pass/sensor/memory count?

## Wall sentence

> **Exact same-medium backpropagation is a reciprocity identity, but useful training in this transient mesh is not knife-edge: the compressed learner survives large fixed return-path mismatch and only degrades strongly once the returned operator has drifted very far from the forward one.**
