# Two-tone cross-frequency beat probe v0.1 — useful null and a group-delay clue

Date: 2026-08-08

## Question

Return to the theta/gamma intuition without using a finite burst whose slow envelope is already obvious in the input spectrum.

For each source location, drive the same two continuous harmonic components

```text
omega_1, omega_2
```

through the frozen arbor.

The linear source spectrum contains only those two frequencies. It contains no Fourier line at

```text
Delta = |omega_2 - omega_1|.
```

At the soma, however, a quadratic / power readout produces the cross-frequency difference phasor

```text
B_j = H_j(omega_1) * conj(H_j(omega_2)).
```

Its phase is

```text
arg B_j = phi_j(omega_1) - phi_j(omega_2),
```

so it measures the phase of the two-tone envelope / beat at the common readout.

The development question was deliberately strong:

> Does the quasi-active material learned only on the low direct frequencies `.03,.04` also make high-carrier difference-frequency / envelope phase more coherent across spatially distributed source locations?

If yes, that would provide one bridge from the confirmed low-band material result toward a gamma/theta-like cross-frequency architecture.

---

## Materials compared

On fresh development bodies 616–621:

```text
zero
uniform material
confirmed hand gradient
self-organized material learned on direct omega=.03,.04
shuffled learned material
```

The learned material was never optimized on the carrier pairs used in this test.

Carrier pairs:

```text
Delta=.03
.10/.13
.16/.19
.24/.27
.36/.39

Delta=.04
.12/.16
.20/.24
.32/.36
.48/.52
```

Frequencies are simulation angular frequencies. The upper pair is not assigned a biological Hz label.

---

# Primary result — the attractive cross-frequency transfer story fails

Across all eight carrier pairs:

```text
mean beat-phase coherence R^2

learned material       .14376
uniform material       .14857
shuffled learned       .14972
```

Therefore

```text
learned - uniform      -.00481
learned - shuffled     -.00597
```

The low-band material learner does **not** generally improve envelope / difference-phase coherence of faster carrier pairs.

Some individual low carrier pairs show small positive learned advantages, but others are negative and the pooled direction is the opposite of the desired story.

So stop the claim:

```text
"theta-optimized / HCN-like material automatically organizes gamma carriers into a coherent theta beat"
```

at this point.

---

## A second result is still interesting

Although the learned material does not improve the beat, the beat phase itself is often much more coherent across source locations than either absolute carrier phase.

Across the six bodies, mean learned-material carrier phase coherence was extremely small for most pairs:

```text
R^2 ~ .001 to .005
```

while the corresponding difference-phase / beat coherence could be an order of magnitude larger.

Examples:

```text
pair .16/.19
learned carrier R^2      .00275
learned beat R^2         .06564

pair .12/.16
learned carrier R^2      .00291
learned beat R^2         .04438

pair .32/.36
learned carrier R^2      .00234
learned beat R^2         .07594
```

This phenomenon also exists in passive/uniform controls, so it is not an HCN-learning result.

---

## Why absolute carrier phase can disappear while envelope phase survives

For nearby frequencies around carrier `omega_c`,

```text
phi(omega_1) - phi(omega_2)
  ~= -Delta * dphi/domega.
```

The absolute propagation phase can wind substantially with source location while the local phase **slope** / group delay varies much less.

Thus a distributed arbor can be poor at preserving absolute carrier phase yet substantially better at preserving nearby-frequency phase difference.

The quadratic readout exposes that phase difference as the beat / envelope phase.

This is standard group-delay algebra in a new place in this project, not a novel theorem.

The useful question for the Geometric Neuron becomes:

```text
which quantities are geometrically stable under projection?

absolute phase?
phase difference?
group delay?
envelope timing?
```

The present result says absolute phase and nearby-frequency phase difference can behave very differently.

---

## Dead-band control caught an attractive artifact

The `.48/.52` pair produced a spectacular apparent beat coherence:

```text
carrier R^2      ~.0054
beat R^2         ~.887
```

But its median beat amplitude was only about

```text
1e-19.
```

So this pair is classified as an effectively untransmitted / dead-band artifact.

Do not use it as evidence for a coherent high-frequency envelope mechanism.

This is the same lesson learned repeatedly in the AIS branch: a normalized timing/phase statistic can look perfect when essentially no physical signal survives.

---

## Direct low-frequency transfer is a different object

For the `.03` difference-frequency pairs, the self-organized material has very high direct `.03` phase coherence across source locations (`R^2 ~ .976` in this development set), but the fast-carrier beat at the same Delta is usually far less coherent.

For example:

```text
pair .10/.13
learned direct .03 R^2    .9763
learned beat R^2          .0336

pair .16/.19
learned direct .03 R^2    .9763
learned beat R^2          .0656
```

So direct low-frequency synchronization and high-carrier group-delay / envelope synchronization are **not interchangeable**.

This explains why the material trained on direct `.03/.04` does not automatically transfer its advantage to high-carrier beats.

---

## Relation to the biological gamma/theta result

Vaidya & Johnston's biological result should not be mapped onto this two-tone mixer one-for-one.

Their gamma-frequency synaptic bursts produce slow components in the synaptic-current waveform, and the HCN gradient affects transfer of those components toward the soma.

The present two-tone experiment instead asks whether a difference-frequency phase relation generated/exposed by a quadratic readout is spatially coherent.

Those are different mechanisms.

The null result is therefore useful: it prevents us from using the musician's beat-frequency analogy as a shortcut around the actual dendritic transfer mechanism.

---

## Next cross-frequency question, if pursued

The interesting surviving object is **group delay**, not gamma-to-theta magic.

A sharper next experiment would optimize material directly for envelope/group-delay coherence of an amplitude-qualified carrier band, then ask whether:

1. the learned map again organizes by distance from consequence;
2. it differs from the direct-low-frequency material map;
3. a single material distribution can trade off direct theta-phase synchrony and fast-envelope synchrony;
4. the result survives an explicit minimum-transfer-amplitude criterion.

Do not preregister this until the carrier band is chosen without using a dead-band maximum.

## Wall sentence

> **The arbor can preserve nearby-frequency phase difference far better than absolute carrier phase, but the material learned for direct low-frequency synchronization does not automatically improve that envelope coordinate. Direct theta-like transfer and fast-carrier beat timing are distinct problems in this model.**
