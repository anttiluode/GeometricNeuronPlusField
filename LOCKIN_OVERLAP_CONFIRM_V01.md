# Rhythmic lock-in gradient readout v0.1 — held-out confirmation

## Result

A balanced local `+/-` phase reference can extract the compressed transient forward/adjoint cross-term from intensity alone, and the exact failure boundary is spectral overlap with the retained waveforms' self-energy difference spectrum.

This was discovered on bodies 472–479, frozen in [`LOCKIN_OVERLAP_CONFIRM_PREREG_V01.md`](LOCKIN_OVERLAP_CONFIRM_PREREG_V01.md), and then tested on fresh bodies **480–491**.

All **6 / 6** registered criteria passed.

## The local measurement

For retained-band local waveforms `u(t)` and `v(t)`, modulate the returned field with a balanced binary reference

```text
s(t) in {+1,-1}
```

and measure only local power

```text
I(t) = |u(t) + s(t) v(t)|^2.
```

A signed accumulator gives

```text
1/2 sum_t s(t) I(t)
 = sum_t Re[conj(u(t)) v(t)]
 + 1/2 sum_t s(t) [|u(t)|^2 + |v(t)|^2].
```

The first term is the desired compressed local gradient. The second is the entire error term.

So the problem becomes a lock-in / heterodyne separation problem rather than a complex multiplication problem.

## Why the modulation rate matters

If `u` and `v` contain only selected Fourier bins `{k_i}`, their self-energy

```text
|u|^2 + |v|^2
```

can contain only pairwise difference frequencies

```text
k_i - k_j  (mod T).
```

The binary reference has its own discrete harmonic spectrum.

Therefore the signed intensity accumulator is exact whenever the reference harmonics are orthogonal to the selected-band difference spectrum.

The development run looked superficially like a simple "fast modulation wins" result. The held-out run shows that **spectral collision**, not speed by itself, is the correct variable.

A particularly useful example is held-out K=8 body seed 488:

```text
half-period 5 frames   2 colliding harmonics   leakage L2 = 0.0429
half-period 7 frames   0 colliding harmonics   leakage ~ machine zero
```

So a slower toggle can be exact while a faster one fails if the slower reference happens to avoid the self-energy difference frequencies.

## Held-out registered results

### C0 — frequency-domain leakage identity

The leakage predicted from the frequency-domain overlap reproduced the directly accumulated time-domain leakage with

```text
max normalized error = 9.71e-15
```

**PASS.**

### C1 — collision-free modulation is exact

Across the fresh bodies there were **80** `(body, K, rate)` points with zero support collision.

```text
max leakage / gradient L2     5.80e-15
median                         1.35e-15
```

Registered maximum: `1e-10`.

**PASS.**

### C2 — collisions mark the failure regime

There were **88** collision-positive points.

```text
median leakage / gradient L2        1.1989
median collision/free leakage ratio 8.86e14
```

**PASS.**

A support collision is not sufficient to specify the final error magnitude because local amplitudes and phases can still cancel. It is the point where leakage becomes spectrally allowed.

### C3 — weighted spectral overlap predicts leakage magnitude

Across every fresh body, K value, and rate:

```text
corr(weighted spectral overlap, measured leakage) = 0.96953
```

Registered requirement: `r > .90`.

**PASS.**

### C4 — a broadband two-state measurement is enough

The polarization identity was first written per retained bin,

```text
(|U_k+V_k|^2 - |U_k-V_k|^2)/4
    = Re[conj(U_k)V_k].
```

That suggests `2K` intensity measurements.

But that count is unnecessarily pessimistic if the retained K components can be replayed together as one band-limited waveform. Parseval gives

```text
1/4 [sum_t |u(t)+v(t)|^2 - sum_t |u(t)-v(t)|^2]
 = sum_k Re[conj(U_k)V_k]
```

for the retained band.

So the complete K-bin cross-term can be measured with **two global phase states**, independent of K, plus a full-window local power integral.

On fresh bodies:

```text
max relative L2 error = 3.34e-15
```

**PASS.**

This does not mean a real device has zero K-dependent cost. It still has to retain/form the selected band, align the two transient waveforms, integrate the measurement window, and implement the common phase state. It means the algebra does not require `2K` separately phase-stepped measurements.

### C5 — modest global phase error is benign in this model

With the nominal `pi` state offset by a fixed `0.10 rad`:

```text
K=8
corr                   0.9999970
relative L2            0.003399
strong-sign agreement  1.000

K=16
corr                   0.9999992
relative L2            0.002968
strong-sign agreement  1.000
```

**PASS.**

This is a model robustness result, not a calibrated tolerance for a particular optical, acoustic, RF, mechanical, or biological device.

## The retained-band approximation itself

The physical readout above reconstructs whichever compressed K-bin gradient it is given. Compression error is separate.

Fresh bodies 480–491:

```text
K=8   mean corr to exact gradient   0.99134
      mean relative L2              0.11783

K=16  mean corr to exact gradient   0.99853
      mean relative L2              0.05210
```

So the hierarchy is now clean:

```text
full transient adjoint
        |
        | spectral compression error
        v
K-bin local gradient
        |
        | intensity-readout error
        v
physical local estimate
```

The second arrow can be essentially exact under the registered measurement protocol.

## Relation to prior art

This is **not** a novelty claim for interference-based physical gradient measurement. In-situ adjoint/backpropagation in wave systems and intensity-based forward/adjoint interference measurements are established prior art, including Hughes et al. (2018) and Pai et al. (2023).

The narrower object studied here is the **finite-time, spectrally compressed broadband local readout** in this transient scattering-mesh testbed, including the explicit spectral-collision rule for a one-run binary lock-in reference.

Novelty relative to the broader time-domain / lock-in / adjoint literature has not been established.

## Biological boundary

This result does **not** establish that gamma oscillations, basket cells, axo-axonic/chandelier cells, the axon initial segment, or theta rhythms implement this protocol.

In particular, inhibitory GABAergic input should not be described as a generic `pi` phase inverter.

What the result does provide is a precise engineering question that can be compared with biology without forcing the answer:

> Can a fast rhythmic gate make a local relation between two temporally structured signals consequential by spectrally rejecting their much larger self-energy terms?

That is now falsifiable in the model. The biological analogy remains a hypothesis.

## Wall sentence

> **A local power detector does not need to multiply the forward and adjoint phasors. A balanced rhythmic reference can demodulate their cross-term directly; it is exact when the reference avoids the retained field's self-energy difference spectrum.**
