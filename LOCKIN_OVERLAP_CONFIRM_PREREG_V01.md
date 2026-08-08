# Lock-in overlap confirmation preregistration v0.1

## Status

Frozen **before** running held-out bodies 480–491.

The development experiment on bodies 472–479 found that a single continuously modulated `+/-` phase reference can recover the compressed transient adjoint gradient from a signed local power accumulator when the reference spectrum does not overlap the self-energy spectrum of the retained K-bin waveforms.

This confirmation tests the mechanism, not a biological claim.

## Algebra frozen before confirmation

For selected-bin local waveforms `u(t)` and `v(t)`, let a balanced reference `s(t) in {+1,-1}` modulate the returned field:

```text
I(t) = |u(t) + s(t) v(t)|^2
```

Then

```text
1/2 sum_t s(t) I(t)
 = sum_t Re[conj(u(t)) v(t)]
 + 1/2 sum_t s(t) [|u(t)|^2 + |v(t)|^2].
```

The first term is the desired compressed local gradient. The second is the only leakage term.

If `u` and `v` contain retained Fourier bins `{k_i}`, the self-energy spectrum can only contain pairwise difference frequencies

```text
k_i - k_j  (mod T).
```

Therefore a sufficient exactness condition is spectral orthogonality between the modulation reference and this difference-frequency set.

## Frozen setup

Use the same implementation and task as `polarization_identity_probe.py` and `lockin_overlap_probe.py`:

- fresh bootstrap bodies: seeds **480–491**;
- lag `20`;
- `T=210` frames;
- boundary-selected `K=8` and `K=16` bins;
- square-wave half-periods `1, 3, 5, 7, 15, 21, 35` frames;
- no change to FunctionalArbor wave physics;
- no learning or structural updates in this mechanism test.

## Registered criteria

### C0 — frequency-domain leakage identity

For every held-out body, K, and tested modulation rate, the frequency-domain Parseval reconstruction of the self-energy leakage must match the directly accumulated time-domain leakage with

```text
normalized L2 error < 1e-10.
```

This is an implementation/algebra positive control.

### C1 — collision-free modulation is exact

Pool all held-out `(body, K, half-period)` points whose discrete reference harmonics have **zero support collision** with the selected-bin difference set.

Require:

```text
at least 24 collision-free points
max leakage / compressed-gradient L2 < 1e-10
```

The count floor prevents passing on a trivial tiny subset.

### C2 — collisions are the failure boundary

Pool all points with at least one support collision.

Require:

```text
at least 24 collision-positive points
median leakage / compressed-gradient L2 > 0.01
median collision-positive leakage
    > 1e6 * median collision-free leakage
```

This does not require every collision to be harmful; amplitudes and phases can cancel. It tests whether support collision marks the regime where measurable leakage becomes possible.

### C3 — weighted spectral overlap predicts leakage magnitude

Across all held-out body/K/rate points, correlate

```text
actual weighted reference x self-energy spectral overlap
```

with measured leakage magnitude.

Require:

```text
Pearson r > 0.90.
```

The development value was roughly `.98-.99`; `.90` leaves substantial room for held-out variation.

### C4 — two-state broadband polarization identity remains exact

The separate two-state broadband measurement from `polarization_identity_probe.py`

```text
1/4 [sum |u+v|^2 - sum |u-v|^2]
```

must reproduce the direct K-bin complex-product map with

```text
max relative L2 < 1e-10
```

for both K=8 and K=16.

This establishes that `2K` independently phase-stepped bin measurements are not algebraically required once the retained band-limited waveforms can be formed/replayed as a common broadband signal.

### C5 — modest phase-setting error remains a small perturbation

For a global nominal-pi phase state with a fixed `0.10 rad` offset, require for both K values:

```text
mean map correlation > 0.999
mean relative L2 < 0.01
mean strong-sign agreement > 0.99.
```

This is a robustness check, not a device-calibrated tolerance claim.

## Kill conditions

The proposed mechanism is rejected or narrowed if:

- collision-free points show appreciable leakage;
- the spectral leakage identity fails;
- support collisions do not separate the low- and high-leakage regimes;
- actual spectral overlap does not predict leakage magnitude;
- or the broadband two-state identity fails on fresh bodies.

## Claim boundary

A pass would support only this engineering statement:

> **In the compressed transient-adjoint model, a balanced rhythmic phase reference can act as a local lock-in readout of the forward/adjoint cross-term. Exactness is controlled by spectral separation between the reference and the retained waveforms' self-energy difference spectrum.**

It would **not** establish that hippocampal theta/gamma rhythms, basket cells, chandelier cells, the AIS, or biological neurons implement adjoint backpropagation.