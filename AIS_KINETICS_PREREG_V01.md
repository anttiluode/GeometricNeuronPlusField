# AIS h/n kinetics v0.1 — do the slow gate time constants set the event-selection band?

The instantaneous-gate ablation established an asymmetry:

- making fast Na activation `m` instantaneous preserves most of the normalized frequency-allocation shape;
- making either Na inactivation `h` or K activation/recovery `n` instantaneous collapses spiking at the unchanged operating point.

That collapse proves those histories matter for maintaining excitability, but it is too strong to say which history shapes the frequency response.

This test perturbs **kinetic speed only**.

## Frozen quantities

Unchanged from AIS active v0.1/v0.2:

- all FunctionalArbor bodies and wave physics;
- soma power input;
- input normalization and gain;
- HH conductances, reversal potentials and capacitance;
- gate steady-state voltage curves;
- frequency battery and burn window.

No gain/conductance is retuned for a kinetic variant.

## Variants

The differential equation for one gate is multiplied by a scalar while leaving `x_inf(V)` unchanged:

```text
dh/dt -> s_h * dh/dt
dn/dt -> s_n * dn/dt
```

Registered variants:

```text
full       h x1,   n x1
h_slow     h x0.5, n x1
h_fast     h x2,   n x1
n_slow     h x1,   n x0.5
n_fast     h x1,   n x2
```

Thus only the memory timescale changes.

## Primary measurement

For each body and variant, form the natural event allocation across modulation frequencies

```text
p_f = N_f / sum_f N_f
```

when at least one event is emitted.

Define the allocation center on a log-frequency axis

```text
mu = sum_f p_f * log2(f)
```

and convert it back to a geometric-center frequency `2^mu` for readability.

The primary mechanistic signature is an **opposite signed paired shift** around the full model:

```text
slow gate -> allocation center shifts one way
fast gate -> allocation center shifts the other way
```

for the same gate across organisms.

Also report total-variation distance from full, mid/high band share (`0.05 <= f < 0.125`) and total event count.

## Interpretation

A gate's kinetics earn a frequency-selection role if:

1. both speed variants remain exposed on most bodies;
2. slow and fast variants move the allocation center in opposite directions on average;
3. the paired ordering is organism-stable rather than a few outliers;
4. normalized allocation shape changes, not only total spike count.

No directional prediction (higher/lower) is registered in advance; only the bidirectional dependence on kinetic speed is.

## Failure conditions

- If either speed perturbation kills spiking in most bodies, report an exposure wall.
- If slow and fast changes affect only total count while normalized frequency allocation stays nearly fixed, the gate controls excitability rather than temporal selection.
- If slow and fast shifts are not opposite or are highly heterogeneous, do not claim a gate-defined passband.
- Do not tune gains after seeing the result.
- This does not revisit the failed v0.2 phase-precision claim.

## Decision

A robust kinetics-dependent allocation shift would finally give the AIS bridge a concrete active temporal coordinate: a gate timescale that selects which upstream field dynamics become events.  Only then is it sensible to ask whether body geometry and AIS geometry/kinetics can compensate one another.
