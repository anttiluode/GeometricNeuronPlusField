# AIS gating-memory ablation v0.1 — fast activation is dispensable; instantaneous h or n collapses spiking

Canonical run: `31243133327`

This experiment kept the v0.1/v0.2 body, normalization, input gain and HH conductances fixed and removed independent gating memory one gate at a time by replacing that gate with its instantaneous steady-state value `x_inf(V)`.

The primary intended metric was the normalized frequency-allocation shape, so a gate would only earn a frequency-selection role if its ablation changed *where* spikes were allocated rather than merely changing total excitability.

## Full boundary

Across 24 bodies the unmodified active boundary emitted on every body:

```text
mean total frequency-battery events   51.54
median total                           47
mean event share in 0.05 <= f < .125  0.194
```

Mean natural event allocation remained dominated by the `0.025` regime, with the smaller `0.05-0.0833` band seen in v0.1.

## Na activation memory (`m`) is not the main ingredient

Making `m` instantaneous did **not** destroy firing:

```text
emitting bodies             24/24
mean total events           79.92
mean TV distance from full  0.102
median TV                   0.085
peak frequency agrees       18/24
mid/high band share delta  -0.019
```

So the finite response time of the fast Na activation gate affects gain, but most of the normalized frequency-allocation shape survives without an independent `m(t)` history.

That is a useful negative result: the observed frequency selection is not simply the finite activation lag of the sodium conductance.

## The registered h/n ablations hit an exposure wall

Making either slower gate instantaneous collapsed the active compartment completely at the frozen operating point:

```text
h_instant      0 events in 24/24 bodies
n_instant      0 events in 24/24 bodies
all_instant    0 events in 24/24 bodies
```

Per the preregistration, this **cannot** be interpreted as “h is the frequency filter” or “n is the frequency filter.”  The shape is undefined when the model emits nothing.

What it does establish is narrower:

> **Independent Na-inactivation and K-activation/recovery dynamics are not decorative state variables in this boundary.  Replacing either by its instantaneous voltage equilibrium destroys the spiking regime reached by the unchanged soma drive, whereas removing fast-m activation memory does not.**

The reason is mechanistically plausible within HH equations: an instantaneous `h_inf(V)` can remove Na availability too aggressively during depolarization, while instantaneous `n_inf(V)` can recruit K current too aggressively.  But that explanation should be tested through kinetics rather than asserted from the collapse.

## Why this blocks geometry co-adaptation for now

We now have two separate facts:

1. v0.2: the full stateful gate does **not** improve rate-matched phase precision beyond its own linearization;
2. this ablation: slow h/n history is necessary to maintain the active spiking regime, but the all-or-none collapse prevents attribution of the frequency-allocation shape to either gate.

So moving or stretching an AIS-like segment now would be premature.  We still do not have a clean active variable for geometry to tune.

## Next clean test

Do not delete h or n.  Change only their **kinetic speed** while preserving their steady-state voltage curves and all conductances:

```text
h kinetics x 0.5, x 1, x 2
n kinetics x 0.5, x 1, x 2
```

No gain retuning.

If the center/shape of the output frequency allocation moves systematically when h or n kinetics are sped up or slowed down, then we have direct evidence that the active state timescale itself is setting a temporal pass region.

That is the quantity worth taking into a later body/AIS co-adaptation experiment.

## Wall sentence

> **Fast Na activation memory is largely unnecessary for the observed frequency-allocation shape, but slow Na inactivation and K recovery cannot be made instantaneous without destroying spiking.  The next question is therefore not whether h/n matter, but whether their time constants set the event-selection band.**
