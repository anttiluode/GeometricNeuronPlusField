# AIS gating-memory ablation v0.1 — historical receipt; h/n arms are now classified as dead-arm controls

> **Interpretation update (Claude audit):** the raw receipt below is retained, but the old wording that independent `h(t)` / `n(t)` history was shown to be dynamically necessary for *this frequency-selection mechanism* is superseded by [`AIS_CLAUDE_AUDIT_V01.md`](AIS_CLAUDE_AUDIT_V01.md).  An independent standard-HH drive sweep reproduced the key asymmetry and showed that making `h` or `n` instantaneous abolishes repetitive firing across the operating sweep.  Those are therefore **dead arms by construction**, not interpretable frequency-mechanism ablations.  The `m`-instantaneous arm remains informative because it stays excitable.  See also [`AIS_N_ISI_V01.md`](AIS_N_ISI_V01.md).

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

## The h/n ablations are exposure failures, not mechanism findings

Making either slower gate instantaneous collapsed the active compartment completely at the frozen operating point:

```text
h_instant      0 events in 24/24 bodies
n_instant      0 events in 24/24 bodies
all_instant    0 events in 24/24 bodies
```

The preregistration already prohibited assigning a frequency-selection shape to a silent arm.  Claude's later standard-HH check strengthens that caution: instantaneous `h` or `n` removes the normal HH timescale separation required for repetitive spiking over a broad drive sweep.

Therefore these rows establish only that **the ablated equations no longer occupy a comparable excitable regime**. They do *not* identify `h` or `n` memory as the source of the observed frequency allocation.

A valid ablation interpretation would first have to restore a matched firing regime without using the frequency profile as the tuning target, and only then compare allocation shapes.

## Why this blocks geometry co-adaptation

The active v0.2 timing test did not earn a precision role, and the h/n instantaneous arms do not isolate a usable active mechanism.  Moving or stretching an AIS-like segment around these results would therefore add degrees of freedom around null/dead-arm observations.

## Follow-up that was actually run

Rather than rescuing the dead arms by arbitrary conductance tuning, the next experiment changed only gate kinetic speed while preserving the steady-state curves and conductances:

```text
h kinetics x 0.5, x 1, x 2
n kinetics x 0.5, x 1, x 2
```

That registered passband-shift prediction also failed; `n` was strongly consequential for excitability but did not translate a clean frequency band.  See [`AIS_KINETICS_V01.md`](AIS_KINETICS_V01.md).

The subsequent preregistered minimum-ISI test also rejected a simple `n`-sets-refractory-period account. See [`AIS_N_ISI_V01.md`](AIS_N_ISI_V01.md).

## Wall sentence

> **The m-instantaneous arm is a useful negative control because it remains excitable. The h- and n-instantaneous arms are dead-arm controls: their silence cannot identify the normal frequency-selection role of those gates.**
