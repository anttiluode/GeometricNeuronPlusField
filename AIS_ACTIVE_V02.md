# AIS active boundary v0.2 — the stateful gate changes frequency allocation but does **not** improve spike-time precision

Canonical GitHub Actions run: `31242925894`

The v0.2 mechanism was frozen from v0.1.  No conductance, gain, upstream wave parameter, normalization rule, modulation frequency, or task lag was retuned after seeing v0.1.

The only change was a stronger analysis/control:

- memoryless and linearized controls were **exactly event-count matched separately at each frequency**;
- the linearized control is the measured small-signal impulse response of the same active membrane, not a hand-written high-pass filter;
- timing was scored only when the active encoder emitted at least four events;
- pairwise phase consistency (PPC) was used alongside vector strength to remove the one-spike `VS=1` trap.

## Registered upper-band result

The registered upper band was `f >= 0.05` cycles / field frame.

There were 21 valid body/frequency pairs from 11 unique frozen bodies.

```text
mean PPC

active       +0.4305
linearized   +0.5816
memoryless   +0.5238

active - linearized   -0.1511
active - memoryless   -0.0933

active vs linearized:   6 wins / 15 losses   two-sided sign p = 0.0784
active vs memoryless:   7 wins / 14 losses   p = 0.1892
```

So the preregistered active-state precision hypothesis **fails**.

The active HH-like boundary does not sharpen event phase timing beyond its own linear small-signal filter.  If anything, over the valid upper-band pairs its spikes are less phase-consistent.

That is the result to keep.  Do not rescue it by tuning the conductances or moving the frequency window.

## Frequency-by-frequency picture

The null is not simply “the active model does nothing.”  It is strongly frequency dependent.

```text
f          valid pairs    PPC memoryless   PPC linearized   PPC active
0.00625        4              0.9383           0.8811          0.2095
0.01250       14              0.7324           0.7356          0.1250
0.02500       24              0.3986           0.4074          0.0403
0.05000       10              0.6777           0.7048          0.6099
0.08333       10              0.4367           0.5308          0.3072
0.12500        1             -0.1429          -0.1429         -0.1306
```

At `0.05` the active gate comes closest: it beats the linearized control in 5/10 valid bodies, but its mean PPC is still lower.  At `0.0833` it loses to the linearized control in all 10 valid bodies.  The highest frequency remains underexposed and is not interpretable as an organism-level result.

## What survives from v0.1

v0.1 used one total-rate-matched threshold over the entire frequency battery.  Under that fixed operating point, the active boundary allocated its finite event budget very differently from the controls:

- it lost badly around `0.025`;
- it produced a strong event band around `0.05-0.0833` where the globally thresholded controls largely stopped firing;
- it largely stopped again at `0.125`.

v0.2 shows that this is **frequency selection / firing allocation**, not superior timing precision.  Once the controls are allowed to emit the same number of events at each frequency, their event phases are at least as precise and usually more precise.

So the clean statement is:

> **The active channel state changes which temporal regimes become output events, but in this first boundary model it does not make those events more phase-precise than a rate-matched linear observation of the same soma signal.**

That is already a useful division of labor.  `h(t)` and `n(t)` can alter event availability without automatically adding information that was absent upstream.

## Temporal-order task

The task-pair controls were also corrected in v0.2: their event count was matched over A->B and B->A only, rather than over the whole frequency battery.

The result exposes another important distinction.

The memoryless and linearized top-k controls usually put all matched events into only one of the two task traces.  Their mean absolute count contrast is therefore high (`~0.792`), but they provide no paired first-event or centroid timing measurement.  The active boundary has a smaller mean absolute count contrast (`~0.266`) but emits events in both task orders often enough to give paired timing in 15/24 bodies (mean first-event separation ~10.28 frames; centroid separation ~12.35 frames).

This is **not evidence that active dynamics improve task classification** — the rate controls actually separate the two conditions more strongly by event count.  It says the active state changes the *code*: from an almost all-or-none count allocation toward two spike trains whose relative timing can be compared.

That coding difference is worth keeping separate from “better.”

## Consequence for the AIS bridge

The experiment was explicitly staged so that position/extent co-adaptation would happen only if the active state first earned a role.

The answer is mixed:

```text
active state as superior timing-precision stage       NO
active state as nonlinear frequency/event selector    YES, descriptively
```

Therefore **do not start body-AIS position/extent co-adaptation yet**.  That would add degrees of freedom immediately after the primary timing claim failed and risk tuning around a null.

The next legitimate question is mechanistic and parameter-free:

> Which internal state is responsible for the frequency redistribution already observed?

Freeze the same HH parameters and ablate state variables one at a time — Na inactivation `h`, K recovery `n`, and delayed activation — without retuning thresholds or gains.  If the mid/high-frequency event band disappears under a specific state ablation, then we can say what the active boundary contributes before asking geometry to co-adapt with it.

## Wall sentence

> **The first AIS-like active boundary does not improve rate-matched spike-time precision beyond its own linearization.  It does, however, nonlinearly redistribute which input frequencies are allowed to become events.  The next step is to identify which gating memory creates that selection, not to tune AIS geometry around the failed precision hypothesis.**
