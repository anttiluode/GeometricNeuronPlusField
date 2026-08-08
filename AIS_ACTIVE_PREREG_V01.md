# AIS active boundary v0.1 — preregistration

This experiment freezes the entire FunctionalArbor body and adds one small active stateful compartment **downstream of the soma trace**.

It is motivated by the operator -> flow -> event picture, but it is designed to be easy to kill.

## Why now

The current body is almost exactly a graph-defined linear resonator bank.  The field can carry temporal structure, but the model presently ends in a passive soma readout.

Claude's `SomaWhyClaude` result also sharpens what the soma means in this construction: the convergence root is close to the amplitude-balance point of the two inputs.  Therefore this experiment does **not** treat the soma as a mysterious privileged biological location.  It uses the soma because it is the existing local convergence signal in the frozen toy.

The question is narrower:

> Once the global analog computation has converged to that local signal, does a genuinely stateful active boundary provide temporal/frequency coding that a rate-matched instantaneous threshold and the active membrane's own linear small-signal filter do not?

## Frozen upstream system

For every seed:

1. bootstrap the v0.9 FunctionalArbor;
2. freeze its morphology;
3. record ordinary soma power traces;
4. do **no growth, credit assignment, or structural adaptation**.

Two drive batteries are used:

- the existing A->B versus B->A temporal-order task at lag 20;
- a continuous single-source carrier with sinusoidal envelope modulation over a fixed frequency sweep.

Envelope frequencies are reported in cycles per field frame.  The downstream ODE maps one field frame to 1 ms only as a simulation scale; this is not a biological calibration.

## Three downstream encoders

### 1. Memoryless control

The normalized soma drive itself is thresholded.  One threshold is selected per body to match the **total event count** of the active boundary over the whole battery.  It is not retuned per frequency or task order.

### 2. Linearized-membrane control

The active membrane is perturbed below threshold at rest and its finite-difference small-signal impulse response is measured.  That measured kernel is convolved with the same soma drive and then thresholded at one total-rate-matched threshold.

This is the important control suggested by Claude.  No hand-written derivative or high-pass filter is allowed.

### 3. Active boundary

A compact Hodgkin-Huxley-like compartment with

- membrane voltage;
- Na activation `m`;
- Na inactivation `h`;
- K activation/recovery `n`;

receives the normalized soma drive as injected current.  Spike events are upward 0-mV crossings.

No derivative term, explicit high-pass term, task label, graph mode, or source identity is available to this compartment.

## Registered primary metric

For every modulation frequency, event phase locking is measured by vector strength

```text
VS = |mean exp(i 2 pi f t_k)|
```

The **primary score** is mean vector strength over the two highest registered modulation frequencies.

The active boundary earns a specifically stateful/high-bandwidth role only if it outperforms **both** controls, especially its own linearization, across organisms rather than merely in the grand mean.

## Secondary task receipt

For A->B and B->A we record, per encoder:

- event-count contrast;
- first-event latency difference;
- event-time centroid difference.

This is secondary because the upstream soma trace already contains task information.  The active boundary is not credited merely for preserving a distinction that a threshold can already read.

## Failure conditions

The AIS bridge does **not** earn the proposed role if any of these occur:

1. the active gate only beats the memoryless threshold but not the membrane's own linearized filter;
2. any apparent advantage disappears when event rate is matched;
3. the active gate improves only low-frequency response while degrading the registered high-frequency score;
4. the result is carried by a few seeds rather than being organism-stable;
5. the active model only works after adding an explicit derivative/high-pass term;
6. numerical instability, zero-event bodies, or threshold saturation make the comparison ill-posed.

A clean null is useful: it would say that the current soma signal does not need an AIS-like nonlinear history state for this task/battery, even if real AIS biology has such state for other reasons.

## What a positive result would and would not mean

A positive result would support this computational statement:

> A local active state can add temporal event coding beyond the frozen body's LTI modal dynamics and beyond the same membrane's linear small-signal filter.

It would **not** show why the biological AIS evolved, nor that this toy reproduces biological AIS frequency scales.
