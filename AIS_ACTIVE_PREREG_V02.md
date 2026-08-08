# AIS active boundary v0.2 — timing precision under per-frequency rate matching

v0.1 exposed a real frequency redistribution but its registered high-frequency vector-strength headline was invalid as a precision claim because some active conditions had only one or two events.

The active membrane parameters are now **frozen exactly as v0.1**.  v0.2 changes the analysis and controls, not the mechanism.

## Frozen mechanism

No changes are allowed to:

- FunctionalArbor body or wave physics;
- soma readout used as boundary input;
- HH-like boundary conductances or gating equations;
- input normalization rule;
- input gain;
- modulation frequencies;
- task lag.

## Stronger frequency control

For each body and each frequency separately:

1. run the active boundary and count its events in the post-burn window;
2. choose a memoryless threshold on that same frequency trace that produces the same number of events;
3. choose a threshold on the active membrane's measured linearized response that produces the same number of events;
4. compare event timing.

This deliberately gives the controls an oracle advantage: they are allowed a separate threshold at every frequency.  They cannot lose merely because the active membrane allocates its global firing budget differently across frequencies.

## Exposure gate

A body/frequency pair is timing-valid only when the active boundary emits at least

```text
N >= 4 events
```

in the analysis window.

Sparse conditions remain reported but do not count as evidence for timing precision.

## Timing metrics

For every rate-matched valid pair report:

- vector strength `VS`;
- pairwise phase consistency `PPC`, the spike-count-unbiased phase-locking estimator;
- events per modulation cycle;
- a coverage-weighted phase score

```text
coverage_score = VS * min(1, N_events / N_cycles)
```

The controls have the same `N_events`, so within a condition the coverage factor is identical.  It exists to stop sparse conditions from dominating aggregate summaries.

## Registered upper-band test

The upper band is

```text
f >= 0.05 cycles / field frame
```

including `0.05`, `0.0833333`, and `0.125`.

Only body/frequency pairs passing `N >= 4` enter the precision comparison.  The v0.1 data already tell us the highest frequency is often underexposed; the gate is therefore part of the correction, not a hidden deletion.

Primary comparison:

```text
mean(active PPC - linearized PPC)
```

across all valid upper-band body/frequency pairs, accompanied by the sign count.

The active-state precision hypothesis earns support only if the active gate beats **its own linearization** on PPC and does so across a majority of valid pairs.  Memoryless results are also reported but are secondary to this harder control.

## Task control

For A->B versus B->A, the memoryless and linearized thresholds are now fitted over the **two task traces only** to match the active total task-event count.  This removes the v0.1 confound where a global frequency-battery threshold often left the controls with no task events.

Report count contrast, first-event latency difference, and event-centroid difference.  The task remains secondary because the soma trace already contains the distinction.

## Failure conditions

The active-boundary precision claim fails if:

1. active PPC does not exceed linearized PPC in the registered upper-band comparison;
2. an apparent VS advantage disappears under PPC;
3. valid-pair count is too small to support an organism-level statement;
4. the active advantage exists only against the memoryless control;
5. the model needs retuned conductances/gain after seeing v0.1;
6. rate matching is not exact or numerically stable.

A null would still leave the v0.1 frequency-allocation result intact: stateful nonlinear dynamics may choose *where* to fire without making individual event timing more precise than the same membrane's linearized filter.
