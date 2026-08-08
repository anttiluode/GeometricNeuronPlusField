# AIS n-kinetics / interspike-interval test v0.1 — preregistration

## Why this test exists

The earlier kinetics sweep changed `n` speed from `0.5x` to `2x` and produced a
large event-count change without the preregistered clean passband translation.
Claude's audit proposed a sharper discriminator:

> if `n` kinetics mainly set a refractory floor, minimum interspike interval
> should move with the kinetic timescale; if event counts change while minimum
> ISI stays roughly fixed, the dominant effect is excitability/threshold.

The existing JSON did not retain spike times, so this test reruns the exact
frozen kinetics battery and records ISIs.

## Frozen pieces

Identical to AIS_KINETICS_V01:

- same 24 body seeds;
- same power interface `|psi_soma|^2`;
- same body-level normalization battery, including the two task traces;
- same HH parameters and input gain;
- same frequency battery and 640-frame traces;
- same burn-in of 160 frames;
- same `n_scale = 0.5, 1, 2`;
- `h_scale = 1` throughout.

No conductance or gain retuning.

## Primary frequency

`f = 0.025 cycles/frame`.

This was the dominant/exposed regime in the original kinetics receipt and is
chosen before examining ISI results.

A body is valid for the primary slow-vs-fast test only if both `n_slow` and
`n_fast` emit at least two post-burn spikes at this frequency.

## Registered refractory prediction

`n_scale` multiplies the gate differential speed, so `n_fast` has the shorter
effective kinetic time constant.

If `n` kinetics set a refractory floor:

```text
minISI(n_slow) > minISI(full) > minISI(n_fast)
log2[minISI(n_fast) / minISI(n_slow)] < 0
```

Primary statistic: paired one-sided Wilcoxon on the per-body
`log2(fast/slow minimum ISI)` values, alternative `< 0`.

Also report the exact sign test and how many valid bodies obey the strict
slow > full > fast ordering.

## Interpretation

- predicted directional minimum-ISI shift -> evidence for refractory-floor
  control;
- strong event-count changes with little/no minimum-ISI shift -> evidence that
  the dominant effect is excitability/effective threshold rather than a simple
  refractory clock;
- too few valid slow/fast bodies -> exposure failure, no mechanism claim.

Secondary receipts include median/p10 ISI and pooled within-condition ISIs, but
they do not replace the registered primary test.
