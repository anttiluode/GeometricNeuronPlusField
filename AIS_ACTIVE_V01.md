# AIS active boundary v0.1 — first run, and why the primary precision score is not yet valid

Canonical run: `31242632344`

The implementation worked end to end on 24 frozen bodies.  The upstream FunctionalArbor was untouched; its soma power fed three downstream encoders:

1. raw soma drive + one total-rate-matched threshold;
2. the active membrane's measured small-signal linear impulse response + one total-rate-matched threshold;
3. a stateful Hodgkin-Huxley-like compartment with `V,m,h,n` and spike threshold crossings.

No derivative or hand-written high-pass term was inserted.

## What the first run appears to say

With one threshold per control matched to the active encoder's **total event count across the whole battery**, event allocation across modulation frequencies differed strongly.

```text
frequency     memoryless VS   linearized VS   active VS
0.00625          0.0000          0.0000         0.0716
0.01250          0.2318          0.2015         0.3297
0.02500          0.5056          0.5108         0.2267
0.05000          0.1161          0.1898         0.5397
0.08333          0.0000          0.0000         0.4953
0.12500          0.0000          0.0000         0.1264
```

The active compartment is **not simply better everywhere**.  At `0.025` cycles/frame it is much worse than both controls, while around `0.05-0.0833` it puts events into a band that the globally thresholded controls mostly abandon.  This already looks more like a nonlinear frequency redistribution than a generic sharpening filter.

The total event counts are exactly matched body by body across all three encoders.

## But the preregistered high-frequency precision score is flawed

The v0.1 preregistration used vector strength as the primary high-frequency score.  In a sparse regime this is not a reliable precision estimator:

```text
one event -> vector strength = 1.0
```

That happened in several bodies.  At `0.125` cycles/frame only 4/24 bodies emitted any active events at all, and only one body emitted at least four.  At `0.0833`, 17/24 emitted at least one event but only 10/24 emitted at least four.

Therefore the v0.1 headline

```text
mean high-frequency VS:
  memoryless  0.000
  linearized  0.000
  active      0.311
```

**must not be interpreted as demonstrated spike-time precision.**  It mixes two effects:

1. frequency-dependent event allocation under one global rate match;
2. phase concentration estimated from sometimes one or two events.

This is the same kind of estimator problem that the earlier lag work exposed with argmax-of-cross-covariance: the interesting-looking number is not yet the right statistic.

## Useful things v0.1 did establish

- The active compartment is numerically stable and produces events on real frozen-arbor soma traces.
- The two controls can be exactly total-rate matched body by body.
- The active system has a strongly non-flat frequency response and is not a monotonic high-pass: it loses badly at `0.025`, gains around `0.05-0.0833`, and largely fails again at `0.125`.
- The globally rate-matched controls and active membrane allocate their finite event budgets to different frequency bands.
- On the A->B/B->A task the active encoder emits usable events in more bodies than the globally thresholded controls, but this too is confounded by the global threshold allocation and is secondary.

## Required v0.2 correction

Before making any AIS-frequency claim, rerun the frozen exact same active model with **stronger oracle controls**:

- rate-match memoryless and linearized controls separately at each frequency to the active event count at that frequency;
- compare event timing only when the active encoder has enough events to estimate timing (`N >= 4` registered minimum);
- report vector strength together with pairwise phase consistency and event-per-cycle coverage;
- treat `f >= 0.05` as the upper-band test, but include only body/frequency pairs passing the exposure gate;
- for the temporal-order task, rate-match controls over the T/D pair rather than over the entire frequency battery.

If the active gate still beats its own linearization under **per-frequency rate matching**, then `h(t)` is buying something beyond merely moving the thresholded event budget to a different frequency band.

If it does not, the correct v0.1 conclusion is simply that the nonlinear membrane reallocates firing across frequency, not that it improves temporal precision.
