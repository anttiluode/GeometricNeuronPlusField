# AIS n-kinetics / ISI v0.1 — n speed does not set a simple refractory clock

Canonical GitHub Actions run: `31244817807`

Claude's HH audit suggested a sharper interpretation test for the earlier
`n_scale` result.  Fast `n` had strongly reduced event count, but that alone
could reflect either a changed refractory period or a changed excitability
window.  The existing kinetics JSON did not retain spike times, so this
experiment reran the exact frozen battery and recorded interspike intervals.

The preregistered primary frequency was `f = 0.025 cycles/frame`, the already
well-exposed dominant regime from AIS_KINETICS_V01.

The prediction for a simple n-set refractory floor was:

```text
minISI(n_slow) > minISI(full) > minISI(n_fast)
log2[minISI(n_fast) / minISI(n_slow)] < 0
```

## Primary result

18 bodies had at least two post-burn spikes in both the slow- and fast-n arms.

```text
mean log2(fast / slow minISI)     +0.00538
median                             +0.13877
fast shorter / fast longer          9 / 9
sign test p                         1.000
one-sided Wilcoxon fast<slow p      0.5169
strict slow > full > fast            6 / 18
```

The preregistered refractory-clock prediction fails completely.

There is no consistent tendency for faster `n` kinetics to shorten the minimum
ISI relative to slow `n`.

## The surprising shape

Across bodies, the native `n_scale = 1` condition actually has the shortest
mean minimum ISI at the primary frequency:

```text
condition      emitting bodies   mean total events   mean minISI@.025   mean medianISI@.025
n slow 0.5x        24 / 24            57.75              15.89                22.02
full 1.0x          24 / 24            51.54              10.47                16.39
n fast 2.0x        19 / 24            14.54              18.45                37.05
```

So both moving `n` slower **and** moving it faster tend to lengthen the short-ISI
edge relative to the native setting, while fast `n` also crushes total firing.

A post-hoc paired look makes that U-shaped tendency concrete:

- slow `n` has a longer minimum ISI than full in all 24 bodies;
- among the 18 fast-n bodies with a measurable primary ISI, fast `n` also tends
  to have a longer minimum ISI than full.

Those comparisons were not the registered primary test, so they are mechanism
clues rather than a new confirmed claim.

## Interpretation

The clean simple stories are now both unsupported:

```text
n timescale translates a frequency passband     NO
n timescale directly sets refractory floor       NO
```

What remains is a dynamical-window picture.  The native balance of delayed K
activation with Na activation/inactivation appears to support the most rapid
repetitive firing in this particular driven boundary.  Moving `n` either way
changes that balance; speeding it strongly can push bodies out of the firing
regime altogether.

That is better described as **state-dependent excitability / regime control**
than as a single tunable frequency or refractory knob.

This also reinforces Claude's criticism of the old instantaneous-n ablation:
when an intervention destroys the firing regime, it cannot by itself identify
the normal computational role of that state variable.

## Consequence for AIS geometry

There is still no justified scalar AIS parameter for body geometry to co-adapt
around.  `n` matters strongly, but its role in this generic HH boundary is
non-monotonic and operating-point dependent.

The next positive question is therefore whether active eventization carries a
useful **selection code** despite sacrificing timing precision — not where to
move or stretch an AIS-like compartment.

## Wall sentence

> **Changing K-gate kinetics strongly changes whether and how often this
> boundary fires, but minimum ISI does not scale with the n time constant in the
> predicted way.  The native n timescale looks like part of an excitability
> optimum, not a simple refractory clock.**
