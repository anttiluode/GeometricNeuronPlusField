# AIS selection-information v0.1 — the "selection at a cost of precision" rescue also fails

Canonical GitHub Actions run: `31245258349`

After the final phase-interface test, the only plausible positive AIS-like role
left in this toy was Claude's sharper tradeoff hypothesis:

> the nonlinear active boundary might sacrifice spike-time precision but make
> each emitted event more informative about which temporal regime is present.

This experiment froze the body, field, HH parameters, interface definitions and
frequency battery.  For each body/feed, the memoryless and own-linearization
controls were given **exactly the same total number of events over the whole
six-frequency battery** as the active boundary.  Unlike v0.2, the controls were
not matched separately at each frequency, so frequency allocation was free to
differ.

The primary feed was `Re(psi_soma)`.

## Registered primary statistic

For six equally weighted frequency conditions and counts `n_f`:

```text
p(f | spike) = n_f / sum_f n_f

I_spike = KL[p(f | spike) || Uniform(6)]
        = log2(6) - H[p(f | spike)]
```

This measures how much **frequency identity is carried by the occurrence of one
spike**, in bits/spike.  All encoders have the same total spike budget within a
body.

The success rule required at least 12 valid Re bodies, positive median
`active-linearized`, and one-sided paired Wilcoxon `active > linearized p < .05`.

## Primary Re(psi) result

All 24 bodies were valid.

```text
mean total events/body        66.67   (identical for all encoders)

mean bits/spike
active                        0.1896
linearized                    1.5962
memoryless                    1.7612

active - linearized mean     -1.4066 bits/spike
active - linearized median   -1.3885
active better / worse          0 / 24
sign-test p                   1.19e-7
Wilcoxon active>linear p      1.000

active - memoryless mean     -1.5717 bits/spike
active better / worse          0 / 24
```

The preregistered positive selection hypothesis fails **unanimously and in the
opposite direction**.

The secondary full binary event-information measure agrees:

```text
I(F; spike/no-spike), mean bits/frame
active        0.00314
linearized    0.03501
memoryless    0.03940
```

So the active boundary is not merely less precise in *when* it emits events. At
the phase-bearing interface, its nonlinear/refractory dynamics also make event
occurrence **less diagnostic of input frequency** than either matched control.

## Why the controls win

The count allocations make the result intuitive.  Mean counts per body under
`Re(psi)`, with identical total event budgets, were approximately:

```text
f            .006   .0125   .025   .050   .083   .125
active       11.08   14.58   9.88   8.17   14.38   8.58
linearized   13.38   31.96   9.08   4.04    7.79   0.42
memoryless   14.38   41.25   5.29   3.62    1.75   0.38
```

The active HH dynamics **flatten** the frequency allocation.  The linear and
memoryless detectors concentrate their matched event budgets much more strongly
in the frequency regimes where their scalar score is largest.  Under this
registered coding metric, that concentration carries substantially more
frequency information per event.

This is almost the reverse of the story we were trying to rescue: refractory /
state dynamics are acting like a normalizer or equalizer across regimes, not an
information-enhancing selector.

## Historical power feed

The original visually striking frequency allocation under `|psi|^2` does not
rescue the claim either.

```text
mean bits/spike
active        1.3833
linearized    2.3173
memoryless    2.3122

active-linear mean   -0.9340
active better/lower    0 / 18   (6 exact ties)
```

Mean power-feed allocations make clear why the old plot looked selective but
was still weaker than the matched controls:

```text
f            .006   .0125   .025   .050   .083   .125
active        .020    .180   .594   .110   .087   .009   mean fraction
linearized    .006    .077   .842   .075   .000   .000
memoryless    .014    .110   .852   .024   .000   .000
```

The active gate is frequency dependent, but the passive scores are **even more
frequency selective** at the same total event count.

Magnitude behaves similarly and loses in all 24 bodies.

## What the AIS-like compartment has and has not earned

The accumulated receipts now separate "consequential" from "computationally
privileged":

```text
changes firing regime / event availability                 YES
sensitive to h/n/K kinetics                                YES
clean tunable passband                                      NO
simple n-set refractory clock                               NO
better rate-matched timing precision                        NO -- worse
special benefit from preserving carrier phase               NO -- worse
more frequency information per matched spike                NO -- much worse
```

So there is currently **no positive performance criterion on which this generic
HH-like downstream compartment beats the simpler observation controls**.

It absolutely transforms the signal.  But transformation alone is not enough
to justify saying that the active state is doing useful computation in this
toy.

## Consequence

This is the stopping rule we wanted when the line got complicated:

> **Do not co-adapt AIS position, extent, channel density, or body geometry around
> this compartment.**

Doing so now would create enough degrees of freedom to tune around a sequence of
registered nulls/negative results.

The upstream geometry/field/soma results remain intact.  What failed is the
specific generic-HH bridge from the soma mixture to a privileged event code.

A future active-boundary branch would need a new reason to exist — for example a
concrete downstream task that requires sparse events, energetic constraints, or
a biologically motivated non-HH channel architecture — and it should start as a
new hypothesis rather than a rescue of this one.

## Wall sentence

> **The HH-like boundary does not buy frequency information at the price of
> timing precision; it loses both.  At a matched total spike budget its events
> are dramatically less frequency-informative than the same membrane's
> linearization.  The AIS co-adaptation branch therefore stops here.**
