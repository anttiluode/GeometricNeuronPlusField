# AIS final phase/interface v0.1 — phase-preserving input does not rescue active timing

Canonical GitHub Actions run: `31244811210`

This experiment changed only the scalar crossing the frozen soma -> active-boundary interface:

```text
power       |psi_soma|^2
magnitude   |psi_soma|
real        Re(psi_soma)
```

The body, wave carrier, HH equations, conductances, gain, clipping magnitude,
frequency battery, burn-in, minimum spike exposure, and exact per-frequency
event-count matching were frozen from AIS active v0.2.

The important conceptual distinction is that `|psi|^2` and `|psi|` are envelope
observables; `Re(psi)` retains the signed carrier waveform and is therefore the
registered final phase-bearing arm.

## Primary result — Re(psi)

The primary statistic first averages the valid upper-band (`f >= 0.05`) active
minus linearized PPC differences **within each body**, so frequencies from one
organism are not counted as independent organisms.

```text
valid bodies                         23 / 24
mean body delta PPC                 -0.16435
median body delta PPC               -0.14955
active better / worse bodies          3 / 20
sign test, two-sided p               0.000488
Wilcoxon, two-sided p                0.0000131
Wilcoxon active > linearized p       0.999995
```

The preregistered success rule required at least 8 valid bodies, positive median
delta, and one-sided Wilcoxon `active > linearized p < .05`.

It fails in the opposite direction.

> **Giving the active boundary the remaining phase-bearing soma scalar does not
> make nonlinear eventization more temporally precise than the same membrane's
> linearized response.  Across bodies, the active spikes are substantially less
> envelope-phase consistent.**

This is much better exposed than the original power test: `Re(psi)` gives 23
valid upper-band bodies rather than 11, so the negative result is not a sparse
spiking artifact.

## Per-frequency Re(psi) receipt

```text
f          valid   mean active events   active - linear PPC
0.00625      22          11.08                 -0.6964
0.01250      24          14.58                 -0.5592
0.02500      23           9.88                 -0.5301
0.05000      18           8.17                 -0.1269
0.08333      23          14.38                 -0.3246
0.12500      20           8.58                 -0.0132
```

The active boundary comes closest at the highest tested modulation frequency,
but there is no registered upper-band timing advantage.

## The carrier itself does not rescue the claim either

Carrier-phase PPC was registered as a secondary descriptive receipt.  `Re(psi)`
is the only feed that actually supplies carrier phase to the HH boundary.

Representative mean carrier PPC values:

```text
f             linearized      active
0.00625          0.795          0.475
0.01250          0.758          0.164
0.02500          0.780          0.321
0.05000          0.941          0.948
0.08333          0.715          0.397
0.12500          0.869          0.781
```

There is a tiny active edge at `f=0.05`, but it does not generalize.  The
body-level upper-band carrier-PPC difference is negative overall as well.

## Envelope feeds also point the same way

The result is not peculiar to signed current:

```text
feed          valid bodies   mean body delta   median delta   active/linear W/L
power              11           -0.1643          -0.1140          2 / 9
magnitude          21           -0.1613          -0.1533          5 / 16
real               23           -0.1644          -0.1495          3 / 20
```

For magnitude, the two-sided body-level Wilcoxon is `p ~= .0043`; for real it is
`p ~= 1.31e-5`.  Power has the same effect direction but weaker exposure.

## What this closes

This is intentionally a lineage-local statement, not a biological universal.

The phase question has now been attacked at the last remaining interface in
this model.  Preserving the signed carrier at eventization **does not earn an
active timing advantage**.  Therefore phase is not currently a justified reason
to introduce AIS position/extent co-adaptation.

The stronger surviving picture is a tradeoff hypothesis:

```text
linearized observation   -> more precise timing
active event boundary    -> nonlinear selection / availability, less precise timing
```

That second half still needs a positive information-coding test before being
called computationally useful.

## Wall sentence

> **Restoring the soma carrier with Re(psi) greatly improves experimental
> exposure but not active timing: the full HH-like event boundary remains less
> phase-precise than its own rate-matched linearization.  Phase therefore fails
> the final eventization test in this lineage.**
