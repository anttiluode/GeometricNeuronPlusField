# AIS h/n kinetics v0.1 — the simple time-constant/passband prediction does not survive

Canonical run: `31243442742`

After the instantaneous-gate ablation, this experiment made a gentler parameter-free mechanistic test.  The active boundary, upstream bodies, normalization, gain, conductances and gate steady-state curves were frozen.  Only the differential speed of `h` or `n` was multiplied by `0.5` or `2.0`.

The preregistered signature was deliberately simple:

> if a gate timescale sets the event-selection band, slowing and speeding that gate should move the normalized frequency-allocation center in opposite directions across organisms.

That signature did **not** appear.

## Full model

```text
emitting bodies             24 / 24
mean total events           51.54
mean allocation center      0.02652 cycles/frame
median allocation center    0.02500
mid/high band share         0.1936
```

## Na inactivation kinetics h

```text
                 h x0.5 (slow)     h x2 (fast)
emitting bodies       24/24             24/24
mean total events      42.83             49.88
mean center             .02731            .02718
mean log2 shift          +.049             +.024 octaves
TV from full             .057              .089
```

Fast versus slow did **not** produce a stable ordering:

```text
fast higher / lower center    8 / 12
mean fast-slow shift          -0.025 octaves
sign p                         0.503
```

So modest changes in `h` speed alter gain and shape a little, but there is no evidence here that the inactivation time constant directly sets a movable passband.

## K activation/recovery kinetics n

`n` speed has a substantially larger effect on the response shape:

```text
                 n x0.5 (slow)     n x2 (fast)
emitting bodies       24/24             19/24
mean total events      57.75             14.54
mean center             .02734            .03077
mean log2 shift          +.038             +.135 octaves
TV from full             .184              .345
```

But the registered bidirectional-timescale signature again fails:

```text
valid slow/fast pairs          19
fast higher / lower center     11 / 7
mean fast-slow shift           +0.084 octaves
sign p                          0.481
```

Faster `n` also pushes five bodies out of the emitting regime entirely and strongly suppresses total firing.  The large TV distance is therefore real sensitivity, but it is not a clean tunable-frequency result.

The mean allocation makes the effect concrete.  With fast `n`, surviving output becomes even more dominated by the `0.025` condition while both the lowest and highest conditions nearly disappear.  Slow `n` spreads more output into both low and high conditions.  That looks more like a change in **excitability window / resonance sharpness** than a simple translation of a passband.

## What the three AIS experiments now say together

### v0.2 precision test

The stateful active boundary does **not** beat its own rate-matched linearization in upper-band phase precision.

### instantaneous gate ablation

Independent `h` and `n` histories are required to maintain the full spiking regime at the frozen operating point; `m` history is largely dispensable.

### kinetics test

Changing the h/n memory speeds certainly changes output, especially for `n`, but it does **not** move a frequency-selection band in the preregistered clean slow/fast manner.

So the current model has earned this statement:

> **Active channel state controls which upstream fluctuations become spikes and how excitable the event boundary is, but this generic HH-like compartment has not earned the stronger claim that its gating time constants implement a clean tunable frequency filter or improve spike-time precision beyond a matched linear filter.**

That is a meaningful stopping point before co-adaptation.

## A newly exposed interface question

There is also a modeling issue that should be confronted before changing AIS geometry: the current active boundary is driven by the historical FunctionalArbor soma **power** readout `|psi_soma|^2`.

That was the correct quantity for the old objective/readout tests, but a biological AIS is driven by a voltage-like membrane signal, not by squared field magnitude.

Therefore the next clean bridge experiment, if this line continues, should change **only the interface variable**, not the HH parameters:

```text
old:     |psi_soma|^2  -> active boundary
compare: |psi_soma|    -> active boundary
and/or:  Re psi_soma   -> active boundary
```

with the same memoryless and linearized controls.

That should be preregistered as a new experiment rather than used to rescue the current null.  If the active-state story depends on giving the boundary a voltage-like waveform rather than an already-rectified power envelope, that would be conceptually important in its own right.

## Wall sentence

> **The AIS-like state is consequential but not yet computationally privileged: it changes event availability, especially through K recovery, while failing both the precision test and the clean kinetic-passband test.  Before co-adapting geometry, test whether we have been feeding the active boundary the wrong physical observable.**
