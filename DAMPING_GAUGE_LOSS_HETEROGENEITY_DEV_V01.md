# Damping-gauge loss calibration and heterogeneity — development v0.1

Date: 2026-08-08

## Result in one sentence

> **The exact conformal echo requires one scalar loss factor, but useful gradient direction is considerably more tolerant: on reused development arbors, a 10% error in calibrating mean residual loss `eps=.005` still gives gradient correlation about `.997`, and substantial cell-to-cell loss variation retains high correlation as long as accumulated mean loss remains moderate.**

This is development-only robustness evidence, not a hardware demonstration or novelty claim.

---

## Setup

The exact uniformly lossy echo compiler from `DAMPING_GAUGE_LOSS_COMPENSATED_DEV_V01.md` assumes

```text
x[n+1] = M x[n] - a x[n-1] + u[n]
```

with one spatially uniform scalar `a=1-eps`.

It uses

```text
terminal reverse scale      a
reverse-source envelope     a^(j+1)
detector-integrator weight  a^(-j)
```

to recover the exact gradient of the actual uniformly lossy device.

`damping_gauge_loss_heterogeneity_probe.py` damaged this in two ways on reused seeds 472–475.

### A. Wrong scalar calibration

The physical core remains uniformly lossy with `eps_true=.005`, but the echo envelopes use

```text
eps_hat = eps_true * (1 + relative_error).
```

### B. Spatially nonuniform loss

Each cell receives a diagonal loss coefficient `eps_i` around a chosen body mean. The actual device recurrence and its exact digital adjoint use that full spatial loss map.

The physical echo is nevertheless allowed only one scalar compensation factor: the known mean loss.

This is deliberately hostile because generic spatial damping does not commute with the wave operator, so a single scalar conformal factor is no longer exact.

---

# A. Calibration error at mean eps=.005

Observed pooled gradient correlations across four bodies:

```text
relative error in eps_hat    mean corr      minimum corr
-10%                         .995783        .992114
 -5%                         .999021        .998173
 -2%                         .999850        .999721
 -1%                         .999963        .999931
-.5%                         .999991        .999983
  0                          1.000000        1.000000
 +.5%                        .999991        .999983
 +1%                         .999964        .999934
 +2%                         .999859        .999738
 +5%                         .999158        .998443
+10%                         .996877        .994259
```

The error primarily changes gradient magnitude before it substantially changes direction.

At +10% calibration error:

```text
mean relative L2 error   .1626
mean norm ratio          1.1401
mean correlation         .9969
```

So exact step-size calibration would be wrong, but the descent direction remains highly aligned.

---

# B. Spatial loss heterogeneity

Loss heterogeneity is reported as the multiplicative coefficient-of-variation parameter used to generate cellwise `eps_i`, renormalized to the requested body mean.

## Mean residual loss eps=.001

```text
loss CV     mean corr       min corr
0           1.000000        1.000000
.01          .9999998        .9999995
.05          .9999943        .9999833
.10          .9999829        .9999566
.20          .9999447        .9998396
.50          .9994536        .9986448
```

Even 50% relative spatial variation is almost irrelevant at this small accumulated loss.

## Mean residual loss eps=.005

```text
loss CV     mean corr       min corr
0           1.000000        1.000000
.01          .9999929        .9999795
.05          .9998644        .9995711
.10          .9994704        .9988205
.20          .9981200        .9917078
.50          .9848318        .9596618
```

At the attractive `eps=.005` operating point, even a very large 50% loss CV leaves the average gradient direction strongly aligned.

## Mean residual loss eps=.01

```text
loss CV     mean corr       min corr
0           1.000000        1.000000
.01          .9999725        .9999363
.05          .9993461        .9988235
.10          .9978307        .9926419
.20          .9884977        .9351399
.50          .9145841        .7861784
```

At larger accumulated loss, nonuniformity begins to matter seriously.

The same nominal fractional heterogeneity is more damaging because the absolute noncommuting loss perturbation is larger and acts over the full transient window.

---

## Engineering interpretation

The clean hierarchy is now:

```text
uniform known loss
    exact scalar conformal compensation

uniform miscalibrated loss
    gradient magnitude biased first
    direction remains excellent over useful errors

moderately nonuniform loss
    no exact scalar compiler
    but direction remains highly aligned while mean loss is small/moderate

large accumulated + strongly nonuniform loss
    gradient direction eventually degrades
```

For the present 210-step task, `eps=.005` is particularly interesting:

```text
end detector gain a^(-T)   ~2.87x
10% mean-loss calibration error -> corr ~.997
20% spatial loss CV             -> corr ~.998
50% spatial loss CV             -> corr ~.985
```

This does not yet combine all errors simultaneously.

---

## Mathematical boundary

For spatially varying damping represented by a diagonal matrix `A`, a natural matrix scaling such as

```text
x_n = A^(n/2) z_n
```

does not in general leave a fixed spatial wave operator because `A` need not commute with the coupling/stiffness operator.

So the exact one-envelope compiler is genuinely a proportional/uniform-loss special case.

The robustness result says only that the error introduced by violating that special case can be small enough to preserve gradient direction in the tested regime.

---

## Current hardware candidate

A plausible development operating region is now:

```text
transient length          T ~210
mean residual loss        eps ~.005 or below
end detector gain         <= ~3x
mean-loss calibration     only percent-to-10%-level accuracy needed for direction
spatial loss disorder     tens of percent tolerable in development
```

Separate development probes also found:

```text
10% momentum-reversal error    gradient corr ~.997
5% reverse-operator drift      gradient corr ~.996
10% RMS gradient readout noise gradient corr ~.995
```

Those errors have not yet been combined in one hostile run.

## Next wall

Combine realistic imperfections rather than varying them one at a time:

```text
residual loss + spatial disorder
loss-calibration error
terminal time-mirror error
reverse-operator drift
local integrator noise
```

Then test whether a gradient update made from that physical surrogate still improves the task, rather than merely correlating with the exact map.

## Wall sentence

> **The conformal echo does not require laboratory-perfect scalar loss: in the present transient window, gradient direction survives surprisingly large loss-calibration errors and substantial spatial loss disorder. The next meaningful test is no longer an identity test—it is closed-loop learning under a realistic combined hardware-error cocktail.**
