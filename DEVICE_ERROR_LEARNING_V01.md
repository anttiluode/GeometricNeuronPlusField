# Device-inspired error learning v0.1 — compressed broadband training survives the registered measurement-error model

The reciprocal broadband learner had already passed two stages:

1. `SPECTRAL_CORRELATION_COMPRESSION_V01.md`: 8–16 boundary-selected spectral channels reconstruct the exact local bond-gradient map with high fidelity;
2. `SPECTRAL_GRADIENT_LEARNING_V01.md`: those compressed gradients preserve useful closed-loop learning when relinearized after every structural update.

The next question was whether that result disappears once the local coherent readout is no longer numerically perfect.

`DEVICE_ERROR_LEARNING_PREREG_V01.md` froze the error model and thresholds before fresh bodies 436–447 were run.

## Error model

This is deliberately called **device-inspired**, not a calibrated model of one particular photonic chip.

For each absolute spectral bin:

```text
phase error       Gaussian phase-setting offset, fixed through training
amplitude error   Gaussian multiplicative calibration error, fixed through training
```

Those two errors are shared across all bonds because they represent port/controller calibration.

Local readout error is different:

```text
tap noise         independent per bond/bin/gradient measurement
                  sigma = stated fraction of that bin's RMS local complex product
```

The registered high-error arm was:

```text
sigma_phi       .050 rad
sigma_amp       .050
sigma_tap       .020
```

with K=8 and K=16 boundary-selected spectral channels.

The moderate descriptive arm used `.025 rad`, `.025`, `.010`.

All arms used the same candidate set, eta=.01, 40 iterations, lag20, and 210-frame transient task.

## C0 — exact learner gate

Registered:

```text
mean exact DeltaC > .015
>=10/12 bodies improve
```

Observed:

```text
mean exact DeltaC             +.026098
improved                        10 / 12
```

Seeds 440 and 442 were essentially flat.

**C0 PASS.**

## C1 — high-error K8 still learns

Registered:

```text
mean K8_high DeltaC > .015
>=10/12 bodies improve
```

Observed:

```text
mean K8_high DeltaC           +.023784
improved                        11 / 12
```

**C1 PASS.**

## C2 — high-error K8 preserves most exact gain

Registered:

```text
mean(K8_high DeltaC) / mean(exact DeltaC) >= .75
mean(K8_high - exact) > -.010
```

Observed:

```text
group gain ratio                .9113
mean K8_high - exact          -.002315
```

**C2 PASS.**

So the high-error K8 arm retained about **91.1%** of the exact learner's group-mean gain on these fresh bodies.

## C3 — high-error K16 preserves most exact gain

Registered:

```text
group gain ratio >= .80
mean difference > -.008
```

Observed:

```text
mean K16_high DeltaC          +.025744
group gain ratio                .9864
mean K16_high - exact        -.000355
improved                        11 / 12
```

**C3 PASS.**

K16 therefore retained about **98.6%** of the exact group-mean gain under this error model.

## C4 — gradient direction remains recognizable throughout learning

Registered:

```text
K8_high  mean map corr > .975
K16_high mean map corr > .990
```

Observed:

```text
K8_high  mean map corr          .99239
         mean relative L2       .12770

K16_high mean map corr          .99556
         mean relative L2       .09373
```

**C4 PASS.**

The noisy compressed map remains directionally close to the exact adjoint while the conductances move.

## C5 — device error versus compression error

Registered:

```text
mean(K8_high DeltaC - K8_ideal DeltaC) > -.007
```

Observed:

```text
K8_ideal mean DeltaC           +.023160
K8_high  mean DeltaC           +.023784
high - ideal                   +.000623
```

**C5 PASS.**

The modeled phase/amplitude/tap errors did not add a measurable group-mean penalty beyond the already-present K=8 spectral truncation in this run.

The fact that the noisy arm is slightly higher is **not** evidence that measurement noise improves the physical gradient. The finite normalized optimizer is path-dependent, and noise can perturb a nearly flat/clipped trajectory in either direction. Seed 440, where exact and ideal K8 are flat while the noisy arms wander to small positive gains, is an obvious example.

## Formal verdict

```text
C0 PASS   exact learner useful on fresh bodies
C1 PASS   high-error K8 learns
C2 PASS   high-error K8 retains 91.1% of exact mean gain
C3 PASS   high-error K16 retains 98.6% of exact mean gain
C4 PASS   map correlations remain .992/.996 through training
C5 PASS   modeled device error does not dominate K8 truncation error
```

**6 / 6 registered criteria pass.**

## Moderate arm, descriptive

```text
K8_mod mean DeltaC             +.023057
improved                        10 / 12
mean map corr                    .99361
mean relative L2                 .10872
```

As expected, it is essentially indistinguishable from the ideal K8 arm at this scale.

## What this earns

The careful statement is:

> **Under this explicit device-inspired error model, the compressed broadband reciprocal gradient remains useful for closed-loop training at 0.05-rad shared spectral phase-setting error, 5% shared spectral-amplitude error, and local tap/readout noise equal to 2% of each selected bin's RMS complex-product magnitude.**

The result does **not** establish that a real silicon, acoustic, mechanical, or transmission-line implementation will have exactly these error statistics or the same tolerance. Mapping a platform's physical phase, detector, loss, thermal, quantization, and nonreciprocity errors into these variables is still necessary.

## Important distinction: reciprocal loss versus pass-to-pass drift

Uniform reciprocal attenuation by itself does not necessarily invalidate the physical adjoint principle; if it is part of the same linear reciprocal operator in both passes, the adjoint includes it.

The more dangerous hardware error is **operator mismatch between forward and backward/adjoint use**:

```text
thermal drift between passes
nonreciprocal components
changing tuner state
frequency-dependent calibration mismatch
asymmetric loss/coupling
```

That should be the next robustness experiment. It attacks the symmetry that makes same-medium adjoint transport exact, rather than merely corrupting the local phasor readout.

## Current hardware status

The line is now:

```text
broadband transient task
      -> reciprocal physical adjoint
      -> K=8/16 port-selected spectral compression
      -> device-inspired coherent-readout errors
      -> local graded update
      -> successful closed-loop training
```

The largest unresolved implementation question is no longer raw 210-frame memory or modest readout noise. It is how to realize the K local phasor products physically and how robust the scheme remains when the **forward and adjoint operators themselves cease to match perfectly**.
