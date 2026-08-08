# Calibrated uniform-loss echo compensation — development v0.1

Date: 2026-08-08

## Result in one sentence

> **A uniformly lossy compiled core does not require anti-damping to recover its exact finite-time gradient: the same lossy dynamics generate an exponentially attenuated reverse-forward echo, and a matching global reverse-input envelope plus global detector-integration envelope cancels that attenuation exactly, recovering the actual device gradient to machine precision without local field-history storage.**

This is an algebraic/compiler development result on reused bodies, not a novelty claim or hardware demonstration.

---

## Why this test followed the first damage sweep

The original ideal damping-gauge compiler transformed the repository's intended uniformly damped recurrence into a conservative reversible core.

A first residual-loss damage model was malformed: it changed the coefficient on the state two steps back without applying the matching onsite shift produced by the same damped integrator. Because the compiled operator contains eigenvalues close to `+2`, that malformed perturbation crossed a discrete stability boundary and generated an artificial catastrophic failure.

`damping_gauge_residual_loss_v02.py` corrected the physical recurrence to

```text
x[n+1] = (Q - eps I) x[n] - (1-eps) x[n-1] + u[n].
```

The corrected un-compensated same-body echo was stable and substantially more robust:

```text
eps       mean gradient correlation
0         1.000000
.0001      .999887
.0005      .996784
.001       .984836
.002       .918659
.005       .359432
.01       -.253987
```

The returned adjoint remained exactly aligned with the digital reverse-mode adjoint; the error came from failing to reconstruct the time-reversed forward trajectory.

That observation suggested a second scalar gauge directly on the physical echo.

---

# Exact lossy reverse identity

Let the actual physical compiled core obey

```text
x[n+1] = M x[n] - a x[n-1] + u[n],
0 < a <= 1,
```

with spatially uniform scalar loss and symmetric `M`.

Define the reverse-indexed forward trajectory

```text
r[j] = x[T-j].
```

Ordinary reverse dynamics would require inverse damping. But define instead

```text
y[j] = a^j r[j]
     = a^j x[T-j].
```

Then exact algebra gives

```text
y[j+2]
  = M y[j+1]
    - a y[j]
    + a^(j+1) u[T-1-j].
```

That is the **same lossy physical recurrence** as the forward device.

Therefore the lossy body can regenerate an attenuated exact reverse-forward trajectory with

```text
y[0] = x[T]
y[1] = a x[T-1]
```

and a globally known reverse-source envelope `a^(j+1)`.

No anti-damped interior and no local waveform tape are required.

---

# Adjoint field

For the symmetric recurrence, the causal returned adjoint already obeys the same lossy operator:

```text
b[j] = p[T-j+1].
```

This was numerically exact in the corrected residual-loss probe.

Thus during the reverse experiment the two local fields are

```text
retraced-forward field    Delta y[j] = a^j Delta x[T-j]
returned adjoint field    Delta b[j] = Delta p[T-j+1].
```

The only mismatch from the desired overlap is the known scalar factor `a^j`.

---

# Global integration-envelope compensation

The exact local bond gradient needs

```text
sum_j Re(conj(Delta x[T-j]) Delta b[j]).
```

The physical `+/-` branch-energy interference gives

```text
Re(conj(Delta y[j]) Delta b[j])
  = a^j Re(conj(Delta x[T-j]) Delta b[j]).
```

Therefore weight the local energy accumulator by the same **global clock envelope**

```text
a^(-j).
```

Then

```text
sum_j a^(-j) Re(conj(Delta y[j]) Delta b[j])
 = exact gradient overlap.
```

Every trainable branch receives the same scalar time-dependent detector weight. No location-specific loss history or field reconstruction is needed as long as loss is spatially uniform and calibrated.

---

# Development audit — reused seeds 472–475

`damping_gauge_loss_compensated_echo.py` compared the physically compensated `+/-` energy readout against the exact reverse-mode gradient of the **actual lossy device**, not merely against the ideal lossless compiler.

For

```text
eps = 0,
      .0001,
      .0005,
      .001,
      .002,
      .005,
      .01,
      .02,
      .05
```

all four bodies gave essentially exact gradients.

Representative pooled result:

```text
eps       mean corr             mean relative L2
0         1.000000000000        5.2e-14
.0001     1.000000000000        7.4e-14
.0005     1.000000000000        1.4e-13
.001      1.000000000000        7.8e-14
.002      1.000000000000        9.7e-14
.005      1.000000000000        1.5e-13
.01       1.000000000000        1.8e-13
.02       1.000000000000        2.3e-13
.05       1.000000000000        8.4e-12
```

The attenuated reverse-forward trajectory itself matched the analytic target to roughly `1e-14`, and the returned adjoint alignment was exact to numerical precision.

---

# Dynamic-range cost

Exactness is not free. The final reverse-time detector gain is

```text
a^(-T)
```

and the final reverse-source envelope is

```text
a^T.
```

For `T=210`:

```text
eps      detector gain at end    reverse-source scale at end
.0001          1.0212                  .9792
.0005          1.1107                  .9003
.001           1.2338                  .8105
.002           1.5226                  .6568
.005           2.8652                  .3490
.01            8.2529                  .1212
.02           69.586                   .01437
.05        47647.8                      2.10e-5
```

So the mathematically exact compiler becomes physically unattractive once `eps*T` is large, even though no instability occurs.

The important engineering quantity is therefore not simply residual loss per step; it is accumulated loss over the transient window.

---

# What this changes

The hardware wall has moved again.

The statement

```text
physical residual damping breaks the echo
```

is too strong for **uniform calibrated loss**.

The more precise statement is:

```text
uniform calibrated loss
    -> exact scalar echo compensation
    -> no local history
    -> dynamic-range cost exp(loss * duration)
```

The genuinely dangerous deviations are now:

```text
wrong estimate of the loss factor
spatially nonuniform damping
frequency-dependent/non-proportional loss
pass-to-pass changes in loss
finite dynamic range / detector floor
imperfect implementation of the terminal conformal time-reversal state
```

---

# Prior-art boundary

The scalar scaling/conjugacy itself belongs to established conformally symplectic / linearly damped Hamiltonian mathematics. Recent work on conformal-symplectic map learning explicitly represents dissipative maps as a symplectic core plus damping/scaling structure.

Hamiltonian Echo Backpropagation and Recurrent Hamiltonian Echo Learning already own the broad use of physical time-reversal echoes for finite-time gradient credit in reversible Hamiltonian systems.

Therefore the candidate research question is not the existence of an exponential integrating factor.

It is the hardware/compiler conjunction:

> **Can a dissipative finite-time scattering computation be represented and trained in a physical echo architecture using only calibrated global conformal envelopes plus local scalar energy accumulation, with acceptable dynamic range and robustness to the ways real loss departs from a scalar uniform factor?**

## Next wall

Damage the scalar assumption:

1. loss-calibration error;
2. spatially heterogeneous loss;
3. time/pass-varying loss;
4. detector-envelope gain error;
5. terminal-state scaling error;
6. longer windows and accumulated dynamic range.

## Wall sentence

> **Uniform loss is not the fundamental memory wall: a lossy reciprocal recurrence can echo an attenuated exact reverse trajectory through the same dynamics, and a shared inverse-loss integration envelope restores the full transient gradient exactly. The real wall is whether physical loss is uniform and calibratable enough—and whether the required exponential dynamic range stays affordable.**
