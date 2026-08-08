# Hardware compiler v0.1 — the toy wave maps directly to a reciprocal LC or mass-spring network

The literature search says the general physical-adjoint idea is established. The useful next step is therefore not to inject an arbitrary “0.025 rad phase error” into the toy. The model parameter must first be mapped to an actual device parameter.

For this linear wave, that mapping is unusually direct.

## 1. Continuous equation behind the discrete update

Ignoring the already-certified negligible saturation term, the mature model is the explicit discretization of

```text
psi_ddot_i
+ gamma psi_dot_i
+ rho psi_i
+ c sum_j k_ij (psi_i - psi_j)
= s_i(t)
```

where `k_ij=k_ji` is a reciprocal local coupling.

The complex `psi` used in the code is an analytic/two-quadrature convenience. The real and imaginary parts obey the same real equation independently.

## 2. Lumped reciprocal LC network

Use node flux `Phi_i(t)` as the physical state. Put at every node:

- shunt capacitance `C0`;
- shunt conductance `G0` for loss;
- optional shunt inductance `L0` for the onsite restoring term.

Connect neighbouring nodes `i,j` with reciprocal coupling inductors `L_ij`.

Kirchhoff's law gives

```text
C0 Phi_ddot_i
+ G0 Phi_dot_i
+ (1/L0) Phi_i
+ sum_j (1/L_ij)(Phi_i-Phi_j)
= I_i(t).
```

Divide by `C0` and identify

```text
psi_i        <-> Phi_i

gamma        = G0 / C0
rho          = 1 / (C0 L0)
c k_ij       = 1 / (C0 L_ij)
s_i          = I_i / C0.
```

So the model's trainable local conductance-like variable is physically an **inverse coupling inductance** (or an equivalent tunable reactive coupling in a transmission-line/metamaterial realization).

This is the right dimensional bridge for component tolerances. A phase-shifter error in radians is not.

## 3. Mechanical mapping

The same equation is also a damped reciprocal mass-spring network:

```text
m x_ddot_i
+ b x_dot_i
+ k0 x_i
+ sum_j kappa_ij (x_i-x_j)
= F_i(t)
```

with

```text
psi_i       <-> x_i
gamma       = b/m
rho         = k0/m
c k_ij      = kappa_ij/m
s_i         = F_i/m.
```

This makes the local forward/backward gradient particularly tangible:

```text
dJ/dkappa_ij
proportional to
integral (x_i-x_j)_forward (x_i-x_j)_adjoint dt.
```

The exact circuit analogue uses the forward and adjoint **branch flux/voltage differences**.

## 4. The soma/objective port

The current readout is local square-law power at one designated node. In hardware, that is a detector at one port/node after the wave medium.

Its task derivative creates the error waveform. Reciprocity then says:

```text
forward trial:       excite input ports -> record local branch histories/output
backward trial:      excite output port with time-reversed derivative waveform
local gradient:      integrate forward-branch difference x backward-branch difference
```

This is the time-domain counterpart of the physical-adjoint / TRIM literature.

## 5. What a local gradient sensor would measure

For coupling edge `e=(i,j)`, the exact sensitivity depends only on the two endpoint differences in the aligned forward and adjoint histories.

A device therefore does **not** need a global field camera to obtain a local gradient. It needs, in principle, access to the branch-local forward and returned signals, or an interference/intensity protocol from which their product can be recovered.

That is the hardware-economics question:

```text
one monitor/multiplier per tunable bond?
time-multiplexed probe access?
local analog accumulation?
three-intensity TRIM-like reconstruction?
```

The mathematics says what must be measured; the device architecture decides its cost.

## 6. Hardware-relevant nonidealities now have model coordinates

For an LC/transmission-line implementation, the first degradation suite should be expressed in quantities such as:

### Coupling calibration / fabrication spread

```text
k_e -> k_e (1 + epsilon_e)
```

or equivalently inverse-inductance error.

### Coupling-control quantization

Quantize `rho_e` or `1/L_e` to a finite number of control levels.

### Loss / Q variation

```text
gamma_i = gamma (1 + epsilon_gamma,i)
```

and possibly branch-dependent loss if the implementation requires it.

### Onsite-frequency disorder

```text
rho_i = rho (1 + epsilon_rho,i)
```

from capacitor/inductor tolerance.

### Error-wave timing fidelity

Jitter, finite waveform bandwidth, delay mismatch, and truncation of the time-reversed derivative waveform.

### Measurement noise

Noise on the local forward/adjoint branch observables before their product is accumulated.

### Thermal / slow drift

A slowly changing multiplicative bias on tunable reactive elements between calibration and update.

These are interpretable in this model. A silicon MZI phase error becomes relevant only after compiling `k_e` into an MZI/scattering parameter.

## 7. A serious practical number already exposed by the compiler

The historical mature parameters are

```text
k_arbor       = 2.5
k_bath        = 0.0002
ratio         = 12,500 : 1.
```

Under the LC mapping,

```text
L_bath / L_arbor = 12,500
```

if both are realized as finite inductive couplers.

That is an extreme continuous dynamic range. A sparse physical network could instead interpret the bath as an **open/absent edge**, with only selected couplers fabricated or switch-enabled, but then topology creation and continuous relaxation need a switchable/tunable coupler architecture.

This is a real hardware issue the abstract model had hidden.

It motivates a future coupling-contrast sweep:

> How much of the temporal computation and local-learning behaviour survives when the strong/weak coupling ratio is reduced from 12,500 toward realistic tunable-device ranges?

That is a more meaningful physical-degradation experiment than transplanting a photonic phase-noise number into an unrelated variable.

## 8. Relation to the newest close prior art

Thakkar & Grbic (2026) describe a tunable 2D transmission-line metamaterial network with local reactive elements and in-situ wave-based training. That is now the most obvious hardware family to compare against directly.

The remaining possible distinction here is not “a tunable 2D wave grid exists.” It is the particular package:

```text
sparse/morphology-like coupling pattern
+ time-domain pulse-history task
+ local square-law/interference mechanism
+ graded structural response
+ exact time-domain reciprocal adjoint
+ compiler/cost comparison against globally optimal spectral coordinates.
```

## Next hardware experiment

Before a fabrication-noise paper, sweep **coupling contrast and component disorder** in this LC-normalized parameterization, while holding the already-preregistered task and tuner budgets fixed.

That will tell us whether the current mechanism lives only in the numerically extreme 12,500:1 topology limit or survives in a realizable reciprocal circuit regime.
