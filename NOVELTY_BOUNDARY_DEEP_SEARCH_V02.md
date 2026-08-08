# Deep prior-art boundary v0.2

Date: 2026-08-08

## Purpose

This note records the result of an intentionally hostile literature search around the mature GeometricNeuronPlusField results.

The goal is not to manufacture novelty by adding qualifiers. The goal is to identify which statements are already owned by prior work, which model-level observations remain useful, and which technical question still appears incompletely assembled in the literature searched so far.

A search that does not find a paper is **not** evidence of novelty.

---

# I. Biological / neuronal novelty has narrowed dramatically

## Local intrinsic conductance self-organization is old

Siegel, Marder & Abbott (1994), *Activity-dependent current distributions in model neurons* (PNAS), already showed that local calcium-dependent regulation in a multicompartment neuron can create nonuniform intrinsic conductance distributions related to morphology and input pattern.

So this repository does not own:

```text
uniform/local intrinsic material
+ local activity-dependent plasticity
-> nonuniform dendritic conductance map
```

## Adjoint gradients for spatial neuronal channel density are old

Steven J. Cox (2006), *An adjoint method for channel localization*, explicitly derives an adjoint analytical gradient for a distributed ion-channel density in a neuronal inverse problem.

So this repository does not own:

```text
adjoint differentiation of spatial ion-channel density
```

## Distance compensation / dendritic democracy is established

Timofeeva, Cox, Coombes & Josić (2008), *Democratization in a passive dendritic tree*, derives distance-dependent scaling needed to normalize somatic impact of distributed dendritic inputs, including a realistic CA1 morphology.

So the large first-order coordinate found by our material learner — graph distance from the readout — is not a new general principle.

## A returned somatic signal can locally reveal distance

Sterratt, Groen, Meredith & van Ooyen (2012) showed that a back-propagating action potential can create local calcium signals predictive of synapse distance from the soma, allowing a purely local plasticity rule to establish synaptic democracy.

So the broad motif

```text
returned consequence signal
implicitly carries a distance coordinate
+
local plasticity
-> distal compensation
```

is prior art.

## Ion-channel distributions have already been inverse-designed for computation

Torben-Nielsen & Stiefel (2009), *Systematic mapping between dendritic function and structure*, optimizes realistic morphology and spatial ion-channel distributions for a chosen single-neuron computation.

So this repository does not own:

```text
optimize a spatial dendritic channel distribution for computation
```

## Modern differentiable neuronal simulators optimize huge branchwise channel fields

Jaxley / Deistler et al. (Nature Methods, 2025) differentiably fit roughly 1,390 branchwise conductance parameters in a realistic morphology and demonstrate training of biophysical neuronal parameters for computational objectives.

So this repository does not own:

```text
large-scale gradient training of branchwise channel densities
```

## Branch-specific channel compensation is already a modeled biological phenomenon

Cirtala & De Schutter (2024) model individual Purkinje branches with distinct channel conductance densities and show that branch-specific channel settings compensate differing branch morphologies and regulate branch responses/spike propagation.

So even

```text
equal/path-related branches can need different electrical material
```

is not a new biological observation.

## HCN gradients and HCN plasticity already alter phase/timing

Vaidya & Johnston (2013) demonstrated HCN-dependent compensation of location-dependent temporal differences and synchronization of rhythmic inputs at the soma.

Sinha & Narayanan (2015) showed that HCN gradients and graded HCN changes alter theta-related spike/LFP phase and phase coherence.

So this repository does not own:

```text
spatial HCN-like material can tune phase coding / synchrony
```

---

# II. What the held-out material experiments still establish

The literature subtraction does not erase the experiments. It changes what they mean.

## First-order result

`MATERIAL_LEARNING_CONFIRM_V02.md` and `MATERIAL_READOUT_RECENTER_CONFIRM_V01.md` establish in this reciprocal wave-arbor model that:

```text
uniform quasi-active material
+ fixed budget
+ exact local material sensitivity
+ phase-coordination objective
-> strong readout-distance organization
```

and moving the consequential readout causes the organization to re-center around the new readout.

That is best interpreted as a clean realization of known adjoint / distance-compensation ideas in this particular transient wave-arbor model, not as a new biological learning principle.

## Second-order result

`MATERIAL_WITHIN_SHELL_CONFIRM_V01.md` now establishes on fresh bodies that after the strongest distance-only material field is optimized and every graph-distance shell's material total is frozen exactly:

```text
mean branch-only gain over radial       +.0041769 R2    12/12 positive
mean gain over same-shell shuffles      +.0056640 R2    12/12 positive
phase-RMS improvement                   +.0076661 rad   12/12 positive
median amplitude ratio                   1.00096
```

The distance-only gain is `+.1188442 R2`, so the branch-only residual is only about 3.5% as large.

This establishes a useful model decomposition:

```text
coarse transfer geometry   graph distance from consequence
fine transfer geometry     equal-distance branch/cell field relation
```

It still does not establish novelty. Branch-specific intrinsic physiology and branch-specific conductance optimization have prior art.

---

# III. Physical-computing prior art also removes most broad claims

## Physical adjoint / in-situ backpropagation is established

Hughes, Minkov, Shi & Fan (2018), *Training of photonic neural networks through in situ backpropagation and gradient measurement*, derives physical adjoint backpropagation and shows that exact gradients can be extracted through internal intensity measurements.

Pai et al. (Science, 2023) experimentally realized in-situ backpropagation in a programmable silicon-photonic network using forward/backward propagation and optical monitoring/interference.

Therefore this repository does not own:

```text
physical backpropagation
physical adjoint training
local forward/adjoint interference as gradient measurement
```

## Physical adjoint optimization of complex multiple-scattering media is established

Guillamon, Wang, Lin & Kottos (Nature Communications, 2025), *In-situ physical adjoint computing in multiple-scattering electromagnetic environments for wave control*, experimentally optimizes complex multi-path electromagnetic environments with local forward/adjoint measurements and two physical field propagations.

Therefore this repository does not own:

```text
in-situ adjoint optimization of a complex multiple-scattering body
```

## A tunable 2-D wave mesh with in-situ physical backpropagation now exists

Thakkar & Grbic (2026), *Wave-based Neuromorphic Circuit Networks: Tunable 2D Transmission-Line Metamaterials*, proposes a 2-D grid of tunable reactive transmission-line unit cells. Computation occurs by wave interference and learned relations are stored in the tunable reactive elements. Gradients are computed from voltage measurements during two **steady-state** excitations: one forward and one adjoint.

This is an extremely close hardware neighbor to GeometricNeuronPlusField.

Its important boundary for our current question is:

```text
input encoding: single-tone sources
training fields: steady-state forward + steady-state adjoint
```

So the generic idea

```text
tunable distributed 2-D wave material
+ in-situ physical adjoint
```

is already prior art.

## Broadband/time-domain adjoints are also established

Park, Boriskina & Chung (2026), *Multi-objective time-domain adjoint via temporal convolution for band-selective electromagnetic topology optimization*, obtains band-selective broadband gradients from one broadband forward simulation and one adjoint simulation by filtering stored time-domain fields and correlating/convolving them.

Park, Miller & Chung (2026), *Nyquist-Sampled Time-Domain Adjoint FDTD for Memory-Efficient Broadband Nanophotonic Inverse Design*, shows that the forward history can be stored only at Nyquist-compliant temporal intervals and used for on-the-fly reverse-time gradient accumulation, reducing field-history memory substantially while preserving broadband gradients.

Therefore this repository does not own:

```text
time-domain adjoint optimization
broadband adjoint gradients
spectral/band decomposition of transient forward-adjoint credit
compressed storage of broadband field history
```

## Broadband physical neural computation is also prior art

Wright et al. (Nature, 2022) train physical neural networks based partly on ultrafast nonlinear optical pulse propagation. Their optical system uses roughly 100-fs broadband pulses whose frequencies mix nonlinearly.

However, their published physics-aware training uses:

```text
physical forward pass
+
differentiable digital model for backward gradients
```

rather than a physically returned adjoint field through the broadband pulse system.

So broadband/transient physical computation itself is not a seam.

---

# IV. The exact conjunction not found in this search

After searching specifically for combinations of

```text
broadband / transient
in-situ backpropagation
physical adjoint
programmable wave mesh
multiple scattering
transmission-line network
photonic pulse network
```

I did **not** find a work that clearly combines all of the following in one trainable system:

1. a **distributed programmable multiple-scattering body** whose local material/coupling parameters are the learned state;
2. computation performed on a **finite-time transient / broadband waveform**, not a single-tone steady state;
3. an adjoint/error waveform **physically propagated through the body**;
4. broadband local credit obtained **in situ from the physical forward and adjoint fields**;
5. local broadband credit accumulated with sufficiently small analog/local memory that a complete full-field digital model or stored waveform history is not required.

This absence is **not a novelty claim**. It is the narrowest technical seam left by the present search.

The closest works divide the pieces:

```text
Hughes / Pai
    physical local adjoint gradient measurement
    programmable photonics
    primarily coherent network / narrowband-style field treatment

Guillamon et al.
    physical adjoint
    complex multiple scattering
    steady-state wave control

Thakkar & Grbic
    tunable distributed 2-D wave mesh
    physical adjoint
    two steady-state passes
    single-tone input

Park et al. 2026
    broadband time-domain adjoint
    band-selective temporal convolution
    digital/numerical field histories

Park/Miller/Chung 2026
    broadband time-domain adjoint
    Nyquist-compressed history
    numerical FDTD gradient accumulation

Wright et al.
    broadband ultrafast pulse physical computation
    physical forward pass
    model-based digital backward pass
```

---

# V. What GeometricNeuronPlusField would have to demonstrate to own a useful technical contribution

The next contribution cannot simply be another numerical adjoint result.

The repository already knows numerically that the exact transient gradient can be represented as delayed local forward/adjoint correlation and compressed into a small set of spectral products.

The hardware question is now:

> **Can a programmable transient wave mesh acquire and accumulate that broadband local gradient physically with a constant or very small number of transient propagations and with bounded local analog state, instead of storing or digitally reconstructing the full field history?**

That question has an operational cost model.

For each local tunable parameter, count:

```text
physical forward propagations
physical adjoint propagations
additional interference/replay propagations
local detectors / voltage taps
local analog filters or resonators
local multiplier / mixer operations
local integrator state
stored temporal samples
ADC samples
phase references
frequency channels K
external digital operations
transpose calibration / reciprocity assumptions
```

The comparison should be made directly against the closest prior systems, not against abstract backpropagation.

---

# VI. Candidate physical gradient protocols

For one complex spectral bin with local forward field `U_k` and adjoint field `V_k`, the required real overlap can be written

```text
Re(conj(U_k) V_k)
  = ( |U_k + V_k|^2 - |U_k - V_k|^2 ) / 4
```

or, if forward and adjoint intensities are also separately available,

```text
2 Re(conj(U_k) V_k)
  = |U_k + V_k|^2 - |U_k|^2 - |V_k|^2.
```

This is the familiar interference route behind physical adjoint measurement; it is not new.

The unresolved broadband engineering question is how to obtain the **weighted sum across the transient spectrum** cheaply.

Three implementation families are worth distinguishing:

### A. Frequency-separated local lock-ins

```text
forward transient
adjoint transient / interference replay
K local resonant/lock-in channels
one accumulator per K or one weighted summer
```

Cost scales with physical K channels but avoids storing the complete waveform.

### B. Local analog time-domain correlation

Band-filter the local forward and adjoint traces with analog kernels, multiply, and integrate:

```text
(h_f * u)(t) × (h_a * v)(t)
        -> local integrator
```

This is physically analogous to the temporal-convolution broadband adjoint formula already known numerically. The research question is whether it can be realized without a digital field history, not whether the formula exists.

### C. Replay/interference protocol

If local storage is expensive, replay the forward excitation while launching the properly time-reversed adjoint so that the overlap is exposed as an intensity difference during a small number of additional physical passes.

The exact number of passes and phase/time alignment requirements must be derived and tested; vague claims of 'the medium computes its own gradient' are no longer enough.

---

# VII. Current research boundary

The biological branch has become a useful interpretation and validation framework, but not the strongest novelty direction.

The physical branch now has the more defensible open question:

> **Can transient/broadband local adjoint credit be measured inside a tunable distributed scattering medium with hardware cost comparable to steady-state in-situ adjoint training?**

A negative result would also be valuable. If broadband credit inevitably requires waveform storage, K parallel channels, or expensive phase-resolved ADC at every trainable site, then the transient GeometricNeuron machine loses much of its claimed physical-learning advantage.

A positive result would have to demonstrate the cost, not merely the mathematics.

## Wall sentence

> **The deep search removes almost every broad novelty claim: intrinsic-material self-organization, distance compensation, branch conductance heterogeneity, neuronal adjoints, physical adjoints, tunable wave meshes, broadband time-domain adjoints and even broadband physical neural computation all have prior art. The narrow seam still not collapsed is an in-situ *broadband transient* adjoint protocol for a distributed programmable scattering body that accumulates local credit physically without reconstructing or storing the full field history digitally.**
