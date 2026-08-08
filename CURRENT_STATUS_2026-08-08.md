# Current status — 2026-08-08

`GeometricNeuronPlusField` has moved substantially beyond the original soma/AIS question.

The shortest current description is:

> **A frozen local geometry defines a broadband moving wave field; a simple local readout converts relations in that field into consequence; the transpose dynamics of the same geometry can route structural consequence back into local couplings; and the return channel can be strongly compressed and temporally multiplexed without losing its useful structural direction.**

The general physical-adjoint principle is established prior art. The active work here is the transient/broadband, locally measurable, geometry-changing version and its possible relation to dendritic electrical organization.

---

## 1. Upstream geometry result

The grown/frozen arbor is well described by its graph-Laplacian modal coordinates.

The geometry determines:

```text
mode eigenvalues / time scales
source-to-mode couplings
mode amplitudes and phases
```

and therefore the possible moving electrical field.

Temporal-order information is not distributed according to electrical energy. The common spatial mode can carry enormous field energy and almost no order information, while much weaker higher modes carry useful order selectivity.

A reduced geometry-defined oscillator bank predicts the full-field modal result at roughly `r=.995` without fitting the answer.

So the earned version of "geometry computes" is:

```text
geometry defines the dynamical coordinates in which the field can compute.
```

---

## 2. Soma/readout result

The soma/root in this constructed arbor is not a magical privileged cell.

Its usefulness is largely explained by being a convergence point where the two source histories become amplitude-balanced enough to interact strongly.

The local square-law readout

```text
|h_A + h_B|^2
```

contains the source-source cross-term

```text
2 Re[h_A conj(h_B)]
```

which carries most of the temporal-order computation.

Mode-pair decomposition showed that the soma is fundamentally a **mode mixer**: the order-sensitive cross-term is sparse over mode pairs and overwhelmingly dominated by cross-mode rather than same-mode interactions.

---

## 3. Generic AIS/HH bridge: consequential, not privileged

A compact HH-like active boundary was tested downstream of the soma field.

It changes firing regime strongly, but under matched controls it did not earn a positive computational advantage:

```text
better spike-time precision                     NO -- worse
clean kinetic passband                          NO
simple n-set refractory clock                   NO
special benefit from carrier phase at eventizer NO -- worse
more frequency information per matched spike    NO -- much worse
```

The AIS branch therefore stops rather than adding position/length/channel-density knobs around negative results.

This does not kill the real AIS. It kills this specific generic-HH rescue story in this model.

---

## 4. Structural credit: the exact local signal is forward field x transpose field

The major turn came when the bond-gradient of the transient task was written explicitly.

For the reciprocal linear model, the exact structural derivative at a bond is proportional to

```text
sum_t Re[(forward field difference)
         * conj(returned adjoint-field difference)]
```

The explicit adjoint map and the physical forward/return overlap match to machine precision.

The same graded geometry is therefore both

```text
the forward computing operator
and
the spatial router of the returned structural derivative.
```

This independently rediscovered a member of the established in-situ/physical-adjoint family used in photonic and other wave systems. Do not claim novelty for physical adjoint backpropagation itself.

---

## 5. What reciprocity actually buys

A genuinely nonsymmetric operator removed an important ambiguity.

Time-reversing the soma derivative and replaying it through the same nonreciprocal operator progressively stops giving the adjoint:

```text
nonreciprocity beta     same-H gradient corr
0                        1.0000
.02                       .9972
.05                       .9802
.10                       .9037
.20                       .7341
.30                       .5389
.40                       .4535
.60                       .2645
```

Using the actual transpose operator restores the exact adjoint at every tested beta to floating-point precision.

The correct wall sentence is:

> **Time reversal is not the magic; transpose propagation is. Reciprocity is valuable because it makes the transpose free.**

---

## 6. Broadband transient gradient compression

The exact local transient correlation uses 210 time samples in the registered task, but its spectral gradient mass is highly concentrated.

Held-out work established:

```text
about 13 / 210 bins carry 95% of absolute gradient mass

K=8 boundary-selected bins
    retain most gradient direction / most learning gain

K=16
    nearly preserve the exact finite-step learner
```

The key point is that a useful common K-bin set can be chosen from **boundary/source-return signals**, not from omniscient inspection of all internal bond gradients.

---

## 7. Local intensity readout: 2K phase steps are not required

For one retained bin,

```text
(|U+V|^2 - |U-V|^2)/4
    = Re[conj(U)V].
```

But the K retained components can be replayed together as one band-limited transient waveform.

Then Parseval gives

```text
1/4 [sum |u+v|^2 - sum |u-v|^2]
    = sum_k Re[conj(U_k)V_k].
```

So the complete compressed cross-term can be recovered with **two global phase states**, independent of K, plus a local full-window intensity integral.

On fresh bodies the broadband two-state identity reproduced the K-bin complex-product map at ~`1e-15` relative numerical precision.

A fixed `0.1 rad` phase-setting error still left map correlation above `.99999` in the registered model.

This does not make the physical implementation cost independent of K: the device still has to retain/form the selected band, align the transient fields, and integrate locally.

---

## 8. One-run rhythmic lock-in readout

A stronger physical protocol uses one balanced binary reference

```text
s(t) in {+1,-1}
I(t) = |u(t) + s(t)v(t)|^2.
```

Then

```text
1/2 sum_t s(t) I(t)
 = sum_t Re[conj(u)v]
 + 1/2 sum_t s(t)(|u|^2+|v|^2).
```

The first term is the desired local gradient.

The second is the entire self-energy leakage.

For retained K-bin waveforms, the self-energy spectrum lies on pairwise difference frequencies. Therefore the lock-in measurement is exact when the reference harmonics avoid that difference-frequency set.

Held-out confirmation:

```text
collision-free points                80
max leakage / gradient L2            5.8e-15

collision-positive points            88
median leakage / gradient L2         1.199

weighted spectral overlap
vs measured leakage                  r=.970
```

A slower rhythm can outperform a faster one if the slower reference happens not to collide with the relevant difference frequencies.

So:

> **speed is not the primitive; spectral separation is.**

---

## 9. Coarse return codes: much of the exact waveform is redundant

The exact soma derivative was damaged toward biologically coarser return codes and tested on fresh bodies.

Formal result: `7/8` preregistered criteria passed; the failed envelope-only criterion remains a failure.

Strong surviving results:

```text
real-valued return Re[g(t)]
    map corr                         .9999986

50%-duty fast periodic gate P14
    median-phase map corr            .998785
    worst-phase mean corr            .997623

32 sparse phase-bearing events
    mean corr                         .8642

32 sparse events retaining only
objective/target-distractor sign
    mean corr                         .8501

same event times forced all-positive
    mean corr                        -.0872
```

So the exact complex return contains more information than the structural direction needs.

Explicit quadrature can be removed almost completely in this model. A great deal of temporal support can be removed if it is removed in the right pattern.

But the **sign of consequence is load-bearing** in the sparse limit.

The envelope-plus-sign arm was often informative but narrowly missed its frozen mean-correlation criterion (`.7952` versus `.8000`).

---

## 10. Equal duty is not equal information

A hostile control compared equal 50%-duty masks.

Fresh held-out bodies:

```text
periodic P2       .999969 mean gradient-map corr
P6                .999886
P10               .999708
P14               .998405
P30               .972457
P42               .851114
P70               .720260

random half       .980956
contiguous half   .652009
```

All masks keep exactly half the return samples and are dose matched.

Therefore the result is not

```text
half the waveform is enough.
```

It is

```text
fast regular temporal interleaving preserves the useful return relation;
slow, random and especially contiguous removal contaminate it increasingly.
```

The hostile control passed `6/6` held-out criteria.

---

## 11. Boundary spectra predict the gate-rate regime

The confirmed mechanism is now measurable before inspecting the gated internal map.

For a 50% mask

```text
m(t)=1/2+r(t)
```

the gated return spectrum is

```text
1/2 G + FFT[r g].
```

The non-DC mask term creates sidebands.

Measure how much of that contamination lands in the K=8/K=16 **boundary-selected important bins**.

Fresh held-out result:

```text
periodic regime damage vs contamination
K8       r=.99856
K16      r=.98181

all periodic + random + block classes
K8       r=.93345
K16      r=.94964
```

P42 produces about `42x` more K8 contamination than P6 and loses `.153` in gradient-map correlation.

A contiguous block is about `5-6x` spectrally dirtier than random half-sampling and loses another `.235` of correlation.

Individual mask realizations are only moderately predicted (`r~.46-.47`), so the score is a **regime/bandwidth predictor**, not a complete map model.

Current wall:

> **The rhythm is not the message. Its useful job can be spectral housekeeping: temporally interleave a return/consequence channel so its modulation sidebands miss the compact boundary spectrum carrying the structural information.**

---

## 12. Biology: what was killed and what became more interesting

### Chandelier = pi phase inverter: rejected

Axo-axonic/chandelier inhibition at the AIS should not be described as a generic `pi` phase inverter. The direct `+V/-V` biological mapping was too strong.

AIS-targeting inhibition is better treated as control over spike initiation, timing, veto and excitability.

### Passive single-arbor reciprocity: real

The earlier objection "biology is nonreciprocal because synapses are one-way" was too broad for this single-neuron geometry problem.

In passive cable theory, dendritic transfer impedance is reciprocal:

```text
Z(i,j,omega) = Z(j,i,omega).
```

So the passive/subthreshold part of one dendritic tree genuinely has the symmetry that makes transpose-like reverse propagation worth asking about.

Real dendrites add active, nonuniform, nonlinear conductances, so an exact biological adjoint is not implied.

### The closest biological return event is the back-propagating consequence

The neuron already sends soma/axon-generated events back into its dendrites. Back-propagating action potentials interact locally with recent dendritic/synaptic depolarization and are involved in plasticity.

That is topologically much closer to

```text
forward local history x returning consequence event
```

than chandelier inhibition is.

### Real theta is spatial inside one cell

In-vivo voltage imaging of CA1 pyramidal cells in 2024 found a systematic intracellular theta phase gradient across dendritic morphology, roughly `-7.9 degrees / 100 micrometers` pooled along the basal-to-tuft axis.

So a real pyramidal neuron does not sit under one spatially uniform theta clock.

However, the first direct model mapping — applying a smooth graph-distance phase mask to the already-formed local plasticity overlap — was a **null**. Smooth phase, shuffled phase and random phase were functionally indistinguishable.

That simple spatial-phase-mask idea stops here.

### HCN is a much better geometric bridge

Vaidya & Johnston (2013) showed that CA1 pyramidal neurons use a **spatial gradient of HCN-mediated inductive membrane properties** to compensate location-dependent input timing and produce temporal synchrony at the soma over theta/gamma ranges.

Critically, local dendritic responses can remain location/phase dependent while the **somatic readout becomes synchronized**.

They also found that gamma-frequency synaptic bursts generate lower theta-frequency components relevant to that synchrony.

That gives a real biological example of

```text
nonuniform / phase-rich dendritic field
        +
spatially graded electrical material
        ->
coordinated consequential readout at soma.
```

This is much closer to the mature Geometric Neuron idea than a dendrite-as-delay-wire or chandelier-as-phase-switch story.

---

## 13. Current biological working picture

Not:

```text
brain runs exact backpropagation
```

and not:

```text
gamma carries the gradient
```

A narrower testable architecture is:

```text
morphology
  + spatially varying electrical material
        |
        v
moving dendritic field / input history
        |
        v
soma/AIS consequential event
        |
        v
returning dendritic event / waveform
        |
        +---- regular temporal multiplexing / gating
        |
        +---- signed consequence / third factor
        v
local coincidence with recent forward state
        |
        v
plastic change
```

The exact transpose remains the engineering positive control.

The biological question is how much useful structural direction survives when the reciprocal passive substrate is replaced by active dendrite and the exact analog return is replaced by real neuronal events.

---

## 14. What to do next

### Engineering next

The temporal gate mechanism is now sufficiently established in the current linear model. Do not keep sweeping gate periods.

The next engineering benchmark should measure whether this compressed/local physical protocol buys anything under a **real implementation cost model**: detector count, local memory, phase/control precision, tuning elements, passes, energy and bandwidth.

### Biology next

Do not add another downstream plasticity phase mask.

Move the spatial biology **upstream into the electrical operator**.

The strongest next candidate is the HCN result:

```text
uniform electrical material
vs
smooth morphology-indexed reactive/impedance gradient
vs
shuffled same material values
```

Ask whether a minimal linearized spatial impedance gradient can compensate geometry-dependent timing at the soma without freezing the local moving field.

Only after that works or fails should a full conductance-based HCN model be justified.

---

## Current one-sentence picture

> **The Geometric Neuron is no longer 'a dendrite whose shape delays signals.' It is becoming a spatially structured electrical material: morphology defines the field's modal possibilities, local material properties reshape their timing, a convergence geometry makes relations consequential, and regular temporal multiplexing can return coarse structural consequence through the same medium without requiring the distributed field itself to settle.**
