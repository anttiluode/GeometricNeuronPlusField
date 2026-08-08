# HCN synchrony bridge v0.1 — real dendrites already compensate geometry with graded impedance

The theta/gamma branch has a much closer biological neighbor than the earlier chandelier-cell phase-inverter speculation.

Sachin P. Vaidya and Daniel Johnston,
**“Temporal synchrony and gamma-to-theta power conversion in the dendrites of CA1 pyramidal neurons,”**
*Nature Neuroscience* 16, 1812–1820 (2013), DOI `10.1038/nn.3562`.

Their central result is directly geometric:

> hippocampal pyramidal neurons use a spatial gradient of HCN-channel-mediated inductive membrane properties to counteract location-dependent temporal differences of dendritic inputs at the soma.

Using simultaneous multisite whole-cell recordings and modeling, they found that this active impedance gradient produces temporal synchrony of rhythmic inputs over wide dendritic regions, especially in theta/gamma ranges.

A particularly important detail for this repository is that the synchrony is a **somatic transfer/readout effect**, not the dendritic voltage field becoming spatially uniform. Their supplementary analysis explicitly reports theta-frequency oscillatory synchrony at the soma while local dendritic responses retain location-dependent phase differences.

Their experiments also report that gamma-frequency synaptic bursts generate lower theta-frequency components that contribute to this oscillatory synchrony.

This is established neuroscience. It is not a result of this repository.

## Why it matters here

Our early FunctionalArbor story was

```text
longer path -> later arrival
```

and the later graph-mode work replaced that with the richer statement

```text
geometry -> distributed poles / modal time scales -> local interference
```

The HCN result says a real CA1 neuron does not passively accept whatever timing its morphology creates.

Its membrane properties vary systematically across that morphology and reshape transfer timing.

A useful abstraction is therefore not

```text
fixed geometry + uniform cable
```

but

```text
geometry G(x)
    +
spatial impedance field Z(x,omega)
    ->
frequency-dependent transfer relation
```

This is extremely close to the mature `GeometricNeuronPlusField` framing: **mass/shape and electrical dynamics co-define the coordinates in which temporal relations arrive at the soma.**

It also gives a biological example of the distinction that kept appearing in the earlier Horizon/Clockfield discussion:

```text
local distributed field remains phase-rich / nonuniform
while
one consequential projection can become coordinated
```

The field does not have to freeze for the readout relation to become useful.

## A subtle connection to the new gating result

The repository has just established on held-out bodies that fast regular temporal interleaving of a return waveform preserves structural direction far better than random or slow equal-duty masks.

The current engineering hypothesis is spectral:

```text
fast regular gate
    -> desired low/baseband copy
       + separated sidebands

slow/random gate
    -> contamination overlaps consequential frequencies
```

Vaidya & Johnston found a biological mechanism with a related but not identical job:

```text
spatially distributed inputs
    -> geometry-dependent delays
    -> graded HCN impedance
    -> frequency-selective phase compensation
    -> synchrony at soma
```

and gamma bursts can contribute theta-frequency components to the synchronized somatic waveform.

The common object is therefore **frequency-dependent compensation of geometry-created timing differences**, not “gamma computes the gradient.”

## This changes the next biological model

The first spatial-theta plasticity-mask experiment in this repo was null. That test only multiplied an already-computed local structural overlap by a phase gate.

The HCN literature suggests that was probably the wrong level to insert spatial biology.

A better intervention is upstream, inside the propagation operator itself.

Instead of a uniform local restoring/damping term, introduce a smooth morphology-indexed reactive term:

```text
q'' + gamma(x) q' + rho(x) q + K L_G q = input
```

or, in a frequency-domain reduced model,

```text
H(omega) = -omega^2 I + i omega Gamma(x) + Rho(x) + K L_G.
```

Then compare:

```text
A. uniform membrane dynamics
B. smooth soma-to-distal impedance gradient
C. shuffled same impedance values
D. independently optimized unconstrained impedance values
```

Primary questions:

1. Can a smooth impedance gradient reduce source-location-dependent phase spread at the soma?
2. Does it improve temporal-order discrimination or merely erase useful delay?
3. Does the best gradient differ for theta-like versus gamma-like driving bands in simulation units?
4. Can high-frequency burst input generate a lower-frequency consequential component through a minimal nonlinearity without destroying the upstream modal relation?
5. Does a smooth gradient outperform a shuffled gradient with identical parameter histogram?

The fifth control is essential. If shuffled impedance performs equally well, morphology-indexed organization has not earned a role.

## Important difference from the paper

HCN channels are active, voltage-dependent conductances. The present FunctionalArbor wave equation is a much simpler linear complex field.

Adding a hand-written frequency filter and calling it HCN would be cosmetic.

The first test should therefore use only the **minimal linearized impedance consequence** of a spatial gradient and ask whether the geometry-dependent synchrony effect appears at all.

Only if that succeeds should a conductance-based HCN current be introduced.

## Relation to the 2024 dendritic theta traveling-wave result

The 2024 in-vivo voltage-imaging result and the 2013 HCN synchrony result are not contradictory.

They address different observables and conditions:

- Vaidya & Johnston: frequency-dependent transfer from distributed dendritic inputs to a common somatic readout can be synchronized by an HCN gradient while local dendritic responses remain phase-dependent.
- Liao et al.: intracellular theta membrane-potential phase itself is spatially organized across the dendritic arbor in vivo.

Together they warn against treating a pyramidal neuron as either a uniform cable or a globally phase-locked lump.

The arbor contains both spatially varying transfer properties and spatially varying oscillatory state.

## Wall sentence

> **A real pyramidal dendrite does something our early delay-line model did not: it spatially grades its own impedance to compensate the timing consequences of its geometry. The local field can remain phase-rich while the soma sees a synchronized relation. The biologically interesting Geometric Neuron is therefore not shape plus signal; it is shape plus a morphology-indexed electrical material.**
