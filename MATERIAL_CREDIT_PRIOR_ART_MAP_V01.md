# Material-credit prior-art map v0.1

Date: 2026-08-08

## Why this note exists

The material-learning branch produced an initially attractive story:

```text
uniform intrinsic electrical material
+ local forward × transpose sensitivity
+ fixed irregular arbor
+ somatic phase-coherence objective
-> a distal / readout-distance material organization emerges
```

That story is real in the present model, but most of its broad ingredients have substantial prior art.

This note deliberately subtracts those ingredients before asking what, if anything, remains distinctive.

---

# Prior-art fence

## 1. Local intrinsic plasticity can self-organize nonuniform channel maps

**Siegel, Marder & Abbott (1994), “Activity-dependent current distributions in model neurons,” PNAS.**

A multicompartment neuron with locally calcium-regulated channel densities develops nonuniform conductance distributions related to dendritic morphology and input pattern. The point is already very close to the broad statement that distributed intrinsic membrane properties can self-organize from local activity.

Therefore this repository does **not** own:

```text
local activity-dependent intrinsic plasticity
-> nonuniform dendritic conductance distribution
```

---

## 2. Adjoint gradients have already been applied to distributed neuronal channel density

**Steven J. Cox (2006), “An adjoint method for channel localization.”**

This work formulates an inverse problem for recovering a nonuniform ion-channel distribution from voltage observations and develops an adjoint analytical gradient for that distributed channel-density parameter.

Therefore this repository does **not** own:

```text
adjoint differentiation of a spatial neuronal channel field
```

The important difference is problem formulation: Cox addresses inverse localization from observations, whereas the present branch reallocates a fixed material budget to improve a computational transfer objective. That is a distinction, not yet a novelty claim.

---

## 3. Distance-to-soma compensation / dendritic democracy is established theory

**Timofeeva, Cox, Coombes & Josić (2008), “Democratization in a passive dendritic tree: an analytical investigation.”**

This work derives distance-dependent synaptic scalings that compensate dendritic attenuation so spatially distributed inputs have normalized impact at the soma, including realistic CA1 morphology.

Therefore the dominant coordinate found by our learner—distance from consequence/readout—cannot by itself be presented as a new principle.

---

## 4. A returning somatic signal can locally reveal distance and establish a distal gradient

**Sterratt, Groen, Meredith & van Ooyen (2012), “Spine calcium transients induced by synaptically-evoked action potentials can predict synapse location and establish synaptic democracy.”**

A back-propagating action potential creates a local calcium signal whose amplitude is informative about distance from the soma. A purely local homeostatic rule can then make distal synapses stronger and reduce location dependence of the somatic EPSP.

Therefore this repository does **not** own the broad idea:

```text
a return signal from consequence
implicitly carries a distance coordinate
and local plasticity can use it to create a distal gradient
```

Our returned field is mathematically different and richer, but the broad credit-coordinate motif is old.

---

## 5. Ion-channel distributions have already been inverse-designed for single-neuron computation

**Torben-Nielsen & Stiefel (2009), “Systematic mapping between dendritic function and structure.”**

The inverse approach jointly explores realistic morphology and spatial ion-channel distributions for a chosen neuronal computation, including structured channel gradients/hotspots.

Therefore this repository does **not** own:

```text
optimize dendritic ion-channel distribution for a computational objective
```

---

## 6. Modern differentiable simulators already train large channel-density parameter sets

**Deistler et al. / Jaxley (Nature Methods, 2025), “Differentiable simulation enables large-scale training of detailed biophysical models of neural dynamics.”**

Automatic differentiation is used to optimize large numbers of branchwise ion-channel conductances in realistic morphologies and to train biophysical neurons on computational objectives.

Therefore this repository does **not** own:

```text
backpropagation / automatic differentiation through a realistic neuron
or
high-dimensional gradient training of branchwise ion-channel densities
```

---

## 7. Local activity × returned/error signal is also an established neural-learning form

Active-dendrite / error-backpropagation models already contain three-factor-like local rules in which local dendritic or presynaptic activity is multiplied by a somatic error/backpropagating signal.

Therefore even the abstract bilinear form

```text
local forward activity × returned consequence signal
```

is not itself a novelty claim.

---

## 8. HCN spatial localization is biologically plastic

CA1 HCN1 distal enrichment is activity dependent and reversible; entorhinal/glutamatergic activity and CaMKII-dependent signaling have been implicated in maintaining the distribution.

This supports the biological plausibility of *spatially regulated electrical material*, but it does not establish that biological HCN trafficking implements our exact adjoint rule or objective.

---

# What the repository has actually added so far

The broad distal-gradient story is now heavily prior-art constrained.

The remaining empirical question is narrower:

> **After the best readout-distance-only solution is allowed to optimize freely, does full local forward×transpose material credit still discover reproducible branch/cell-specific structure that improves the consequential readout?**

This is where the present model can still say something concrete.

---

# Development diagnostic 1 — scalar returned-field cues are insufficient

`material_credit_decomposition.py` was run on development bodies 622–627.

The exact learned material histogram was held fixed and the same density values were reassigned according to scalar coordinates measured in the uniform state.

Mean phase-coherence `R^2`:

```text
exact full placement       .60023
explicit distance rank     .58406
return-delay rank          .25607
return-amplitude rank      .44920
forward-amplitude rank     .31393
random shuffle             .45395
uniform                    .49165
```

The exact map beat:

```text
distance rank       5 / 6 bodies
return delay        6 / 6
return amplitude    6 / 6
forward amplitude   6 / 6
shuffle             6 / 6
```

So a simple scalar returned-field delay or attenuation does not reproduce the learned placement in this model.

However, explicit graph distance already captures most of the useful placement.

Do **not** overinterpret this as proving that no biologically plausible return-only signal could work. The transpose field is complex and frequency structured; these controls test only simple scalar reductions.

---

# Development diagnostic 2 — optimize the distance-only solution, not a hand gradient

`material_radial_vs_full.py` gives the distance-only control every advantage:

```text
same frozen arbor
same quasi-active material physics
same .03/.04 objective
same uniform initialization
same total material budget
same per-cell cap
exact analytic gradient
```

but all cells at the same integer graph distance from the readout must share one density.

The radial learner therefore has one optimized degree of freedom per distance shell rather than receiving a hand-drawn linear gradient.

Development bodies 622–627:

```text
mean coherence R^2

uniform                 .491650
hand linear gradient    .545380
optimized radial        .587817
full cell-by-cell       .600233

radial - uniform       +.096167
full - uniform         +.108582
full - radial          +.012416
full beats radial       6 / 6
```

The optimized radial solution captures approximately **88.6%** of the full learner's improvement over uniform material.

So the correct decomposition at this stage is roughly:

```text
large effect:
    readout-distance compensation
    already close in spirit to dendritic-democracy literature

smaller effect:
    branch/cell-specific correction beyond distance
    +.0124 R^2 in this development set
    full wins 6/6
```

The branch-specific residual is **development only** until a frozen held-out test is run.

---

# Current novelty boundary

Do not currently claim novelty for:

```text
self-organized intrinsic conductance gradients
adjoint channel-density gradients
somatic distance compensation
return signals that encode soma distance
optimized dendritic channel maps
AD-trained conductance distributions
local activity × returned error
physical adjoint training in reciprocal media
```

A defensible present description is narrower:

> **GeometricNeuronPlusField asks whether an exact local forward×transpose sensitivity acting on distributed intrinsic electrical material produces useful spatial structure beyond the best readout-distance-only compensation on a fixed irregular wave arbor, and whether that residual can ultimately be implemented as a local physical measurement rather than a digitally assembled gradient.**

I have not found, in the current search, a single prior work that exactly combines all of those conditions. That is not yet evidence of novelty.

---

# Next wall

Freeze the full-vs-optimized-radial comparison on untouched bodies.

The branch-specific seam survives only if:

1. optimized radial material again explains most of the improvement over uniform;
2. the full local learner nevertheless beats optimized radial material with a reproducible positive margin;
3. that margin is not purchased by pathological attenuation;
4. the result remains present at both frozen low frequencies.

If the full learner does not beat the optimized radial learner on fresh bodies, then the strongest material result should be described primarily as a new implementation/testbed of established distance-compensation ideas rather than as a distinct geometric mechanism.

## Wall sentence

> **The emergent distal gradient is not the seam: prior work already owns local intrinsic self-organization, dendritic democracy, return-derived distance cues, channel optimization, and neuronal adjoints. The only material-learning seam still standing is the reproducible branch-specific correction beyond an optimally learned readout-distance profile—and that now has to earn itself on fresh arbors.**
