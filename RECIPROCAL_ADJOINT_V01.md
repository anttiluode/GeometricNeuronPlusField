# Reciprocal adjoint v0.1 — the exact credit field is the same wave run backward in task time

This is the cleanest closure so far between the forward field computation and the old retrograde-credit idea.

`ANALOG_FRONTIER_LEARNING_V01.md` showed that repeated exact-adjoint updates can learn continuous frontier conductances. The remaining objection was important: an algorithmic reverse pass is not itself a physical mechanism.

For the reciprocal linear FunctionalArbor wave, that distinction collapses farther than expected.

> **The bond-credit field used by the exact discrete adjoint can be generated exactly by the same damped wave medium, started from zero and driven from the soma by the time-reversed objective-derivative waveform.**

No separate credit graph, anti-damping medium, or learned backward operator is required in this model.

## Why the identity exists

Write the mature linear wave as

```text
a = 1 - gamma dt
H = c L(k) - rho I

v[n+1]   = a v[n] + dt H psi[n] + dt s[n]
psi[n+1] = psi[n] + dt v[n+1].
```

Eliminating velocity gives

```text
psi[n+1] = (1+a) psi[n] - a psi[n-1]
           + dt^2 H psi[n] + dt^2 s[n].
```

The spatial operator `H` is symmetric because the bond medium is reciprocal.

Let `p_psi,p_v` be the exact discrete adjoint and define the combination that actually appears in the bond gradient

```text
mu[n] = dt p_psi[n] + p_v[n].
```

Its backward recurrence is

```text
mu[n] = (1+a) mu[n+1] - a mu[n+2]
        + dt^2 H mu[n+1] + dt g[n]
```

where `g[n]` is the derivative of the scalar objective with respect to the soma field at that time.

Reverse the time index. This becomes the **same forward damped recurrence** as the original wave if the new soma source is

```text
s_retro[r] = g[T-r] / dt.
```

So time reversal moves the apparent adjoint anti-causality into the source ordering. The medium itself does not have to run with negative damping.

## Development check — reused seeds 240-241

Before the fresh run, the indexing/source scaling was checked on two already-used bodies.

```text
bond-map correlation               1.000000000000
relative L2 error                   ~2e-15
frontier-score correlation          1.000000000000
one-step rho difference             <5e-17
```

That established the implementation before fresh bodies were touched.

## Fresh preregistered run — seeds 264-275

### R1 — complete bond map

Registered:

```text
mean map correlation > 0.999999
mean relative L2 error < 1e-8
```

Observed:

```text
mean bond-map correlation           0.9999999999999999
minimum body correlation            0.9999999999999993
mean relative L2 error              2.73e-15
maximum relative L2 error           3.57e-15
```

**R1 PASS.**

The physical replay and algorithmic adjoint are the same bond-sensitivity map to floating-point precision.

### R2 — frontier structural scores

Registered pooled frontier correlation `>0.999999` and max normalized error `<1e-7`.

Observed:

```text
pooled frontier correlation         1.000000000000
max normalized absolute error       4.85e-15
```

**R2 PASS.**

### R3 — actual projected learning step

Registered max updated-`rho` discrepancy `<1e-8`.

Observed:

```text
max rho difference                  3.99e-17
```

**R3 PASS.**

So all registered identity tests pass at machine precision.

## Time reversal is not cosmetic

As a preregistered descriptive ablation, the exact same soma derivative waveform was launched through the same medium **without reversing its temporal order**.

```text
mean bond-map correlation           0.84787
mean relative L2 error              0.65777
```

The correlation remains moderately high because these task waveforms have substantial slow/shared structure, but the map is no longer remotely an identity. Individual body correlations ranged from about `0.47` to `0.98`, with relative errors as high as `1.42`.

The temporal ordering of the returned waveform therefore carries real credit information.

## The local rule

The exact conductance sensitivity can now be viewed physically as

```text
forward field difference across bond e
                  ×
retrograde credit-wave difference across bond e
                  ↓
          local sensitivity dJ/dk_e
```

summed over the aligned forward and reversed-credit histories.

That is a genuine three-part structure:

```text
local forward state
× local backward/task state
× local structural degree of freedom.
```

The global objective does not need to be delivered as one scalar label to every recent event. It is encoded into the **waveform emitted by the soma**, and reciprocity lets the existing geometry spatially distribute that waveform into the correct sensitivity field.

## What this says about v0.8-v0.9

The old retrograde carrier was not wrong because "backward credit transport" was impossible. v0.8 already showed transport itself worked.

What was missing was the content of the carrier.

A scalar reward, broad activity tag, or newborn identity cannot generally tell each bond whether its change helped the soma's temporal interference. The exact returned object is richer:

> **a time-structured, complex task derivative waveform.**

The arbor then performs the spatial credit assignment by propagating that waveform back through the same reciprocal coupling geometry.

So the chain becomes

```text
forward anatomy
    -> forward task field
    -> soma interferometric readout
    -> task-derivative waveform
    -> time reverse / replay from soma
    -> same anatomy used as reciprocal credit medium
    -> forward × backward local overlap
    -> graded conductance update.
```

## A stronger symmetry

Earlier we wrote:

> Observation geometry and credit geometry are the same task boundary viewed forward and backward.

Now there is a more literal statement for this model:

> **Propagation geometry and credit geometry are the same physical operator. The distinction is the direction and waveform of use.**

The body is not merely the object being optimized. It is simultaneously:

- the forward computational medium;
- the mechanism that maps source histories into the soma interference;
- the reciprocal medium that maps soma consequence back into local structural sensitivity.

That is unusually economical.

## Relation to the interferometer / Moire intuition

The forward square-law soma readout and the backward reciprocal wave are two different uses of interference.

Forward, quadratic readout exposes pairwise phase relationships:

```text
|sum_n u_n|^2
 = sum_n |u_n|^2
 + sum_{n != m} u_n conj(u_m).
```

Those cross terms are mathematically the same kind of **difference-coordinate extraction** underlying beating and Moire-like amplification: relational phase differences become a lower-dimensional observable. The current result does not prove a literal spatial Moire pattern, and mode 0 is better described as a common/self-reference than a classical frequency-offset heterodyne local oscillator.

Backward, reciprocity reuses the same wave geometry to distribute the derivative of that interferometric consequence.

So the recurring motif is more precise than "we keep seeing waves":

> **Geometry creates relational interference in the forward direction, and the same geometry converts relational consequence into local sensitivity in the reverse direction.**

## What remains physically nontrivial

This identity does **not** mean a biological neuron automatically performs exact adjoint learning.

The exact construction requires the soma to provide the appropriate objective-derivative waveform in reverse temporal order. That implies some mechanism for storing, delaying, reversing, or otherwise generating the conjugate task history.

That is now the real physical wall.

## Next clean experiment

Destroy the exact reversal in controlled ways and ask how much learning survives:

```text
exact reversal
coarse temporal bins
low-pass reversed waveform
phase-only / envelope-only return
finite delay window
noisy reversal
non-reversed return
scalar reward pulse
```

Compare each returned waveform first against the exact adjoint bond map, then in the analog frontier learner.

If a coarse/local approximation preserves most of the learning advantage, the mechanism begins to look physically plausible rather than merely mathematically implementable.
