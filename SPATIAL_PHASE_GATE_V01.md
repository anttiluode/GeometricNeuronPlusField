# Spatial dendritic phase gate v0.1 — smooth morphology-indexed phase is null in the first direct mapping

## Result

The first direct attempt to turn the in-vivo dendritic theta phase gradient into a useful geometry-indexed plasticity gate is a **null**.

The biological motivation remains real: Liao et al. (Nature Communications, 2024) reported that intracellular theta phase varies systematically across the basal-to-tuft axis of single CA1 pyramidal neurons, with a pooled phase gradient of about `-7.9 degrees / 100 micrometers`, and described the subthreshold theta pattern as a traveling wave across the dendritic tree.

But simply assigning a smooth phase offset by graph distance and using it to gate the already-computed local forward x return structural overlap does **not** help this model.

Development bodies: seeds **528-533**.

No held-out positive confirmation is justified.

## What was tested

The exact reciprocal local structural density was first computed:

```text
c_e(t) = 2 dt K Re[conj(dmu_e(t)) dpsi_e(t)]
```

for each real bond `e`.

The exact bond gradient is

```text
g_e = sum_t c_e(t).
```

Only after that computation, a 50%-duty local gate was applied:

```text
g_e(gated) = 2 sum_t m_e(t) c_e(t).
```

The factor 2 compensates the common 50% duty.

The gate period was `P=42`, chosen because the global return-gating experiments had already shown that this temporal scale is phase-sensitive.

Four spatial phase organizations were compared:

```text
global
    same phase at every bond

smooth
    phase_e = phi0 + alpha * graph_distance_from_soma

shuffled
    exact same set of smooth phase values,
    randomly permuted across real bonds

random
    independent local phases
```

For smooth gates,

```text
alpha = .25, .5, 1, 2, 3, 4
```

simulation frames of phase shift per graph edge.

Every global offset `phi0` was scanned.

## Development result

The global gate itself gave

```text
mean map correlation      0.64703588
mean phase-min correlation 0.56052090
```

Adding a smooth spatial phase gradient changed essentially nothing.

Representative group means:

```text
alpha=.25
smooth       0.64703922
shuffled     0.64704041
random       0.64703131
smooth-shuffled   -0.00000119
smooth-global     +0.00000334

alpha=1.0
smooth       0.64703316
shuffled     0.64703421
random       0.64703470

alpha=2.0
smooth       0.64703334
shuffled     0.64703250
random       0.64703665

alpha=3.0
smooth       0.64702544
shuffled     0.64703531
random       0.64703308
```

The differences are at roughly `1e-6` to `1e-5` scale.

There is no evidence that smooth graph-distance organization is privileged.

## What this kills

This particular mapping is not supported:

```text
real dendritic theta has a spatial phase gradient
        therefore
smooth morphology-indexed theta phase directly improves
local forward x return plasticity gating
```

It does not, in this implementation.

The smooth phase field is effectively interchangeable with shuffled or independent phase at the level of the measured gradient-map fidelity.

## Why the null is informative

There are at least three reasons this simple mapping may have been too downstream.

### 1. The geometry had already done the computation

The gate was applied only after the exact forward and adjoint fields had propagated and after their local overlap density `c_e(t)` had been formed.

So the spatial phase field could only decide **when to keep pieces of an already-computed local answer**.

It did not alter:

- source arrival phase;
- dendritic transfer;
- mode excitation;
- local membrane impedance;
- inhibition/excitation balance;
- or the return propagation itself.

A real intracellular theta phase gradient may affect some of those upstream quantities.

### 2. Graph distance is not the biological coordinate

The Liao et al. result is an empirical basal-to-tuft phase organization in a laminar hippocampal circuit. It should not be reduced automatically to

```text
phase proportional to path length from soma.
```

Layer-specific excitation and inhibition, active conductances, and morphology all contribute to intracellular theta.

### 3. Averaging across every global phase is hostile to a phase-specific role

The present analysis scans every global offset and summarizes the resulting map correlations.

If the biological function occurs at a specific theta phase — for example, a preferred plasticity or excitability window — averaging the entire cycle can deliberately wash out that function.

That is not a reason to rescue this result post hoc. It is a reason that a future test must freeze a different, independently motivated observable before running new bodies.

## Relation to the confirmed global-gating result

The global return-gating branch is much stronger.

On held-out bodies, equal-duty masks were emphatically not equivalent:

```text
fast regular interleaving     ~ exact structural direction
random half sampling          worse
one contiguous half-window    much worse
slow periodic gating          progressively worse
```

That result concerns **temporal spectral organization** and survived hostile controls.

The present spatial-phase test says that simply painting a smooth phase gradient over the local plasticity accumulator does not add a geometric benefit.

These results are compatible:

```text
temporal regularity matters
spatial phase ordering, in this downstream mapping, does not
```

## Correct next move

Do not tune `alpha`, period, duty cycle, or coordinate system until a positive mechanism is specified independently.

The next load-bearing experiment remains the frequency-domain explanation of the confirmed periodic return-gating effect:

```text
return code g(t)
        x
periodic mask m(t)
        |
        v
spectral replicas / sidebands
        |
        v
return-to-gradient sensitivity
        |
        v
gradient-map error
```

If that mechanism is established, a spatial phase field can later be reintroduced only where it changes an identified spectral or propagation interaction.

## Biological claim boundary

The null does **not** dispute the in-vivo spatial theta result.

It says only:

> **A real morphology-indexed theta phase field exists in CA1 pyramidal neurons, but the first naive way of using such a field in GeometricNeuronPlusField — as a smooth graph-distance phase mask on an already-formed local plasticity signal — is functionally indistinguishable from shuffled phase.**

## Wall sentence

> **The dendritic theta traveling wave is real; our first attempt to make its spatial phase itself a plasticity coordinate was not. Temporal interleaving survives. The simple spatial-phase mask does not.**
