# Operator -> flow -> event

A working synthesis after the graph-mode, modal-locality, soma-tap and AIS comparisons.

This is a hypothesis note, not a biological result.

## 1. The slow object is the operator

In the current toy the frozen arbor is not best understood as a stored image of the electrical field.  It is an operator that determines which motions the field can support.

```text
body graph G
   -> Laplacian L(G)
   -> modes {lambda_n, phi_n}
   -> transfer function available to the field
```

Morphology is therefore slow memory in the form of **constraints on dynamics**.

## 2. The fast object is the flow

The instantaneous electrical state is a moving superposition of those allowed modes:

```text
psi(x,t) = sum_n q_n(t) phi_n(x)
```

The graph-mode mechanism probe shows that the mature FunctionalArbor is almost exactly reducible to these geometry-defined oscillators.  The field need not settle.  Its job is to move through the dynamical coordinates supplied by anatomy.

## 3. A missing object is the event boundary

A passive LTI resonator bank plus a quadratic readout has a hard ceiling.  It can filter temporal structure but it does not by itself supply a localized, history-dependent, directional output event.

The AIS suggests a third object:

```text
operator / anatomy G          slow
        |
        v
flow / field psi(t)           fast continuous
        |
        v
AIS active state h(t)         local history + nonlinear transfer
        |
        v
spike event times {t_k}       sparse output code
```

The conceptual job of an AIS-like boundary would therefore not be to "understand" the whole modal state.  It would **eventize** a consequential projection of it.

That word is useful here:

> **eventization = converting a distributed continuous field trajectory into sparse, time-addressable output events.**

A spike has an origin and a time.  Those are precisely the properties a continuously distributed reciprocal field does not naturally provide as a compact output symbol.

## 4. Frequency filtering becomes central, not decorative

If the AIS is the event boundary, its frequency response determines which temporal components of the analog field are allowed to become spike timing.

This is not merely analogy.  Experiments on disrupted AIS architecture / Nav density report reduced dynamic-gain bandwidth and reduced action-potential timing precision (Lazarov et al., *Science Advances*, 2018, DOI `10.1126/sciadv.aau8621`).

Auditory neurons also tune AIS geometry with characteristic frequency: AIS length and location vary across frequency-specialized cells (reviewed in Kuba, 2012).

A later modeling literature warns that geometry/channel count alone is not enough: active channel voltage sensitivity and kinetics can dominate the high-frequency boost.  So the relevant object is not `G_AIS` alone but an active dynamical operator:

```text
A_AIS = A(position, length, channel density, channel kinetics, state)
```

## 5. The soma-local observability result changes the picture

A small soma ball is not expected to reconstruct the global mode coordinates one-by-one.  That is probably the wrong requirement.

The interface can instead be deliberately compressive:

```text
high-dimensional global mode state
              |
              v
      low-dimensional local mixture
              |
              v
       active AIS temporal filter
              |
              v
            event
```

The next empirical question in this repo is whether the soma is unusually good at the **task scalar** even when it is poor at full modal reconstruction.

## 6. This also reframes credit assignment

Claude's modal-locality audit shows that one local anatomical edit can perturb many global modes.  That makes per-cell consequence assignment badly aligned with the computational coordinates.

But global modal impact does **not** mathematically prove that local eligibility is impossible.  A local structural action can carry an event ID and later receive a scalar global consequence, exactly as policy-gradient / three-factor rules can assign credit to an action whose environmental effect is global.

The sharper distinction is:

```text
bad coordinate:
    reward every cell in proportion to recent activity/birth

better coordinate:
    remember STRUCTURAL EVENT i
    -> let global dynamics change
    -> later assign consequence to EVENT i as one causal action
```

v0.5 accidentally had this property because an entire proposed topology change was kept or reverted atomically.  v0.7-v0.9 lost the event identity while making development more realistic.

This suggests two parallel routes rather than one:

1. **AIS route:** eventize fast electrical computation into spikes.
2. **morphogenesis route:** eventize slow structural experiments so delayed global reward can address an anatomical action rather than a diffuse set of cells.

That symmetry may be important.

## 7. Three time scales

The old frozen/moving split may actually need a middle layer.

```text
slow:         morphology / scaffold / position
              G

intermediate: active channel state, inactivation, modulation
              h(t)

fast:         electrical field / modal motion
              psi(t)
```

Biological AIS plasticity spans these scales: channel gating is fast, channel modulation can be intermediate, and AIS position/length/scaffold change slowly.

A system with only `G + psi` is a geometry-defined filter bank.
A system with `G + h + psi` can become a history-dependent event generator.

## Current wall sentence

> **Geometry supplies the operator, the field supplies the moving computation, and an AIS-like active boundary may supply the missing event: a frequency-selective, history-dependent conversion from distributed analog dynamics to sparse directional spike timing.**
