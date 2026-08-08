# Reciprocity is not time reversal, source-exchange symmetry, or absence of temporal-order information

The photonic in-situ-backpropagation connection makes reciprocity central, but it also creates a tempting overstatement:

> “Reciprocity costs the arrow and buys the gradient.”

The second half is useful shorthand for this model. The first half is too broad.

## Three different symmetries

### 1. Spatial reciprocity

For the linear coupling operator,

```text
H = H^T
```

or equivalently for the Green/transfer function between two spatial points,

```text
G_ij(omega) = G_ji(omega).
```

This is what lets the same physical medium implement the transpose/adjoint spatial propagation.

### 2. Time-reversal invariance

A damped wave equation is **not** invariant under naive `t -> -t`:

```text
psi_tt + gamma psi_t + K psi = s(t)
```

becomes

```text
psi_tt - gamma psi_t + K psi = s(-t).
```

Our adjoint identity works because reversing the *adjoint recurrence/source ordering* maps the backward adjoint problem onto another causal forward run of the same damped medium. That is not the same claim as saying the original dissipative dynamics has no arrow of time.

### 3. Source-exchange symmetry

Reciprocity says the transfer from spatial point `i` to `j` equals transfer from `j` to `i`. It does **not** say two distinct source locations A and B have identical transfer histories to a third point S:

```text
h_A->S(t) need not equal h_B->S(t).
```

Therefore A-then-B and B-then-A can produce different local waveforms and different square-law interference even in a reciprocal medium.

## Our own model is already a counterexample to the universal slogan

The mature FunctionalArbor operator is reciprocal, yet the point soma square-law readout has substantial A/B temporal-order contrast. The computation comes from unequal geometry-shaped source transfer histories plus a nonlinear/quadratic observation.

What failed in the earlier skew/quadrature probes was a more specific proposal: extracting a robust order arrow from a passive local antisymmetric phase statistic. That null should not be promoted into the theorem that reciprocal media cannot encode or discriminate temporal order at all.

## The actual trade

A better engineering sentence is:

> **Reciprocity constrains directional spatial transport, but it buys direct physical reuse of the forward operator for the adjoint.**

Nonreciprocal devices do not make gradients impossible; they mean the adjoint generally requires the physical transpose/conjugate operator rather than simply sending the error field backward through an unchanged device.

## Biological consequence

Chemical synapses between neurons are directional, so an entire synaptic network cannot generally implement its exact transpose just by reversing signal propagation.

But a passive dendritic cable inside one neuron is a different physical object. Linear cable transfer impedance is reciprocal. Thus the photonic result does not prove a biological same-medium learning mechanism, but neither does it invalidate the useful comparison between a branching passive geometry and a reciprocal wave/cable operator.

The biological wall remains the active, nonlinear, stateful system:

```text
passive reciprocal cable component
+ spatially distributed active conductances
+ directional synaptic boundary inputs
+ AIS/spike eventization
+ biochemical/plasticity machinery.
```

## Wall sentence

> **Reciprocity is a spatial operator symmetry, not a synonym for “no time arrow.” In this project it is exactly what makes the physical adjoint cheap; temporal-order information can still arise from geometry-shaped histories and nonlinear observation.**
