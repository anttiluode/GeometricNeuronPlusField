# Modal mechanism probe v0.1

Canonical GitHub Actions run: `31240432612`

After the graph-mode band survived held-out bodies, the next question was why those modes see temporal order at all.

The probe used all graph modes of 12 frozen FunctionalArbor bodies and compared the measured full-field modal contrast with a parameter-free reduced model:

```text
q_n'' + damping q_n' + (restoring + stiffness * K * lambda_n) q_n
    = phi_n(A) s_A(t) + phi_n(B) s_B(t)
```

No parameter was fitted to the modal contrasts. `lambda_n`, `phi_n(A)`, and `phi_n(B)` came from each body's graph; `dt`, damping, restoring, stiffness, K, pulse envelope, carrier phase and source amplitude came directly from FunctionalArbor.

## Result

The reduction is almost exact.

```text
12 bodies
all non-constant modes

predicted vs measured signed contrast
    mean per-body r       0.99484
    median per-body r     0.999994
    pooled r              0.99518

predicted vs measured |contrast|
    mean per-body r       0.99307
    median per-body r     0.999992
    pooled r              0.99361

top-10 selective-mode overlap
    mean                  9.67 / 10
```

Two bodies are somewhat less exact (`r ~ .99` and `.95`); most are visually almost on the identity line. The small departures are where the full lattice, weak mature bath coupling, pointwise saturation, and the reduced isolated-body linear model cease to be identical.

## Important deflation

This near-perfect result is **not** a mysterious empirical miracle.

The mature FunctionalArbor is almost exactly a fixed-K wave equation on the body graph. The graph Laplacian eigenvectors are therefore nearly the dynamical normal modes of the thing we built. Diagonalizing that operator is mathematically expected to turn one distributed wave equation into many scalar oscillators.

That is precisely why this result matters conceptually, but it should not be sold as an unexpected prediction.

The useful statement is:

> **The computation that looked like a travelling electrical field over a grown dendritic body can be reduced, almost without loss and without fitted parameters, to geometry-defined modal oscillators driven according to where the input terminals sit in those modes.**

## What does NOT explain selectivity by itself

Simple static geometry scores do poorly across the complete mode set:

```text
pooled corr(|C|, lambda)                  -0.038
pooled corr(|C|, |phi(A)-phi(B)|)         -0.012
pooled corr(|C|, |phi(A) phi(B)|)         +0.176
pooled corr(|C|, A-vs-B path alignment)   -0.030
```

So there is no useful story such as "higher modes are better" or "the most A-vs-B-shaped mode wins."

Temporal-order selectivity comes from the **combination**:

```text
mode geometry / eigenvalue
        x
coupling of A into that mode
        x
coupling of B into that mode
        x
pulse lag and mode's temporal response
```

That is why the full reduced transfer function succeeds while any one static scalar mostly fails.

## Where the Geometric Neuron picture has arrived

We started with a vague statement that dendritic/soma geometry might compute.

The current toy now has a much more exact description:

```text
grown anatomy G
      |
      +--> graph Laplacian spectrum {lambda_n, phi_n}
      |
input terminal A ---- phi_n(A) ---\
                                  +--> modal oscillator q_n(t)
input terminal B ---- phi_n(B) ---/
                                           |
                                           v
                                distributed field psi(x,t)
                                           |
                                  observation / projection
                                           |
                                           v
                                      consequence
```

The geometry is doing two concrete jobs:

1. it defines the normal coordinates and time scales in which the electrical field can move;
2. it determines how each physical input couples into those coordinates.

The field does not need to freeze. It can remain a superposition of moving modal states. Computation can live in how those states are driven and which projection is consequential.

## The next real difficulty

The graph eigenbasis is a microscope, not yet a soma.

A biological or physically local soma cannot simply reach across the entire arbor and multiply every cell by an arbitrary signed eigenvector weight. The next question is therefore an **observability/localization** question:

> Can the confirmed task-relevant modal band be observed by a small, local, physically plausible soma aperture because the field itself brings those modes together there?

That reconnects this result to the soma tap experiment instead of replacing it.

A good next probe is to measure how observable the confirmed band is from graph-distance balls around the soma, and whether the coherent Gaussian tap succeeded because it couples to that band better than uniform or misplaced apertures.
