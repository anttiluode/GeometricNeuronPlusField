# Soma cross-term mode-pair discovery v0.1 — preregistration

## Question

The transfer decomposition reduced most soma temporal-order selectivity to the coherent cross-source term

```text
2 Re[h_A(t) conj(h_B(t-tau))]
```

while the graph-mode branch showed that anatomy supplies a natural Laplacian modal basis. The next question is therefore:

> **Which graph-mode pairs create the order-sensitive cross term at the soma, and is that contribution concentrated in a small modal subspace or distributed broadly across the arbor spectrum?**

This is a discovery experiment. It is not allowed to claim held-out confirmation from the same bodies. Any fixed modal subset or quantitative concentration prediction suggested by this run must be frozen in a separate preregistration before testing new bodies.

## Frozen setup

Use previously unused FunctionalArbor bodies:

```text
seeds       72-83
lag         20
steps       210
source gain A/B = 1/1
```

No growth, learning, readout, wave, pulse, or morphology parameters are changed.

For each frozen body, record the two complex single-source field histories

```text
h_A(x,t)
h_B(x,t)
```

and diagonalize the full unweighted 4-neighbour body graph Laplacian

```text
L phi_n = lambda_n phi_n.
```

The basis is a diagnostic microscope only; no claim is made that a biological soma explicitly computes eigenvectors.

## Exact modal decomposition

For occupied-cell field vector `z_A(t)`, define

```text
q_A,n(t) = <phi_n, z_A(t)>
q_B,n(t) = <phi_n, z_B(t)>.
```

At the soma `s`, each mode contributes

```text
u_A,n(t) = q_A,n(t) phi_n(s)
u_B,n(t) = q_B,n(t) phi_n(s).
```

The reconstructed soma histories are

```text
h_A(s,t) = sum_n u_A,n(t)
h_B(s,t) = sum_n u_B,n(t).
```

The implementation must verify near-machine-precision reconstruction before interpreting pair contributions.

## Order-sensitive pair matrix

Using the coherent single-source reconstruction, find the peak times of the target A->B and distractor B->A soma power traces. At those fixed peak times define the ordered mode-pair contribution

```text
M_nm = 2 Re[u_A,n(t_T) conj(u_B,m(t_T-tau))]
       -
       2 Re[u_B,n(t_D) conj(u_A,m(t_D-tau))].
```

Thus

```text
sum_nm M_nm
```

must exactly equal the target-minus-distractor difference of the soma cross term at the same peak times.

For concentration analysis convert this to an unordered pair matrix without changing the sum:

```text
U_nn = M_nn
U_nm = M_nm + M_mn,  n < m.
```

All pair-concentration metrics use `|U_nm|`; signed cancellation is reported separately.

## Discovery metrics

For each body report:

1. soma single-source modal reconstruction error;
2. pair-matrix reconstruction error;
3. number and fraction of unordered mode pairs needed to contain 50% and 80% of absolute pair mass;
4. participation-ratio effective pair count;
5. cancellation ratio `|sum U| / sum |U|`;
6. fraction of absolute pair mass on diagonal pairs `n=n` versus cross-mode pairs;
7. absolute involvement of each mode, assigning each off-diagonal pair half to each member;
8. involvement and enrichment of the previously confirmed graph-mode band `{18,19,20}` relative to its dimensional share;
9. the highest-involvement mode indices and strongest unordered mode pairs;
10. a coarse spectral-block interaction matrix.

## Interpretation rules

No arbitrary cutoff will be invented after the run and called a preregistered definition of "sparse." The discovery receipt will state the observed concentration directly.

Evidence for a compact modal mechanism would look like a small pair fraction carrying most absolute mass and/or a small fixed set of mode indices repeatedly carrying high involvement across bodies.

Evidence for a distributed mechanism would look like a large participation ratio, many pairs required for 50/80% mass, unstable top modes, and weak enrichment of the old 18-20 band.

A strong cancellation ratio effect is itself informative: the computation may depend on many large modal interactions whose signed contributions cancel rather than on a few dominant positive pairs.

## Held-out rule

After inspecting seeds 72-83, freeze one or more concrete predictions — for example a fixed mode subset, a concentration range, or a spectral block — and test only those predictions on new bodies in a separate confirmation run.
