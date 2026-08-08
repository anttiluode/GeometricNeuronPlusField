# Pole containment v0.1 — the full reciprocal graph wave is exactly a pole bank

`BENCHMARK_SCOPE_V01.md` noted a necessary sanity check before any graph-vs-filter benchmark: a frozen finite reciprocal linear wave operator diagonalizes into independent damped modes, so a sufficiently free pole/residue bank must contain the graph model.

`full_pole_containment.py` diagonalized the **exact 31x31 mature weighted operator** used by the adjoint surrogate, including weak mature bath coupling and leakage to the zero exterior boundary. Reused seed 264 was sufficient because this is an identity check, not a statistical discovery.

## Result

```text
state dimension                              961
operator construction relative error       1.04e-16

soma target-trace relative L2 error         1.85e-14
soma distractor-trace relative L2 error     2.19e-14

energy target
 direct                                     0.3941558964131438
 modal                                      0.3941558964131442

energy distractor
 direct                                     0.2998329200205082
 modal                                      0.2998329200205105

contrast direct                             0.1359142599406042
contrast modal                              0.1359142599406008
absolute contrast error                     3.44e-15
```

The exact spatial simulation and the unconstrained modal/pole reconstruction are identical to floating-point precision.

## Consequence

This closes off one tempting but wrong benchmark claim:

> A linear reciprocal graph cannot have more abstract input-output expressivity than a fully free pole/residue bank that contains all of its modes.

Any “graph wins” result must therefore be about something else:

```text
fewer trainable physical tuners
structured parameter sharing
locality
robustness
training/measurement cost
fabrication constraints
```

not about escaping the pole-bank function class.

## The interesting property of a local conductance tuner

A direct pole coordinate changes one chosen spectral degree of freedom.

A local conductance change instead perturbs the operator itself:

```text
L -> L + DeltaL_e
```

and therefore moves many poles **and** residues coherently. One physical scalar can create a high-dimensional but constrained transfer-function deformation.

That is the parameter-efficiency hypothesis worth benchmarking.

## Next benchmark

The next comparison should start from the same exact base transfer function and give both arms the same number of additional real tuners:

- local conductance coordinates in the spatial operator;
- direct modal coordinates (pole frequency and input-specific residue coordinates).

Both should optimize the same multi-lag temporal objective. The free-modal arm should be allowed to choose its most sensitive coordinates, making it a strong baseline. The held-out question is then whether the **structured deformation induced by local geometry** provides better generalization per tuner than direct spectral tuning.
