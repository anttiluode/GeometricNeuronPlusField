# Free spectral coordinate locality — development audit

The held-out matched-tuner benchmark counts `G8` and `F8` as eight mathematical scalar coordinates each. That is the correct abstract parameter-count comparison, but it is **not automatically a matched physical-hardware-cost comparison**.

To quantify the distinction, `spectral_coordinate_locality.py` reconstructed the 32 direct modal coordinates selected on reused development bodies 240-243 and expressed the relevant mode vector `phi_n` back in the 31x31 physical grid basis.

## Pole-frequency coordinate

A pure direct eigenvalue / modal-stiffness shift that preserves the mode basis has physical matrix perturbation proportional to

```text
Delta K_n = phi_n phi_n^T.
```

That is rank one but generally dense.

A nearest-neighbour bond network can directly change only diagonal and adjacent off-diagonal entries. Therefore the Frobenius norm of `Delta K_n` outside that sparsity pattern is an unavoidable lower bound on the error of realizing the pure modal coordinate with one-step local bond hardware.

For the 18 selected `F` coordinates:

```text
mean mode participation ratio                         29.95 nodes
mean local-bond matrix residual lower bound           0.9518
minimum residual lower bound                          0.9186
```

So even the most physically local selected pole coordinate has **at least 91.9% relative Frobenius residual** if one insists on the nearest-neighbour matrix sparsity pattern while trying to realize the pure rank-one modal perturbation directly.

That is not an optimization failure; it is a support mismatch. Most of the pure modal matrix lives in non-neighbour couplings.

## Source/output residue coordinate

A pure modal input or output residue perturbation corresponds in the physical basis to a vector proportional to `phi_n`.

For the 14 selected `A/B/C` coordinates:

```text
mean mode participation ratio                         38.82 nodes
mean best-one-node relative vector residual           0.9741
minimum best-one-node residual                        0.9538
```

A single local input/output site therefore cannot realize these direct residue coordinates closely. Their useful degree of freedom is spatially distributed.

## What this changes

It does **not** undo `MATCHED_TUNER_CONFIRM_V01.md`.

At equal abstract scalar count, the free spectral coordinates are better. That result stands.

But the scalar coordinate has different physical meaning:

```text
G8 coordinate
    one local bond conductance

F8 coordinate
    oracle global spectral deformation
    (pure pole or modal residue coordinate)
```

The free coordinate is exactly the sort of representation one would choose in software if arbitrary spectral control were free. It is not equivalent to one nearest-neighbour physical actuator.

## Important caveat

Do not convert the residual numbers directly into “number of required phase shifters.” A programmable interferometric mesh can synthesize global transforms through many local components, and the compiler/circuit architecture determines the actual actuator cost.

The earned statement is narrower:

> **The free baseline that wins per mathematical scalar is strongly nonlocal in the spatial basis of this nearest-neighbour medium. Therefore the next fair hardware comparison must count the physical resources required to synthesize those spectral coordinates, rather than identifying one global spectral coordinate with one local bond tuner.**

That is now the clearest hardware benchmark direction.
