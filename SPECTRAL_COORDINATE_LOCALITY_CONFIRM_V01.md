# Spectral coordinate locality on the held-out benchmark — the winning free knobs are strongly nonlocal

After `MATCHED_TUNER_CONFIRM_V01.md` was complete, the exact eight winning free-modal coordinates from the same fresh bodies were audited in the physical 31x31 grid basis. This is a **descriptive hardware-cost audit**, not a second confirmatory test on the already-opened seeds.

Across seeds 288-299 there were 96 selected free coordinates:

```text
pole/stiffness F        51
soma/output residue C   26
source-B residue B      11
source-A residue A       8
```

## Pure pole coordinates

For a normalized mode `phi_n`, a pure direct modal-stiffness shift that preserves the eigenvectors corresponds to

```text
Delta K_n ∝ phi_n phi_n^T.
```

A nearest-neighbour grid can directly alter only diagonal and adjacent off-diagonal entries. The Frobenius mass of `phi_n phi_n^T` outside that pattern is therefore an unavoidable lower bound on the relative operator error of realizing the pure spectral coordinate inside the same local-bond sparsity class.

For the 51 selected `F` coordinates:

```text
mean mode participation ratio                         29.89 nodes
mean local-bond matrix residual lower bound           0.95334
minimum residual lower bound                          0.90620
maximum residual lower bound                          0.96997
```

Even the easiest selected pure pole coordinate leaves at least about **90.6% relative matrix residual** under the local nearest-neighbour sparsity pattern.

## Pure residue coordinates

A direct `A/B/C` modal residue coordinate corresponds in the physical node basis to a distributed vector proportional to `phi_n`.

For the 45 selected residue coordinates:

```text
mean mode participation ratio                         36.93 nodes
mean best-one-node vector residual                    0.97455
minimum best-one-node residual                        0.95722
maximum best-one-node residual                        0.98363
```

So the winning residue knobs are likewise not equivalent to one local physical port/readout adjustment.

## Combined reading with the benchmark

Both facts must be retained simultaneously:

1. **Abstract parameter count:** F8 beats G8 decisively — 12/12 bodies, mean held-out advantage `+0.1400`.
2. **Physical coordinate locality:** the winning F8 variables are strongly nonlocal when expressed in the nearest-neighbour grid basis.

Therefore it is wrong to say either:

```text
“geometry wins”
```

or

```text
“eight free modal knobs are physically equivalent to eight local bond tuners.”
```

The actual result is:

> **Direct spectral coordinates are much better optimization coordinates for this temporal task, but they are oracle/global coordinates relative to the local scattering medium.**

The next hardware benchmark must price the compiler that turns a desired spectral change into realizable local device controls.

## What should be counted next

A useful physical-resource comparison should count, at minimum:

```text
number of independently controlled local elements
number/location of monitors
forward/backward physical passes per gradient step
precision per control element
loss and thermal/crosstalk budget
energy per update
latency
calibration/compiler overhead
```

The result here motivates that benchmark; it does not supply those costs by itself.
