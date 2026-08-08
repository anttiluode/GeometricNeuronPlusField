# Hardware pivot v0.1 — from neuron analogy to differentiable reciprocal scattering mesh

The reciprocal-adjoint result changes the center of gravity of this repo.

The mechanism is not a new category of learning rule. It belongs to an established physical-computing literature:

- Tyler W. Hughes, Momchil Minkov, Yu Shi, and Shanhui Fan, **“Training of photonic neural networks through in situ backpropagation and gradient measurement,”** *Optica* 5, 864–871 (2018), DOI `10.1364/OPTICA.5.000864`.
- Sunil Pai et al., **“Experimentally realized in situ backpropagation for deep learning in photonic neural networks,”** *Science* 380, 398–404 (2023), DOI `10.1126/science.ade8450`.
- Earlier physical-dynamical backpropagation work includes Michiel Hermans et al., **“Trainable and Dynamic Computing: Error Backpropagation through Physical Media,”** arXiv:1407.6637 (2014).

Hughes et al. derive physical adjoint/backpropagation and local gradient measurement from interference of forward and adjoint fields. Pai et al. experimentally trained a silicon photonic mesh by measuring backpropagated gradients from forward/backward optical interference.

So the correct interpretation of our machine-precision reciprocal replay is:

> **an independent derivation, in this specific damped transient scattering model, of the same physical-adjoint principle used by in-situ-trained wave hardware.**

That is a promotion in physical plausibility and a deflation in novelty of the general gradient mechanism.

## Reciprocity: what it buys and what it costs

For the mature linear medium,

```text
H = H^T
```

because the bond couplings are reciprocal. That symmetry is exactly what lets the adjoint spatial operator be implemented by the same physical body.

The repo's older arrow-of-time failures and the exact reciprocal-gradient result are therefore two sides of one property:

> **Reciprocity does not supply an intrinsic directional computational arrow, but it makes physical adjoint transport exceptionally cheap.**

That is now a design trade, not a bug to be patched away.

## Biology: narrow the claim rather than making an absolute slogan

The reciprocal-gradient result should **not** be promoted as evidence that biological neural networks perform exact backpropagation this way.

Chemical synaptic transmission between neurons is strongly directional, active neuronal membranes are nonlinear/stateful, and no mechanism has been shown here that stores and emits the exact time-reversed objective derivative required by the mathematical adjoint.

However, “therefore every biological analogy is dead” is also broader than the result supports. Passive subthreshold cable propagation within a dendritic tree can be approximately reciprocal in a linearized regime, while the network-level synaptic graph is not. The safe split is:

```text
forward dendritic/field analogies                 still a separate empirical question
exact same-medium adjoint as biological credit    NOT established here
reciprocal engineered wave hardware               directly relevant
```

The hardware pivot is therefore about where the strongest evidence now points, not a theorem that every upstream neuronal comparison is impossible.

## What the missing selection control settled

The old graded-coupling result selected bonds because the initial adjoint derivative liked them. `RANDOM_FRONTIER_CONTROL_V01.md` removed that conditioning.

Fresh random legal frontier additions:

```text
72 total
53 / 72 = 73.6% interior optimum
median interior alpha = .15
87.5% slope-sign reversal
mean binary regret = .01459
```

The positive-gradient subset was more interior (`86.1%`), proving that selection enriches the effect, but it does not create it.

So continuous coupling strength remains a genuine computational coordinate in this model.

## The abstract scalar benchmark also settled

The differentiable topology finally enabled the benchmark that had been waiting in the background.

`MATCHED_TUNER_CONFIRM_V01.md` compared eight local bond coordinates with eight freely selected modal pole/residue coordinates on a multi-lag held-out task.

Fresh 12 bodies:

```text
G8 local bond tuners mean test improvement       +.04915
F8 free spectral tuners mean test improvement    +.18918
F8 - G8                                           +.14004
F8 beats G8                                       12 / 12
```

So the local spatial mesh **does not** win on performance per unconstrained abstract scalar. That claim is closed negatively.

The reason to build the mesh must be physical:

```text
local actuation
physical parallel propagation
in-situ gradient acquisition
fabrication / wiring economics
robustness / self-calibration
latency / bandwidth / energy
multi-port spatial tasks
```

not superior free-pole expressivity.

## The broadband wrinkle

Standard steady/coherent photonic in-situ backpropagation has a very natural local interference readout. Our toy's task is instead a damped finite-time transient.

The exact causal reciprocal replay exposed a real systems wall:

```text
adjoint transport     cheap: same medium
local gradient        requires forward-history x returned-history correlation
```

Naive simultaneous products failed.

`SPECTRAL_CORRELATION_COMPRESSION_V01.md` then showed that the anti-diagonal time correlation has an exact DFT decomposition and is strongly sparse for this task. On fresh bodies:

```text
8 boundary-selected bins     mean gradient-map corr .9919
16 bins                      mean corr .9973, mean rel-L2 .0672
mean bins for 95% absolute spectral gradient mass = 13 / 210
```

All five held-out criteria passed.

`SPECTRAL_GRADIENT_LEARNING_V01.md` then put the compressed gradient inside the optimizer. On fresh bodies:

```text
exact mean DeltaC       +.02879
K8 mean DeltaC          +.02469   (85.8% of exact group gain; 12/12 improve)
K16 mean DeltaC         +.03025
```

Again all registered criteria passed.

That gives this repo a more specific hardware contribution to investigate:

> **How should broadband/time-domain in-situ backpropagation be implemented when the exact local correlation is delayed in time but compressible into a small common set of coherent spectral channels?**

## Current compiler picture

```text
TASK / PORT WAVEFORMS
        |
        v
RECIPROCAL GRADED SCATTERING MESH
local coupling field rho_e
        |
        +------ forward physical computation ------+
        |                                           |
        v                                           v
source-specific transfer histories          output / objective
                                                    |
                                                    v
                                         derivative waveform
                                                    |
                                             reverse in time
                                                    |
                                                    v
                                  same mesh carries adjoint field
                                                    |
            local K-bin forward phasors             |
                       \                            /
                        \                          /
                         +-- local coherent product
                                  |
                                  v
                           dJ / d rho_e
                                  |
                                  v
                         local graded update
```

The software role is no longer best described as “simulate a neuron.”

It is becoming a **differentiable compiler / training harness for reciprocal analog wave meshes**, with FunctionalArbor morphology providing one structured family of scattering geometries.

## What is still not solved

1. **Physical phasor readout.** The current compressed learner numerically multiplies stored complex phasors. A chip-level intensity/phase-stepping or analog multiplier implementation remains to be specified.
2. **Device-calibrated error.** Phase/amplitude calibration error, tap/readout noise, loss, drift, and reciprocity mismatch should be tested at ranges corresponding to actual hardware rather than arbitrary corruption.
3. **Physical cost benchmark.** `8 mathematical scalars vs 8 mathematical scalars` favors the free spectral bank. The next fair comparison must count actuators, detectors/taps, calibration operations, passes, memory, and energy.
4. **Task class.** The present task is small and highly structured. A spatial mesh needs a task where its locality/multi-port structure can plausibly matter.

## Current wall sentence

> **The reciprocal mesh does not beat a free spectral model per abstract scalar; its opportunity is that the same physical wave operator computes forward, carries the exact adjoint backward, and now appears able to compress broadband local gradient readout to a small number of coherent channels.**
