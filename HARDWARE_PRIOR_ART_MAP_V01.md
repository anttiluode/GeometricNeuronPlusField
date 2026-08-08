# Hardware prior-art map v0.1 — where this repo sits after the reciprocal-adjoint result

The exact same-medium adjoint should not be presented as a new learning principle. There is already a substantial physical-backpropagation literature.

This note records the nearest known neighbors so later claims stay narrow.

## 1. Hughes et al. 2018 — photonic in-situ backpropagation

Tyler W. Hughes, Momchil Minkov, Yu Shi, and Shanhui Fan,
**“Training of photonic neural networks through in situ backpropagation and gradient measurement,”** *Optica* 5, 864–871 (2018), DOI `10.1364/OPTICA.5.000864`.

They derive the photonic analogue of backpropagation with the adjoint-variable method and show that parameter gradients can be obtained physically from interference/intensity measurements inside the device.

Our reciprocal-wave derivation belongs to that same mathematical family.

## 2. Pai et al. 2023 — experimentally realized silicon photonic training

Sunil Pai et al., **“Experimentally realized in situ backpropagation for deep learning in photonic neural networks,”** *Science* 380, 398–404 (2023), DOI `10.1126/science.ade8450`.

Their packaged silicon photonic network used bidirectional propagation, internal grating-tap monitors, universal optical amplitude/phase control, and interference of forward and adjoint optical fields to measure gradients. They also explored an analog-domain gradient/update protocol.

This establishes that the broad physical idea is not merely theoretical: in-situ backpropagation through a reciprocal wave network has been built and used to train a physical system.

## 3. Thakkar & Grbic 2026 — the closest architectural cousin found so far

Shrey Thakkar and Anthony Grbic,
**“Wave-based neuromorphic circuit networks: tunable 2D transmission-line metamaterials,”** *Optical Materials Express* 16(8), 2542–2559 (2026), DOI `10.1364/OME.599576`; arXiv:2606.00194.

This is especially relevant because it moves beyond an MZI matrix mesh toward a **2D grid of subwavelength transmission-line unit cells with tunable reactive elements**. Computation occurs through wave propagation and interference across the grid, and training is derived from in-situ adjoint backpropagation using forward and adjoint excitations.

That is very close to the hardware class suggested by `GeometricNeuronPlusField`: a locally tunable scattering medium whose state is stored in distributed physical couplings.

The paper describes/proposes and numerically trains the transmission-line networks; unlike Pai et al., it should not be cited here as an experimental fabrication of this exact 2D-grid architecture unless a later hardware paper establishes that.

## 4. What is therefore NOT novel here

Do not claim novelty for:

```text
physical adjoint backpropagation in reciprocal wave media
same-device forward and backward propagation
local forward/adjoint overlap as a parameter gradient
wave-interference computing in a tunable 2D medium
continuous local reactive/coupling parameters as learned state
```

There is already strong prior art for all of those themes.

## 5. What this repo may still contribute

The current branch has several narrower wrinkles that are not erased by that prior art.

### A. Broadband finite-time rather than single-tone / steady-state training

The FunctionalArbor task is a damped transient temporal-order problem. The exact adjoint is therefore a **time-reversed derivative waveform**, not merely a second same-frequency phasor excitation.

That exposed a delayed local-correlation problem that is largely absent in a single-tone steady-state derivation.

### B. Spectral compression of the delayed local correlation

`SPECTRAL_CORRELATION_COMPRESSION_V01.md` shows that the exact 210-sample local time correlation can be decomposed into frequency-bin products, and on held-out bodies about 13/210 bins carry 95% of the absolute gradient mass. A common K=8 or K=16 frequency set can be selected from boundary signals alone.

### C. Closed-loop learning from the compressed broadband gradient

`SPECTRAL_GRADIENT_LEARNING_V01.md` shows that the K=8 compressed gradient retains 85.8% of the exact learner's group gain on fresh bodies, while K=16 essentially preserves it under the registered finite-step optimizer.

### D. Continuous graded spatial coupling is empirically important in this particular transient mesh

The random-frontier control confirmed that intermediate bond strengths are commonly optimal even without conditioning on gradient sign.

### E. Negative benchmark against free spectral coordinates

The repo already knows where **not** to claim superiority: direct free modal/pole-residue coordinates beat local bonds decisively per abstract scalar on the current temporal benchmark.

So any hardware advantage must be measured in physical cost rather than mathematical parameter count.

## 6. The useful claim boundary

A defensible current description is:

> **GeometricNeuronPlusField is now a transient/broadband differentiable reciprocal-scattering-mesh testbed. Its general adjoint-training principle is established prior art; the active research question is how cheaply a distributed physical mesh can measure and apply broadband local gradients, and whether local physical implementation wins under a real hardware cost model.**

That is a much stronger place to work from than pretending the entire mechanism was newly invented here.
