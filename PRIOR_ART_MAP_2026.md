# Prior-art map — where this project now sits

The reciprocal-adjoint result triggered the right literature search. Several pieces of the current project have close and important precedents. That is useful: it separates established physical-computing machinery from the parts that still need to earn novelty.

## 1. Hughes et al. 2018 — physical adjoint / TRIM

Tyler W. Hughes, Momchil Minkov, Yu Shi, Shanhui Fan, **“Training of photonic neural networks through in situ backpropagation and gradient measurement,”** *Optica* 5(7), 864–871 (2018), DOI `10.1364/OPTICA.5.000864`.

Core overlap:
- adjoint-variable derivation;
- reciprocal photonic propagation used for backward/error propagation;
- **time-reversal interference method (TRIM)**;
- local gradient obtained from forward/adjoint field interference/intensity measurements;
- gradient cost independent of number of tunable physical parameters.

Our exact soma-launched time-reversed credit wave therefore belongs to an established physical-adjoint class rather than being a new general training principle.

## 2. Hughes et al. 2019 — wave physics as an analog RNN

Tyler W. Hughes, Ian A. D. Williamson, Momchil Minkov, Shanhui Fan, **“Wave physics as an analog recurrent neural network,”** *Science Advances* 5, eaay6946 (2019), DOI `10.1126/sciadv.aay6946`.

This is even closer to the forward-computation side of the project:
- second-order wave dynamics are mapped directly to an RNN recurrence;
- spatial material parameters become trainable recurrent dynamics;
- raw temporal waveforms are processed by scattering through an inhomogeneous medium;
- they demonstrate inverse-designed vowel classification with wave propagation providing memory/computation.

So the broad claim

> “a spatial wave medium can act as a trainable recurrent temporal computer”

is established prior art.

## 3. Pai et al. 2023 — experimental silicon in-situ backpropagation

Sunil Pai et al., **“Experimentally realized in situ backpropagation for deep learning in photonic neural networks,”** *Science* 380(6643), 398–404 (2023), DOI `10.1126/science.ade8450`.

Core overlap:
- physical forward and backward/adjoint propagation;
- interference-based gradient measurement;
- programmable silicon photonic mesh;
- demonstrated multilayer learning on fabricated hardware;
- analog-domain gradient/update route explored.

This establishes that physical adjoint learning is not only a simulation method.

## 4. Thakkar & Grbic 2026 — the closest device-class match found so far

Shrey Thakkar and Anthony Grbic, **“Wave-based Neuromorphic Circuit Networks: Tunable 2D Transmission-Line Metamaterials,”** *Optical Materials Express* 16(8), 2542–2559 (2026), DOI `10.1364/OME.599576`; arXiv `2606.00194`.

This paper is strikingly close to the hardware interpretation that emerged here:
- a **2D grid** of interconnected subwavelength transmission-line unit cells;
- tunable reactive local elements;
- computation by wave propagation and interference across the grid;
- learned input-output relation stored in local tunable physical parameters;
- in-situ backpropagation from a forward and adjoint excitation;
- training demonstrated for allostery and classification;
- robustness to damage considered.

This is a much closer comparison than saying only “photonic MZI meshes exist.” A programmable 2D reactive wave network with local tuners and physical adjoint training is now explicit 2026 prior art.

Important difference: Thakkar/Grbic use **single-tone steady-state** excitation in the description/abstract, whereas this repo's core task is **time-domain pulse history / temporal order** in a damped second-order medium. That temporal distinction may matter.

## 5. Time-domain acoustic/metamaterial adjoints are also established

For example, Lin et al., **“Topology and shape optimization of broadband acoustic metamaterials and phononic crystals,”** *Acoustical Science and Technology* 38(5), 254–260 (2017), DOI `10.1250/ast.38.254`, uses a **time-dependent adjoint** for topology/shape optimization of acoustic wave media.

So “time-domain wave topology + adjoint sensitivity” is not itself new either.

## 6. What the repo should no longer present as a novelty

Do not claim novelty for these broad ideas:

```text
wave propagation as recurrent computation
continuous material/coupling optimization by adjoints
reciprocal backward physical propagation for gradients
forward/backward interference as local gradient measurement
2D tunable wave grids as neuromorphic hardware
```

There is strong prior art for all of them.

## 7. What remains specific and potentially interesting here

The distinctive package is narrower:

### A. Sparse morphology-like coupling geometry inside a weak background

The current object is not a generic dense trainable sheet. It begins as a field-grown branching structure embedded in weak coupling, then admits graded local strengthening/weakening.

### B. A temporal-order task analyzed mechanistically rather than only by accuracy

The repo has explicitly decomposed the task into:

```text
geometry-defined modes
source-specific transfer histories
sparse off-diagonal mode mixing
common/reference mode contribution
local square-law observation
lagged complex compatibility
```

That mechanistic chain is more specific than “the wave device classifies a signal.”

### C. The soma/root as a designated local task bottleneck

Given a body organized around a convergence root, task contrast is strongly concentrated there. This is a property of this morphology/task/readout construction, not a general wave-computing theorem.

### D. Binary structural events vs graded structural response

The repo quantified the radius in which the adjoint is useful for a structural conductance perturbation, then showed on a preregistered random-bond control that many response curves have genuinely interior optima and strong nonmonotonicity.

### E. Exact discrete time-domain same-medium credit identity for this model

The general adjoint principle is established, but the repo has an exact discrete derivation and machine-precision check for its specific damped update law, objective, and bond-local sensitivity.

### F. The possible compiler problem

The most useful open engineering question may be:

> **Can a sparse/locally tunable morphology-like wave network deliver a useful temporal transformation with materially cheaper control, sensing, fabrication, or adaptation than a generic programmable mesh or direct spectral filter?**

That is a hardware-architecture question, not a claim to have invented wave-based neural computation.

## 8. The bar just got higher — use that

The existence of very close prior work is not a reason to abandon the project. It gives us mature baselines and a vocabulary.

The repo should now behave like an architecture paper:

```text
known physical-computing principle
        +
our constrained geometry / temporal task / readout
        +
matched strong baselines
        +
measured hardware-relevant tradeoffs
```

The fastest way to find whether there is something new is to keep trying to kill the supposed advantage with those baselines.
