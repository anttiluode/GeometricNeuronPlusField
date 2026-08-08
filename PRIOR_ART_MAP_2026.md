# Prior-art map — where this project now sits

The reciprocal-adjoint result triggered the right literature search. Several pieces of the current project have close and important precedents. That is useful: it separates established physical-computing machinery from the parts that still need to earn novelty.

## 1. Hermans et al. 2015 — reciprocal dynamical media + physical BPTT

Michiel Hermans, Michaël Burm, Thomas Van Vaerenbergh, Joni Dambre, Peter Bienstman et al., **“Trainable hardware for dynamical computing using error backpropagation through physical media,”** *Nature Communications* 6, 6729 (2015), DOI `10.1038/ncomms7729`.

This is the earliest close temporal precedent found in this search and is extremely relevant to the exact time-reversal result here.

Their theory starts from a reciprocal **linear dynamical system** with time-dependent impulse responses. During the backward phase, source and receiver roles are exchanged; reciprocity supplies the transposed impulse-response operator, while physical time is run forward after **time-reversing the external error/Jacobian signals**. They explicitly state that this physically implements error backpropagation through the same reciprocal system.

They experimentally demonstrate the principle with an **acoustic wave medium** (speaker + 6 m tube + microphone) whose travelling/reflected waves provide the task memory, and they discuss electro-optical extensions.

Thus the broad idea

```text
time-dependent reciprocal physical medium
+ time-reversed error waveform
+ same hardware used backward
= physical backpropagation through temporal dynamics
```

predates the photonic TRIM work. Their experiment trains masks rather than the internal acoustic medium itself, but the temporal same-medium backpropagation principle is already explicit.

## 2. Hughes et al. 2018 — physical adjoint / TRIM for internal photonic parameters

Tyler W. Hughes, Momchil Minkov, Yu Shi, Shanhui Fan, **“Training of photonic neural networks through in situ backpropagation and gradient measurement,”** *Optica* 5(7), 864–871 (2018), DOI `10.1364/OPTICA.5.000864`.

Core overlap:
- adjoint-variable derivation;
- reciprocal photonic propagation used for backward/error propagation;
- **time-reversal interference method (TRIM)**;
- local gradient obtained from forward/adjoint field interference/intensity measurements;
- gradient cost independent of number of tunable physical parameters.

This is especially close to the **internal bond-gradient** side of our result: the forward and adjoint fields physically meet at the parameter location, so the local overlap exposes the derivative.

Our exact soma-launched time-reversed credit wave therefore belongs to an established physical-adjoint class rather than being a new general training principle.

## 3. Hughes et al. 2019 — wave physics as an analog RNN

Tyler W. Hughes, Ian A. D. Williamson, Momchil Minkov, Shanhui Fan, **“Wave physics as an analog recurrent neural network,”** *Science Advances* 5, eaay6946 (2019), DOI `10.1126/sciadv.aay6946`.

This is even closer to the forward-computation side of the project:
- second-order wave dynamics are mapped directly to an RNN recurrence;
- spatial material parameters become trainable recurrent dynamics;
- raw temporal waveforms are processed by scattering through an inhomogeneous medium;
- they demonstrate inverse-designed vowel classification with wave propagation providing memory/computation.

So the broad claim

> “a spatial wave medium can act as a trainable recurrent temporal computer”

is established prior art.

## 4. Pai et al. 2023 — experimental silicon in-situ backpropagation

Sunil Pai et al., **“Experimentally realized in situ backpropagation for deep learning in photonic neural networks,”** *Science* 380(6643), 398–404 (2023), DOI `10.1126/science.ade8450`.

Core overlap:
- physical forward and backward/adjoint propagation;
- interference-based gradient measurement;
- programmable silicon photonic mesh;
- demonstrated multilayer learning on fabricated hardware;
- analog-domain gradient/update route explored.

This establishes that physical adjoint learning is not only a simulation method.

## 5. Li & Mao 2024 — fabricated mechanical networks and the same bond-local product

Shuaifeng Li and Xiaoming Mao, **“Training all-mechanical neural networks for task learning through in situ backpropagation,”** *Nature Communications* 15, 10528 (2024), DOI `10.1038/s41467-024-54849-z`.

This is the clearest mechanical analogue of our local bond-gradient formula found so far.

For a linear spring network with symmetric stiffness matrix `D`, they derive the exact local stiffness gradient as

```text
grad_i = e_forward,i * e_adjoint,i
```

where `e_i` is the elongation of bond `i`. They then **fabricate 2D mechanical networks**, physically measure forward and adjoint bond elongations, and recover gradients locally with high experimental accuracy.

The paper also states an important limitation that directly echoes our alpha-sweep result: the learning rule is exact in the **linear / small-deformation regime**; nonlinear deformation degrades the gradient correspondence.

So our dynamic result

```text
integral over time of
forward bond difference * adjoint bond difference
```

is not merely analogous to photonic TRIM. It has a very direct static mechanical sibling that has been physically built and measured.

## 6. Thakkar & Grbic 2026 — the closest device-class match found so far

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

## 7. Time-domain acoustic/metamaterial adjoints are also established

For example, Lin et al., **“Topology and shape optimization of broadband acoustic metamaterials and phononic crystals,”** *Acoustical Science and Technology* 38(5), 254–260 (2017), DOI `10.1250/ast.38.254`, uses a **time-dependent adjoint** for topology/shape optimization of acoustic wave media.

So “time-domain wave topology + adjoint sensitivity” is not itself new either.

## 8. What the repo should no longer present as a novelty

Do not claim novelty for these broad ideas:

```text
wave propagation as recurrent computation
continuous material/coupling optimization by adjoints
reciprocal backward physical propagation for temporal gradients
time-reversing an error waveform to run BPTT through reciprocal dynamics
forward/backward interference as local gradient measurement
bond-local forward x adjoint stiffness learning
2D tunable wave grids as neuromorphic hardware
```

There is strong prior art for all of them.

## 9. What remains specific and potentially interesting here

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

## 10. The bar just got higher — use that

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
