# GeometricNeuronPlusField

A small experimental branch of the Geometric Neuron / Functional Arbor line.

The question is now more precise than "does geometry compute?"

> **What does frozen anatomy contribute to a moving electrical field, where does task information become locally available, and what active boundary is needed to turn that analog field into a temporally precise output event?**

The source soma-tap experiment is in [`FunctionalArbors/SomaTapTestsClaude`](https://github.com/anttiluode/FunctionalArbors/tree/main/SomaTapTestsClaude).

## Current working picture

```text
slow anatomy / operator G
        |
        v
geometry-defined modal filter bank
        |
        v
moving field psi(x,t)
        |
        v
convergence / soma task bottleneck
        |
        v
[active AIS-like state h(t) -- not built yet]
        |
        v
spike event timing
```

The graph basis is a **microscope**, not a claim that biological somata calculate eigenvectors.

## What survived

### 1. Temporal-order information is spectrally structured

On frozen FunctionalArbor bodies the energy-dominant spatial common mode is almost perfectly blind to A/B order, while a higher graph-Laplacian band (modes 18–20 in the registered setup) carries appreciable order information and gains selectivity from coherent projection.

The discovery-set band survived a held-out 12-body confirmation. See [`DISCOVERY_V01.md`](DISCOVERY_V01.md) and [`CONFIRMATION_V01.md`](CONFIRMATION_V01.md).

### 2. The distributed wave is almost exactly a geometry-defined resonator bank

A parameter-free reduced model

```text
q_n'' + damping q_n' + (restoring + stiffness*K*lambda_n) q_n
    = phi_n(A) s_A(t) + phi_n(B) s_B(t)
```

predicts the full-field modal contrast almost exactly (`pooled r ~= .995`). That is mathematically expected for the fixed-K graph wave, but conceptually useful: anatomy selects the poles/time scales and input couplings of an LTI resonator bank. See [`MODAL_MECHANISM_V01.md`](MODAL_MECHANISM_V01.md).

### 3. A one-cell anatomical edit is global in modal coordinates

The modal-locality audit perturbs real frozen bodies one cell at a time. One added cell substantially disturbs about 24% of the modal identities on average; the confirmed task band loses about 8% identity, and the eigenvector change is only weakly localized around the edit. See [`MODAL_LOCALITY_V01.md`](MODAL_LOCALITY_V01.md).

This makes local per-cell credit a badly aligned coordinate for the global computation. It does **not** prove that event-level local eligibility is impossible.

### 4. Locality compresses the global modes, but the soma is a favorable local mixing point

A compact patch anywhere on the tree is poor at reconstructing three global graph modes compared with scattered graph-wide taps. With the fairer control — same-radius compact balls around every other body cell — the soma neighborhood is consistently better-conditioned than most equally local apertures, without containing unusually high task-band energy. See [`LOCAL_OBSERVABILITY_V02.md`](LOCAL_OBSERVABILITY_V02.md).

### 5. The actual task scalar is strongly concentrated at the soma/root

The stronger functional result is simpler. Across 24 frozen bodies:

```text
mean soma |temporal-order contrast|      0.2283
mean over all occupied cells             0.0541
median soma percentile among body cells  0.8786
soma in top quartile                     24 / 24
```

Coherent graph balls around the soma remain far more task-selective than same-radius local balls elsewhere from radius 1 through 6. The soma is therefore better described as a **task bottleneck** than as a miniature reconstruction of the global field. See [`TASK_BOTTLENECK_V02.md`](TASK_BOTTLENECK_V02.md).

Construction caveat: FunctionalArbor begins with a designated soma/root and grows source connectivity around it. This result shows what convergence geometry does; it does not show that an unbiased developmental process discovered the output site.

## The AIS bridge

The current toy stops at a passive quadratic readout. Biology does not.

The axon initial segment sits immediately downstream of the somatodendritic convergence region and is an active, stateful, plastic compartment enriched in voltage-gated channels. Its length, location, channel composition and channel kinetics affect excitability and temporal encoding.

That suggests a narrower and testable architectural hypothesis:

> **Distributed geometry supplies an analog modal computation; convergence makes a task-relevant mixture locally available; an AIS-like active boundary converts that moving mixture into sparse, history-dependent, temporally precise output events.**

See [`AIS_BRIDGE.md`](AIS_BRIDGE.md) and [`OPERATOR_FLOW_EVENT.md`](OPERATOR_FLOW_EVENT.md).

The next model should not ask an AIS-like compartment to reconstruct the global modes. It should operate only on the local soma mixture and be tested against rate-matched passive/memoryless controls.

## Experiments

- [`graph_mode_probe.py`](graph_mode_probe.py) — graph spectral microscope + original live-field/settled-readout test.
- [`confirm_graph_band.py`](confirm_graph_band.py) — held-out confirmation of the discovery band.
- [`modal_mechanism_probe.py`](modal_mechanism_probe.py) — reduced modal oscillator mechanism.
- [`modal_locality.py`](modal_locality.py) — one-cell structural perturbation audit.
- [`local_observability_probe.py`](local_observability_probe.py) — soma-local observability versus scattered and contiguous controls.
- [`task_bottleneck_probe.py`](task_bottleneck_probe.py) — actual task selectivity over every cell and local aperture.

## Live-field / settled-readout result

The original strict criterion remains a null. The field stayed live, but zero modes passed the prewritten combination of stability, energy and selectivity thresholds. Those thresholds were not relaxed after seeing the data.

That null is retained because "the field can keep moving while the computation settles" still needs a better operational demonstration than a stable near-zero observable.

## Run

Clone this repo beside `FunctionalArbors`:

```text
parent/
  FunctionalArbors/
  GeometricNeuronPlusField/
```

Then, for example:

```bash
pip install -r requirements.txt
python graph_mode_probe.py --functional-arbors ../FunctionalArbors --seeds 12 --modes 24
python modal_locality.py --functional-arbors ../FunctionalArbors --seeds 12
python local_observability_probe.py --functional-arbors ../FunctionalArbors --seeds 24
python task_bottleneck_probe.py --functional-arbors ../FunctionalArbors --seeds 24
```

GitHub Actions runs the canonical receipts and uploads their JSON/figure artifacts.

## What this is not

This is an experimental computational model, not a claim that real dendrites are literal graph-Laplacian computers or that the AIS exists for one single evolutionary reason. Real AIS biology also includes maintenance of axonal polarity and trafficking barriers, among other functions.

The useful target is narrower: identify which pieces of a geometry/field/event architecture are actually necessary in the toy, then compare that division of labor with known neuronal compartments without retrofitting biology to the model.
