# GeometricNeuronPlusField

Experimental branch of the Geometric Neuron / Functional Arbor line.

The question is now more precise than "does geometry compute?"

> **What does frozen anatomy contribute to a moving electrical field, where does task information become locally available, and what active boundary is needed to turn that analog field into output events?**

The source soma-tap experiment is in [`FunctionalArbors/SomaTapTestsClaude`](https://github.com/anttiluode/FunctionalArbors/tree/main/SomaTapTestsClaude).  The amplitude-balance follow-up is in this repo under [`SomaWhyClaude/`](SomaWhyClaude/).

## Current working picture

```text
slow anatomy / operator G
        |
        v
geometry-defined modal resonator bank
        |
        v
moving field psi(x,t)
        |
        v
convergence root / amplitude-balance point
        |
        v
active state V,m,h,n
        |
        v
spike events
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

The modal-locality audit perturbs real frozen bodies one cell at a time. One added cell substantially disturbs about 24% of modal identities on average; the confirmed task band loses about 8% identity, and the eigenvector change is only weakly localized around the edit. See [`MODAL_LOCALITY_V01.md`](MODAL_LOCALITY_V01.md).

This makes local per-cell credit a badly aligned coordinate for the global computation. It does **not** prove that event-level local eligibility is impossible.

### 4. Locality compresses the global modes, but the soma region is a favorable local mixture

A compact patch anywhere on the tree is poor at reconstructing three global graph modes compared with scattered graph-wide taps. With the fair control — same-radius compact balls around every other body cell — the soma neighborhood is consistently better-conditioned than most equally local apertures, without containing unusually high task-band energy. See [`LOCAL_OBSERVABILITY_V02.md`](LOCAL_OBSERVABILITY_V02.md).

### 5. The soma/root task effect is largely an amplitude-balance effect

The original task-bottleneck result was strong: soma temporal-order contrast is much larger than the average body cell.  Claude's `SomaWhyClaude` follow-up reproduced it and supplied a simpler mechanism.

Across its 12 frozen bodies, per-cell order selectivity correlated with single-source amplitude balance

```text
b(x) = min(p_A,p_B) / max(p_A,p_B)
```

at about `r = +0.70`, while total energy had essentially no explanatory power.  The soma/root is therefore better described here as the **A/B amplitude-balance / coincidence point produced by the convergence construction**, not as an independently discovered magical soma location.

This yields a clean falsifiable prediction: move one source and the most selective point should track the new balance point rather than remain fixed at the designated soma.

## The AIS bridge was built

The old repo stopped at a passive power readout.  We have now frozen the entire upstream system and added a compact HH-like active boundary immediately downstream of the soma readout.

[`ais_active_probe.py`](ais_active_probe.py) contains state

```text
V, m, h, n
```

with Na activation/inactivation and K activation/recovery.  No derivative or hand-written high-pass term is inserted.

It is compared with:

1. the raw soma signal under a rate-matched memoryless threshold;
2. the **measured small-signal linearization of the same active membrane**, convolved with the same soma signal and rate matched;
3. the full active stateful boundary.

See [`AIS_ACTIVE_PREREG_V01.md`](AIS_ACTIVE_PREREG_V01.md), [`AIS_ACTIVE_V01.md`](AIS_ACTIVE_V01.md), [`AIS_ACTIVE_PREREG_V02.md`](AIS_ACTIVE_PREREG_V02.md), and [`AIS_ACTIVE_V02.md`](AIS_ACTIVE_V02.md).

## The important AIS result is a null in the right place

The first run made the active boundary look spectacular at high frequency, but that was an estimator trap: sparse conditions can give `vector strength = 1` from one spike.  We did not keep that result.

v0.2 froze the mechanism and strengthened the control:

- exact **per-frequency** event-count matching;
- timing scored only for active conditions with at least four events;
- PPC (pairwise phase consistency) alongside vector strength;
- the active gate compared primarily against its own linearized membrane response.

Registered upper-band result (`f >= 0.05`, 21 valid body/frequency pairs from 11 bodies):

```text
mean PPC
active        0.4305
linearized    0.5816
memoryless    0.5238

active - linearized   -0.1511
wins / losses          6 / 15
```

So the first active boundary **does not improve rate-matched spike-time precision beyond its own linearization**.  See [`AIS_ACTIVE_V02.md`](AIS_ACTIVE_V02.md).

That result blocks the tempting next move: we are **not** co-adapting AIS position/extent around a failed precision claim.

## What the active state *does* do

Under one fixed operating point, the active membrane allocates spikes very differently across modulation frequencies.  It is non-monotonic rather than a generic high-pass: one regime is suppressed, another opens, and the highest tested regime closes again.

So the remaining question became mechanistic: which gating memory makes that selection?

### Gate-memory ablation

[`ais_gate_ablation_probe.py`](ais_gate_ablation_probe.py) keeps all gains/conductances fixed and removes independent gate history by making one gate instantaneous.

Results across 24 bodies:

```text
m instantaneous      fires 24/24; normalized response shape mostly survives
h instantaneous      fires  0/24
n instantaneous      fires  0/24
all instantaneous    fires  0/24
```

Fast Na activation memory therefore is not the main ingredient.  Independent h/n history is required to remain in the spiking regime, but the complete collapse prevents assigning the frequency-selection shape to either gate. See [`AIS_GATE_ABLATION_V01.md`](AIS_GATE_ABLATION_V01.md).

### Kinetic-speed test

[`ais_kinetics_probe.py`](ais_kinetics_probe.py) then changed only `dh/dt` or `dn/dt` speed by `0.5x` or `2x`, preserving steady-state curves, conductances and gain.

The preregistered prediction was that a true gate-set passband should move in opposite directions when the relevant gate is slowed versus sped up.

It did **not**.

`n` kinetics strongly alter excitability and normalized frequency allocation (`TV ~= .18` slow, `.35` fast), but faster/slower `n` do not produce a reliable bidirectional shift of the allocation center.  `h` effects are smaller and likewise fail the registered ordering. See [`AIS_KINETICS_V01.md`](AIS_KINETICS_V01.md).

Current earned statement:

> **Active channel state changes which upstream fluctuations become spikes and strongly controls event availability, especially through K recovery, but this generic HH-like boundary has not earned the stronger claim that its gating time constants implement a clean tunable frequency filter or improve spike-time precision beyond a matched linear filter.**

## A newly exposed interface problem

The active boundary experiments deliberately inherited the historical FunctionalArbor soma output:

```text
|psi_soma|^2
```

That is appropriate for the old power/readout objective, but a biological AIS receives a voltage-like membrane waveform, not a squared field magnitude.

Before changing AIS geometry, the next clean bridge test should therefore alter **only the interface observable**, with the active parameters frozen:

```text
power:      |psi_soma|^2
amplitude:  |psi_soma|
signed:     Re psi_soma
```

and repeat the memoryless + own-linearization controls.  This must be treated as a new preregistered experiment, not a rescue of the current null.

## Operator -> flow -> event

The current synthesis remains useful even with the AIS nulls:

```text
operator / morphology G      slow constraints on possible dynamics
        |
        v
flow / field psi(t)          continuous execution of those dynamics
        |
        v
local convergence            amplitude-balanced consequential mixture
        |
        v
active boundary state        nonlinear event availability/history
        |
        v
events                       sparse addressable output
```

See [`AIS_BRIDGE.md`](AIS_BRIDGE.md) and [`OPERATOR_FLOW_EVENT.md`](OPERATOR_FLOW_EVENT.md).

## Experiments

- [`graph_mode_probe.py`](graph_mode_probe.py) — graph spectral microscope + original live-field/settled-readout test.
- [`confirm_graph_band.py`](confirm_graph_band.py) — held-out confirmation of the discovery band.
- [`modal_mechanism_probe.py`](modal_mechanism_probe.py) — reduced modal oscillator mechanism.
- [`modal_locality.py`](modal_locality.py) — one-cell structural perturbation audit.
- [`local_observability_probe.py`](local_observability_probe.py) — soma-local observability versus scattered and contiguous controls.
- [`task_bottleneck_probe.py`](task_bottleneck_probe.py) — actual task selectivity over every cell and local aperture.
- [`ais_active_probe.py`](ais_active_probe.py) — first stateful active boundary and matched linearization.
- [`ais_active_probe_v02.py`](ais_active_probe_v02.py) — strict per-frequency rate-matched timing test.
- [`ais_gate_ablation_probe.py`](ais_gate_ablation_probe.py) — gating-memory ablation.
- [`ais_kinetics_probe.py`](ais_kinetics_probe.py) — h/n kinetic-speed test.

## Live-field / settled-readout result

The original strict criterion remains a null. The field stayed live, but zero modes passed the prewritten combination of stability, energy and selectivity thresholds. Those thresholds were not relaxed after seeing the data.

## What this is not

This is an experimental computational model, not a claim that real dendrites are literal graph-Laplacian computers or that the AIS exists for one single evolutionary reason. Real AIS biology also includes maintenance of axonal polarity and trafficking barriers, among other functions.

The useful target is narrower: identify which pieces of a geometry/field/event architecture are actually necessary in the model, then compare that division of labor with known neuronal compartments without retrofitting biology to the model.
