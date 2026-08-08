# GeometricNeuronPlusField

Experimental branch of the Geometric Neuron / Functional Arbor line.

The question is now more precise than "does geometry compute?"

> **What does frozen anatomy contribute to a moving electrical field, where does task information become locally available, and what active boundary turns that analog field into output events?**

The source soma-tap experiment is in [`FunctionalArbors/SomaTapTestsClaude`](https://github.com/anttiluode/FunctionalArbors/tree/main/SomaTapTestsClaude). The amplitude-balance follow-up is in [`SomaWhyClaude/`](SomaWhyClaude/). Claude's independent HH dead-arm audit is preserved in [`HHablationTestClaude/`](HHablationTestClaude/) and summarized in [`AIS_CLAUDE_AUDIT_V01.md`](AIS_CLAUDE_AUDIT_V01.md).

## Current working picture

```text
slow anatomy / operator G
        |
        v
geometry-defined modal resonator bank
        |
        v
moving complex field psi(x,t)
        |
        v
convergence root / A-B amplitude-balance point
        |
        v
active state V,m,h,n
(history-dependent event availability)
        |
        v
spike events
```

The graph basis is a **microscope**, not a claim that biological somata calculate eigenvectors.

## What survived upstream

### 1. Temporal-order information is spectrally structured

On frozen FunctionalArbor bodies the energy-dominant spatial common mode is almost perfectly blind to A/B order, while a higher graph-Laplacian band (modes 18–20 in the registered setup) carries appreciable order information and gains selectivity from coherent projection.

The discovery-set band survived a held-out 12-body confirmation. See [`DISCOVERY_V01.md`](DISCOVERY_V01.md) and [`CONFIRMATION_V01.md`](CONFIRMATION_V01.md).

### 2. The distributed wave is almost exactly a geometry-defined resonator bank

A parameter-free reduced model

```text
q_n'' + damping q_n' + (restoring + stiffness*K*lambda_n) q_n
    = phi_n(A) s_A(t) + phi_n(B) s_B(t)
```

predicts the full-field modal contrast almost exactly (`pooled r ~= .995`). This is mathematically expected for the fixed-K graph wave, but conceptually useful: anatomy selects poles/time scales and input couplings of an LTI resonator bank. See [`MODAL_MECHANISM_V01.md`](MODAL_MECHANISM_V01.md).

### 3. A one-cell anatomical edit is global in modal coordinates

The modal-locality audit perturbs real frozen bodies one cell at a time. One added cell substantially disturbs about 24% of modal identities on average; the confirmed task band loses about 8% identity, and the eigenvector change is only weakly localized around the edit. See [`MODAL_LOCALITY_V01.md`](MODAL_LOCALITY_V01.md).

This makes local per-cell credit a badly aligned coordinate for the global computation. It does **not** prove that event-level local eligibility is impossible.

### 4. Locality compresses the global modes, but the soma region is a favorable local mixture

A compact patch anywhere on the tree is poor at reconstructing three global graph modes compared with scattered graph-wide taps. With the fair control — same-radius compact balls around every other body cell — the soma neighborhood is consistently better-conditioned than most equally local apertures, without containing unusually high task-band energy. See [`LOCAL_OBSERVABILITY_V02.md`](LOCAL_OBSERVABILITY_V02.md).

### 5. The soma/root task effect is largely an amplitude-balance effect

The original task-bottleneck result was strong: soma temporal-order contrast is much larger than the average body cell. Claude's `SomaWhyClaude` follow-up reproduced it and supplied a simpler mechanism.

Across its 12 frozen bodies, per-cell order selectivity correlated with single-source amplitude balance

```text
b(x) = min(p_A,p_B) / max(p_A,p_B)
```

at about `r = +0.70`, while total energy had essentially no explanatory power. The soma/root is therefore better described here as the **A/B amplitude-balance / coincidence point produced by the convergence construction**, not as an independently discovered magical soma location.

This yields a clean falsifiable prediction: move one source and the most selective point should track the new balance point rather than remain fixed at the designated soma.

## The AIS bridge was built

The old repo stopped at a passive power readout. We froze the upstream system and added a compact HH-like active boundary immediately downstream of the soma readout.

[`ais_active_probe.py`](ais_active_probe.py) contains state

```text
V, m, h, n
```

with Na activation/inactivation and K activation/recovery. No derivative or hand-written high-pass term is inserted.

It is compared with:

1. the same soma scalar under a rate-matched memoryless detector;
2. the **measured small-signal linearization of the same active membrane**, convolved with the same input and rate matched;
3. the full active stateful boundary.

See [`AIS_ACTIVE_PREREG_V01.md`](AIS_ACTIVE_PREREG_V01.md), [`AIS_ACTIVE_V01.md`](AIS_ACTIVE_V01.md), [`AIS_ACTIVE_PREREG_V02.md`](AIS_ACTIVE_PREREG_V02.md), and [`AIS_ACTIVE_V02.md`](AIS_ACTIVE_V02.md).

## Timing: the nonlinear boundary costs precision

The first run made the active boundary look spectacular at high frequency, but that was an estimator trap: sparse conditions can give `vector strength = 1` from one spike. We discarded that interpretation.

v0.2 froze the mechanism and strengthened the control with exact per-frequency event-count matching, a minimum of four active spikes, and PPC.

The original pooled upper-band receipt had 21 body/frequency pairs from 11 bodies:

```text
mean PPC
active        0.4305
linearized    0.5816
memoryless    0.5238

active - linearized   -0.1511
wins / losses          6 / 15
sign-test p             0.078
```

Claude correctly flagged that the 6/21 sign test is directional but not formally below `.05`. We therefore re-reduced the stored data to **one mean upper-band delta per body** instead of treating repeated frequencies from the same body as independent organisms:

```text
11 bodies
mean body delta active-linear   -0.1348
median                            -0.0849
active better / worse bodies       2 / 9
paired Wilcoxon, two-sided p       0.0322
```

The earned statement is stronger and cleaner than "no advantage": in this registered sample, nonlinear HH eventization **degrades timing precision relative to its own linearized response**. See [`AIS_CLAUDE_AUDIT_V01.md`](AIS_CLAUDE_AUDIT_V01.md).

## The instantaneous h/n ablation was a dead-arm experiment

The historical gate-memory experiment found:

```text
m instantaneous      fires 24/24
h instantaneous      fires  0/24
n instantaneous      fires  0/24
```

Claude independently swept standard HH constant drive over orders of magnitude and reproduced the qualitative asymmetry: `m`-instantaneous remains excitable, while instantaneous `h`/`n` destroys the repetitive-spiking regime.

That means the old h/n 0/24 result **cannot identify a frequency-selection mechanism**. The intervention destroyed the operating regime. The historical receipt is retained but explicitly superseded on this point. See [`HHablationTestClaude/`](HHablationTestClaude/), [`AIS_CLAUDE_AUDIT_V01.md`](AIS_CLAUDE_AUDIT_V01.md), and the corrected [`AIS_GATE_ABLATION_V01.md`](AIS_GATE_ABLATION_V01.md).

The m-instantaneous arm remains useful because it stays alive: most of the normalized frequency-allocation shape survives without independent fast-m history.

## Kinetics: n matters strongly, but not as a simple passband knob

[`ais_kinetics_probe.py`](ais_kinetics_probe.py) changed only `dh/dt` or `dn/dt` speed by `0.5x` or `2x`, preserving steady-state curves, conductances and gain.

A true gate-set passband was preregistered to move in opposite directions when the relevant gate was slowed versus sped up. It did **not**.

`n` kinetics strongly alter excitability and response shape (`TV ~= .18` slow, `.35` fast), but faster/slower `n` do not reliably translate the allocation center. `h` effects are smaller and likewise fail the registered ordering. See [`AIS_KINETICS_V01.md`](AIS_KINETICS_V01.md).

### Minimum-ISI test

Claude then proposed the clean refractory discriminator: if `n` kinetics set a refractory floor, minimum ISI should scale with the n time constant.

[`ais_n_isi_probe.py`](ais_n_isi_probe.py) reran the frozen kinetics battery and preregistered `f=0.025` as the primary condition.

```text
18 valid slow/fast bodies
mean log2(fast/slow minISI)   +0.0054
fast shorter / longer            9 / 9
sign p                           1.000
Wilcoxon fast<slow p             0.517
strict slow > full > fast         6 / 18
```

So the simple refractory-clock prediction also fails. Curiously, the native `n_scale=1` condition has the shortest mean minimum ISI; moving n either slower or faster tends to lengthen the short-ISI edge, while fast n also suppresses firing severely.

That looks more like an **excitability / dynamical-regime optimum** than either a translated passband or a single refractory knob. See [`AIS_N_ISI_PREREG_V01.md`](AIS_N_ISI_PREREG_V01.md) and [`AIS_N_ISI_V01.md`](AIS_N_ISI_V01.md).

## Final phase/interface test

The early AIS tests inherited the historical FunctionalArbor soma scalar:

```text
|psi_soma|^2
```

But the complex field still exists internally. `|psi|^2` and `|psi|` remove the carrier phase; only `Re(psi)` preserves a signed carrier waveform.

So before touching AIS geometry we preregistered a final interface test with all active parameters frozen:

```text
power       |psi_soma|^2
magnitude   |psi_soma|
real        Re(psi_soma)   <-- primary phase-bearing arm
```

The result is unusually clean. `Re(psi)` increases upper-band exposure to **23/24 bodies**, but the active boundary still loses to its own exact-rate-matched linearization:

```text
Re(psi), body-level upper-band delta PPC
valid bodies                 23
mean active-linear          -0.16435
median                      -0.14955
active better / worse         3 / 20
sign p                        0.000488
Wilcoxon two-sided p          0.0000131
```

Magnitude and power point in the same direction. Carrier-phase PPC itself also fails to reveal a general active advantage.

Therefore restoring the last phase-bearing scalar **does not rescue active timing**. In this lineage, phase is not currently justified as a special resource at eventization. See [`AIS_FINAL_PHASE_PREREG_V01.md`](AIS_FINAL_PHASE_PREREG_V01.md) and [`AIS_FINAL_PHASE_V01.md`](AIS_FINAL_PHASE_V01.md).

## What is actually left of the AIS story?

A narrower hypothesis survives:

> **The active boundary may buy nonlinear event selection at a cost in timing precision.**

That is not yet a positive computational result. The next fair test is to ask whether, under a matched total spike budget, the active spike pattern carries more information about the input frequency/regime than the memoryless and own-linearization controls.

Until that is demonstrated, AIS position/length/extent co-adaptation remains blocked.

## Operator -> flow -> event

The current synthesis remains useful even with the active-boundary nulls:

```text
operator / morphology G      slow constraints on possible dynamics
        |
        v
flow / complex field psi(t)  continuous execution of those dynamics
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
- [`ais_gate_ablation_probe.py`](ais_gate_ablation_probe.py) — historical instantaneous-gate experiment; h/n arms now classified as dead arms.
- [`ais_kinetics_probe.py`](ais_kinetics_probe.py) — h/n kinetic-speed test.
- [`ais_n_isi_probe.py`](ais_n_isi_probe.py) — preregistered n-kinetics minimum-ISI mechanism test.
- [`ais_interface_phase_probe.py`](ais_interface_phase_probe.py) — final `|psi|^2` / `|psi|` / `Re(psi)` phase-interface test.

## Live-field / settled-readout result

The original strict criterion remains a null. The field stayed live, but zero modes passed the prewritten combination of stability, energy and selectivity thresholds. Those thresholds were not relaxed after seeing the data.

## What this is not

This is an experimental computational model, not a claim that real dendrites are literal graph-Laplacian computers or that the AIS exists for one single evolutionary reason. Real AIS biology also includes maintenance of axonal polarity and trafficking barriers, among other functions.

The useful target is narrower: identify which pieces of a geometry/field/event architecture are actually necessary in the model, then compare that division of labor with known neuronal compartments without retrofitting biology to the model.
