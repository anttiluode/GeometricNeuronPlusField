# Biological rhythm boundary v0.1 — what the lock-in result does and does not suggest

The engineering branch has now produced a precise result: a fast balanced phase reference can demodulate a local forward/adjoint cross-term from intensity, and its exactness is controlled by spectral overlap with the retained field's self-energy difference spectrum.

That result makes theta/gamma/AIS biology tempting to reinterpret. The temptation needs a hard boundary.

## What is established in biology

### 1. The AIS is a consequence gate, not a generic phase inverter

Axo-axonic/chandelier cells target the axon initial segment and regulate action-potential initiation. CA3 axo-axonic cells fire rhythmically at specific theta phases and are gamma-coupled in vivo (Viney et al., *Nature Neuroscience*, 2013, DOI `10.1038/nn.3550`).

But GABA at the AIS should **not** be modeled as an automatic `pi` phase flip. Modern measurements support a strong inhibitory/shunting control over spike initiation, threshold, and onset timing, even in developmental conditions where the local GABA reversal potential can be relatively depolarized (Lipkin & Bender, *Journal of Neuroscience*, 2023; PMCID `PMC10500977`).

So the direct analogy

```text
chandelier firing = multiply returned field by -1
```

is rejected.

A more defensible engineering analogue is

```text
AIS inhibition = gate / delay / veto / retime the consequential output event.
```

### 2. Theta phase separates physiological operations

In CA1, spike excitability and dendritic LTP are maximal at different phases of theta, and basal versus apical dendritic plasticity has different preferred phase structure (Leung & colleagues, *eNeuro*, 2018, `ENEURO.0236-18.2018`).

This supports theta as a **temporal context / multiplexing variable**. It does not establish forward-pass and backward-pass phases in the backpropagation sense.

### 3. Medial septum coordinates multiple timescales

The medial septum is not merely a slow 4–12 Hz metronome. Király et al. (*Nature Communications*, 2023, DOI `10.1038/s41467-023-41746-0`) found medial-septal firing phase-coupled to theta-nested beta/gamma components and showed that optogenetic activation of medial-septal PV neurons can elicit theta-nested beta/gamma activity in CA1.

So the biological system genuinely has a hierarchy

```text
slow phase context
    containing / organizing
faster rhythmic packets
```

which is at least architecturally compatible with temporal multiplexing and demodulation ideas.

### 4. A neuron already has a natural 'return event'

The closest biological object to the model's returned credit field is not chandelier inhibition. It is the neuron's own somatic/axonal consequence propagating back into the dendrites.

Back-propagating action potentials interact locally with EPSPs and dendritic oscillations in a narrow timing window and can strongly amplify dendritic depolarization; this coincidence has long been linked to synaptic plasticity (Stuart & Häusser, *Nature Neuroscience*, 2001, DOI `10.1038/82910`).

Biology therefore already contains the structural motif

```text
local forward synaptic history
            x
returning consequence-related dendritic event
            |
            v
local plasticity
```

That is genuinely reminiscent of the model's local forward-field x return-field overlap.

It is **not** evidence that the returning biological event is the mathematical transpose/adjoint field.

### 5. The single-arbor reciprocity issue is subtler than the network reciprocity issue

An earlier version of this note leaned too hard on the fact that chemical synapses are nonreciprocal. That is true for a network, but it is not the decisive objection for this repository because the optimized geometry is a **single dendritic arbor**.

For a linear passive dendritic cable, transfer impedance is reciprocal:

```text
Z(i,j,omega) = Z(j,i,omega).
```

This is standard cable theory and remains true for arbitrarily branching passive models; passive dendritic reductions explicitly use this symmetry. See, for example, Major & Evans (*Biophysical Journal*, 1994, DOI `10.1016/S0006-3495(94)80836-7`) and the later Neuron_Reduce formulation (Amsalem et al., *Nature Communications*, 2020, DOI `10.1038/s41467-019-13932-6`).

That means the spatial statement

```text
forward subthreshold transfer through a passive tree
and
reverse subthreshold transfer through the same tree
```

really can share a reciprocal operator.

This is unexpectedly close to the engineering condition that made the physical adjoint cheap.

The mismatch moves elsewhere:

- real dendrites are not purely passive;
- channel densities are highly nonuniform;
- back-propagating action potentials are regenerative nonlinear events, not small reciprocal test waves;
- dendritic spikes can arise locally;
- plasticity modifies synapses/channels and morphology, not merely one symmetric linear coupling coefficient;
- network-level chemical synapses are directional.

So passive cable reciprocity **reopens the single-neuron analogy**, but does not give the biological neuron an exact adjoint for free.

## The key mismatch remains — but it is now correctly located

The exact physical gradient in this repository depends on transpose propagation. Reciprocity is useful because it lets the same physical spatial operator realize that transpose.

For the **passive/subthreshold part of one dendritic tree**, biology actually possesses the relevant reciprocal transfer symmetry. For the active return event and the larger circuit, it generally does not preserve the same linear operator.

So the unsupported statement remains

```text
brain implements exact in-situ adjoint backpropagation
```

but the reason is no longer simply "synapses are one-way."

A better question is:

> **How much of the useful transpose-like spatial credit pattern survives when a reciprocal passive dendritic substrate is driven backward by a biologically coarse, active return event rather than by the exact analog adjoint waveform?**

That is directly testable in this repository.

## What the new result does suggest

The confirmed lock-in experiment changes the biological question.

We no longer need to ask whether gamma literally changes `V` into `-V`.

The narrower question is:

> **Can nested rhythmic gating make a weak local relation between a recent input history and a returning consequence event observable against much larger self-energy/background activity?**

That is a demodulation question.

The model gives a concrete design principle:

```text
use a fast reference whose spectral support avoids
background/self-energy difference frequencies
```

and the desired cross-term survives while background cancels.

The hippocampal literature independently shows:

- theta-phase-dependent plasticity;
- theta-phase-dependent spike excitability;
- theta-nested gamma bands;
- AIS-targeting interneurons with theta/gamma phase structure;
- local dendritic coincidence between synaptic input and back-propagating spikes;
- reciprocal passive transfer inside a dendritic cable approximation.

Those facts make **rhythmic coincidence gating** worth testing. They do not identify the exact signal being demodulated.

## A non-inverting protocol is the better bridge

The engineering machine does not actually require a `pi` inverter if separate measurement windows are allowed.

For local forward and return signals `u` and `v`, three non-inverting intensity states are enough:

```text
P_u   = |u|^2
P_v   = |v|^2
P_uv  = |u+v|^2
```

Then

```text
(P_uv - P_u - P_v) / 2 = Re[conj(u)v].
```

This is established interference algebra, not a biological claim.

But it suggests a much better rhythm analogy than the earlier `+V/-V` story:

```text
one phase/window: local forward activity
one phase/window: return/consequence activity
one phase/window: overlap permitted
slow local chemistry: subtract / normalize baselines over time
```

The exact biological implementation, if any, is unknown.

## The experiment this commits us to

Do **not** build a biological model that is guaranteed to reproduce the adjoint.

Instead compare increasingly biological return codes against the exact gradient:

```text
A. exact transpose waveform                         positive control
B. three-state non-inverting intensity protocol    physics control
C. phase-windowed / duty-gated return waveform     rhythm test
D. fixed-amplitude sparse return events            bAP-like test
E. AIS-delayed / vetoed sparse return events       axo-axonic timing test
F. local forward x return coincidence only         biological extreme
```

For C–F, score:

- bond-gradient correlation with the exact adjoint;
- strong-sign agreement;
- closed-loop learning gain;
- robustness to timing jitter;
- dependence on slow/fast scale separation.

The important possible outcome is a failure. If a bAP-like event code destroys the gradient, the analogy stops at 'local coincidence'. If a surprisingly coarse phase-windowed code retains useful structural direction, then the next question becomes what information that coarse return is preserving.

## Current wall sentence

> **The engineering result does not predict that chandelier cells invert a credit wave. Passive dendritic geometry really is reciprocal enough to make the single-neuron analogy worth reopening; the hard question is whether a coarse biological return event plus oscillatory gating preserves any useful fraction of the transpose-like local credit pattern.**
