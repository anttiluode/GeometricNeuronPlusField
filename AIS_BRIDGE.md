# AIS bridge — from a passive geometry/field bank to an active output boundary

The graph-mode work has now made one limitation of the current toy explicit.

A mature FunctionalArbor is almost exactly a fixed-K linear wave equation on a graph.  In the graph Laplacian basis it becomes a bank of independent damped oscillators, followed by a readout.  That is useful temporal filtering, but it is still reciprocal LTI dynamics until the readout/nonlinearity is added.

Claude's modal-locality audit adds a second constraint: a one-cell structural change is not a one-mode/local event.  On the 12-body x 25-perturbation receipt, one added cell disturbs about a quarter of the spectrum and the eigenvector change is only weakly enriched near the changed cell.  So a local structural tag is a poor coordinate for a consequence that lives in global modes.

This makes the axon initial segment (AIS) a particularly interesting biological comparison, because it sits exactly at the boundary our toy currently lacks:

```text
many distributed somatodendritic degrees of freedom
                    |
                    v
        local axon initial segment
      active channels + internal state
                    |
                    v
       temporally precise spike output
                    |
                    v
             directional axon
```

## What biology already says

This is **not** a claim that the toy explains the AIS.  The comparison is useful because several experimentally established AIS properties line up with the computational role that is missing from the toy.

- Leterrier's review describes the AIS as the compartment that generates and shapes the action potential before axonal propagation, with very high Nav channel concentration and multiple Kv channel complexes.  Its composition, length and position are regulated and plastic.
- Hamada et al. (PNAS 2016, DOI 10.1073/pnas.1607548113) showed that AIS distance covaries with dendritic morphology and that this covariation can normalize the somatic action potential.  The output boundary is therefore not independent of the input geometry.
- Lazarov et al. (Sci Adv 2018, DOI 10.1126/sciadv.aau8621) found that high AIS Nav density is not simply required to make *any* axonal spike; its loss strongly reduces the bandwidth and temporal precision of spike encoding.  This is especially relevant here: the AIS is not merely a threshold location, it helps determine the frequency transfer from analog input fluctuations to spike timing.
- In auditory neurons, AIS length and/or location varies systematically with characteristic sound frequency (reviewed by Kuba, 2012).  That is a direct biological example of output-segment geometry being tuned to a temporal/frequency regime.

A later modeling literature also cautions against a one-variable story: high-frequency encoding depends strongly on active channel kinetics/voltage sensitivity, and not only on AIS distance or channel count.  So the hypothesis here is **geometry + active state**, not geometry alone.

## The working hypothesis

The dendritic/somatic body and the AIS may solve different halves of the same problem.

```text
frozen morphology G_D
       |
       +--> defines modal coordinates and time constants
       |
       v
moving analog field psi(x,t)
       |
       v
local AIS boundary G_AIS + active channel state h(t)
       |
       +--> frequency/temporal selection
       +--> thresholding
       +--> history dependence / refractory state
       |
       v
discrete directional spike train
```

The important distinction is that the AIS is not another passive soma tap.

A passive aperture can project the existing field.  An AIS-like compartment can also change state.  That is the minimal ingredient needed to escape the LTI-filter-bank ceiling of the current FunctionalArbor.

In that sense the architecture may be:

> **geometry defines what analog motions are possible; the AIS is an active spectral boundary that decides which of those motions become irreversible output events.**

"Irreversible" here means computationally/history-dependent, not a thermodynamic claim.

## Why put the active machinery there?

The topology offers a plausible economy argument.

The somatodendritic field is distributed and its modal coordinates are global.  Making every dendritic compartment equally active/history-dependent would be expensive and would make credit/state assignment even harder.  A single active segment placed **after convergence of the dendritic computation but before the long output cable** can exert global control with local machinery.

That is exactly the architectural location of the AIS.

This repo can test the computational side of that statement without pretending to prove the evolutionary reason for the AIS.

## Experiments now in line

### A. Local observability of the confirmed task band

`local_observability_probe.py` asks whether the confirmed global modes 18-20 are already recoverable from small graph-distance balls around the soma, compared with equal-size random apertures on the same body.

If the soma region is unusually good, passive geometry itself funnels the task band into a local interface.

If it is ordinary, that is equally useful: the AIS-like active boundary must create/select the useful temporal code rather than merely sitting at a privileged passive projection.

### B. Active-boundary probe

The next model should then keep the entire arbor frozen and add only a reduced AIS-like stateful encoder downstream of a local voltage trace.  It should compare at least:

1. passive quadratic/point readout;
2. memoryless threshold control;
3. active stateful spike-initiation gate with fast activation and recovery/inactivation;
4. matched low-bandwidth/low-Nav-like control.

The questions are not "can a threshold make spikes?" but:

- Does active state broaden or reshape the frequency transfer from the arbor field to spike timing?
- Does it create history dependence that cannot be diagonalized away with the body modes?
- Can AIS position/extent or active gain compensate differences between frozen body geometries, analogous to the dendrite/AIS covariation reported by Hamada et al.?
- Does a parameter set that improves bandwidth also sharpen temporal-order reporting, or are those separate functions?

### C. Body–AIS co-adaptation

Only after A and B should we let the two geometries move.

```text
body spectrum / modal couplings  <---->  AIS location, extent, kinetics
                 |                         |
                 +------ spike code -------+
```

The clean target would be output normalization across heterogeneous bodies, not simply maximizing excitability.

## Failure conditions

Keep this branch killable.

- If small soma-local apertures are no better than matched random apertures, do not call the soma a passive modal focus.
- If a stateful active gate does no more than a rate-matched memoryless threshold, the proposed AIS bridge has not earned a role.
- If apparent bandwidth changes are entirely firing-rate confounds, calibrate/match rate before interpreting them.
- If changing AIS-like position or extent has the same effect on every body, there is no body/readout co-adaptation story in this toy.
- If an active gate only works because a derivative/high-pass term was inserted by hand, that does not count as an emergent frequency filter.

## Wall sentence

The current hypothesis is narrower than "the AIS exists because of geometry":

> **A distributed passive geometry can implement an analog resonator bank, but a local active, plastic boundary is a natural place to turn that global moving field into a high-bandwidth, history-dependent, directional output code.  The biological AIS has exactly those active and plastic properties; the next toy experiments ask whether that architectural division of labor is computationally necessary here.**
