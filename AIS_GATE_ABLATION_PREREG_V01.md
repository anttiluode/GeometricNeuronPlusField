# AIS gating-memory ablation v0.1 — what creates the nonlinear frequency selection?

`AIS_ACTIVE_V02` killed the stronger claim that the first active boundary improves rate-matched spike-time precision beyond its own linearization.

One result remains real and unexplained: under a fixed operating point the full active boundary redistributes its output events strongly across input modulation frequency.  It is not a monotonic high-pass; it suppresses one regime, opens another, and then closes again.

Before changing AIS geometry, identify the internal state responsible for that frequency-selection shape.

## Frozen upstream and active parameters

This experiment reuses **exactly** the v0.1/v0.2:

- FunctionalArbor bootstrap and mature field physics;
- soma power readout;
- normalization battery/rule;
- input gain;
- HH conductances and reversal potentials;
- modulation frequencies and burn window.

No parameter is retuned for an ablation.

## Ablations

The full model evolves

```text
V, m, h, n
```

with differential gating state.

The following variants remove *memory* from one gate by replacing that gate after each voltage update with its instantaneous steady-state value at the current voltage:

```text
h_instant     Na inactivation has no independent history
n_instant     K activation/recovery has no independent history
m_instant     Na activation has no independent history
all_instant   m,h,n all have no independent gating history
```

The voltage/capacitance state remains in every model.  This is therefore an ablation of **gating memory**, not a claim to make the membrane completely memoryless.

Using `x_inf(V)` rather than deleting a conductance keeps the gate's voltage dependence and avoids replacing the model with a different ion-current inventory.

## Primary measurement: frequency-allocation shape

For each body and model, count natural spike events after burn at every registered modulation frequency.

Turn that count vector into a normalized allocation profile:

```text
p_f = N_f / sum_f N_f
```

provided the model emitted at least one event.

Compare each ablation with the full model by total-variation distance

```text
TV = 0.5 * sum_f |p_f(ablation) - p_f(full)|
```

This removes the first-order confound that an ablation may simply make the membrane more or less excitable overall.

Also report:

- total event count;
- event share in the previously interesting `0.05 <= f < 0.125` band;
- peak-allocation frequency;
- peak-frequency agreement with the full model.

## Interpretation rules

A gate earns a role in the observed frequency selection if removing its independent history changes the **normalized frequency-allocation shape**, not merely the total spike count.

The strongest candidate is the single-gate ablation with the largest organism-stable TV distance from the full model.

`all_instant` is a positive-control-style test for whether independent gating histories matter collectively.

## Failure / caution conditions

- If all single-gate TV distances are small and `all_instant` also preserves the full profile, the v0.1 frequency selection does not depend materially on gating memory; it is mainly voltage nonlinearity/capacitance or upstream drive.
- If an ablation produces zero events in most bodies, its shape is not identifiable; report the exposure failure instead of assigning that gate a frequency role.
- If a large TV distance is driven only by a few bodies, report heterogeneity.
- Do not tune input gain separately for ablations after seeing the result.
- This experiment does not resurrect the failed spike-time-precision claim.

## Decision after this test

Only if a specific active state has an identifiable, robust frequency-selection role do we have a well-defined active quantity to later co-adapt with AIS position/extent.  Otherwise AIS geometry co-adaptation remains premature.
