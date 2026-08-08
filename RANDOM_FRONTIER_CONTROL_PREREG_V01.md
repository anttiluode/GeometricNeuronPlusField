# Random frontier control preregistration v0.1

`BOND_RESPONSE_V01.md` found that 80.5% of **gradient-favored** held-out additions had an interior optimum on the bath-to-arbor conductance path, with median interior `alpha_best=0.10`.

That conditioning can itself enrich interior optima: selecting an event by a positive derivative at `alpha=0` guarantees that the curve initially rises, while the derivative is only local. This control removes that conditioning.

## Frozen protocol

Fresh FunctionalArbor bodies: seeds **276-287**.

For each body, use the existing structural-event candidate sampler to choose up to 6 legal tip-like additions before inspecting any gradient. Then evaluate every sampled addition along exactly the same conductance grid used in `BOND_RESPONSE_V01.md`:

```text
alpha = 0,
.005, .01, .02, .03, .05, .075,
.10, .15, .20, .30, .40, .50, .60,
.75, 1.0
```

The base adjoint derivative is recorded only after random selection, for descriptive positive/negative-gradient stratification. No gradient is recomputed along a response curve.

## Primary estimands

For all randomly selected additions pooled across bodies:

```text
interior fraction = P(0 < alpha_best < 1)
alpha_best = 0 fraction
alpha_best = 1 fraction
median interior alpha_best
slope-sign-reversal fraction
mean binary regret = best partial gain - full-binary gain
```

## Fixed classification

### BROAD_GRADED

If random interior fraction is at least 0.70, the broad graded-coupling interpretation survives removal of gradient-sign conditioning.

### SELECTION_DOMINATED

If random interior fraction is at most 0.55 **and** the `alpha_best=0` fraction is at least 0.30, the earlier ~80% interior enrichment is substantially selection-conditioned. The valid statement becomes narrower: among additions locally favored at the base state, partial strengthening is often preferable to forcing the binary endpoint.

### MIXED

Anything between these regions is reported as mixed evidence. Thresholds are fixed before fresh seeds 276-287.

For context only, the same statistics are reported for the post-selection positive- and negative-derivative subsets, including `favored interior fraction - random interior fraction`.

This control does not test whether continuous conductance can be optimized. It tests only whether the generalized impedance-field interpretation survives the missing selection control.
