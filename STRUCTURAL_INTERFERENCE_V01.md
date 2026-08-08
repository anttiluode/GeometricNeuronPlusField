# Structural interference v0.1 — one-cell edits act mainly through transfer-history compatibility

The confirmed local computation from `INTERFERENCE_FACTORIZATION_V01.md` is

```text
geometry
  -> h_A, h_B
  -> interference visibility V
  -> lag-directional complex compatibility
  -> coherent temporal-order contrast.
```

This experiment asks what a **single local structural edit** actually changes in those coordinates.

## Counterfactual event design

On frozen mature v0.9 bodies, with soma and source terminals fixed, two legal one-cell event classes were tested:

```text
add     empty cell with exactly one occupied neighbor
        (a tip-like one-bond extension)

delete  occupied non-source/non-soma cell whose removal
        preserves connectivity and the exact source terminals
```

Up to six additions and six deletions were sampled per body. Every event was evaluated as a before/after frozen-body counterfactual; no growth or credit dynamics were active.

Discovery used fresh seeds 156-167. Held-out confirmation used fresh seeds 168-179.

## Two-factor exact decomposition

At the soma, the integrated interference statistic is

```text
C_int = f(V,R)
R = (Re rho_plus, Re rho_minus).
```

For each event, the change was split by a two-factor Shapley decomposition:

```text
contrib_V + contrib_R = Delta C_int
```

where `contrib_V` is the part of the change mediated through amplitude visibility and `contrib_R` is the part mediated through the lagged complex relation of the source histories.

The accounting identity held to machine precision (`max error < 1e-17`).

## Discovery — useful mixed result

The original broad prediction that `Delta C_int` would precisely track the historical peak-power task did not meet its preregistered body-level threshold:

```text
mean body corr(dC_int, dC_peak)             0.4276
positive bodies                              9 / 12
D1                                           FAIL
```

But event direction was already useful:

```text
pooled sign agreement                        0.7048
D2                                           PASS
```

The preregistered claim that full interference change would beat visibility change in most bodies also failed its count/sign criterion even though the average improvement was large:

```text
mean corr improvement over dV               +0.3451
improved / worse bodies                       7 / 5
sign p                                        0.774
D3                                           FAIL
```

The unexpected mechanism result was much cleaner:

```text
mean |contrib_V|                             0.00128
mean |contrib_R|                             0.01506
absolute compatibility / visibility ratio    11.76
body-mean compatibility share                 0.894
bodies with compatibility share > .80         10 / 12
events with |R| > |V|                        134 / 144
```

That dominance was not part of the original discovery success criteria, so it was preregistered separately before new bodies were run.

## Held-out confirmation — seeds 168-179

### C1 — compatibility mediation dominates

```text
mean body compatibility share                0.93052
median                                        0.97733
bodies with share > .75                       11 / 12
registered threshold mean > .82, >=9 bodies
```

**C1 PASS.**

### C2 — dominance is event-by-event

```text
events with |contrib_R| > |contrib_V|        138 / 144
fraction                                      0.95833
registered threshold                         > .85
```

**C2 PASS.**

### C3 — the integrated change retains task direction

Among 108 events with `|dC_peak| > 1e-5`:

```text
sign(dC_int) == sign(dC_peak)                 0.80556
registered threshold                         > .65
```

**C3 PASS.**

### C4 — Shapley accounting

```text
max absolute reconstruction error            6.94e-18
registered threshold                         < 1e-12
```

**C4 PASS.**

The held-out set was stronger than discovery even on the descriptive magnitude relation:

```text
mean body corr(dC_int, dC_peak)               0.65219
median                                        0.68431
positive bodies                               11 / 12
pooled corr(dC_int, dC_peak)                  0.72648

mean improvement over dV                     +0.65625
improved / worse                              10 / 2
sign p                                        0.03857
```

These latter values are descriptive with respect to the confirmation preregistration; the registered confirmation centered on mediation dominance and sign relevance.

## Additions and deletions agree

The result is not carried by only one event class.

```text
                       add            delete
n                       72                72
mean |dC_peak|        .04852            .02966
mean |dC_int|         .02004            .01201
mean |contrib_R|      .02048            .01223
mean |contrib_V|      .00114            .00062
sign agreement         .817              .792
pooled corr             .726              .729
```

So both adding and removing one local piece mainly alter the order computation by reshaping the **relative transfer history**, not by simply making the two sources more or less equally strong at the root.

## Distance matters descriptively

On held-out events, farther edits tended to have smaller effects:

```text
corr(distance from soma, |dC_peak|)          -0.574
corr(distance from soma, |dC_int|)           -0.492
corr(distance from soma, |contrib_R|)        -0.508
corr(distance from soma, |contrib_V|)        -0.103
```

This is consistent with a local edit perturbing the globally shaped transfer relation with decreasing leverage at the readout, but it is not by itself a credit rule.

## What this resolves

Earlier experiments showed:

```text
WHY the soma/root is a good location
    -> convergence gives high A/B amplitude visibility.

WHAT carries temporal-order value there
    -> the lag-directional coherent relation of h_A and h_B.

WHAT one-cell structural edits mainly change
    -> that directional compatibility, not visibility.
```

That distinction is useful. Amplitude balance is primarily a **location/opportunity** variable. Once a convergence point already exists, structural tuning of the computation occurs mainly by changing the detailed timing/complex relationship between the two transfer histories.

This also explains why a growth cue based only on field magnitude or recent activity can miss the real causal variable: most of the value of a structural event is not a large change in local amplitude opportunity.

## Important limitation

This is still an omniscient before/after counterfactual. To know `contrib_R` exactly, we recompute the whole transfer histories after the edit.

So it identifies the right **causal target** more clearly, but it does not solve biological/local credit assignment.

## Wall sentence

> **Convergence determines where the two inputs can meet, but one-cell anatomical edits tune the computation mainly by changing how their geometry-shaped histories meet in time and complex phase. In held-out bodies, roughly 93% of the absolute change in the integrated interference computation was mediated through transfer-history compatibility rather than amplitude visibility.**

## Next clean question

The next wall is now precise:

> **Can the effect of a local structural event on the soma interference objective be predicted from a local sensitivity signal, without actually building the event and recomputing the whole arbor?**

For a reciprocal linear wave operator this is an adjoint problem. A task-conditioned backward field from the soma should, in principle, convert the global consequence into a local sensitivity through its overlap with the forward field. The next experiment should test that prediction against the exact one-cell counterfactuals before treating any biological analogue seriously.
