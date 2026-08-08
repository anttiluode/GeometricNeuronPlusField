# Structural-event interference held-out confirmation preregistration v0.1

Discovery on fresh seeds 156-167 gave a mixed result for predicting the historical peak objective but a strong mechanistic result inside the confirmed interference statistic.

Discovery:

```text
mean body corr(dC_int, dC_peak)             0.4276   original D1 FAIL
positive bodies                              9/12
pooled sign agreement                       0.7048   D2 PASS

mean improvement over dV                    0.3451
improved bodies                              7/12   original D3 FAIL

mean |Shapley visibility contribution|      0.00128
mean |Shapley compatibility contribution|   0.01506
compatibility / visibility absolute ratio   11.76
compatibility share, body mean              0.894
bodies with compatibility share > .80       10/12
individual events |R| > |V|                134/144
```

The refined hypothesis is therefore not that `dC_int` is already a complete structural credit signal for the peak objective. It is narrower:

> **A one-cell structural edit changes the integrated interference computation primarily by changing the lag-directional complex relationship of the two source histories, not by changing their amplitude visibility.**

No discovery body is reused below.

## Held-out bodies

```text
seeds 168-179
```

Same event generator, max six legal additions + six legal deletions per body, source-terminal protection, lag 20, 210 steps, and exact two-factor Shapley decomposition from `structural_interference_probe.py`.

## Registered confirmation criteria

### C1 — compatibility mediation dominates visibility mediation

For each body define

```text
Rshare = sum_e |contrib_R(e)| /
         sum_e (|contrib_R(e)| + |contrib_V(e)|).
```

PASS if:

```text
mean Rshare > 0.82
and at least 9/12 bodies have Rshare > 0.75.
```

### C2 — dominance holds event by event, not only after averaging

PASS if more than `85%` of all evaluated held-out events satisfy

```text
|contrib_R| > |contrib_V|.
```

### C3 — the integrated interference change retains directional relevance to the peak task

Among events with `|dC_peak| > 1e-5`, PASS if pooled sign agreement between `dC_int` and `dC_peak` is `> 0.65`.

This is intentionally weaker than claiming precise magnitude prediction; discovery did not earn that claim.

### C4 — exact Shapley accounting remains numerical identity

PASS if maximum absolute reconstruction error is `<1e-12`.

## Descriptive only

Report without additional pass/fail thresholds:

- mean/median body correlation `corr(dC_int,dC_peak)` and number positive;
- `corr(contrib_R,dC_peak)`;
- additions versus deletions;
- effect magnitude versus graph distance from soma;
- mean absolute `dC_peak`, `dC_int`, `contrib_R`, and `contrib_V`.

## Interpretation fixed in advance

- C1+C2 pass: structural edits alter the confirmed interference statistic mainly through **transfer-history compatibility** rather than amplitude opportunity.
- C3 pass as well: that compatibility-mediated quantity retains useful directional relation to the historical task and is a plausible target for future eligibility approximations.
- C3 fail: compatibility is mechanistically dominant for `C_int` but not sufficient as a credit coordinate for the peak objective.

No claim of biological locality follows from this counterfactual probe; it explicitly compares before and after whole-arbor transfer functions.
