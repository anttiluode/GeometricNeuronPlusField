# Source-gain balance causality v0.1 — preregistration

## Why this test follows the source-move result

The source-move experiment gave a replicated positive result for map relocation: when one injection site was moved inward on a frozen body, changes in the single-source amplitude-balance map predicted changes in the temporal-order selectivity map. But the strongest root-displacement criterion failed; the designated soma often remained highly balanced and highly selective even when another cell had the absolute maximum balance ratio.

Source relocation changes **two things at once**:

1. relative amplitudes;
2. path/transfer timing.

The next experiment therefore changes amplitude only.

## Intervention

Freeze fresh v0.9 bootstrap bodies. Keep anatomy, source locations, carrier, lag, wave parameters, and soma fixed. Multiply only one source's pulse amplitude.

Five conditions:

```text
baseline      A=1.0  B=1.0
A_half        A=0.5  B=1.0
A_double      A=2.0  B=1.0
B_half        A=1.0  B=0.5
B_double      A=1.0  B=2.0
```

No parameter is retuned after observing the result.

For every condition and every body cell measure single-source peak powers `p_A`, `p_B`, amplitude balance

```text
b(x) = min(p_A,p_B) / max(p_A,p_B)
```

and temporal-order selectivity `|C(x)|` for A->B versus B->A at lag 20.

## Primary P1 — same-location causal balance test at the soma

For each body, subtract the baseline soma values from the four gain-perturbed conditions and compute

```text
r_soma = corr( delta b_soma, delta |C_soma| ).
```

This holds geometry and location fixed and asks whether experimentally making the two source contributions more or less balanced at the *same soma* changes its order selectivity in the same direction.

Registered prediction:

```text
mean r_soma > 0
and body-level two-sided sign test p < .05.
```

## Primary P2 — whole-map gain causality

For every non-baseline gain condition compute across body cells

```text
r_delta_map = corr( b_gain - b_baseline,
                    |C|_gain - |C|_baseline ).
```

Average the four perturbation values within each body.

Registered prediction:

```text
mean body r_delta_map > 0
and body-level two-sided sign test p < .05.
```

This is the gain-only counterpart of the already positive source-move relocation test.

## Secondary receipts

Report without promoting them to primary criteria:

- within-body Spearman/Pearson relation between absolute soma balance and soma selectivity across all five gains;
- which gain maximizes soma balance and which maximizes soma selectivity;
- soma balance/selectivity percentiles across the body;
- per-cell correlation `corr(b, |C|)` for each gain;
- top-decile overlap of the two maps.

## Fresh discovery / confirmation sets

The earlier source-move test used seeds 0-23. This experiment uses unseen bodies:

```text
discovery      seeds 24-35
confirmation   seeds 36-47
```

The same frozen protocol and criteria are used for both sets.

## Interpretation

If P1 and P2 pass on both sets, amplitude balance has earned a causal role rather than merely being a spatial correlate of the designated root. If map relocation passes but P1 fails, balance is probably a proxy for a transfer-function property that source movement changes. If both fail, the source-move result should be reinterpreted as path/timing rather than amplitude balance.
