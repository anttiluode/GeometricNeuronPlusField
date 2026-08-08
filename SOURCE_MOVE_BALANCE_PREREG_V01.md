# Source-move balance-point test v0.1 — preregistration

## Question

`SomaWhyClaude` found that the designated soma/root is highly order-selective but that per-cell selectivity is strongly associated with single-source amplitude balance

```text
b(x) = min(p_A(x), p_B(x)) / max(p_A(x), p_B(x)).
```

The current explanation is therefore not "the soma is intrinsically special" but:

> the grown convergence geometry places the soma near an A/B amplitude-balance point, and temporal-order selectivity is large where the two source contributions are comparable.

That explanation has not yet been tested causally.

## Intervention

Freeze each v0.9 bootstrap body completely. Do not regrow, prune, retune wave physics, or change the designated soma.

Keep one source at its historical terminal and move the other source *inward along its already-existing source-to-soma path*. Use the same pulse shape, amplitude, carrier, task lag, and wave dynamics. Test the following seven source configurations:

```text
original
B at 0.75, 0.50, 0.25 of its original graph distance from soma
A at 0.75, 0.50, 0.25 of its original graph distance from soma
```

Fraction 1.0 is the historical terminal and fraction 0.0 would be the soma. Only the injection site changes; anatomy does not.

For every condition and every occupied cell measure:

```text
p_A(x)       max single-source power from A
p_B(x)       max single-source power from B
b(x)         amplitude balance min/max
|C(x)|       absolute A->B versus B->A order contrast at lag 20
```

The original designated soma remains physically present and unchanged throughout.

## Primary causal tests

### P1 — map relocation

For each moved-source condition compute, over all body cells,

```text
r_delta = corr( b_moved - b_original,
                |C|_moved - |C|_original ).
```

Average the six `r_delta` values within each body before any organism-level test.

Registered prediction:

```text
mean body relocation correlation > 0
and two-sided sign test of body means has p < .05.
```

This asks whether *changes in the balance landscape* predict *changes in the selectivity landscape* within the exact same frozen anatomy.

### P2 — balance point beats fixed soma when balance actually moves

For each moved condition, identify the cell with maximal `b(x)`. If that cell is at least 3 graph edges from the designated soma, call the condition `displaced`.

On displaced conditions compare

```text
|C|(best-balance cell)  versus  |C|(designated soma).
```

Reduce to one mean difference per body before significance testing.

Registered prediction:

```text
mean difference > 0
and two-sided paired/sign test across bodies p < .05.
```

If this passes, the fixed root loses its special status when the amplitude-balance point is experimentally moved away from it.

## Secondary receipts

For every condition report:

- correlation `corr(b, |C|)` across cells;
- percentile of the designated soma in `b` and `|C|`;
- graph distance from soma to the best-balance and best-selectivity cells;
- graph distance between those two maxima;
- top-decile overlap between balance and selectivity maps.

These are descriptive unless explicitly promoted in a later preregistration.

## Discovery / confirmation split

Run the exact same frozen protocol on:

```text
discovery     seeds 0-11
confirmation  seeds 12-23
```

The code and criteria are frozen before either canonical run. No parameters or source fractions may be changed between sets.

## Failure conditions

The balance-point explanation is weakened if source relocation changes the balance landscape but selectivity stays pinned to the designated soma, if `r_delta` is not reliably positive, or if displaced balance maxima do not outperform the soma.

A null is not to be rescued by changing fractions, lag, pulse amplitude, carrier, body size, or correlation metric after seeing the data.

## Why this is next

The generic HH/AIS branch is frozen after failing timing, phase, refractory-clock, and matched-budget frequency-information tests. This source-move experiment returns to the strongest surviving upstream claim and directly distinguishes:

```text
fixed anatomical root is special
```

from

```text
order information follows a geometry-created amplitude-balance point.
```
