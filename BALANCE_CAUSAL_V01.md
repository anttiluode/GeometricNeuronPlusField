# Causal amplitude-balance experiments v0.1

Two interventions now test the `SomaWhyClaude` amplitude-balance explanation without changing the growth rule.

## 1. Move one source on a frozen body

Preregistered in `SOURCE_MOVE_BALANCE_PREREG_V01.md` and run on seeds 0-11 (discovery) and 12-23 (confirmation).

One source was moved inward along its existing source-to-soma path while anatomy, wave physics, pulse, lag and designated soma remained fixed.

Primary P1 asked whether the *change* in the amplitude-balance map predicts the *change* in the order-selectivity map:

```text
r_delta = corr( b_moved - b_original,
                |C|_moved - |C|_original )
```

Body-level result:

```text
                         discovery       confirmation
mean relocation r          +0.3363          +0.3484
median                     +0.3475          +0.4288
positive / negative          12 / 0           10 / 2
sign p                      0.000488         0.03857
registered P1                  PASS             PASS
```

So the spatial association found by Claude is not only static: experimentally changing where the two source amplitudes balance moves the selectivity landscape in the same direction within the same frozen anatomy.

### But the fixed soma was not overthrown

P2 asked whether, when the absolute best-balance cell moved >=3 graph edges away from soma, that cell would become more selective than the designated soma.

```text
                         discovery       confirmation
valid bodies                    12                11
mean best-balance - soma C    +.0256            +.0203
positive / negative             9 / 3             5 / 6
sign p                         .146              1.000
Wilcoxon p                     .301              .831
registered P2                  FAIL              FAIL
```

The reason is visible in the descriptive receipts: source relocation often moved the absolute balance maximum several edges away while the soma itself remained a high-percentile balance point and a high-percentile selectivity point. The intervention changed path timing as well as amplitude and did not cleanly isolate balance at one fixed location.

## 2. Change source gain only on fresh bodies

Preregistered in `SOURCE_GAIN_BALANCE_PREREG_V01.md` and run on previously unused bodies:

```text
discovery      seeds 24-35
confirmation   seeds 36-47
```

Source positions and all geometry were fixed. Five gain conditions were used:

```text
A/B = 1/1, 0.5/1, 2/1, 1/0.5, 1/2
```

### Whole-map causal result is strong and replicated

P2 again asked whether gain-induced changes in balance predict gain-induced changes in selectivity across body cells.

```text
                         discovery       confirmation
mean body delta-map r       +0.5928          +0.6066
median                      +0.5763          +0.6546
positive / negative          12 / 0           12 / 0
sign p                      0.000488         0.000488
registered P2                  PASS             PASS
```

This is the cleanest causal result so far for the balance variable because source gain changes amplitude without moving the source or changing path length.

### Same-soma dose response is not universal

P1 held location fixed and asked whether, within each soma, gain-induced changes in soma balance predict changes in soma selectivity.

```text
                         discovery       confirmation
mean soma delta r           +0.5005          +0.8333
median                      +0.9999          +1.0000
positive / negative           9 / 3           11 / 1
sign p                         .146             .00635
registered P1                  FAIL             PASS
```

The near +/-1 values arise because the symmetric half/double manipulations generate a very constrained four-point perturbation geometry. Three discovery bodies show the opposite relation strongly enough that the preregistered discovery criterion fails.

Therefore amplitude balance is **causal but not sufficient**. Across the body, changing balance robustly changes where order selectivity lives. At one fixed soma, however, the exact effect of balancing the two inputs depends on the body's source-to-soma transfer functions.

## Current interpretation

The most economical picture is now:

```text
geometry
   -> two source-specific transfer functions h_A(x,t), h_B(x,t)
   -> amplitude balance sets whether both sources can matter locally
   -> their relative temporal/complex structure sets whether swapping order matters
```

Balance is therefore an **opportunity variable** for coincidence/order computation, not the whole computation.

This also explains why the soma/root can remain privileged after the absolute balance maximum moves: it is not enough for two peak amplitudes to be equal. Their full transfer histories must interact with the task lag in an order-sensitive way.

## Next question

The next clean experiment should use only the two single-source impulse responses at each cell and decompose the pair response into:

```text
self/envelope terms:
    |h_A|^2 + |shift(h_B)|^2

cross/interference term:
    2 Re[h_A * conj(shift(h_B))]
```

Then ask whether the full A->B/B->A selectivity is already predicted by the incoherent self/envelope terms, or whether the cross-source complex term is required. That identifies the second factor without adding another biological compartment.
