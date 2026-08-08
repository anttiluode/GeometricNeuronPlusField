# Single-source transfer decomposition v0.1 — preregistration

## Motivation

Two causal balance experiments now show that changing the relative A/B amplitudes changes the spatial order-selectivity map, but amplitude balance is not sufficient to explain the exact soma response in every body.

For a nearly linear complex field, the missing information is contained in the full source-specific transfer histories.

At one cell `x`, let the complex single-source responses be

```text
h_A(x,t)
h_B(x,t)
```

for the exact frozen body and exact gain condition.

At task lag `tau`, the coherent synthetic pair responses are

```text
psi_T = h_A(t)       + h_B(t-tau)
psi_D = h_B(t)       + h_A(t-tau)
```

and power expands as

```text
|psi|^2 = |h_1|^2 + |h_2|^2 + 2 Re[h_1 conj(h_2)].
```

The first two terms are an incoherent/envelope-only model. The last term is the cross-source complex interference term.

The question is:

> **Is temporal-order selectivity already explained by the two source power envelopes, or is the cross-source complex term required?**

## Frozen protocol

Use fresh v0.9 bootstrap bodies. No growth or parameter changes.

Use the same five source-gain conditions as `SOURCE_GAIN_BALANCE_PREREG_V01.md`:

```text
baseline      A=1.0  B=1.0
A_half        A=0.5  B=1.0
A_double      A=2.0  B=1.0
B_half        A=1.0  B=0.5
B_double      A=1.0  B=2.0
```

For each condition:

1. simulate A alone and record the full complex field history `h_A(x,t)`;
2. simulate B alone and record `h_B(x,t)`;
3. simulate the real paired A->B and B->A drives;
4. synthesize A->B and B->A from the two single-source histories without running the pair dynamics.

Compute three signed contrast maps using the same peak-power contrast as the historical task:

```text
C_actual     real pair simulation
C_coherent   |h_A + shifted h_B|^2
C_incoherent |h_A|^2 + |shifted h_B|^2
```

`C_incoherent` explicitly deletes the cross-source term while preserving each source's complete power-envelope history.

## P1 — decomposition validity

Across all occupied cells and all five gain conditions within each body, compare `C_coherent` with `C_actual`.

Registered prediction:

```text
mean body signed correlation > 0.95
```

This is expected if pair dynamics are sufficiently close to superposition at the tested gains. If it fails, the coherent/incoherent decomposition is not a valid microscope for this regime and P2/P3 are not interpreted mechanistically.

## P2 — is the complex cross term load-bearing across the body?

For each body compute

```text
MAE_coherent   = mean |C_actual - C_coherent|
MAE_incoherent = mean |C_actual - C_incoherent|
D_body         = MAE_incoherent - MAE_coherent
```

over all cells and five gain conditions.

Registered cross-term prediction:

```text
mean D_body > 0
and two-sided body-level sign test p < .05.
```

If P1 is valid and P2 passes, envelope timing alone is insufficient: the cross-source complex term is required to reproduce where/order selectivity appears.

## P3 — does the same conclusion hold specifically at the soma?

At the designated soma, average absolute prediction error over the five gain conditions and form

```text
D_soma = error_incoherent - error_coherent.
```

Registered prediction:

```text
mean D_soma > 0
and two-sided body-level sign test p < .05.
```

This directly tests the missing factor exposed by the non-universal same-soma gain response.

## Secondary receipts

Report:

- signed and absolute correlations for coherent and incoherent predictions;
- coherent/incoherent MAE by gain condition;
- mean `|C|` at soma and over the body for all three models;
- the fraction of actual selectivity magnitude retained by the incoherent model;
- mean magnitude of the explicit cross term at the actual target/distractor peak times.

These are descriptive.

## Fresh discovery / confirmation sets

Earlier tests used seeds 0-47. Use unseen bodies:

```text
discovery      seeds 48-59
confirmation   seeds 60-71
```

The code and criteria are frozen before either canonical run.

## Interpretation

Possible outcomes:

```text
P1 pass, P2/P3 pass:
    amplitude balance is an opportunity variable; order computation additionally
    requires the relative complex transfer histories / interference term.

P1 pass, P2/P3 fail:
    source power-envelope timing is sufficient; complex phase is not needed here.

P1 fail:
    pair nonlinearity is too strong for this decomposition and a nonlinear
    transfer analysis is required before claiming either mechanism.
```
