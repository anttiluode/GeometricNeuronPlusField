# Interference factorization discovery preregistration v0.1

## Question

`KERNEL_ALIGNMENT_V01.md` showed that the soma is strongly privileged in amplitude balance but only moderately special in normalized task-kernel orientation. The fixed descriptive product `balance * alignment` was useful but only partly confirmed.

The next step is to replace that heuristic product with the exact algebra of coherent superposition.

> **Can temporal-order selectivity be predicted by an explicit interference-visibility × lagged-complex-compatibility decomposition of the two single-source transfer histories, without fitted coefficients?**

## Bodies

Fresh, previously unused FunctionalArbor bodies:

```text
seeds 132-143
```

Same frozen v0.9 bootstrap, lag `tau=20`, and 210-step single-source traces.

## Local factorization

For each occupied cell `x`, record

```text
h_A(x,t), h_B(x,t).
```

Let

```text
E_A = sum_t |h_A|^2
E_B = sum_t |h_B|^2.
```

The amplitude/opportunity factor is the standard coherent-interference visibility

```text
V = 2 sqrt(E_A E_B) / (E_A + E_B),        0 <= V <= 1.
```

Define the normalized lagged complex overlaps

```text
rho_plus  = sum_{t=tau}^{T-1} h_A(t) conj(h_B(t-tau)) / sqrt(E_A E_B)
rho_minus = sum_{t=tau}^{T-1} h_B(t) conj(h_A(t-tau)) / sqrt(E_A E_B).
```

The directional complex-compatibility factor is

```text
Delta_rho = Re(rho_plus) - Re(rho_minus).
```

If each single-source response is embedded in a zero-padded `T+tau` window before the second response is shifted, then both orders have the same self-energy `E_A+E_B`. Their **exact integrated-energy temporal-order contrast** is therefore

```text
C_int = V * Delta_rho /
        (2 + V * (Re(rho_plus) + Re(rho_minus))).
```

This identity uses no fitted coefficients.

For comparison, historical peak-power contrast remains

```text
C_peak = (max |h_A + shift(h_B)|^2 - max |h_B + shift(h_A)|^2) /
         (max |h_A + shift(h_B)|^2 + max |h_B + shift(h_A)|^2).
```

## Registered discovery predictions

### D1 — integrated interference sign tracks peak computation

Within each body, across occupied cells:

```text
r_signed = corr(C_int, C_peak).
```

Prediction: mean `r_signed > 0.65`, with at least `10/12` positive bodies.

### D2 — integrated interference magnitude tracks peak selectivity

```text
r_abs = corr(|C_int|, |C_peak|).
```

Prediction: mean `r_abs > 0.60`.

### D3 — visibility × directional compatibility improves over visibility alone

Define the numerator score

```text
N = V * |Delta_rho|.
```

Compare across-cell correlations with historical absolute selectivity:

```text
r_V = corr(V, |C_peak|)
r_N = corr(N, |C_peak|).
```

Prediction: mean `(r_N-r_V) > 0.05` and at least `9/12` bodies improve.

### D4 — soma location privilege should lie mainly in visibility

Prediction:

```text
median soma V percentile > 0.85.
```

No pass/fail claim is preregistered for soma `|Delta_rho|` percentile; it is diagnostic.

## Controls

- Directly construct the zero-padded target/distractor traces and verify that their integrated-energy contrast equals the closed-form `C_int` to numerical precision.
- Report peak-balance `B=min(PA,PB)/max(PA,PB)` and its correlation with visibility `V`.
- Report soma percentiles for `V`, `|Delta_rho|`, `V|Delta_rho|`, `|C_int|`, and `|C_peak|`.
- No parameter fitting, per-body weights, or lag tuning.

## Interpretation

If this succeeds, the previous amplitude-balance and coherent-cross-term results reduce to a familiar local statement:

```text
order computation = amplitude visibility × directional complex compatibility
                     (plus the exact normalization in C_int).
```

Geometry remains essential because it creates `h_A` and `h_B`; the local readout then exposes their relation through coherent square-law superposition.
