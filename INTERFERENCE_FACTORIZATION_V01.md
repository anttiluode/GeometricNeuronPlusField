# Interference factorization v0.1 — amplitude visibility × directional complex compatibility

This experiment replaces the heuristic `balance * kernel-alignment` product with the exact algebra of coherent superposition.

The question was:

> **Can temporal-order selectivity be predicted by an explicit interference-visibility × lagged-complex-compatibility decomposition of the two single-source transfer histories, without fitted coefficients?**

The answer is **yes, strongly and on held-out bodies**.

## Local algebra

At each occupied cell `x`, the frozen geometry produces two single-source complex transfer histories

```text
h_A(x,t), h_B(x,t).
```

Define total source energies

```text
E_A = sum |h_A|^2
E_B = sum |h_B|^2.
```

The coherent interference-visibility factor is

```text
V = 2 sqrt(E_A E_B) / (E_A + E_B),   0 <= V <= 1.
```

Define normalized lagged overlaps at the registered lag `tau=20`:

```text
rho_plus  = sum h_A(t) conj(h_B(t-tau)) / sqrt(E_A E_B)
rho_minus = sum h_B(t) conj(h_A(t-tau)) / sqrt(E_A E_B).
```

The directional compatibility term is

```text
Delta_rho = Re(rho_plus) - Re(rho_minus).
```

With zero-padding so delaying either source preserves its total self-energy, the integrated temporal-order contrast is exactly

```text
C_int = V * Delta_rho /
        [2 + V * (Re(rho_plus) + Re(rho_minus))].
```

No fitted coefficient appears anywhere in this expression.

The direct zero-padded energy calculation and the closed-form expression agreed to relative error about `5e-15` in both discovery and confirmation.

## Discovery — fresh seeds 132-143

All four registered discovery predictions passed.

```text
mean corr(C_int, C_peak)                 0.78259
positive signed-correlation bodies       12 / 12

mean corr(|C_int|, |C_peak|)             0.78415

mean corr(V, |C_peak|)                   0.62972
mean corr(V*|Delta_rho|, |C_peak|)       0.83883
mean improvement                         0.20911
improved / worse                         11 / 1
sign p                                   0.00635

median soma V percentile                 0.92857
mean soma |Delta_rho| percentile         0.61310
```

So visibility alone is useful, but adding the directional complex relationship between the two transfer histories improves the body-wide prediction markedly.

## Held-out confirmation — fresh seeds 144-155

The confirmation criteria were frozen before these bodies were run.

### C1 — signed integrated interference tracks signed peak computation

```text
mean corr(C_int, C_peak)                 0.75178
positive bodies                          12 / 12
threshold                                >0.72 and >=11 positive
```

**C1 PASS.**

### C2 — magnitude relation replicates

```text
mean corr(|C_int|, |C_peak|)             0.73040
threshold                                >0.72
```

**C2 PASS.**

### C3 — compatibility adds information beyond visibility

```text
mean corr(V, |C_peak|)                   0.62108
mean corr(V*|Delta_rho|, |C_peak|)       0.80693
mean improvement                         +0.18584
positive / negative bodies               12 / 0
sign p                                   0.000488
```

**C3 PASS.**

This is much cleaner than the earlier `balance * normalized-kernel-alignment` result. The standard interference factorization improves over visibility in every held-out body.

### C4 — soma privilege is mainly amplitude opportunity

```text
median soma V percentile                 0.90714
median soma |Delta_rho| percentile       0.77857
registered requirements                  >0.88 and <0.80
```

**C4 PASS.**

The soma/root is therefore unusually favorable in amplitude visibility, while directional complex compatibility is less consistently extreme and remains body-specific.

## Relation to SomaWhyClaude

The old amplitude-balance variable was

```text
B = min(P_A,P_B) / max(P_A,P_B).
```

Energy visibility `V` is a smoother and physically standard version of the same opportunity concept. Across held-out bodies, `corr(V,B)` averaged about `0.83`.

The important improvement is that the second factor is no longer an abstract normalized kernel score. It is the actual lag-directional complex overlap of the two source transfer histories.

## Relation to the mode-pair result

`MODE_PAIR_V01.md` showed that soma point power is a sparse, mostly off-diagonal modal mixer. The present result says what that mixing means in local transfer-function language.

The two descriptions are equivalent views of the same coherent interaction:

```text
modal view
---------
q_A,n, q_B,m
  -> sparse pair interactions through K_s

local transfer view
-------------------
h_A, h_B
  -> visibility V
  -> lagged complex compatibility Delta_rho
  -> coherent order contrast
```

The graph basis explains how geometry creates the histories. The interference factorization explains what relation between those histories matters at a local readout.

## Current mechanism

The strongest current reduction is now:

```text
frozen anatomy G
   -> graph-defined resonator / transfer structure
   -> source histories h_A(x,t), h_B(x,t)
   -> convergence/root gives high interference visibility V
   -> geometry-dependent lagged complex relation Delta_rho
   -> coherent source-source cross term
   -> local square-law readout
   -> temporal-order selectivity
```

This is more precise than saying "amplitude balance plus phase."

`V` answers:

> are both histories locally present with enough comparable strength to interfere?

`Delta_rho` answers:

> does their geometry-shaped temporal/complex relationship distinguish +tau from -tau?

Both are needed.

## Important scope

`C_int` is an **integrated-energy** statistic, while the historical FunctionalArbor objective uses **peak power**. They are not numerically identical. What replicated is the spatial and signed structure: a coefficient-free integrated interference statistic predicts the peak-order map with correlations around `0.73-0.75` on held-out bodies.

That is strong evidence that the peak task is not an unrelated artifact of one sample time. It is substantially organized by the same directional coherent interaction present in the complete source histories.

## Wall sentence

> **The root is privileged because convergence makes the two geometry-shaped transfer histories highly visible to one another; temporal order is then supplied by the asymmetry of their lagged complex overlap. Geometry creates both the amplitude opportunity and the directional compatibility, and the local square-law readout exposes their coherent product.**

## Next clean question

The unresolved learning problem can now be stated in these coordinates:

> **When one structural event changes the arbor, is its change in task value predictable from how it changes interference visibility `V` and directional compatibility `Delta_rho` at the soma?**

If so, structural eligibility can be defined against a physically meaningful intermediate quantity rather than broad activity tags or cell identity.
