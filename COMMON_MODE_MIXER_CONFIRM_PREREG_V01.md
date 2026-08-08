# Common-mode mixer confirmation v0.1 — preregistration

This confirmation is frozen after the exploratory `common_mode_mixer_probe.py` run on reused seeds 72-95 and after the separately held-out mode-pair structure test on seeds 84-95.

## Hypothesis

The graph's unique constant mode (mode 0) is almost blind to temporal order when measured by its own power, but quadratic power readout can use it as a reference component that mixes with the nonconstant residual field:

```text
psi_soma = c + r
|psi_soma|^2 = |c|^2 + 2 Re[c conj(r)] + |r|^2
```

where `c` is the mode-0 soma contribution and `r` is the sum of all nonconstant modal contributions.

The exploratory 24-body result suggested that the first-order approximation

```text
P_ref = |c|^2 + 2 Re[c conj(r)]
```

tracks the signed full temporal-order contrast much better than the residual-only power `|r|^2`, despite mode 0 alone being essentially order blind.

This is mathematically analogous to reference-assisted square-law detection. The term "mixer" is used in that signal-processing sense only.

## Held-out bodies

Use previously unseen FunctionalArbor bodies:

```text
seeds       96-107
lag         20
steps       210
source gain A/B = 1/1
```

No growth, wave, source, soma, lag, gain, or readout parameter is changed.

## Registered predictions

### P1 — mode 0 alone is order blind

Mean absolute temporal-order contrast of `|c|^2` must be

```text
< 0.001.
```

### P2 — first-order common-mode mixer tracks the full soma computation

Across the 12 bodies all three must hold:

```text
corr(C_ref, C_full)          > 0.95
mean |C_ref - C_full|        < 0.08
same signed preference       >= 10 / 12 bodies
```

### P3 — the reference term is more informative about the full result than residual power alone

Define per-body absolute errors

```text
e_ref = |C_ref - C_full|
e_res = |C_residual - C_full|.
```

Require:

```text
mean e_ref < mean e_res
and e_ref < e_res in at least 10 / 12 bodies.
```

Also report `corr(C_residual, C_full)` descriptively.

## Interpretation

If P1-P3 pass together, the earned statement is:

> **The constant mode is not itself the temporal-order code. Instead, under the soma's quadratic readout it acts as a reference component whose cross term with nonconstant geometry-shaped modes recovers much of the signed order computation.**

That would explain how the globally energy-dominant mode can be nearly blind when read alone yet still be load-bearing in the local modal interaction structure.

If P2 or P3 fails, the stronger reference/mixer interpretation is rejected even though the earlier absolute pair-mass hub result remains valid.
