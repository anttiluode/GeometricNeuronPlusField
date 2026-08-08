# Material branch-residual preregistration v0.1 — superseded unexecuted

Date: 2026-08-08

`MATERIAL_BRANCH_RESIDUAL_CONFIRM_PREREG_V01.md` was committed before seeds 628–639 were examined.

It is now **superseded without execution**. Seeds 628–639 have not been run for that confirmation and remain available as untouched bodies.

## Why it was superseded

The preregistration was based on `material_radial_vs_full.py`, whose distance-shell learner used an inconsistent combination:

```text
exact shell derivative = sum of cell derivatives in the shell
+
projection metric corresponding to expanded per-cell coordinates
```

For a shell parameter `x_s` with shell population `n_s` and material-budget constraint

```text
sum_s n_s x_s = B,
```

the ordinary Euclidean shell-parameter projection instead has KKT form

```text
x_s = clip(v_s - lambda n_s, 0, cap).
```

`material_radial_vs_full_v02.py` corrected this before any held-out branch-residual body was examined.

## Effect of the correction on already-seen development bodies 622–627

Original weaker radial control:

```text
uniform R2          .491650
radial R2           .587817
full R2             .600233
full - radial      +.012416
full wins            6 / 6
```

Corrected hostile radial control:

```text
uniform R2          .491650
radial R2           .594962
full R2             .600233
full - radial      +.005271
full wins            5 / 6
```

The corrected radial learner captures about 95% of the full learner's gain over uniform material.

The per-frequency comparison also invalidated the old preregistered two-frequency story:

```text
omega=.03   mean(full-radial)  -.001911   full wins 2/6
omega=.04   mean(full-radial)  +.012453   full wins 4/6
```

Therefore it would be inappropriate to execute the old thresholds on fresh bodies.

## Current replacement mechanism test

Rather than comparing two separately optimized solutions, the new development test is nested:

1. optimize the corrected radial material field;
2. freeze the total material in every graph-distance shell exactly;
3. release only cell-to-cell redistribution within each shell;
4. optimize the same `.03/.04` objective with the exact material gradient;
5. compare the released placement with within-shell shuffles preserving every shell's total and value histogram.

Only after that development mechanism is clean should any new held-out branch-specific claim be preregistered.

## Historical rule

The original preregistration is retained unchanged. This note does not rewrite it; it records why it was never executed.
