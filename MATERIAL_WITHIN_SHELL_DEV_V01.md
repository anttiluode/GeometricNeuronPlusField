# Within-shell material release — development v0.1

Date: 2026-08-08

## Question

After the strongest distance-only material field is optimized, is there any useful material organization left that depends on *which cell or branch* occupies a given graph-distance shell?

This is a development result on already-seen bodies 622–627. It is not a held-out confirmation and not a novelty claim.

---

## Why this control was necessary

The first full-vs-radial comparison used an inconsistent shell-projection metric and overstated the non-radial residual. The corrected hostile radial learner in `material_radial_vs_full_v02.py` captured about 95% of the full learner's gain over uniform material:

```text
uniform R2          .491650
radial R2           .594962
full R2             .600233
full - radial      +.005271
full wins            5 / 6
```

Moreover the independent full-vs-radial difference split by frequency:

```text
omega=.03   -.001911
omega=.04   +.012453
```

So the original idea of a large broadband branch-specific residual was rejected.

The cleaner nested test starts from the corrected radial optimum and does not allow the radial profile to move at all.

---

# Nested release protocol

For each frozen arbor:

1. optimize the distance-shell material field for the joint `.03/.04` objective;
2. record the total material in every integer graph-distance shell;
3. freeze every one of those shell totals exactly;
4. release only cell-to-cell redistribution *within each shell*;
5. optimize the same objective using the exact cellwise material gradient projected to zero sum inside every shell.

Thus the branch-only learner cannot change:

```text
total material budget
material in any graph-distance shell
the radial material profile
```

It can change only which equal-distance cells receive more or less material.

Numerical shell-total conservation was audited after learning.

---

## Development result — seeds 622–627

Pure within-shell release:

```text
mean joint R2 gain over optimized radial     +.00455534
positive bodies                               6 / 6

mean circular phase-RMS improvement          +.01217083 rad
positive bodies                               6 / 6

median released/radial amplitude ratio        .996067
max shell-total numerical error               2.78e-17
```

Per frequency, even though the learner was trained jointly:

```text
omega=.03
  mean R2 gain     +.00397850
  positive          4 / 6

omega=.04
  mean R2 gain     +.00513218
  positive          6 / 6
```

So branch-only redistribution can improve the joint objective without changing the distance profile.

---

# Within-shell placement shuffle

A stronger control then preserved, separately in every graph-distance shell:

```text
the shell total
+
the exact multiset / histogram of released density values
```

but randomly permuted those values among cells in the same shell.

Twelve shuffles were averaged per body.

Observed:

```text
released - within-shell shuffled R2     +.00570322
positive bodies                           6 / 6

phase-RMS gain vs within-shell shuffle   +.01518680 rad
median amplitude ratio                    .987006
```

Per frequency:

```text
omega=.03
  released - shell-shuffle R2    +.00652201
  positive                         5 / 6

omega=.04
  released - shell-shuffle R2    +.00488443
  positive                         6 / 6
```

Therefore the development residual is not explained merely by adding heterogeneity within a distance shell. The *placement* of those values among equal-distance cells matters.

---

## Current decomposition

The material-learning effect now has a quantitative hierarchy:

```text
large component:
    distance from consequential readout
    approximately 95% of the full gain in the corrected development comparison

small component:
    within-distance cell / branch placement
    about +.0046 R2 over optimized radial in the nested release
    about +.0057 R2 over shell-histogram-matched shuffles
```

This small component is the only material-geometry residual currently worth a held-out test.

---

## Prior-art boundary

Do not describe the residual as a new biological phenomenon merely because it is branch specific.

Prior work already includes:

- branch-specific ion-channel conductance densities related to branch morphology;
- dendritic democracy / distance compensation;
- local activity-dependent intrinsic plasticity;
- adjoint gradients for distributed neuronal channel density;
- gradient-based optimization of large branchwise conductance parameter sets.

The present model-level question is narrower:

> With the best readout-distance material profile held exactly fixed, does consequence-driven exact local material credit reproducibly find additional equal-distance branch/cell placement that improves a temporal phase-coordination objective?

## Development wall sentence

> **Distance explains almost all of the learned electrical-material organization, but not all: with every distance shell's material total frozen, exact local credit can still redistribute material among equal-distance cells to gain a small amount of phase coordination, and shuffling that within-shell placement destroys the gain.**
