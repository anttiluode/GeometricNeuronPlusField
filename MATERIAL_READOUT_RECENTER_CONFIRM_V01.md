# Material readout re-centering confirmation v0.1

Date: 2026-08-08

## Result in one sentence

> **Moving the consequential readout across the same frozen arbor re-centers the learned quasi-active material organization around the new readout rather than the anatomical soma; all 6/6 preregistered held-out criteria passed on eight fresh bodies.**

This is a model-level mechanism result. It does not claim that biological neurons relocate a soma or that HCN trafficking implements this exact optimizer.

---

## Why this control matters

The material-learning development experiment starts from uniform electrical material and optimizes only somatic phase coherence. After learning, density strongly increases with graph distance from the soma even though graph distance never enters the objective.

That left a serious ambiguity:

```text
is distance-from-soma genuinely the useful coordinate?

or

is the soma simply hard-wired into the model in some hidden way?
```

The clean hostile intervention is to move the place where consequence is read out while leaving the same frozen morphology and learning rule intact.

If the material map remained soma-centered, the geometric interpretation would weaken sharply.

If it re-centered on the new readout, the useful coordinate would be **consequence-relative**.

The candidate, procedure and thresholds were frozen in `MATERIAL_READOUT_RECENTER_CONFIRM_PREREG_V01.md` before seeds 596–603 were examined.

---

## Frozen intervention

For each fresh frozen FunctionalArbor body:

1. run the original material learner with the anatomical soma as readout;
2. choose the occupied cell farthest in graph distance from the soma;
3. restart from the same uniform material state;
4. optimize the identical two-frequency phase-coherence objective around that moved readout;
5. only after learning, measure density correlation with distance from:
   - the moved readout;
   - the anatomical soma.

Frozen learning parameters:

```text
omega          .03, .04
tau_h           2
mu              .5
budget ref g0   .005
budget ratio    10
local cap       .05
steps           50
step_fraction   .10
```

No graph-distance coordinate enters training.

---

# Held-out confirmation — seeds 596–603

All eight requested bodies were available.

## R0 — moved-readout learning remains useful

Registered:

```text
mean moved-readout coherence gain > .05
positive gain on >= 75% of bodies
```

Observed:

```text
mean gain      +.1179686
positive        8 / 8
```

**R0 PASS.**

---

## R1 — density organizes around the moved readout

Registered:

```text
mean Spearman rho(density, distance-to-new-readout) > .40
positive rho on >= 75% of bodies
```

Observed:

```text
mean rho       +.6684531
positive        8 / 8
```

Per body:

```text
596   +.691
597   +.507
598   +.853
599   +.682
600   +.643
601   +.650
602   +.613
603   +.708
```

**R1 PASS.**

---

## R2 — the new readout explains density better than the old soma

Registered:

```text
delta_rho = rho_to_new_readout - rho_to_old_soma

mean delta_rho > .40
positive delta on >= 75% of bodies
```

Observed:

```text
mean delta_rho   +1.316491
positive delta    8 / 8
```

**R2 PASS.**

---

## R3 — the learned map is no longer soma-centered

Registered:

```text
mean rho(density, distance-to-anatomical-soma) < .15
```

Observed:

```text
mean rho to old soma   -.6480376
```

Every fresh body was negatively correlated with old-soma distance after moving the readout.

**R3 PASS.**

---

## R4 — the original soma-centered result reproduces on the same bodies

Registered:

```text
mean soma-centered coherence gain > .05
mean soma-centered rho(distance-to-soma) > .40
```

Observed:

```text
mean soma gain    +.1389899
mean soma rho     +.7761055
```

**R4 PASS.**

So the re-centering result is not obtained by choosing a peculiar held-out set where the original material-learning effect disappeared.

---

## R5 — moving the readout does not trivialize learning

Registered:

```text
mean moved gain / mean soma gain > .45
```

Observed:

```text
ratio   .8487563
```

**R5 PASS.**

---

# Confirmation verdict

```text
R0 PASS   moved-readout learning remains useful
R1 PASS   density organizes by distance from new readout
R2 PASS   new-readout coordinate beats old-soma coordinate
R3 PASS   population is no longer soma-centered
R4 PASS   original soma result reproduces on same bodies
R5 PASS   moved learning gain stays same order as soma learning

6 / 6 PASS
```

---

## What this changes conceptually

The earlier shorthand

```text
distal material enrichment from the soma
```

is now too specific.

The stronger model-level abstraction is

```text
material organization by transfer distance from consequence.
```

The soma happens to be one consequential readout in the original model. When consequence is moved elsewhere, the transpose field changes automatically and the material distribution follows it.

This turns a superficially anatomical result into a more general geometric one.

---

## Local view of the mechanism

For local material density `d_i`, the harmonic sensitivity is

```text
dH/dd_i = -c(omega) y_i x_i
```

where

```text
x_i   local forward field
y_i   local transpose field launched from the readout.
```

Moving the readout changes `y` without giving each local tuner a new global coordinate system.

A useful interpretation is:

```text
forward field
    says what is happening locally

transpose / return field
    says how that location matters to consequence elsewhere

local overlap
    changes electrical material
```

The emergent distance coordinate is therefore not represented explicitly at the local tuner. It is implicit in the geometry of the two fields.

---

## Current wall sentence

> **The learned electrical-material coordinate follows consequence: moving the readout across a fixed arbor causes the same local material learner to re-center its spatial organization around the new readout and anti-correlate with the old soma coordinate.**

This result now has both a development observation and an independent held-out confirmation.
