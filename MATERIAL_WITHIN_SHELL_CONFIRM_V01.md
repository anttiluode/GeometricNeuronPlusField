# Within-shell material placement confirmation v0.1 — held-out 8/8 pass

Date: 2026-08-08

## Result in one sentence

> **After the strongest graph-distance-only material profile is optimized and every distance shell's total material is frozen exactly, local material credit still finds a small but reproducible placement among equal-distance cells that improves somatic phase coordination and beats within-shell histogram-matched shuffles; all 8/8 preregistered criteria passed on 12 fresh arbors.**

This is a model-level residual, not a novelty claim. The dominant effect remains readout-distance compensation.

---

## Why this test was narrow

Prior-art review and corrected hostile controls eliminated the larger interpretations.

The repository no longer treats any of the following as candidate novelties:

```text
self-organized intrinsic conductance gradients
adjoint gradients for distributed channel density
dendritic distance compensation / democracy
return signals that encode soma distance
optimized branchwise channel maps
high-dimensional differentiable conductance training
local forward activity × returned error
physical adjoint training in reciprocal media
```

The corrected radial learner captured roughly 95% of the full learner's development gain over uniform material.

The remaining question was therefore deliberately small:

```text
with the best distance profile frozen exactly,
does exact local credit care which equal-distance cell gets which material?
```

The protocol and thresholds were frozen in `MATERIAL_WITHIN_SHELL_CONFIRM_PREREG_V01.md` before seeds 628–639 were examined.

---

# Frozen protocol

## Stage 1 — optimized radial material

One shared material density per integer graph-distance shell from the soma/readout.

```text
omega               .03 and .04
tau_h                2
mu                   .5
budget reference g0  .005
budget ratio          10
local density cap     .05
uniform start
fixed total material budget
160 radial projected-gradient steps
```

The radial learner uses the corrected ordinary Euclidean shell-parameter projection.

## Stage 2 — release only within shells

Starting from the radial optimum:

```text
freeze total material in every distance shell
release independent cell densities within each shell
80 projected-gradient steps
project every update separately inside each shell
```

No material can move between graph-distance shells.

## Stage 3 — same-shell shuffle

For each learned released map, 12 controls per body preserve, independently in every distance shell:

```text
shell total
exact multiset / histogram of learned density values
```

but randomly permute those values among equal-distance cells.

This destroys only the exact branch/cell placement.

---

# Held-out confirmation — seeds 628–639

All 12 requested bodies bootstrapped.

Per-body joint gains:

```text
seed   radial-uniform   release-radial   release-shellshuffle
628      +.109721          +.000587           +.000649
629      +.183742          +.000541           +.001086
630      +.156922          +.009351           +.010787
631      +.132056          +.003539           +.006590
632      +.107595          +.011764           +.016773
633      +.104097          +.009771           +.007923
634      +.044363          +.001588           +.003812
635      +.117274          +.005268           +.005757
636      +.177007          +.002741           +.006041
637      +.118267          +.001390           +.001047
638      +.073922          +.003100           +.006573
639      +.101164          +.000483           +.000929
```

Every fresh body improved after within-shell release, and every fresh body beat its within-shell placement shuffles.

---

## W0 — usable population

Registered at least 10 usable bodies.

Observed:

```text
12 / 12 usable
```

**W0 PASS.**

---

## W1 — distance-only material remains the dominant useful component

Registered:

```text
mean(radial - uniform) > .06 R2
positive on >=75% bodies
```

Observed:

```text
mean gain     +.1188442
positive       12 / 12
```

**W1 PASS.**

The large material-learning effect remains primarily a distance-from-readout compensation effect.

---

## W2 — branch-only release improves optimized radial material

Registered:

```text
mean(released - radial) > .002 R2
positive on >=75% bodies
```

Observed:

```text
mean gain     +.0041769
positive       12 / 12
```

**W2 PASS.**

This gain is obtained without moving any material between distance shells.

---

## W3 — exact within-shell placement matters

Registered:

```text
mean(released - within-shell shuffle) > .0025 R2
positive on >=75% bodies
```

Observed:

```text
mean gain     +.0056640
positive       12 / 12
```

**W3 PASS.**

Thus the effect is not merely that heterogeneity inside a distance shell is useful. The mapping of those same values onto equal-distance cells matters.

---

## W4 — direct circular phase spread agrees

Registered:

```text
mean(radial phase RMS - released phase RMS) > .005 rad
positive on >=2/3 bodies
```

Observed:

```text
mean gain     +.0076661 rad
positive       12 / 12
```

**W4 PASS.**

---

## W5 — amplitude remains ordinary

Registered pooled median released/radial soma-amplitude ratio between `.75` and `1.25`, with at least 75% of bodies between `.60` and `1.40`.

Observed:

```text
median ratio      1.000959
in range           12 / 12
```

**W5 PASS.**

The small coherence residual is not purchased by suppressing or amplifying the transfer channel.

---

## W6 — both frozen frequencies retain a non-radial placement signal

Registered positive mean released-minus-radial and released-minus-shuffle effects at both frequencies.

Observed:

```text
omega=.03
  radial - uniform              +.0903751
  released - radial             +.0023968
  positive vs radial             7 / 12
  released - shell shuffle      +.0052373
  positive vs shell shuffle     10 / 12

omega=.04
  radial - uniform              +.1473134
  released - radial             +.0059570
  positive vs radial            10 / 12
  released - shell shuffle      +.0060906
  positive vs shell shuffle     10 / 12
```

Both frozen frequencies have positive pooled non-radial effects.

**W6 PASS.**

The residual is stronger and more body-consistent at `.04`, but it is not exclusive to `.04`.

---

## W7 — distance remains quantitatively dominant

Registered:

```text
G_radial  > 0
G_release / G_radial < .15
```

Observed:

```text
G_radial                  .1188442
G_release                 .0041769
G_release / G_radial      .0351461
```

So the within-shell residual is only about **3.5% as large as the radial gain** in this fresh population.

**W7 PASS.**

This criterion is important: a confirmed branch residual does not overturn the hierarchy.

---

# Confirmation verdict

```text
W0 PASS   usable population
W1 PASS   optimized distance profile is strongly useful
W2 PASS   within-shell release improves optimized radial material
W3 PASS   exact equal-distance placement beats shell-matched shuffles
W4 PASS   direct circular phase spread agrees
W5 PASS   transfer amplitude remains ordinary
W6 PASS   both .03 and .04 retain positive pooled residuals
W7 PASS   distance remains overwhelmingly dominant

8 / 8 PASS
```

---

## What has actually been earned

The mature material result is now hierarchical rather than monolithic:

```text
FIRST ORDER
readout-distance material organization
    large effect
    robust
    close in spirit to established dendritic-democracy / HCN compensation ideas

SECOND ORDER
within-distance branch/cell placement
    small effect
    +.00418 R2 over optimized radial
    +.00566 R2 over same-shell histogram-matched shuffles
    12/12 positive on both pooled comparisons
    amplitude neutral
```

The second-order result says that graph distance is not a sufficient statistic for the present irregular wave arbor, even though it explains almost all of the useful material organization.

---

## Local interpretation

For local material density `d_i`, exact sensitivity is

```text
dH/dd_i = -c(omega) y_i x_i
```

where `x_i` is the local forward field and `y_i` is the readout-launched transpose field.

Cells in the same graph-distance shell can therefore receive different credit because their forward and returned complex fields differ through branch connectivity, interference and frequency-dependent loading.

The held-out result does not prove that biological intrinsic plasticity computes this quantity. It shows only that the quantity contains a reproducible second-order coordinate beyond graph distance in this model.

---

## Prior-art boundary

This result should not be advertised as the discovery of branch-specific channel physiology.

Prior work already includes:

- morphology-dependent branch conductance heterogeneity;
- dendritic democracy and distance compensation;
- activity-dependent intrinsic plasticity;
- adjoint channel-density gradients;
- inverse-designed and differentiably trained channel distributions.

The specific model-level contribution is the controlled decomposition:

```text
optimize the radial coordinate
freeze it exactly
release only equal-distance placement
show a small gain
shuffle only within those same shells
show the exact placement matters
confirm on fresh bodies
```

That decomposition is useful whether or not it ultimately proves novel.

---

## Relation to the broader Geometric Neuron idea

The result gives a precise meaning to the statement that geometry is more than path length.

For this task:

```text
path length / graph distance
    is the dominant coarse coordinate

branch topology + wave interference
    provide a smaller correction among locations with the same coarse coordinate
```

So a useful hierarchy is:

```text
consequence / readout
      ↓
coarse transfer geometry: distance
      ↓
fine transfer geometry: branch-specific field relation
      ↓
local electrical material
```

## Wall sentence

> **Graph distance explains almost all of the self-organized electrical-material solution, but it is not sufficient: after every distance shell's material budget is frozen, exact local credit reproducibly finds a small equal-distance placement correction, and shuffling only that within-shell placement removes the gain on held-out arbors.**
