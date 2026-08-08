# Material-adjoint learning discovery v0.1

Date: 2026-08-08

## Question

After confirming that a hand-imposed morphology-indexed quasi-active material gradient can compensate location-dependent transfer phase at the soma, remove graph distance from the objective entirely.

Start every occupied arbor cell with the **same material density**, keep the total amount of material fixed, and ask whether the exact local material-adjoint sensitivity can discover a useful spatial distribution on its own.

Only after learning inspect whether the learned density has any relation to graph distance or branch geometry.

This is a development/discovery result. It has not yet been held out.

---

## Frozen material physics reused from the confirmed HCN-like bridge

```text
v' = K L psi - damping*v - restoring*psi
     - d(x)*psi
     - mu*d(x)*z
     + source

z' = (psi-z)/tau_h
```

with

```text
tau_h = 2
mu    = .5
omega = .03 and .04
```

The total density budget equals the amount contained in the previously confirmed hand profile with

```text
g0                  .005
soma->distal ratio  10x
```

but the learner is never shown that profile.

Constraints:

```text
d_i >= 0
sum_i d_i = fixed total material budget
d_i <= .05
```

---

## Objective contains no morphology coordinate

For each frequency, independently inject the same unit harmonic source at many occupied locations and collect the soma transfers

```text
H_j.
```

Normalize every transfer to unit phase phasor

```text
u_j = H_j / |H_j|
```

and maximize squared circular coherence

```text
R^2 = |mean_j u_j|^2.
```

The training objective is the mean `R^2` across `omega=.03,.04`.

Important consequences:

- graph distance is absent;
- distal enrichment is absent;
- branch identity is absent;
- absolute transfer amplitude does not weight one source location more than another after phase normalization.

The optimizer only sees whether the phases observed at the soma line up.

---

## Exact local material gradient

For harmonic operator

```text
A(d,omega) x_j = b_j
```

and soma transpose field

```text
A^T y = e_s,
```

local density sensitivity factorizes as

```text
dH_j/dd_i
  = -c(omega) y_i x_{j,i}.
```

`hcn_material_learning.py` differentiates the circular-coherence objective through this local factorization.

Finite-difference audit over the four development bodies:

```text
maximum relative derivative error   4.22e-8
```

so the learned redistribution is not being driven by a numerically guessed material gradient.

See `MATERIAL_ADJOINT_DERIVATION_V01.md` for the algebra.

---

# Development learning — seeds 580–583

All four bodies were trained for 50 projected-gradient steps from their uniform-density state.

## Objective result

```text
mean phase coherence R^2

uniform start          .48438
learned                .61750
hand linear gradient   .55241
shuffled learned map   .45339
```

Body by body:

```text
seed 580   uniform .5019 -> learned .6248   hand .5471   shuffled .4558
seed 581   uniform .5066 -> learned .6187   hand .5618   shuffled .4943
seed 582   uniform .5306 -> learned .6411   hand .5830   shuffled .4995
seed 583   uniform .3985 -> learned .5855   hand .5177   shuffled .3640
```

Therefore, on the development set:

```text
learned beats uniform       4 / 4
learned beats shuffled      4 / 4
learned beats hand gradient 4 / 4
```

The corresponding mean circular soma phase RMS was

```text
uniform             .87720
learned             .68199
hand gradient       .75907
shuffled learned    .91798
```

---

## A distance coordinate emerged without being supplied

Only after learning, compare density with graph distance from the soma.

```text
Spearman(density, graph distance)

seed 580    +.907
seed 581    +.676
seed 582    +.901
seed 583    +.609

mean        +.773
```

Pearson correlation averaged `+.786`.

All four learned maps therefore independently put more material farther from the readout despite having no distance term in the objective or update.

This is the strongest discovery in this branch so far.

---

## But the learned object is not simply the hand-drawn linear gradient

The optimization is constrained by a fixed material budget and a hard `.05` per-cell cap.

The learned solutions are often strongly polarized:

- a substantial proximal region is driven to approximately zero density;
- many mid/distal cells reach the cap;
- in some arbors the farthest tips then fall below the cap again.

Examples from the development maps:

- seeds 580 and 582 are close to a proximal-zero / distal-saturated allocation;
- seed 581 peaks through a broad mid-distal region and declines at the far end;
- seed 583 peaks around the mid/distal arbor and declines substantially toward its most extreme tips.

A linear graph-distance coordinate explains roughly 35%–82% of density variance depending on body.

So the result should **not** be summarized as

```text
learning rediscovered the exact biological HCN gradient.
```

A better description is

```text
learning discovered a strong distance-from-readout backbone
plus morphology-specific local corrections.
```

---

## Amplitude caveat

The objective phase-normalizes every injection location, so weak sites do not disappear from the coherence score.

Nevertheless the learned material does alter transfer amplitude.

Median learned/uniform soma-amplitude ratios on the four bodies were approximately

```text
.780
.595
.842
.602
```

This is not catastrophic suppression, but it is stronger attenuation than the confirmed hand-gradient result, whose smooth/uniform and smooth/shuffled amplitude ratios remained near one.

A later confirmation should retain an explicit amplitude sanity bound rather than silently optimizing phase at any electrical cost.

---

# Is the learned coordinate only graph distance?

`hcn_material_geometry_controls.py` reran the same learning and then attacked the learned map with controls that preserve distance organization.

### Radialized control

For every graph-distance shell, replace each local density by the shell mean.

This keeps the learned density-versus-distance profile but removes within-shell branch differences.

### Within-distance shuffle

Within every graph-distance shell, preserve the exact learned density multiset but randomly permute which branch/cell receives each value.

This keeps **both** the full distance profile and every shell's density histogram while destroying branch-specific assignment.

### Global shuffle

Shuffle the learned values over the whole arbor, destroying both distance and branch relation.

Observed mean coherence:

```text
full learned             .61750
radialized               .60584
within-distance shuffle  .60615
global shuffle           .45291
```

Body-by-body full learned versus radialized:

```text
580   .6248 vs .6211
581   .6187 vs .5945
582   .6411 vs .6390
583   .5855 vs .5688
```

Full learned beat radialized and within-distance shuffled on **4/4** bodies.

Mean advantage:

```text
full - radialized               +.01166 R^2
full - within-distance shuffle  +.01135 R^2
```

The RMS magnitude of branch-specific residual density after subtracting the radial profile was about 23% of total learned-density standard deviation.

---

## Geometry hierarchy implied by the controls

The hostile control gives a useful decomposition:

```text
readout-relative graph distance
        dominant learned coordinate

branch/topological placement at fixed distance
        smaller but consistent correction
```

Destroying the distance relation globally costs about `.165 R^2` relative to the learned map.

Removing only branch-specific detail while preserving distance costs about `.012 R^2`.

So this discovery does **not** support the dramatic claim that the learner builds an arbitrary branch code. Most of the useful material organization is radial in graph distance from the consequential readout.

But path distance is not the whole answer either: local branch placement adds a reproducible smaller gain on every development body tested.

---

## Current interpretation

This connects several formerly separate results:

```text
frozen morphology
    defines the wave-transfer problem

local quasi-active material
    changes the frequency-dependent transfer operator

somatic objective
    defines which transfer relation is consequential

transpose field
    returns sensitivity through the same reciprocal operator

local forward x transpose product
    redistributes electrical material

emergent result
    strong material organization by distance from the readout
    + smaller branch-specific correction
```

The useful coordinate was not explicitly coded into the learning objective. It emerged because cells at different locations have different leverage over phase arrival at the readout.

---

## Biological boundary

CA1 HCN1 distribution is known to be spatially graded and activity dependent, but this experiment does not establish that biological HCN trafficking follows this material-adjoint rule.

The present learned profile is also much more strongly polarized than a realistic continuously regulated channel-density field because the optimizer has no trafficking, synthesis, smoothness, homeostatic or energetic penalty beyond a fixed total budget and a hard local cap.

The biological value of the comparison is therefore structural, not literal:

> **a real neuron has an activity-regulated spatial electrical material; in the model, a local sensitivity rule can make such a material self-organize around the timing requirements imposed by the morphology and its readout.**

---

## Next hostile tests

Before confirmation, attack the emergent distance result rather than tuning it.

Two particularly strong controls are now motivated:

1. **constraint robustness:** vary the local density cap / redistribution regularity while keeping total material budget and phase objective fixed; the distance tendency should not depend on one bang-bang boundary condition;
2. **move the consequential readout:** optimize the identical objective around a different occupied readout location. If the learned material coordinate re-centers on the new readout rather than remaining soma-centered, then the result is genuinely readout-relative geometry rather than an accidental anatomical prior.

After those controls, freeze the learning algorithm and take it to untouched bodies.

## Wall sentence

> **When local quasi-active material is allowed to adapt under a fixed budget, the reciprocal material adjoint discovers a strong coordinate by distance from the consequential readout without being told that coordinate; most of the gain is radial, while branch geometry contributes a smaller second-order correction.**
