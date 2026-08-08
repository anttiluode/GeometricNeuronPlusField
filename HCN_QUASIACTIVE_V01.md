# Quasi-active morphology-indexed material v0.1

## Result in one sentence

> **On frozen FunctionalArbor geometries, a soma-to-distal distribution of a complete minimal quasi-active restorative material synchronizes the phase of spatially distributed harmonic inputs at the soma better than uniform, shuffled, or reversed placement, while preserving substantial local dendritic phase structure; the frozen candidate passed all 7/7 held-out criteria on 12 fresh bodies.**

This is a computational result in the present wave-arbor model. It is not a claim that the numerical material parameters are biological HCN conductances or that CA1 neurons implement this exact equation.

---

## Why this branch was started

Vaidya & Johnston (2013) showed experimentally and with biophysical modeling that CA1 pyramidal neurons use a spatial gradient of HCN-channel-mediated inductive membrane properties to counteract location-dependent temporal differences of dendritic inputs at the soma.

That result is a close biological neighbor of the mature Geometric Neuron question:

```text
morphology creates location-dependent transfer relations
+
spatial electrical material reshapes those relations
->
a coordinated consequential readout can emerge at the soma
```

without requiring every local dendritic voltage to become globally phase locked.

---

## The first HCN-like proxy failed its mechanism test

The first development model in `hcn_impedance_probe.py` added only a delayed restorative state:

```text
v'   = ... - g_h(x) z
z'   = (psi-z)/tau_h
```

A smooth soma-to-distal gradient strongly reduced somatic phase spread in a low-frequency band. At one attractive development point the smooth profile gave approximately

```text
soma phase RMS
smooth      .815
uniform    1.107
shuffle    1.131
zero       1.650
```

But `hcn_inductive_diagnostic.py` showed the wrong biological mechanism sign. Relative to uniform material, the smooth-gradient phase shift became *less* advancing with graph distance (`r ~= -.81`), and local phase spread collapsed to about 55% of passive at that operating point.

Therefore the first result was explicitly rejected as an HCN mechanism despite its attractive synchronization metric.

This failure mattered because it exposed a missing term rather than merely a bad parameter choice.

---

## The missing quasi-active term

Linearized voltage-dependent-current theory separates a channel contribution into both:

```text
resting / static membrane conductance
+
delayed voltage-dependent feedback
```

For an HCN-like restorative current, the delayed component provides negative feedback.

The first proxy had modeled only the delayed branch.

`hcn_quasiactive_probe.py` therefore ties both terms to the same local density field `d(x)`:

```text
v' = K L psi - damping*v - restoring*psi
     - d(x)*psi
     - mu*d(x)*z
     + source

z' = (psi-z)/tau_h
```

with `mu > 0` restorative.

In the harmonic domain,

```text
A(d,omega)
  = A0(omega)
    + diag[d_i * (1 + mu*a_h(omega))].
```

This is still deliberately minimal: a second-order reciprocal wave model plus a one-state quasi-active local material, not a conductance-based neuronal membrane.

---

## Development sweep

The complete quasi-active model was explored on seeds 560–563 over:

```text
g0        .005, .01, .02
ratio     5, 7, 10
channel tau  2, 4, 6
mu/static ratio  .5, 1, 2
omega     .02, .03, .04, .05, .06, .08
```

Profiles were compared with identical frozen anatomy:

```text
zero       no added material
uniform    same mean density everywhere
smooth     density increasing with graph distance from soma
shuffle    exact smooth-profile density values shuffled across occupied cells
reverse    exact same density histogram in reversed distance order
```

Several parameter regions showed the desired combination of:

1. reduced somatic phase spread;
2. positive distal-vs-proximal phase advance;
3. substantial retained local phase structure.

The highest aggregate development score was **not** selected automatically.

Body-by-body robustness instead froze:

```text
g0                  .005
soma->distal ratio  10x
tau_h                2
mu                    .5
omega                 .03 and .04
```

Across the four development bodies × two frequencies:

```text
smooth beats shuffled      8 / 8
smooth beats uniform       8 / 8
required distal lead > 0   8 / 8

mean gain vs shuffled      +.1757 rad
mean gain vs uniform       +.1559 rad
mean distal lead           +1.5399 rad
mean local retention        1.170
mean amplitude ratio
  smooth/shuffle             .961
```

The candidate and thresholds were then frozen in `HCN_QUASIACTIVE_CONFIRM_PREREG_V01.md` before seeds 568–579 were examined.

---

# Held-out confirmation — seeds 568–579

All 12 requested bodies bootstrapped, giving 24 body-frequency observations.

## Q0 — smooth placement beats shuffled and uniform material

Registered:

```text
mean gain vs shuffle > .05 rad
mean gain vs uniform > .05 rad
```

Observed:

```text
mean gain vs shuffle    +.12216 rad
mean gain vs uniform    +.10505 rad
```

**Q0 PASS.**

---

## Q1 — body-frequency robustness

Registered at least 75% positive observations for both controls.

Observed:

```text
smooth better than shuffle    20 / 24  (.833)
smooth better than uniform    19 / 24  (.792)
```

**Q1 PASS.**

The few losses are retained. The result is a strong population tendency, not a universal per-body theorem.

---

## Q2 — HCN-like distance sign

Registered:

```text
mean distal-minus-proximal phase advance > .30 rad
positive in at least 75% of observations
```

Observed:

```text
mean distal-minus-prox phase advance   +1.87770 rad
positive observations                    24 / 24
```

**Q2 PASS.**

This is the crucial difference from the rejected delayed-only proxy.

The morphology-indexed quasi-active material does not merely reduce an RMS metric. Distal input locations acquire systematically more phase advance than proximal locations relative to the no-material field.

---

## Q3 — local dendritic phase structure survives

Registered:

```text
mean local phase retention > .70
```

Observed:

```text
mean local phase retention      1.17553
median                           1.16927
```

**Q3 PASS.**

The local field is not globally flattened to obtain the somatic result. On average the smooth quasi-active profile actually leaves **more** local phase spread than the passive baseline while making the common somatic readout more coordinated.

This is the qualitative separation that motivated the experiment:

```text
local field can remain phase-rich
while
somatic projection becomes more phase-coherent.
```

---

## Q4 — amplitude control

Registered both pooled median soma-amplitude ratios between `.5` and `2`.

Observed:

```text
median smooth/shuffle amplitude ratio   .97743
median smooth/uniform amplitude ratio   .96149
```

**Q4 PASS.**

The phase result is not explained by simply shutting down or explosively amplifying the smooth-gradient condition.

---

## Q5 — both frozen frequencies contribute

Registered positive synchronization and phase-lead margins independently at `.03` and `.04`.

Observed:

```text
omega=.03
  gain vs shuffle     +.16990 rad
  gain vs uniform     +.14229 rad
  distal phase lead   +2.53002 rad

omega=.04
  gain vs shuffle     +.07441 rad
  gain vs uniform     +.06781 rad
  distal phase lead   +1.22538 rad
```

**Q5 PASS.**

The effect is stronger at `.03`, but `.04` independently clears the frozen criteria.

---

## Q6 — gradient direction matters

The reversed profile contains the same density values as the smooth profile, assigned in the opposite graph-distance order.

Registered:

```text
mean(reverse soma phase RMS - smooth soma phase RMS) > .05 rad
```

Observed:

```text
+.38310 rad
```

**Q6 PASS.**

The density histogram alone is insufficient. Its relation to morphology matters strongly.

---

# Confirmation verdict

```text
Q0 PASS   smooth beats shuffled and uniform
Q1 PASS   effect is robust across body-frequency observations
Q2 PASS   distal phase-advance sign is correct
Q3 PASS   local phase field remains structured
Q4 PASS   amplitude remains ordinary
Q5 PASS   both frozen low-band frequencies contribute
Q6 PASS   reversed gradient is substantially worse

7 / 7 PASS
```

The held-out result generalized with useful margin rather than barely crossing the thresholds.

---

## What has been earned

Within this model, we can now state:

> **Electrical material is part of the effective geometry. A morphology-indexed quasi-active restorative field can compensate location-dependent transfer timing at a common readout even while the distributed local field remains nonuniform and phase-rich.**

That is a materially stronger Geometric Neuron statement than the old delay-line picture.

The computational object is not merely

```text
shape G
+
signal psi
```

but

```text
morphology G
+
spatial electrical material M(x,omega)
+
moving field psi(x,t).
```

Morphology and material jointly determine the transfer coordinates seen by the soma.

---

## Connection to the physical-adjoint branch

The local material density is not just another hand-tuned knob.

For

```text
A(d,omega) x_j = b_j
H_j = e_s^T x_j,
```

`MATERIAL_ADJOINT_DERIVATION_V01.md` shows

```text
dH_j / dd_i
  = - c(omega) * y_i * x_{j,i}
```

where

```text
A^T y = e_s.
```

So the sensitivity of local membrane/material density is again a **local forward-field × returned-transpose-field product**.

The same reciprocal-medium logic that routed credit to structural bonds can therefore route credit to the electrical material laid over the morphology.

---

## Biological boundary

The correspondence should remain narrow.

Established biology motivates the *type* of intervention:

- CA1 HCN properties are spatially nonuniform;
- the HCN gradient can compensate location-dependent timing at the soma;
- HCN1 spatial distribution is itself activity dependent and reversible.

This repository has **not** shown:

- that the numerical `10x` density ratio maps to a measured CA1 ratio;
- that real HCN channels are optimized for this exact objective;
- that biological dendrites execute exact physical-adjoint learning;
- that the present one-state quasi-active model captures conductance-based HCN kinetics in detail.

The value of the biological comparison is that it supplied a falsifiable architecture. The first incomplete proxy failed that architecture; the corrected complete quasi-active version survived a held-out test.

---

## Next experiment

Stop hand-drawing channel gradients.

Start all occupied arbor cells with a **uniform material density under a fixed total-density budget**.

Then optimize a differentiable somatic phase-coordination objective using the exact material-adjoint sensitivity.

Only after learning should we inspect the resulting profile and ask:

```text
does learned density correlate with graph distance?
does it become distal enriched?
does it outperform shuffled learned-density controls?
does the same solution transfer across nearby frequencies?
```

If distal enrichment emerges from the objective, it becomes an output of the geometry/field/material interaction rather than an assumption copied from biology.

That is the next wall.
