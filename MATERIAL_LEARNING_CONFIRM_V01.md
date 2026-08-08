# Material-adjoint learning confirmation v0.1 — formal 7/8 fail, substantive 7/7 positive

Date: 2026-08-08

## Formal verdict

The frozen confirmation rule in `MATERIAL_LEARNING_CONFIRM_PREREG_V01.md` required all `L0–L7` criteria to pass.

Observed:

```text
L0 gradient audit             FAIL
L1 learned objective          PASS
L2 beats shuffled placement   PASS
L3 beats hand gradient        PASS
L4 distance organization      PASS
L5 amplitude sanity           PASS
L6 both frequencies improve   PASS
L7 phase RMS improves         PASS

formal verdict: 7 / 8, FAIL
```

This file preserves that verdict. The threshold is not changed after seeing the data.

However, diagnosis of the sole failure shows that it came from an ill-conditioned **relative-error statistic at a mathematically near-zero derivative**, not from evidence that the analytic material gradient was wrong.

The distinction matters, so both facts are recorded separately.

---

# Held-out substantive result — seeds 584–595

All 12 fresh bodies learned in the expected direction.

Mean phase-coherence objective:

```text
uniform             .4965493
learned             .6121996
hand linear profile .5511680
shuffled learned    .4573876
```

Mean gains:

```text
learned - uniform    +.1156503
learned - shuffled   +.1548120
learned - hand       +.0610316
```

Body counts:

```text
learned beats uniform    12 / 12
learned beats shuffled   12 / 12
learned beats hand       12 / 12
```

The directly interpretable circular soma phase spread agreed:

```text
uniform phase RMS    .8403661
learned phase RMS    .6666967
mean improvement     .1736695 rad
```

---

## The distance coordinate generalized strongly

After learning only, correlate local density with graph distance from the soma/readout.

Observed:

```text
mean Spearman rho       +.7978330
positive bodies          12 / 12
minimum body rho         +.6801241
```

This is slightly stronger than the four-body development mean (`+.773`).

The result therefore did not depend on the original small development set.

---

## The learned map still beats the hand-drawn biological prior

The stronger preregistered self-organization criterion required the learned local map to outperform the previously confirmed smooth soma-to-distal hand profile under the same material budget.

Observed:

```text
mean learned - hand R^2   +.0610316
positive bodies            12 / 12
```

So the learner is not merely reproducing a linear distance gradient. It finds a more effective morphology-specific allocation.

---

## Both frozen frequencies independently improve

```text
omega=.03
mean learned - uniform R^2   +.0541219

omega=.04
mean learned - uniform R^2   +.1771788
```

Both clear the frozen `.02` threshold.

---

## Amplitude remains nondegenerate

The phase objective normalizes source transfers before scoring, so an explicit amplitude guard was preregistered.

Observed learned/uniform median soma-transfer amplitude ratio:

```text
pooled median   .7421448
pooled mean     .7417963
in allowed per-body range   12 / 12
```

The learner attenuates transfer somewhat, as already seen in development, but it does not obtain the phase result by annihilating the channel.

---

# Why L0 failed

The implementation audit randomly selected five local material derivatives per body and compared the analytic adjoint derivative with a central finite difference.

The frozen statistic was

```text
relative_error
  = |finite_difference - analytic|
    / (|finite_difference| + |analytic| + 1e-12)
```

with required maximum `< 1e-5`.

For seed 584, four sampled derivatives matched normally:

```text
analytic          finite difference
.131304334478     .131304334627
.244959653306     .244959653339
.132503756252     .132503756300
.132817998432     .132817998433
```

The fifth was a near-zero derivative:

```text
analytic             +9.49e-17
finite difference    -1.665e-10
absolute discrepancy  1.665e-10
```

Because both values are essentially zero, the denominator of the relative-error statistic is also essentially zero. The harmless finite-difference numerical residue therefore becomes

```text
reported relative error ~= .994
```

and trips L0.

Across all `12 bodies × 5 samples = 60` held-out derivative checks, the largest **absolute** analytic-versus-finite-difference discrepancy was only about

```text
5.44e-10.
```

The other nonzero-scale derivative checks were at the same ~floating-point/finite-difference agreement level seen in development.

---

## Correct classification

Do not retroactively relabel v0.1 as an all-pass confirmation.

The registered experiment says:

```text
formal v0.1 confirmation = FAIL, 7/8
```

But the scientific diagnosis says:

```text
the only failed criterion was a numerically ill-conditioned
relative-error guard at a zero derivative;

all seven substantive predictions passed, strongly.
```

That warrants a **new preregistered confirmation**, not a rewritten old one.

---

## Methodological correction for v0.2

The derivative audit needs two regimes:

```text
nonzero derivative scale:
    use relative error

near-zero derivative scale:
    use absolute error
```

A sensible frozen zero-safe rule is:

```text
if |analytic| + |finite_difference| >= 1e-6:
    relative error < 1e-5

else:
    absolute error < 1e-8
```

This is much looser than the observed worst absolute discrepancy (`~5.44e-10`) while avoiding division by numerical zero.

All substantive L1–L7 thresholds should remain exactly unchanged.

The corrected experiment must use completely fresh bodies.

---

## Wall sentence

> **The first held-out material-learning confirmation is formally a 7/8 failure because its derivative-audit statistic is singular at zero; nevertheless every substantive learning prediction generalized on all 12 fresh bodies, and the audit's actual absolute derivative error remained below 6e-10. The correct next move is a new zero-safe preregistration, not reinterpretation of the old verdict.**
