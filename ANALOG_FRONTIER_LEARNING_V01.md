# Analog frontier learning v0.1 — the continuous credit loop learns, but some mechanism claims soften on holdout

This is the first experiment in the line that uses the soma-conditioned adjoint as an **optimizer**, not merely as a sensitivity audit.

The design variable is a fixed set of frontier bond maturities

```text
rho_e in [0,1]
k_e = k_bath + rho_e (k_arbor-k_bath).
```

All candidate frontier bonds start at `rho=0`. Three arms see the same body and candidate set:

```text
A  relinearized adjoint     recompute the true soma adjoint every iteration
B  frozen adjoint           reuse the initial gradient for all iterations
C  shuffled adjoint         recompute current gradients but permute them across bonds
```

Fixed schedule:

```text
eta = 0.01
40 iterations
up to 8 frontier candidates
```

No line search, optimizer tuning, momentum, material penalty, or post-hoc stopping.

## Discovery — seeds 240-251

```text
mean final Delta C
relinearized                  +0.031996     12/12 improved
frozen                        -0.055036
shuffled                      +0.007010

mean paired advantage
relinearized - frozen         +0.087032     10/12 bodies
relinearized - shuffled       +0.024986     10/12 bodies

median nondecreasing-step fraction
relinearized                   1.000

mean frontier material sum rho
relinearized                   0.12976
frozen                         0.99313
shuffled                       0.24422

initially positive candidates           44
later acquired negative true gradient   35 / 44 = 79.5%
```

All registered discovery D1-D4 criteria passed.

The discovery therefore suggested a strong picture: repeated relinearization learns a sparse/low-material analog frontier; stale gradients overshoot badly; spatially shuffled credit is much weaker; and many initially favorable bonds later turn unfavorable as the field reorganizes.

## Held-out confirmation — seeds 252-263

The confirmation criteria were deliberately tightened before these bodies were run.

```text
mean final Delta C
relinearized                  +0.023927     12/12 improved
frozen                        -0.023325
shuffled                      +0.009319

mean paired advantage
relinearized - frozen         +0.047252     10/12 bodies
relinearized - shuffled       +0.014608     11/12 bodies

median nondecreasing-step fraction
relinearized                   1.000

mean frontier material sum rho
relinearized                   0.10787
frozen                         0.55532
shuffled                       0.16447

initially positive candidates           35
later acquired negative true gradient   18 / 35 = 51.4%
```

### C1 — does the relinearized learner actually improve the objective?

Registered:

```text
mean Delta C > 0.020
at least 11/12 bodies improve
```

Observed:

```text
mean Delta C = +0.023927
12/12 improve
```

**C1 PASS.**

This is the central result. The local adjoint loop itself works on held-out bodies.

### C2 — is stale credit decisively worse?

Registered:

```text
mean relinearized-frozen > 0.050
at least 9/12 bodies favor relinearization
```

Observed:

```text
mean difference = +0.047252
10/12 favor relinearization
```

**C2 FAIL, narrowly on the mean-effect threshold.**

The direction replicated and frozen credit was negative on average, but the preregistered `0.050` margin was missed by `0.00275`.

### C3 — is correct spatial credit decisively better than shuffled credit?

Registered:

```text
mean relinearized-shuffled > 0.015
at least 9/12 bodies favor relinearization
```

Observed:

```text
mean difference = +0.014608
11/12 favor relinearization
```

**C3 FAIL, extremely narrowly on the mean-effect threshold.**

Again the direction is highly consistent, but the fixed threshold is the threshold.

### C4 — is small-step ascent stable?

Registered median nondecreasing-step fraction `>0.90`.

Observed:

```text
median = 1.000
```

**C4 PASS.**

### C5 — does relinearization usually reverse initially positive bond gradients?

Registered pooled reversal fraction `>0.60`.

Observed:

```text
18 / 35 = 0.5143
```

**C5 FAIL.**

So the discovery's particularly strong self-correction story did not replicate at the preregistered frequency. About half of initially favorable candidate bonds do reverse sign, which is still substantial, but not enough to support "usually" on the holdout.

## Formal confirmation score

```text
C1  PASS    the learner improves held-out objectives
C2  FAIL    stale-credit advantage misses mean threshold narrowly
C3  FAIL    spatial-credit advantage misses mean threshold extremely narrowly
C4  PASS    ascent is stable
C5  FAIL    gradient sign reversal is not frequent enough
```

This is **2/5 strict confirmation**, and should be reported exactly that way.

But the failures do not erase the central positive result. C1 and C4 were the direct tests of whether repeated local adjoint updates can optimize the analog frontier at all. They passed cleanly. C2 and C3 both replicated their predicted direction and body counts but missed deliberately stronger magnitude thresholds. C5 is the substantive mechanism softening.

## Combined 24-body description — not a new confirmatory test

Pooling discovery and confirmation only for descriptive context:

```text
relinearized mean Delta C              +0.02796     24/24 improve
frozen mean Delta C                    -0.03918
shuffled mean Delta C                  +0.00816

relinearized - frozen                  +0.06714     20/24 bodies
relinearized - shuffled                +0.01980     21/24 bodies

initially positive -> later negative    53/79 = 67.1%

mean sum rho
relinearized                             0.1188
frozen                                   0.7742
shuffled                                 0.2043
```

These pooled values were not preregistered and are not substitutes for the held-out verdict. They do show why the discovery and confirmation feel qualitatively consistent despite the stricter confirmation misses.

## What is now actually established

The strongest statement is narrower than the discovery suggested:

> **In this model's continuous frontier state space, repeated soma-conditioned adjoint updates reliably increase the smooth temporal-order objective on held-out bodies using small analog conductance changes.**

It is also fair to say:

- a frozen initial gradient is poor over the full 40-step path and is negative on average;
- spatially shuffled current gradients are weaker than the correctly assigned gradients in 11/12 held-out bodies;
- relinearized learning uses much less total frontier conductance than the frozen arm;
- gradient sign reversal is common but not held-out-confirmed as the dominant self-limiting mechanism.

## Why this closes an old loop

The v0.7-v0.9 line was trying to discover an eligibility tag for **binary structural events**. We now know why that was such a hard target:

```text
binary event
  -> huge conductance jump
  -> leaves local differential regime
  -> first-order credit can reverse before endpoint
```

The continuous formulation asks a different question:

```text
current field + soma objective
  -> adjoint sensitivity
  -> small local conductance update
  -> recompute field
  -> recompute sensitivity
```

That problem is tractable, and the held-out learner works.

## Interferometer + adjoint = one forward/backward geometry

The emerging architecture can now be written compactly:

```text
analog coupling geometry rho
        ↓
forward propagation operator A(rho)
        ↓
source histories h_A, h_B
        ↓
interference visibility + directional complex compatibility
        ↓
local square-law soma readout
        ↓
scalar task objective J
        ↓
adjoint seeded by that readout/objective
        ↓
local forward × backward bond overlap
        ↓
small update to rho
```

The soma is therefore doing two logically linked jobs:

1. **forward:** its local quadratic observation decides which interference is consequential;
2. **backward:** that same consequence supplies the terminal condition for structural sensitivity.

A useful sentence is:

> **Observation geometry and credit geometry are the same task boundary viewed in the forward and adjoint directions.**

## Next wall

The exact adjoint is still an algorithmic reverse pass. It is not yet a physical or biological mechanism.

The next clean question is therefore:

> **Can a reciprocal soma-launched retrograde wave approximate the exact adjoint bond-sensitivity map well enough to preserve analog frontier learning?**

That is where the old "retrograde carrier" idea can return in a much more precise form: not as scalar reward transported backward, but as a task-conditioned field whose local overlap with the forward field approximates the adjoint sensitivity.
