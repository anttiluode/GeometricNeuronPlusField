# Combined-error conformal echo closed-loop learning — development v0.1

Date: 2026-08-08

## Result in one sentence

> **With residual loss, spatial loss disorder, loss-calibration error, imperfect terminal time reversal, reverse-operator mismatch, and local gradient-readout noise present simultaneously, the conformal echo gradient trained graded frontier bonds successfully on all eight reused development arbors and retained 99.81% of the exact physical-gradient learner's mean task gain; permuting the same echo credit among candidate bonds removed the gain on most bodies.**

This is a development result on reused bodies 472–479. It is not a held-out confirmation and it still assumes that the `+` and `-` interference phase states experience the same reverse operator during one gradient acquisition. Differential pass-to-pass drift is tested separately.

---

# Why this test matters

The preceding damping-gauge branch established increasingly realistic but mostly static statements:

```text
exact damped task -> exact reversible gauge compiler
uniform residual physical loss -> exact conformal echo compensation
loss calibration error -> gradient direction remains aligned
spatial loss disorder -> gradient direction remains aligned
terminal time-mirror error -> gradient direction remains aligned
reverse-operator mismatch -> gradient direction remains aligned
local gradient-map noise -> gradient direction remains aligned
```

A high gradient-map correlation is not sufficient for a physical-learning claim.

The next criterion is functional:

```text
does the damaged physical-gradient surrogate actually improve the task
when it is used repeatedly to change the body?
```

`damping_gauge_combined_closed_loop.py` tests that question.

---

# Frozen development error cocktail

The echo learner used all of the following simultaneously:

```text
mean residual physical loss eps          .005
cellwise loss coefficient of variation   .20
scalar loss calibration error            +.05
terminal time-mirror alpha               .95
reverse-operator coupling drift sigma    .02
local gradient-readout noise RMS         .05
```

Interpretation:

```text
20% spatial variation around residual loss
5% error in the common loss factor used for conformal envelopes
5% terminal reverse-state error
2% multiplicative reverse-operator mismatch
5% RMS additive error in the final local candidate-gradient values
```

The residual-loss disorder was fixed for each body across learning iterations, representing fabrication/device variation.

The reverse-operator perturbation and readout noise were redrawn during learning.

---

# Learning problem

For each mature FunctionalArbor body:

1. select up to eight frontier bonds adjacent to the frozen arbor;
2. initialize their graded participation `rho=0`;
3. map `rho` continuously between mature-bath and arbor coupling strengths;
4. optimize the original target-vs-distractor finite-time contrast for 24 updates;
5. relinearize after every update.

Step size:

```text
eta = .01
```

All learning arms receive the same candidate set and physical loss realization.

---

# Four arms

## 1. nominal_exact

The repository's ordinary exact adjoint gradient for the intended damped numerical task.

This is a reference for the original model, not the actual imperfect physical implementation.

## 2. physical_exact

A digital reverse-mode gradient of the **actual compiled physical device**, including its fixed spatial residual-loss disorder.

This is the upper-bound reference for training that physical implementation.

## 3. echo_combined

The memory-free conformal echo estimate under the full error cocktail.

The forward physical task is evaluated in the actual lossy compiled core. The reverse trajectory uses the calibrated scalar conformal envelope, imperfect terminal time reversal, reverse-operator mismatch, and noisy local gradient readout.

## 4. shuffled_echo

The same echo candidate-gradient values are randomly permuted among the frontier bonds before every update.

Thus its gradient histogram/norm is inherited from the echo calculation, but the local credit assignment is destroyed.

---

# Development result — seeds 472–479

Per body, physical-task contrast improvement after 24 updates:

```text
seed    physical exact    combined echo    shuffled echo    mean echo corr
472       +.0274             +.0268           +.0263          .9973
473       +.0035             +.0035           -.0039          .9562
474       +.0188             +.0183           -.0100          .8948
475       +.0225             +.0224           -.0230          .9837
476       +.0639             +.0642           -.0098          .9922
477       +.0358             +.0366           -.0021          .9973
478       +.0131             +.0129           -.0026          .9974
479       +.0199             +.0199           +.0099          .9935
```

The static echo-vs-exact candidate-gradient correlation is useful diagnostically, but the learning outcome is the primary result.

---

## Physical task gain

```text
arm                 mean delta C     improved bodies
nominal exact        +.0186780          8 / 8
physical exact       +.0256259          8 / 8
combined echo        +.0255771          8 / 8
shuffled echo        -.0018858          2 / 8
```

The corrupted echo learner therefore retained

```text
mean echo gain / mean exact-physical gain
    = .998098
```

or about **99.81%** of the exact physical-gradient learner's mean improvement.

This is not a statement that every intermediate gradient was exact. One body had mean echo correlation below `.90`, yet the closed loop still improved nearly as much as the exact physical-gradient arm.

---

## Echo placement matters

Comparing body-by-body physical-task gains:

```text
mean(echo gain - shuffled-echo gain)   +.0274629
combined echo beats shuffled echo       8 / 8
```

The shuffled-credit control therefore rejects the interpretation that any similarly sized frontier perturbation would improve the body.

Seed 472 is retained as an instructive exception to the stronger statement “shuffling always hurts”: its particular shuffled trajectory also improved substantially. The population-level and 8/8 paired comparison remain strongly in favor of the correctly placed echo credit.

---

# Transfer back to the intended nominal task

Although the echo learner was trained using the actual imperfect physical forward objective, its learned geometry also improved the original nominal damped-model contrast:

```text
arm                 mean nominal delta C   improved bodies
physical exact          +.0292948              8 / 8
combined echo           +.0292841              8 / 8
shuffled echo           -.0012331              2 / 8
```

So in this development regime the physical implementation errors did not drive the learner toward a hardware-specific solution that damages the intended task.

---

# What has been earned

Before this experiment the strongest statement was

```text
the damaged echo gradient points approximately the right way.
```

The stronger development statement is now

> **Under a combined non-pristine hardware-error model, repeated conformal-echo credit updates can train a changing distributed wave body, and their closed-loop gain is essentially the same as a digital exact gradient of that physical device.**

That is the first functional learning result in the damping-gauge branch.

---

# What has NOT been paid yet

The experiment still gives the `+` and `-` energy-interference states the same reverse operator during one gradient acquisition.

This means a static 2% reverse mismatch is present, but not independent drift between the two phase states.

That distinction matters because

```text
|w+a|^2 - |w-a|^2
```

cancels large self-energy terms only when the two measurements refer to sufficiently identical underlying fields.

If the mesh changes between separate `+` and `-` trials, self-energy leakage can dominate the much smaller cross term.

`damping_gauge_pass_drift_probe.py` therefore attacks differential `+/-` pass drift before any fresh held-out confirmation is frozen.

Other unpaid hardware assumptions remain:

```text
physical implementation of the terminal conformal time mirror
synchronization of reverse source/error envelopes
finite detector bandwidth
finite energy-integrator dynamic range
component nonlinearities
frequency-dependent loss not captured by one recurrence coefficient
large-scale calibration of local tunable couplings
```

---

# Prior-art boundary

Hamiltonian Echo Backpropagation and Recurrent Hamiltonian Echo Learning already establish physical time-reversal/echo approaches for finite-time temporal credit in non-dissipative Hamiltonian systems.

Scattering Backpropagation addresses physical training in driven-dissipative systems through scattering experiments, but around steady scattering responses rather than this finite transient-history problem.

The current GeometricNeuronPlusField result should therefore not be described as inventing physical echo learning.

The narrower technical question remains:

> **Can a finite-time dissipative scattering task be conformally compiled into a physically echo-trainable representation whose local broadband credit remains useful under realistic non-Hamiltonian implementation errors and without storing the local trajectory?**

The present development result is affirmative under one combined error cocktail, subject to the differential-pass wall above.

## Wall sentence

> **The conformal echo is no longer merely an exact algebraic gradient identity: with multiple hardware imperfections present simultaneously it trains the distributed body in closed loop almost as well as the exact physical gradient. The next wall is whether the local cross term can be measured without requiring two unrealistically identical reverse passes.**
