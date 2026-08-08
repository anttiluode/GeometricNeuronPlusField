# Analog frontier learning held-out confirmation preregistration v0.1

Discovery on fresh seeds 240-251 was the first direct optimization test of the soma-conditioned adjoint in the continuous structural state space.

Discovery:

```text
relinearized mean Delta C             +0.031996   (12/12 improved)
frozen-gradient mean Delta C          -0.055036
shuffled-gradient mean Delta C        +0.007010

relinearized - frozen                 +0.087032   (10/12 bodies)
relinearized - shuffled               +0.024986   (10/12 bodies)
median relinearized monotone fraction  1.000

mean frontier material sum rho
relinearized                           0.1298
frozen                                 0.9931
shuffled                               0.2442

initially positive candidates          44
later acquired negative true gradient  35
```

All discovery D1-D4 criteria passed. This document freezes held-out criteria before any new body is run.

## Held-out bodies

```text
seeds 252-263
```

No algorithmic changes:

```text
max frontier candidates = 8
eta                     = 0.01
iterations              = 40
lag                     = 20
steps                   = 210
```

Same three arms: relinearized true adjoint, frozen initial adjoint, and relinearized-but-spatially-shuffled adjoint.

## Registered confirmation criteria

### C1 — repeated local adjoint updates learn the frontier

PASS if:

```text
mean final Delta C_relinearized > 0.020
and at least 11/12 bodies improve.
```

### C2 — stale credit is worse

PASS if:

```text
mean (Delta C_relinearized - Delta C_frozen) > 0.050
and at least 9/12 bodies favor relinearization.
```

### C3 — spatial assignment is necessary

PASS if:

```text
mean (Delta C_relinearized - Delta C_shuffled) > 0.015
and at least 9/12 bodies favor the true spatial adjoint.
```

### C4 — small-step ascent remains stable

PASS if the median body fraction of nondecreasing relinearized iterations is `>0.90`.

### C5 — relinearization is genuinely self-correcting

Among candidate bonds whose **initial** directional gradient is positive, count those that later acquire a negative true gradient during the relinearized trajectory.

PASS if the pooled reversal fraction is `>0.60`.

This tests the mechanism suggested by the one-bond response curves: the optimizer should not merely keep pushing an initially favored bond; it should discover its turning point and back away.

## Descriptive controls

Also report:

- final material `sum rho` in all arms;
- objective gain per unit frontier material;
- final `rho` distribution;
- initial-gradient/final-rho correlation;
- per-body objective trajectories.

No material-efficiency threshold is preregistered because the three projected arms can end with different total `rho` after clipping.

## Interpretation fixed in advance

If C1-C5 pass, the continuous-state credit loop is no longer only a derivative identity. It is a working learner:

```text
soma objective
  -> adjoint field
  -> local bond sensitivities
  -> small analog updates
  -> recomputed field
  -> recomputed credit
  -> self-limiting conductance tuning.
```

The next unresolved step would then be physical plausibility: can a reciprocal/retrograde wave approximate the exact algorithmic adjoint well enough to preserve learning?
