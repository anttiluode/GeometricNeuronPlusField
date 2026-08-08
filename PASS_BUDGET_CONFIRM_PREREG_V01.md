# Local relinearization pass-budget confirmation preregistration v0.1

`COORDINATE_CURVATURE_CONFIRM_V01.md` established that local bond coordinates are much more curved than direct spectral coordinates. `PASS_BUDGET_DEV_V01.md` then found on reused bodies that halving the local task-space trust radius from `.01` to `.005` and doubling the number of relinearized iterations from 40 to 80 substantially improved held-out performance, although it still did not catch the 40-step spectral reference. A further halving to `.0025 x160` did not help.

This fresh run freezes that pass-budget comparison.

## Protocol

Fresh bodies: seeds **336-347**.

Common task:

```text
P = 8 selected coordinates
train lags = 16,20,24
test lags  = 14,18,22,26
```

The graph coordinate set is selected **once** at the common base state using the `.01` scale-invariant selection rule and then reused for all graph schedules.

Schedules:

```text
G40     local bonds       delta=.0100   40 iterations
G80     local bonds       delta=.0050   80 iterations
G160    local bonds       delta=.0025  160 iterations

F40     free spectral     delta=.0100   40 iterations
```

Every schedule uses the scale-invariant update from `matched_tuner_trust.py`: each step is normalized by the current task-space Jacobian so its pre-clipping predicted RMS change has the declared `delta`.

Nominal cumulative trust budget is `.4` in all schedules.

## Registered tests

### C1 — extra local relinearization helps

Primary development-selected comparison:

```text
mean [C_test(G80) - C_test(G40)] > .020
and G80 > G40 in at least 8/12 completed bodies.
```

### C2 — twofold extra local pass budget does not fully erase the spectral advantage

```text
mean [C_test(F40) - C_test(G80)] > .030
and F40 > G80 in at least 9/12 completed bodies.
```

This criterion can fail in either direction; if local G80 catches the spectral reference, that is the more interesting result and will be reported as such.

### C3 — over-relinearization is not monotonically better

Development found G160 slightly worse than G80. Confirm the weaker statement:

```text
mean [C_test(G160) - C_test(G80)] < .020.
```

This does not require G160 to be worse, only that another doubling of local passes does not produce a large additional gain at fixed nominal trust budget.

### C4 — all schedules remain useful

Mean held-out improvement over the common base must be positive for G40, G80, G160, and F40. This is a sanity condition, not a superiority claim.

## Physical pass interpretation

One relinearized update of a reciprocal physical medium requires at least a task forward measurement and an adjoint/backward measurement; exact implementation details can add extra interference/calibration measurements. Therefore iteration count is a direct proxy for repeated physical gradient acquisition, even though it is not yet a complete energy/time cost model.

If C1-C3 pass, the earned hardware statement is:

> **Curvature imposes a training-pass tax on local physical coordinates: more frequent reciprocal gradient refresh improves them, but simply taking proportionally more tiny local steps does not make the curved coordinate chart equivalent to direct spectral control.**

No schedule or threshold will be changed on seeds 336-347 after inspection.
