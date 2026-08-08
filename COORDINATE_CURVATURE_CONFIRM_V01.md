# Coordinate-curvature confirmation v0.1 — local physical coordinates bend far faster than spectral coordinates

`MATCHED_TUNER_TRUST_CONFIRM_V01.md` established small-step parity between local bond and direct spectral coordinates once task-space step size is normalized, while the relative performance shifts toward spectral coordinates as the allowed trust radius grows.

`COORDINATE_CURVATURE_CONFIRM_PREREG_V01.md` then preregistered the direct mechanism test: move one selected coordinate at a time by an amount whose **initial training-lag tangent** predicts a fixed RMS task change, and measure how far the actual seven-lag response departs from that first-order prediction.

Fresh bodies: seeds **324-335**. Eight graph and eight spectral coordinates per body, 96 probes per arm.

Primary curvature metric:

```text
E = ||Delta C_actual - Delta C_pred|| / ||Delta C_pred||
```

with

```text
C = [C14,C16,C18,C20,C22,C24,C26].
```

## Fresh result

Median relative linearization error:

```text
delta       local graph        free spectral      graph/free
.001          0.08019             0.00406            19.8x
.0025         0.22026             0.01001            22.0x
.005          0.46363             0.01957            23.7x
.010          0.79974             0.03626            22.1x
```

No tested probe in either arm was clipped by a parameter bound.

## C1 — curvature gap already at the smallest step

Registered:

```text
median E_graph / median E_free > 5
median E_graph > .015
```

Observed at `.001`:

```text
19.76x
E_graph = .08019
```

**C1 PASS.**

## C2 — large-step curvature gap

Registered at `.010`:

```text
ratio > 5
median E_graph > .15
```

Observed:

```text
22.06x
E_graph = .79974
```

**C2 PASS.**

## C3 — graph curvature grows strongly with trust radius

Registered:

```text
E_graph(.010) - E_graph(.001) > .10
```

Observed:

```text
+0.71955
```

**C3 PASS.**

## C4 — spectral coordinates remain comparatively straight

Registered:

```text
median E_free(.010) < .10
```

Observed:

```text
0.03626
```

**C4 PASS.**

## C5 — no clipping confound

Registered fewer than 5% clipped probes.

Observed:

```text
0 / 384 graph probe-steps clipped
0 / 384 free probe-steps clipped
```

**C5 PASS.**

All five registered mechanism tests pass.

## An important nuance: direction often survives longer than amplitude

Median tangent cosine remains high even as relative linearization error grows:

```text
at delta=.010
median cosine graph       0.99434
median cosine free        0.99998
```

So for the typical selected local coordinate, the finite response is not immediately pointing somewhere completely unrelated. Much of the curvature initially appears as **wrong gain / wrong distance along the intended functional direction**, with a minority of local coordinates eventually suffering severe directional failures that pull down the mean cosine.

That is exactly why relinearization works: the local direction can remain useful while its first-order step size becomes unreliable.

## The optimization geometry is now explicit

Three independent results now line up:

```text
1. one-cell / large bond edits leave the adjoint's linear regime;
2. random bond response curves are strongly nonmonotone and often have interior optima;
3. matched local and spectral task tangents are similarly useful at small trust radius,
   but local finite steps depart from their tangent ~20x faster in median relative error.
```

This earns a sharper statement:

> **The main disadvantage of local physical bond coordinates on this task is not first-order direction. It is curvature. Direct modal coordinates are close to the natural coordinates of the linear dynamics; local bond coordinates are a highly curved physical chart over those modal dynamics.**

## Hardware meaning

This is not necessarily bad news for local physical hardware. The exact in-situ adjoint is cheap in a reciprocal medium, so a local implementation can compensate for curvature by using **smaller updates and more frequent physical relinearization**.

That makes the next practical quantity clear:

> **How many extra forward/adjoint passes does local physical tuning need to match a larger-step spectral optimizer?**

That is a training-cost benchmark rather than an expressivity benchmark.
