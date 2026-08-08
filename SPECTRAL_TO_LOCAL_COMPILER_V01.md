# Spectral-to-local compiler v0.1 — the winning spectral directions are functionally local-looking, but much larger per unit parameter

`SPECTRAL_COORDINATE_LOCALITY_CONFIRM_V01.md` showed that a pure modal pole/residue coordinate is strongly nonlocal in the **operator basis** of the nearest-neighbour grid.

That does not tell us how hard it is to reproduce the same *task effect*. `spectral_to_local_compiler.py` therefore compared tangent vectors in a 25-dimensional functional space:

```text
C(lag), lag = 8,9,...,32.
```

For each of the 96 winning free-modal coordinates from the opened benchmark bodies 288-299:

1. orient the free coordinate in its training-objective improving direction;
2. compute its derivative vector across all 25 lags;
3. compute the derivative vector of every one of the 1860 local bonds, oriented into that bond's feasible one-sided direction at the binary base state;
4. greedily fit the free tangent using nonnegative local-bond directions.

This is descriptive, not a new held-out test.

## Surprise: one local bond often matches the *shape* extremely well

With only one selected local bond:

```text
mean relative tangent residual              0.04269
median relative tangent residual            0.02762
fraction residual < 0.10                    0.90625
fraction residual < 0.20                    0.97917
```

So 87/96 winning free-coordinate effects across the entire lag family can be matched in *shape* to within 10% by the tangent of a single feasible local bond.

This sharply changes the interpretation of the earlier operator-support audit. A pure modal operator perturbation is globally dense, but **its effect on this low-dimensional task manifold is usually almost collinear with some local physical perturbation**.

With more local directions, shape residuals become tiny:

```text
local directions     median relative residual
1                     0.02762
2                     0.01755
4                     0.00590
8                     0.000618
16                    0.000146
24                    0.000078
```

The multi-column fits become highly ill-conditioned, so coefficient sums at large `k` should not be interpreted as hardware costs. The one-bond result is much cleaner.

## But amplitude exposes a benchmark-normalization problem

The fitted coefficient for one local bond answers:

```text
one unit of free-coordinate tangent
≈ coefficient × one unit of local rho tangent.
```

Across all 96 selected coordinates:

```text
median one-bond coefficient                 14.90
```

Two thirds have coefficient `<=100`, meaning that a free-coordinate step of `0.01` could, at the tangent level, be matched in magnitude by at most a full `[0,1]` local-rho move using the best single bond.

This matters because `MATCHED_TUNER_CONFIRM_V01.md` used the same numerical normalized step `eta=0.01` in **different coordinate systems**:

```text
local arm      max ~0.01 rho change per iteration
free arm       max ~0.01 pole/residue-coordinate change per iteration.
```

Those numbers are not physically or functionally commensurate. A median free `0.01` tangent corresponds to roughly `0.149` local-rho units under the best one-bond tangent fit — about fifteen times the local arm's permitted per-iteration move.

The exact factor varies enormously across coordinates, and some fits are ill-conditioned, so `15x` is not a universal conversion constant. But the existence of the scale mismatch is unambiguous.

## Canonical correction to the matched-tuner result

The held-out result that the chosen free-coordinate optimizer achieved a larger objective remains numerically true.

What must be demoted is the phrase:

> “free spectral coordinates decisively beat local bonds **per scalar**.”

The experiment matched scalar count and iteration count, but it did **not** match functional trust-region size, physical parameter range, or induced transfer-function step size. Therefore it cannot by itself establish a parameter-efficiency theorem.

The more accurate statement is:

> **Under the specific coordinate normalization and `eta=0.01` optimizer used in v0.1, direct spectral coordinates train much faster/better. Yet the winning free tangents are usually almost reproducible by a single local-bond tangent across the lag task manifold, at a different parameter scale.**

## Next fair benchmark

The next comparison must be invariant to arbitrary coordinate scaling.

A clean protocol is:

1. select coordinates by **normalized task-space alignment**, not raw derivative magnitude;
2. at every iteration normalize each coordinate by the norm of its current multi-lag tangent;
3. give both arms the same predicted RMS change in the training-output vector per iteration (a task-space trust radius);
4. relinearize after every step;
5. enforce each arm's actual physical/declared parameter bounds;
6. test the same held-out lags.

Only after that should we make a claim about geometry vs spectral coordinates as optimization parameterizations.

## Why this is useful

This is exactly what a compiler analysis is supposed to reveal. The spectral coordinate is globally nonlocal as an operator, yet the task only observes a very small slice of the transfer function. On that slice, local physical directions can shadow the oracle spectral directions surprisingly well.

That reopens the hardware question in a much better form:

> **How much local actuation range and how many local tuners are required to reproduce a useful spectral task-space move?**

That is measurable without pretending the coordinate systems have the same units.
