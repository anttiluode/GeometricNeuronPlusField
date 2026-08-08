# True nonreciprocity v0.1 — same-medium physical backprop fails, transpose propagation restores the exact adjoint

The reciprocal-adjoint result had an obvious logical loophole: perhaps time-reversing the soma derivative and sending *some* wave backward through the body was enough, with reciprocity merely convenient.

This experiment removes that loophole.

A genuinely nonsymmetric local wave operator was introduced, and the return waveform was propagated either through the **same nonsymmetric operator** or through its **actual transpose**.

The result is clean:

> **same-H replay progressively stops being the adjoint as nonreciprocity grows; transpose replay remains the exact adjoint to floating-point precision.**

`NONRECIPROCAL_ADJOINT_PREREG_V01.md` froze all thresholds before fresh bodies 460–471 were run.

## Operator

Start from the ordinary symmetric weighted Laplacian `L` and add a local real skew-symmetric nearest-neighbour operator `A`:

```text
H_beta = L + beta A
A^T = -A
```

Therefore

```text
H_beta^T = L - beta A = H_-beta.
```

The local skew pattern is random but fixed per body.

The skew background is **not** itself optimized. The audited structural coordinates remain the usual symmetric bond conductances. This isolates the effect of using the wrong versus correct adjoint propagation operator.

This is an abstract nonreciprocity model, not a calibrated optical/microwave isolator or circulator model.

## Three return paths

For each forward task under `H_beta`:

```text
exact
  explicit discrete adjoint using H_beta^T

same-H physical replay
  reverse the soma derivative in time
  and send it through H_beta again

transpose replay
  send the same reversed waveform through H_beta^T = H_-beta
```

The primary observable is the gradient map with respect to the symmetric bond conductances.

## C0 — reciprocal baseline

At beta=0:

```text
mean same-H gradient corr       1.000000000000
mean same-H relative L2         2.09e-15
```

Registered requirements were `corr > .999999` and `rel-L2 < 1e-10`.

**C0 PASS.**

## C1 — transpose replay remains exact after nonreciprocity is added

At beta=.10:

```text
mean transpose gradient corr    1.000000000000
mean relative L2                2.32e-15
```

At beta=.20:

```text
mean transpose gradient corr    1.000000000000
mean relative L2                2.48e-15
```

No body approached the registered `1e-10` error ceiling.

**C1 PASS.**

The transpose positive control is therefore not approximate: it reproduces the exact discrete adjoint at machine precision.

## C2 — same-H replay is already wrong at beta=.10

Registered:

```text
mean same-H corr < .95
>= 10/12 bodies corr < .98
```

Observed:

```text
mean same-H corr                .903679
bodies corr < .98               12 / 12
mean relative L2                .693863
```

Individual correlations ranged from about `.7338` to `.9645`.

**C2 PASS.**

## C3 — failure strengthens at beta=.20

Registered:

```text
mean same-H corr < .88
mean relative L2 > .50
```

Observed:

```text
mean same-H corr                .734103
mean relative L2               3.316589
```

**C3 PASS.**

The same-return operator is now badly wrong even though the transpose replay is still exact.

## C4 — transpose specifically rescues the gradient

Registered correlation rescue:

```text
beta=.10    transpose - same-H > .03
beta=.20    transpose - same-H > .10
```

Observed:

```text
beta=.10    +.096321
beta=.20    +.265897
```

**C4 PASS.**

## Formal verdict

```text
C0 PASS   reciprocal same-medium identity recovered
C1 PASS   transpose replay exact under nonreciprocity
C2 PASS   same-H replay wrong at beta=.10
C3 PASS   same-H failure much larger at beta=.20
C4 PASS   using the transpose specifically restores the gradient
```

**5 / 5 registered criteria pass.**

## Descriptive failure curve

Fresh 12-body mean same-H bond-gradient correlations:

```text
beta      same-H corr      transpose corr
0          1.0000            1.0000
.02         .9972            1.0000
.05         .9802            1.0000
.10         .9037            1.0000
.20         .7341            1.0000
.30         .5389            1.0000
.40         .4535            1.0000
.60         .2645            1.0000
```

So the effect is not a binary switch at infinitesimal asymmetry. Small nonreciprocity gives a useful approximation, just as the pass-mismatch experiment suggested. But the causal direction is unmistakable: the farther `H` moves from `H^T`, the less the same-device return resembles the adjoint.

## This sharpens the reciprocity trade

The correct engineering statement is now:

> **Reciprocity buys a zero-extra-operator physical adjoint because the forward operator is already its own transpose.**

That is stronger and more precise than “backward waves somehow carry credit.”

For a nonreciprocal system, backpropagation itself is not impossible. The mathematical adjoint still exists. What is lost is the cheap identity

```text
physical forward operator = physical adjoint spatial operator.
```

A nonreciprocal machine therefore needs some way to realize `H^T`/`H†`:

- a separately configured reverse device;
- reversible/nonreciprocal elements whose control can be transposed;
- a calibrated digital model;
- a second physical mesh implementing the adjoint;
- or an approximation good enough for the optimizer.

## Relation to the earlier “arrow vs gradient” sentence

This experiment earns the **gradient half** directly.

The older passive-medium tests showed that the reciprocal FunctionalArbor did not spontaneously create the desired temporal arrow. The current result shows that the same transpose symmetry that limits directional asymmetry also makes the adjoint physically cheap.

So a careful wall sentence is:

> **Reciprocity removes one source of intrinsic directionality and, in exchange, makes forward and adjoint propagation the same spatial operation.**

That is not a universal theorem that every nonreciprocal system has a useful computational arrow, nor that every reciprocal system lacks all temporal asymmetry. It is the precise trade demonstrated by this model family.

## Hardware consequence

The hardware pivot is now more constrained:

```text
if the mesh is reciprocal:
    one calibrated operator can serve computation + adjoint transport

if the mesh is nonreciprocal:
    the return path must implement the transpose, not merely reuse the forward operator
```

That places reciprocity in the **cost model**. It is not merely a material property.

A reciprocal local scattering mesh potentially saves an entire separately calibrated adjoint network.

## Next wall

The remaining local-readout issue is now the sharper one.

The broadband learner uses `K=8/16` complex local spectral products. Numerically they are simply multiplied. A physical protocol must implement

```text
Re[conj(U_k) V_k]
```

from measurable intensities/phase steps, count how many passes or detector samples are required, and state what local analog memory is actually needed.

That is the next place where the existing Hughes/Pai in-situ-interference machinery should be adapted explicitly to the transient K-bin setting rather than reinvented abstractly.

## Wall sentence

> **Time reversal is not the magic; transpose propagation is. Reciprocity is valuable because it makes the transpose free.**
