# Single-source transfer decomposition v0.1 — the missing factor is the complex cross-source term

Canonical workflow: `transfer-decomposition`, run `31246490320`.

Fresh bodies were used:

```text
discovery      seeds 48-59
confirmation   seeds 60-71
```

No growth, wave, gain-set, lag, source position, or readout parameters were changed after preregistration.

## What was decomposed

For every body, every occupied cell, and each of five source-gain conditions, the two single-source complex responses were recorded:

```text
h_A(x,t)
h_B(x,t)
```

The actual pair simulation was then compared with two reconstructions.

Coherent reconstruction:

```text
psi_T = h_A(t) + h_B(t-tau)
psi_D = h_B(t) + h_A(t-tau)
```

Incoherent/envelope-only reconstruction:

```text
P_T = |h_A(t)|^2 + |h_B(t-tau)|^2
P_D = |h_B(t)|^2 + |h_A(t-tau)|^2
```

The only thing removed in the incoherent reconstruction is the cross-source complex term

```text
2 Re[h_1 conj(h_2)].
```

The historical peak-power order contrast was computed for actual, coherent and incoherent traces.

## P1 — coherent synthesis reproduces the real pair almost exactly

```text
                         discovery          confirmation
mean signed corr         0.999999840        0.999999840
median                   0.999999845        0.999999852
mean MAE                  0.00002224         0.00002262
registered P1                 PASS               PASS
```

This is effectively an exact reduction of the pair task to the two independently measured source transfer histories in this operating regime.

The mild field saturation contributes only a tiny residual at these gains.

## P2 — deleting the complex cross term loses most of the computation

Body-wide contrast prediction error:

```text
                         discovery          confirmation
MAE coherent              0.0000222          0.0000226
MAE incoherent            0.0350019          0.0346121
incoh - coh MAE          +0.0349796         +0.0345895
positive / negative        12 / 0             12 / 0
sign p                     0.000488           0.000488
registered P2                 PASS               PASS
```

Mean absolute order selectivity over body cells:

```text
                         discovery          confirmation
actual                     0.04803            0.04772
coherent                   0.04804            0.04773
incoherent                 0.01345            0.01362
```

The envelope-only model retains only about 28% of the observed absolute contrast magnitude on average.

So amplitude balance and source power-envelope timing are not the whole computation. The cross-source complex interaction is load-bearing.

## P3 — the conclusion is even stronger at the soma

Soma prediction error averaged over the five gain conditions:

```text
                         discovery          confirmation
MAE coherent              0.0000626          0.0000771
MAE incoherent            0.118587           0.088091
incoh - coh error        +0.118524          +0.088014
positive / negative        12 / 0             12 / 0
sign p                     0.000488           0.000488
registered P3                 PASS               PASS
```

Mean soma absolute selectivity:

```text
                         discovery          confirmation
actual                     0.19338            0.14366
coherent                   0.19337            0.14367
incoherent                 0.07626            0.05923
```

Thus the non-universal same-soma amplitude-balance dose response is not a mystery left for an AIS-like compartment. The missing factor already exists upstream: the two geometry-shaped complex transfer histories must interact coherently at the readout.

## What this changes

The current causal/mechanistic chain can now be stated more precisely:

```text
anatomy G
   -> graph-defined transfer dynamics
   -> h_A(x,t), h_B(x,t)
   -> amplitude balance determines whether both inputs have local leverage
   -> coherent cross-source term 2 Re[h_A conj(h_B shifted)]
   -> temporal-order contrast
```

Amplitude balance is an **opportunity variable**. It says whether both inputs are present strongly enough/comparably enough to interact.

The actual order computation is in their **relative transfer history**, expressed in this complex wave model by the cross-source term.

This also resolves the apparently mixed source-gain result: changing balance moves the selectivity landscape reliably across the body, but at a fixed soma the outcome can differ by body because equal amplitudes can still meet with different temporal/complex relationships.

## Important distinction from the failed AIS phase tests

The AIS/eventizer branch found no special benefit from restoring a phase-bearing `Re(psi)` signal to a generic HH boundary. That remains a valid null.

The present result is upstream and different. The historical soma power readout itself is generated from coherent field superposition. Removing the complex cross term *before* that power measurement destroys most of the temporal-order signal.

So for this model:

```text
phase/complex relation as a growth cue       not established
phase/complex relation as AIS advantage      failed
complex cross-source relation in field readout REQUIRED
```

There is no contradiction: the same mathematical resource can be load-bearing in one stage and useless in another.

## Wall sentence

> **Frozen geometry defines two source-specific complex transfer histories. Amplitude balance determines where both histories can matter, but temporal-order computation is carried largely by their coherent cross term. The soma/root is a privileged coincidence location because geometry brings two globally shaped wave histories together there, not because it measures total energy or because a downstream HH eventizer adds the missing computation.**

## Next clean question

Now that the computation has been reduced to the cross-source term, the next question is no longer "does phase matter?" in the vague sense.

It is:

> **Which geometry-defined modes contribute the cross term at the soma, and is the order-sensitive interference concentrated in a small modal subspace or distributed broadly across the arbor spectrum?**

That is the direct bridge back to the confirmed graph-mode result and should be tested before adding new biological machinery.
