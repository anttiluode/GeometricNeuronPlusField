# Held-out graph-band confirmation v0.1

Canonical GitHub Actions run: `31240309039`

The discovery set was seeds 0–11. Before looking at seeds 12–23, `DISCOVERY_V01.md` froze three claims and `confirm_graph_band.py` froze the band to modes 18–20.

## Result

All three registered tests passed on the 12 held-out bodies.

```text
held-out seeds completed                 12 / 12

1. common-mode blindness
   coherent mode-0 mean |C|              0.00008415
   registered limit                      < 0.005
   PASS

2. modes 18–20 coherence advantage
   mean coherent-minus-incoherent |C|   +0.071796
   positive bodies                       10 / 12
   exact two-sided sign test             p = 0.038574
   registered rule                       >=9/12 and p<.05
   PASS

3. modes 18–20 are informative
   mean coherent band |C|                0.115138
   registered limit                      > 0.05
   PASS

context only:
   point-soma mean |C|                   0.235888
```

There was explicitly **no** registered claim that the graph-mode band would beat the point soma detector. It does not.

## What is now established in this toy

The first spectral bump was not just a seed-0-to-11 accident.

On new morphologies, the spatial common mode is again almost perfectly blind to A/B temporal order, while a higher graph-Laplacian band again gains temporal selectivity specifically from coherent projection relative to its matched incoherent energy readout.

So we can now say more cleanly:

> **Temporal-order information in the frozen arbor is spectrally structured. It is nearly absent from the spatial common mode and is exposed in a higher geometry-defined band by preserving cross-location phase/sign relationships.**

That is stronger than the hand-designed soma tap result because the readout family was not point/cross/ring/Gaussian geometry chosen by us. It came from the body's own graph.

It is still not a biological claim that a soma computes graph eigenvectors. The graph basis is a microscope.

## What is not rescued

The v0.1 live-field / settled-readout criterion remains a null. We do not relax its thresholds after the fact.

The point soma readout also remains a strong detector, so this does not explain away the FunctionalArbor v0.9 credit-assignment null.

## The next mechanistic question

Why modes 18–20?

The graph basis suggests an actual physics reduction. If `phi_n` is a body mode, its modal amplitude should approximately obey a forced damped-oscillator equation:

```text
q_n'' + damping q_n' + (restoring + K * lambda_n) q_n
    ~= phi_n(A) s_A(t) + phi_n(B) s_B(t)
```

That gives two geometry-derived ingredients for temporal-order sensitivity:

1. `lambda_n` sets the modal time scale / phase response;
2. `phi_n(A)` and `phi_n(B)` set how strongly the two terminals drive that mode.

The next experiment should therefore stop treating the spectral bump as a mysterious index range and ask whether measured modal contrast is predicted by **eigenvalue + source coupling**, and whether the modes align with the A-vs-B unique-path geometry.

If that works, the Geometric Neuron picture becomes much sharper:

```text
anatomy
  -> graph spectrum
  -> source coupling into modes
  -> travelling/oscillating modal state
  -> geometric projection
  -> consequence
```

The geometry would not merely carry the field. It would define the field's computational coordinates.
