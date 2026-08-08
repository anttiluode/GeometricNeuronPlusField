# Graph-mode discovery v0.1

Canonical GitHub Actions run: `31240168834`

12 frozen FunctionalArbor v0.9 bootstrap bodies, seeds 0–11. No growth, credit, learning, or wave-physics changes. Every graph-mode readout watched the same `psi(x,t)` as the point soma readout.

This is the **discovery set**. Anything selected after looking at these 12 bodies is treated as exploratory until it survives held-out seeds.

## 1. The constant graph mode is almost perfectly order-blind

At lag 20:

```text
point readout mean |C|             0.22077
mode 0 coherent mean |C|           0.000098
mode 0 incoherent mean |C|         0.02262
```

The coherent mode-0 result is not a weak-signal accident in the repeated-drive run. Mode 0 carries, on average, about **52.1% of total field energy** (range about 36.5%–64.5%) while carrying essentially no A/B order contrast.

That is the cleanest result of the first run:

> **The energy-dominant spatial common mode is almost blind to temporal order.**

This is exactly the graph-basis version of the earlier whole-body coherent tap result. The information is not simply where most of the field energy is.

## 2. A higher graph-mode band exposes temporal order through coherence

The strongest post-hoc cluster is modes 18–20 (roughly `lambda ~ 0.55–0.75` across these bodies).

Across seeds 0–11:

```text
modes 18–20 mean coherent |C|                 0.12533
modes 18–20 mean coherent-minus-incoherent    0.09466
seeds with positive band coherence gain       12 / 12
```

The three modes together carry only about **0.34% of total repeated-drive field energy on average**.

So the first graph microscope shows an intriguing separation:

```text
common mode:
    huge energy
    almost no temporal-order information

higher spatial band:
    tiny energy
    appreciable temporal-order selectivity
    selectivity strengthened by coherent projection
```

Individual mode 18 and mode 20 each had positive coherent-over-incoherent gain in 11/12 discovery bodies. But those are post-hoc per-mode observations and there are 24 modes, so they are not treated as confirmed single-mode findings.

The **band** is the object to confirm, because it is the structure that stood out after looking at the spectrum.

## 3. The first live-field / settled-readout criterion is a null

The field remained live in every body:

```text
median normalized field motion across seeds   ~0.104
field-live bodies                              12 / 12
```

But **zero** modes passed the prewritten combined criterion:

```text
cycle CV <= 0.05
modal energy fraction >= 0.005
|C| >= 0.05
field motion >= 0.01
```

That result stays a null. We do not loosen the threshold after seeing it just to manufacture a success.

Mode 19 came closest in one sense: strong selectivity and median cycle CV ~0.097, but its energy fraction was only ~0.00045. That is interesting for designing a different question later, not a pass under the registered v0.1 definition.

## 4. Interpretation earned so far

The graph basis is already useful as a microscope.

It says the body's electrical activity is not organized as:

```text
more field energy = more task information
```

Instead, task information can live in weak spatial deviations from an energetically dominant common field. Coherent geometric projection can expose some of those deviations far more strongly than a phase-destroyed energy measurement over the same mode support.

That is a concrete version of the phrase:

> **geometry as the transformation between dynamics and observables**

It is not yet evidence that biological somata implement graph eigenmodes.

## 5. Held-out confirmation registered before looking at seeds 12–23

`confirm_graph_band.py` freezes the choices made from discovery and tests new bodies only.

Registered claims:

1. **Common-mode blindness:** coherent mode 0 mean `|C| < 0.005` at lag 20.
2. **Band coherence advantage:** for modes 18–20, the per-seed mean `(coherent |C| - incoherent |C|)` is positive in at least 9/12 held-out bodies and the exact two-sided sign test is `< 0.05`.
3. **Band is informative:** held-out mean coherent `|C|` across modes 18–20 is `> 0.05`.
4. No registered claim says the band beats the point soma detector.
5. The live-field/settled-readout null is not retested by changing thresholds in this confirmation.

If the band fails on held-out seeds, the spectral bump was a discovery-set feature and we keep the failure.
