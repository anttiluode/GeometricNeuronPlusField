# Local task-band observability v0.2 — locality-matched controls

Canonical GitHub Actions run: `31241432164`

v0.1 compared a soma-centered graph ball with **equal-size scattered random cell sets**.  That was a useful first screen but an unfair geometry control: scattered taps naturally span global modes better than any compact local patch.

v0.2 fixes that.  For every frozen 70-cell body and every radius `R`, the soma ball is now compared with **same-radius graph balls centered on every other occupied cell**.  The scattered equal-count control remains in the receipt only as context.

The confirmed task band is still modes 18–20.  No wave simulation, learning, fitted readout, or active AIS model is used.

## Result

The conclusion changes in an interesting way.

Compact local apertures everywhere are poor at independently reconstructing three global modes.  The near-singularity seen in v0.1 is therefore mostly a consequence of **locality**, not a pathology of the soma.

Within that physically relevant class of local apertures, the soma is actually a comparatively favorable location.

```text
R   mean soma cells   capture percentile*   s_min percentile*   isotropy percentile*
1       3.21               0.643                 0.786                0.750
2       5.58               0.557                 0.714                0.693
3       7.75               0.557                 0.679                0.686
4      10.08               0.579                 0.729                0.764
5      12.62               0.521                 0.693                0.707
6      15.08               0.493                 0.629                0.657

* median percentile of the soma-centered aperture among same-radius local balls
  on the same body.
```

For the smallest-singular-value criterion — the hardest-to-see direction in the 3-D task band — the soma ball is above the median local aperture in:

```text
R1   22 / 24 bodies   exact sign p = 0.000036
R2   21 / 24          p = 0.000277
R3   21 / 24          p = 0.000277
R4   22 / 24          p = 0.000036
R5   22 / 24          p = 0.000036
R6   20 / 24          p = 0.001544
```

The isotropy metric gives the same qualitative result.  Energy capture itself is much less special: by radius 5–6 the soma is around the middle of the local-ball distribution.

So the interesting property is **not that more task-band energy is concentrated at the soma**.  It is that, among compact local windows, the task-band components tend to be somewhat less collapsed into the same direction there.

## The corrected picture

v0.1 tempted the sentence:

> the soma is a lossy bottleneck.

That remains true in the absolute sense, but it was incomplete.  Every compact local aperture is a lossy bottleneck for global modes.

The better sentence is:

> **The soma is a lossy local bottleneck that is nevertheless unusually well-conditioned relative to other equally local bottlenecks.**

That is much more interesting for the AIS comparison.

```text
global modal field
       |
       |  locality necessarily compresses
       v
soma / axon-hillock neighborhood
       |  comparatively favorable local mixture
       v
active AIS boundary ?
       |
       v
spike timing
```

The result does not establish that biology deliberately positions the soma/AIS to observe graph eigenmodes.  These graph modes remain a microscope for this toy.

But it kills one alternative: the soma is not merely an arbitrary local point that happens to have been selected by the code.  In these grown bodies, its neighborhood is repeatedly above-average among local neighborhoods at preserving independent access to the confirmed task-mode subspace.

## Why the scattered control looked disastrous

At radius 3, for example, a soma ball has mean `s_min ~= 0.0094`, whereas an equal-count **scattered** random aperture has `s_min ~= 0.131` — more than an order of magnitude larger.

That does not mean the soma is uniquely blind.  Distributed samples get to reach across the graph and naturally separate global eigenvectors.  Same-radius local balls have mean `s_min ~= 0.0071`, putting the soma above the local mean.

This distinction matters for any biological interpretation: a real soma or AIS is compact.  It cannot implement a scattered graph-wide measurement for free.

## Next question

Full modal observability is stronger than neuronal function requires.  A neuron need not reconstruct modes 18, 19 and 20 individually.  It needs a useful low-dimensional task decision.

`task_bottleneck_probe.py` therefore measures the actual A-then-B versus B-then-A wave contrast at every occupied cell and at coherent local balls, and asks whether the soma is also privileged for **the task scalar itself**, not merely for modal conditioning.

If both survive, the architecture becomes unusually specific:

1. global geometry creates task-bearing modes;
2. locality compresses those modes;
3. the soma region is a comparatively good local compression point;
4. an AIS-like active boundary can then be tested as the history-dependent frequency-selective event encoder.

## Wall sentence

> **Global task modes are intrinsically difficult to separate from any compact patch, but the soma neighborhood is consistently better-conditioned than other equally local neighborhoods while not containing unusually high modal energy.  The soma looks less like an energy focus and more like a favorable local mixing point.**
