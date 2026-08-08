# Local task-band observability v0.1

Canonical GitHub Actions run: `31241224985`

The graph-mode microscope uses global eigenvectors.  This probe asked the stricter physical question: **can a small graph-distance ball around the soma independently observe the confirmed task band (modes 18-20)?**

For each of 24 frozen 70-cell FunctionalArbor bodies, the three confirmed mode vectors were restricted to cells within graph radius `R` of the soma.  The singular values of that local mode matrix measure whether the three global coordinates remain independently visible.  Every soma aperture was compared with 250 equal-size random cell sets on the same body.

No wave simulation, learning, or fitted readout was used.

## Result: the soma is not a passive modal focus

The result is a strong null for the idea that the soma neighborhood is a particularly good place to reconstruct the three global task modes.

```text
R   mean cells   task-band energy capture / random   median s_min percentile   median isotropy percentile
0      1.0                 0.653x                          0.500                     0.500
1      3.2                 0.836x                          0.110                     0.118
2      5.6                 0.924x                          0.032                     0.034
3      7.8                 0.919x                          0.006                     0.008
4     10.1                 0.942x                          0.002                     0.000
5     12.6                 0.969x                          0.000                     0.000
6     15.1                 0.975x                          0.000                     0.000
```

By radius 3-6, a soma-centered aperture contains roughly the amount of task-band modal energy expected from its cell count, **but the three restricted modes are extraordinarily poorly conditioned compared with random apertures of the same size**.

At radius 3, for example:

```text
mean aperture size                      7.75 cells
mean task-band energy capture           0.1024
mean matched-random capture             0.1114
capture enrichment                      0.919x
mean smallest singular value            0.00937
mean matched-random smallest SV         0.13117
median soma-aperture s_min percentile   0.006
full algebraic rank                     24 / 24 bodies
```

So the matrix is technically rank 3, but it is nearly singular in the soma neighborhood: the three global task modes look much more alike there than they do at random collections of cells.

## What this means

This does **not** say the soma cannot report the task.  The earlier tap experiment already showed that the ordinary point soma is a strong temporal-order detector.

It says something more specific:

> **The soma is a compression point, not a faithful local reconstruction of the global mode bank.**

That is arguably more interesting for the AIS bridge.

A biological output compartment does not need to recover `q18`, `q19`, and `q20` separately.  It needs a consequential combination of whatever the global field presents locally.  The present result says a soma-local patch has already collapsed several distinct global modal coordinates into a nearly low-dimensional mixture.

In other words:

```text
global modal state
      |
      v
soma-local field
  (lossy compression)
      |
      v
??? active readout / AIS
      |
      v
spike decision
```

The missing question is therefore no longer "can the soma see all of the modes?"  It cannot, robustly.

The useful question is:

> **Is the compressed soma mixture nevertheless unusually good for the task-relevant scalar decision?**

That is now tested separately by `task_bottleneck_probe.py`, which maps the actual A->B versus B->A point-power contrast over every occupied cell and compares soma-centered coherent apertures with matched random apertures.

## Why this matters for the AIS idea

This result argues against a passive "soma as local Fourier analyzer" picture.

If the AIS bridge survives, it should look instead like:

1. dendritic geometry creates a high-dimensional global analog state;
2. convergence near the soma compresses it;
3. an **active, history-dependent** AIS-like boundary acts on that compressed signal;
4. spike initiation converts the analog mixture into a temporally precise directional output.

That architecture fits the biological fact that the AIS is an active channel-rich compartment rather than merely a passive spatial aperture.

## Wall sentence

> **The confirmed task modes are global and independently observable from many distributed samples, but they become almost collinear in small soma-centered neighborhoods.  The soma is therefore not a passive modal microscope; it is a lossy bottleneck.  Any AIS-like computation must exploit the task-relevant mixture available at that bottleneck rather than reconstruct the global field.**
