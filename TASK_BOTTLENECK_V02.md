# Soma task bottleneck v0.2 — the local task projection is concentrated at the root

Canonical GitHub Actions run: `31241537322`

The local-observability experiment asked whether a soma patch can reconstruct the confirmed global graph modes individually.  That is a demanding state-estimation problem.

This experiment asks the simpler question a neuron actually needs:

> **Is the soma neighborhood a good place to read the task-relevant scalar distinction, even if it cannot reconstruct the global modal state?**

The same 24 frozen 70-cell bodies were driven with A->B and B->A at lag 20.  No growth, credit or learning occurred.  Every occupied cell saw the identical field evolution.

For each body we measured peak-power temporal-order contrast at every cell, then coherent uniform graph balls around the soma.  Soma balls were compared both with equal-count scattered controls and with the more important **same-radius contiguous balls centered at every other occupied cell**.

## Result 1 — the soma point is already a strong task projection

```text
24 / 24 bodies: soma lies in the top quartile of body cells
10 / 24 bodies: soma lies in the top decile

mean soma |C|                    0.22833
mean |C| over all body cells     0.05413
median soma percentile           0.87857
```

Across the 24 bodies the soma percentile never falls below about `0.779`.

So the soma is not where field energy in general is maximal; it is repeatedly near the high end of **temporal-order selectivity**.

## Result 2 — the neighborhood result survives locality-matched controls

```text
R   mean soma cells   soma-ball |C|   mean other local-ball |C|   median local percentile
1       3.21              0.2363              0.0578                    0.879
2       5.58              0.2338              0.0594                    0.879
3       7.75              0.2246              0.0600                    0.864
4      10.08              0.1859              0.0607                    0.850
5      12.62              0.1539              0.0599                    0.821
6      15.08              0.1174              0.0578                    0.757
```

At **every radius 1–6, all 24/24 soma-centered balls are above the median same-radius local aperture** on their own body.

Under a simple two-sided sign test against 0.5, `24/24` has `p ~= 1.19e-7`.  The exact p-value is not the most important part because the soma is a designated root in this model, not an exchangeable random cell.  The organism was built around convergence there.  The important result is the geometry of that convergence: the task-relevant scalar becomes strongly concentrated at the root even though the confirmed modal coordinates remain global.

## Put this beside local observability

The pair of results is more informative than either alone.

`LOCAL_OBSERVABILITY_V02`:

```text
small local patch
   -> poor absolute reconstruction of global modes
   -> soma somewhat better-conditioned than other local patches
```

`TASK_BOTTLENECK_V02`:

```text
same small local patch
   -> excellent task scalar at soma
   -> far better than other local patches
```

So the soma does **not** behave like a local Fourier microscope that recovers the global state.

It behaves much more like a **task bottleneck**:

```text
many global modal coordinates
            |
            v
     convergence geometry
            |
            v
      low-dimensional mixture
            |
            +---- contains the task distinction strongly
            |
            v
           SOMA
            |
            v
     [missing active boundary]
```

That missing boundary is exactly where the AIS comparison becomes interesting.

## Why this is relevant to an AIS-like stage

A downstream active compartment does not need access to every `q_n(t)`.  It only needs a local signal in which the globally computed distinction has become available.

The current toy now has precisely that decomposition:

1. **body geometry** defines a global analog resonator bank;
2. **moving field** carries temporal computation through those modes;
3. **convergence at the soma** compresses the global state into a locally strong task variable;
4. the model presently stops at a passive power readout;
5. a biological neuron instead places the **active axon initial segment** immediately downstream of this convergence region.

The next experimental branch should therefore not ask an AIS-like model to recover the global modes.  It should ask whether active channel state can convert the existing soma mixture into a temporally precise event code and whether AIS parameters compensate differences between body geometries.

## Important construction caveat

This result does **not** show that an unbiased morphogenetic process spontaneously chose the soma as the optimal output site.  `FunctionalArbor` starts with a designated soma/root and grows source connectivity around it.

What is earned is narrower:

> **Given a branching body organized around a convergence root, the distributed wave computation produces a task-selective bottleneck at that root.**

That is enough to motivate the next architecture without turning a construction choice into a biological discovery.

## Wall sentence

> **The global modal field is not locally reconstructable, but its task-relevant distinction is strongly compressed at the soma: the root is a task bottleneck, not a miniature copy of the field.  This gives an AIS-like active boundary something concrete to operate on.**
