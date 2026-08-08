# GeometricNeuronPlusField

A small experimental branch of the Geometric Neuron / Functional Arbor line.

The question here is narrower than "does geometry compute?" The Functional Arbor soma tap test already showed that the same frozen wave-carrying body can report a different temporal preference when the spatial readout is changed. The readout is therefore not a neutral voltmeter; it is part of the computation.

This repo takes the next step:

> **Stop choosing soma shapes by hand. Use the arbor's own graph modes as the readout basis, then ask whether a non-settled electrical field can support a settled, informative observable.**

The source experiment is in [`FunctionalArbors/SomaTapTestsClaude`](https://github.com/anttiluode/FunctionalArbors/tree/main/SomaTapTestsClaude).

## Working picture

There are now two geometries around one moving field:

```text
input history
    |
    v
arbor / propagation geometry G_D
    |
    v
moving complex field psi(x,t)
    |
    v
readout geometry G_S
    |
    v
scalar consequence
```

The old implementation fixed `G_S` to one point:

```text
y(t) = |psi(soma,t)|^2
```

The tap test showed that this choice is load-bearing. Near-soma apertures mostly rescale the same signal, distal apertures can reverse preference, and a uniform whole-body coherent aperture becomes order-blind. A smooth coherent Gaussian aperture was the one extended readout that showed a robust coherent-over-incoherent advantage.

The new question is whether the arbor itself supplies a more natural family of readouts.

## Experiment 1 — graph-mode microscope

For a frozen body, build the unweighted body graph and its combinatorial Laplacian

```text
L = D - A
L phi_n = lambda_n phi_n
```

The `phi_n` are not hand-drawn soma masks. They are the body's own spatial modes.

For every mode we measure a coherent projection

```text
a_n(t) = sum_x phi_n(x) psi(x,t)
y_n(t) = |a_n(t)|^2
```

and a phase-destroyed matched control

```text
y_n,incoh(t) = sum_x phi_n(x)^2 |psi(x,t)|^2
```

Because the eigenvectors are orthonormal, the modal powers also give a clean decomposition of where field energy sits in the geometry.

The probe sweeps A/B pulse lag exactly as the soma tap test did and asks:

1. Which graph modes carry temporal-order selectivity?
2. Is selectivity concentrated in a small spectral band or spread across the body?
3. Does coherent projection buy information beyond the matched incoherent modal-energy control?
4. Does the mode that best reports path asymmetry vary with morphology?
5. Is the spatially constant mode (`lambda ~= 0`) blind to order, as the whole-body tap result suggests?

## Experiment 2 — live field / settled readout

A field does not need to freeze for an observation of it to become usable.

We repeatedly drive a frozen body with the same A-then-B pulse pair. The instantaneous field is required to remain active:

```text
D_field(t) = ||psi(t) - psi(t-1)|| / (||psi(t)|| + eps)
```

For each graph mode, we integrate its coherent power over one drive cycle and ask whether that cycle-level observable converges even while `D_field` stays nonzero.

A candidate "settled readout of a live field" must pass three gates:

- the field is still moving;
- the modal observable is stable from cycle to cycle;
- the mode carries non-trivial energy **and** non-trivial temporal selectivity.

This last gate matters. A perfectly stable zero readout is not computation.

## Registered failure conditions

This branch should be easy to kill.

- If graph modes merely reproduce the point detector with no spectral structure, the graph-basis idea buys nothing.
- If coherent and incoherent mode readouts are indistinguishable, the interference claim does not survive this basis.
- If only near-zero-energy or task-blind modes look "settled," the live-field/settled-readout idea fails.
- If a result appears only in the seed average and is not organism-stable, report the heterogeneity rather than the mean story.
- If the constant mode carries the same order information as nonzero modes, the simple "common mode is order-blind" interpretation from the tap test was too strong.

## Run

Clone this repo beside `FunctionalArbors`:

```text
parent/
  FunctionalArbors/
  GeometricNeuronPlusField/
```

Then:

```bash
pip install -r requirements.txt
python graph_mode_probe.py --functional-arbors ../FunctionalArbors --seeds 12 --modes 24
```

A quick internal check that does not require FunctionalArbors:

```bash
python graph_mode_probe.py --selftest
```

Outputs are written under `runs/graph_modes/` by default:

- `graph_mode_results.json` — full seed-level receipts
- `graph_mode_summary.png` — selectivity and live/settled summary by mode index
- `mode_maps_seed0.png` — low graph modes painted back onto the first frozen arbor

## What this is not

This is an experimental computational model, not a claim that biological somata literally diagonalize graph Laplacians. The graph basis is a microscope: a geometry-derived coordinate system for asking where the information already present in the simulated field lives.

If something survives, the later biological/computational question is more interesting: whether a physical soma/dendritic interface could approximate a useful low-dimensional projection of those modes, and whether body geometry and readout geometry can co-adapt.
