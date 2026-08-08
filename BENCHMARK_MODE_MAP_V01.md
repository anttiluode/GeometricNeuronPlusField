# Benchmark mode map v0.1 — the winning free tuners target the same low body modes used by soma mixing

The matched-tuner benchmark produced an initially odd-looking detail: the winning full-grid spectral coordinates had indices around `892-907` in a 961-state operator.

`full_body_mode_map.py` resolves what those indices are.

The exact full 31x31 operator contains a weak bath surrounding the 70-cell strong arbor. For each selected full-grid eigenmode, the mode was restricted to the occupied body and matched against the original unweighted 70-cell body graph basis used by the graph-mode microscope.

## The mapping is essentially exact

Across the 12 held-out benchmark bodies:

```text
selected free coordinates                         96
unique selected full-grid modes                   71
mean |overlap| with best body graph mode          0.999999827
minimum |overlap|                                 0.999996823
mean full-mode power located on the body           0.999742842
minimum body power fraction                        0.996960375
```

So indices `~892-907` are not mysterious bath modes. They are almost perfectly embedded **body resonances**.

## Which body modes does the free optimizer choose?

Every one of the 96 selected coordinates maps to body modes **1-16**.

```text
body mode : number of selected coordinates

1   8
2  20
3  14
4  10
5   6
6   5
7   5
8   2
9   5
10  5
11  3
12  5
13  2
14  3
16  3
```

Registered earlier descriptive fractions:

```text
selected in modes 0-17       96 / 96 = 1.000
selected in modes 18-20       0 / 96 = 0.000
```

## This closes a loop with the soma mode-pair result

This is exactly what `MODE_PAIR_V01.md` would lead us to expect.

The old direct graph-mode microscope found modes 18-20 informative **when each mode was observed alone**.

But the later exact soma pair decomposition showed that the actual local square-law soma computation is different:

- ~97-98% of absolute pair mass is off-diagonal;
- the common mode 0 is the largest interaction hub;
- most total mode involvement lives in modes 0-17;
- the old 18-20 direct-readout band is not enriched in soma pair mixing.

The free benchmark optimizer now independently selects **only modes 1-16** as the best direct spectral knobs.

That is not an arbitrary optimizer artifact. It says the best place to edit the transfer function is the same low-mode subsystem that the soma's interference decomposition identified as consequential.

## Why mode 0 itself is not selected

Mode 0 is nearly order-blind on its own and acts as a common/reference component in the local square-law mixer. The useful tunable degrees of freedom are therefore naturally its **nonconstant partners**: change their pole/residue histories while leaving the large common reference available for mixing.

In schematic form:

```text
large common/reference component (mode 0)
                   x
      tunable low nonconstant modes 1-16
                   ↓
          soma cross-mode interference
                   ↓
             order contrast
```

This is a much cleaner version of the earlier “heterodyne” intuition. The data support a self-reference / reference-assisted mixer, not a literal external local oscillator.

## Why the free baseline wins

The graph bond coordinates and the free modal coordinates are now easy to distinguish mechanistically.

A local bond edit:

```text
one physical coupling
    -> correlated perturbation of many body modes/residues
```

A free spectral edit:

```text
choose one consequential low body mode
    -> move its pole or one residue directly
    -> leave the others independently controlled
```

For this task, the second coordinate system is simply better aligned with the soma computation. That explains the decisive held-out F8 win without invoking greater state dimension or a different simulator.

The cost is physical nonlocality: `SPECTRAL_COORDINATE_LOCALITY_CONFIRM_V01.md` shows that those direct spectral coordinates are not equivalent to one local bond actuator.

## Current reduction

The benchmark has therefore produced a mechanism result in addition to a leaderboard result:

> **The temporal-order task is controlled most efficiently by direct edits to low nonconstant body resonances that participate in soma cross-mode mixing. Local geometry reaches those same resonances only through coupled, constrained perturbations.**

That is the precise trade now worth carrying into a hardware compiler analysis.
