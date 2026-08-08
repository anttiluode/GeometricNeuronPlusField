# Free spectral coordinate ablation v0.1 — the benchmark win is not an artifact of one oracle coordinate type

This is a descriptive ablation on the already-opened matched-tuner benchmark bodies 288-299. It asks which direct spectral coordinate classes drive the strong `F8` advantage over eight local physical bond tuners.

All arms receive 8 trainable mathematical scalars, the same 40 normalized gradient steps, training lags `16,20,24`, and evaluation lags `14,18,22,26`.

Coordinate classes:

```text
F      direct modal stiffness / pole-frequency coordinate
A      source-A modal residue
B      source-B modal residue
C      soma/output modal residue
```

A restricted arm selects its eight strongest base gradients only from the allowed coordinate classes, then trains those fixed coordinates exactly as in the original benchmark.

## Result

```text
mean base test C                    -0.020657
mean graph G8 test C                +0.028492

free F-only test C                  +0.111105   Delta +0.131762
free ABC-only test C                +0.120278   Delta +0.140935
free FC test C                      +0.154691   Delta +0.175348
free FABC test C                    +0.168527   Delta +0.189184
```

Every restricted free arm beats the graph arm in **12/12 bodies**.

Mean test advantage over G8:

```text
F only          +0.08261   minimum body advantage +0.00830
ABC only        +0.09179   minimum body advantage +0.00102
F + C           +0.12620   minimum body advantage +0.03839
F + A+B+C       +0.14004   minimum body advantage +0.05393
```

The selected coordinate counts were:

```text
F-only          96 F
ABC-only        46 C, 25 B, 25 A
FC              66 F, 30 C
FABC            51 F, 26 C, 11 B, 8 A
```

## Interpretation

The free spectral win cannot be dismissed as being caused only by source-specific residue knobs.

Even **pole-frequency coordinates alone**, with no direct source or output residue editing, beat eight local bonds in every body and more than double the graph arm's mean held-out improvement.

Likewise, residues alone beat the graph in every body.

The strongest compact combination is `F+C`: directly tune consequential body resonances and their soma participation. Adding source-specific residues gives a further but smaller gain.

This fits the mechanism exposed by `BENCHMARK_MODE_MAP_V01.md`:

```text
low body resonances 1-16
        +
soma/output participation
        ↓
low-mode cross-interference at the square-law readout.
```

## What remains physical

This strengthens rather than weakens the hardware-cost distinction.

The `F-only` result shows that the abstract advantage survives after removing source-specific input knobs, but `SPECTRAL_COORDINATE_LOCALITY_CONFIRM_V01.md` shows that a pure pole coordinate is itself a strongly nonlocal rank-one operator deformation in the nearest-neighbour spatial basis.

So the current split is clean:

> **Spectral coordinates are the better optimization coordinates for this task. Local bonds are the cheaper/natural physical coordinates.**

The next useful question is the compiler question: how many feasible local bond directions are needed to approximate one winning spectral direction over a family of temporal conditions? `spectral_to_local_compiler.py` is the first task-space audit of that trade.
