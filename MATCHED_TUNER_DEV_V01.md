# Matched-tuner benchmark development v0.1 — direct spectral coordinates win the first comparison

This is a **development result on reused seeds 240-243**, not a held-out confirmation.

`BENCHMARK_SCOPE_V01.md` first established that a fully free pole/residue bank contains the frozen reciprocal graph transfer, so the meaningful contest is parameterization efficiency per added trainable scalar, not raw expressivity.

`matched_tuner_benchmark.py` therefore started both arms from the same exact mature transfer and gave each **8 additional trainable real coordinates**, 40 normalized gradient steps, and the same train/test lag sets.

Train lags:

```text
16, 20, 24
```

Development test lags:

```text
14, 18, 22, 26
```

Graph arm:
- chose the 8 most favorable legal frontier conductance additions from ~44-48 candidates;
- trained their normalized conductances in `[0,1]` with the relinearized exact adjoint.

Free-modal arm:
- diagonalized the same exact 961-state base operator;
- searched all direct modal pole-frequency and source-A/source-B residue coordinates;
- selected the 8 largest initial absolute gradients from 2883 candidate scalar coordinates;
- retrained those same coordinates for 40 steps.

## Development result

```text
mean base development-test C                 0.14850
mean graph development-test C                0.20750
mean free-modal development-test C           0.29688

mean graph Delta test                       +0.05899
mean free-modal Delta test                  +0.14837

mean graph - free-modal test                -0.08938
graph beats free-modal bodies                 0 / 4
```

Per body graph-minus-free test:

```text
seed 240   -0.01898
seed 241   -0.09081
seed 242   -0.09271
seed 243   -0.15502
```

Both arms generalized from the three training lags to the four interleaved lags, but the direct spectral coordinates improved much more.

Base graph/modal identity error was at most `1.57e-11`, so this is not a different-simulator artifact.

## Do not rescue the graph by weakening the baseline

This first result is already an important correction to the hoped-for story. There is no evidence here that local geometric tuning beats direct spectral tuning per added scalar. On this simple temporal-order objective, the opposite is true in all four development bodies.

However, the comparison still has one asymmetry worth auditing before freezing a held-out conclusion:

- the graph compiler searched only legal **frontier additions** (~48 candidate locations), while the free compiler searched 2883 direct spectral coordinates;
- a local conductance change also changes soma/output residues through eigenvector motion, but the free baseline did not yet include an explicit output-residue coordinate.

The second point strengthens the free baseline; the first strengthens the graph arm if corrected.

The next development version should therefore be deliberately harder to explain away:

1. allow the graph compiler to search **every physically feasible local bond tuner** in the exact weighted medium: weak bonds may strengthen and strong arbor bonds may weaken;
2. allow the free baseline to search pole frequency, source-A residue, source-B residue, **and soma/output residue** coordinates;
3. retain 8 trainable scalars, the same objective, the same iteration budget, and the same reused development seeds.

If the free spectral parameterization still wins after that symmetry audit, there is little reason to spend fresh held-out seeds testing a graph-performance advantage on this task. The hardware case would then rest on locality, measurement/training implementation, fabrication, robustness, or energy — not superior abstract task performance per scalar.
