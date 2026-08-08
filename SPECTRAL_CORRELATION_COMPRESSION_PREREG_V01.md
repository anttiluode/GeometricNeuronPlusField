# Spectral correlation compression — held-out preregistration v0.1

## Question

The exact time-domain internal bond gradient needs the local anti-diagonal correlation

```text
g_e ~ sum_r forward_e[T-1-r] * conj(retro_e[r]).
```

The development probe showed that this correlation has an exact DFT representation and that its gradient map is strongly concentrated in a small number of frequency bins. The important hardware-facing observation was that bins ranked using only boundary-accessible spectra (external source spectrum times soma-return source spectrum) were nearly as effective as an oracle ranking that inspects internal gradient contributions.

This confirmation freezes that implementation and tests fresh bodies.

## Frozen implementation

Script: `spectral_correlation_compression_probe.py`

No code or thresholds in this document may be changed after the fresh run begins.

Task:

```text
lag        20 frames
steps      210
fresh seeds 400-411
```

The same mature reciprocal FunctionalArbor wave and the same contrast-energy objective used by the reciprocal-adjoint line are retained.

For every bond, the exact physical replay correlation can be written

```text
sum_n reverse(f_e)[n] conj(b_e[n])
 = (1/T) sum_k exp(+i 2 pi k/T) F_e[-k] conj(B_e[k]).
```

The all-bin reconstruction is the algebraic positive control.

Two bin rankings are computed exactly as in development:

```text
oracle
  total absolute internal gradient contribution of each bin

boundary
  external/source spectral magnitude x soma-return-source spectral magnitude
```

Only `boundary` is relevant to the proposed port-selected compression mechanism. `oracle` is retained as a lower-bound/context control.

## Registered criteria

### C0 — all-bin identity

The DFT representation must reproduce the exact adjoint bond map:

```text
mean all-bin correlation > 0.999999
maximum body relative L2 error < 1e-10
```

If this fails, stop: the implementation is wrong.

### C1 — eight boundary-selected bins preserve map direction

At `K=8` boundary-selected bins:

```text
mean bond-map correlation > 0.985
at least 10 / 12 bodies have correlation > 0.970
```

Both conditions are required.

### C2 — sixteen boundary-selected bins preserve the map quantitatively

At `K=16`:

```text
mean bond-map correlation > 0.995
mean relative L2 error < 0.10
```

Both conditions are required.

### C3 — gradient spectral mass is sparse

For each body, sort bins by the oracle absolute gradient-contribution mass and count the bins required to reach 95% cumulative mass.

Registered group criterion:

```text
mean K95 <= 20 bins
```

This is descriptive of the task/medium, not a hardware-available selection rule by itself.

### C4 — boundary ranking is close to the internal oracle at K=8

Let the body-level K8 metrics be averaged across the 12 bodies.

Registered:

```text
mean_corr(boundary K8) - mean_corr(oracle K8) > -0.020
mean_relL2(boundary K8) - mean_relL2(oracle K8) < +0.050
```

Both conditions are required.

## Interpretation rules

If C0-C4 all pass, the earned statement is:

> **For this registered temporal task, the exact local forward×adjoint history correlation is spectrally sparse enough that a small bank of port-selected phasor/lock-in accumulators can approximate the full bond-gradient map without retaining all 210 local time samples.**

Do **not** claim from this experiment alone that:

- a two-pass intensity-only TRIM circuit has been constructed;
- the required local memory is O(1) independent of task bandwidth/duration;
- eight frequencies are universal across tasks or bodies;
- this is a biological learning mechanism;
- real-valued hardware needs exactly the same one-sided complex-bin convention as the analytic toy field.

The hardware claim remains one step weaker: spectral compression changes the candidate local storage/readout burden from a full time trace toward a modest number of coherent spectral accumulators. A physical intensity/interference realization and learning-through-the-compressed-gradient are separate tests.
