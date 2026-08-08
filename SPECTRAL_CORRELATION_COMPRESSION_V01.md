# Spectral correlation compression v0.1 — the time-domain credit history compresses to a small phasor bank

The reciprocal-adjoint result solved **transport** of the task derivative through the same reciprocal wave medium, but `TIME_DOMAIN_IN_SITU_WALL.md` exposed the remaining hardware cost: the exact internal bond gradient is a time-aligned local correlation between the forward bond history and the causally replayed adjoint history.

The naive memory-free fixes failed. Ordinary simultaneous forward×retro correlation was poor, and simply reversing the external input did not reconstruct the damped internal forward history.

This experiment asked whether the required local temporal history can instead be represented spectrally.

## Exact identity

For one local bond let

```text
f[n] = forward edge difference
b[n] = causal retro/adjoint-replay edge difference
```

for `n=0..T-1`. The physical replay needs

```text
sum_n f[T-1-n] conj(b[n]).
```

With DFTs `F`, `B`,

```text
FFT(reverse(f))[k]
  = exp(+i 2 pi k/T) F[-k]
```

and therefore

```text
sum_n f[T-1-n] conj(b[n])
 = (1/T) sum_k exp(+i 2 pi k/T) F[-k] conj(B[k]).
```

So the exact long time-domain correlation can be accumulated as independent frequency-bin products. The question becomes whether many of those bins can be discarded.

## Development

Reused seeds 240-243 suggested strong concentration:

```text
mean bins for 50% absolute gradient mass      3.00
mean bins for 80%                             5.75
mean bins for 95%                            13.00
```

A ranking based only on boundary-accessible spectra — external source spectrum × soma-return source spectrum — approximated the internal oracle ranking unusually well:

```text
boundary-selected bins     corr      relative L2
K=4                        .9818       .2015
K=8                        .9940       .1016
K=16                       .9981       .0576
K=24                       .9997       .0227
```

This motivated the frozen held-out test in `SPECTRAL_CORRELATION_COMPRESSION_PREREG_V01.md`.

## Held-out confirmation — fresh seeds 400-411

### C0 — all-bin identity

Registered:

```text
mean correlation > .999999
max relative L2 < 1e-10
```

Observed:

```text
mean correlation             1.000000000000
max relative L2              2.73e-15
```

**C0 PASS.**

The DFT representation is the exact same bond gradient to floating-point precision.

### C1 — eight port-selected bins preserve map direction

Registered:

```text
mean K8 boundary corr > .985
>= 10/12 bodies corr > .970
```

Observed:

```text
mean K8 boundary corr          .9919
bodies > .970                  11 / 12
```

The one weaker body was seed 404 at about `.958`; this is retained rather than hidden.

**C1 PASS.**

### C2 — sixteen bins preserve the map quantitatively

Registered:

```text
mean K16 corr > .995
mean relative L2 < .10
```

Observed:

```text
mean K16 boundary corr         .9973
mean relative L2               .0672
```

**C2 PASS.**

### C3 — spectral gradient mass is sparse

Registered mean `K95 <= 20`.

Observed:

```text
mean K50                       2.50 bins
mean K80                       5.92 bins
mean K95                      13.00 bins
```

Individual K95 values were between 10 and 15 bins.

**C3 PASS.**

### C4 — the boundary ranking stays close to the internal oracle

At K=8:

```text
boundary mean corr             .9919
oracle mean corr               .9954
difference                    -.0035

boundary mean relative L2      .1120
oracle mean relative L2        .0909
difference                    +.0211
```

Registered tolerances were `corr difference > -.020` and `relative-L2 difference < +.050`.

**C4 PASS.**

## Confirmation verdict

```text
C0 PASS   spectral identity is exact
C1 PASS   8 boundary-selected bins preserve gradient direction
C2 PASS   16 bins give <10% mean relative-L2 error
C3 PASS   95% of absolute spectral gradient mass uses ~13/210 bins
C4 PASS   boundary-only bin selection remains close to internal oracle
```

**5 / 5 registered criteria pass.**

## What changed

The previous hardware wall was stated as:

```text
same medium gives adjoint transport in O(1) wave passes
BUT
local gradient readout may require a 210-sample forward history at every tuner.
```

That is now too pessimistic for this task.

A more accurate candidate cost is:

```text
same reciprocal medium
    -> physical adjoint transport

local tuner
    -> small coherent spectral accumulator bank
       rather than full raw temporal trace
```

At the registered operating point, 8 port-selected bins recover a bond-gradient map with mean correlation `.992`; 16 recover `.997` with about `6.7%` relative L2 error.

This is a **~210 samples -> O(10) complex phasor-state compression** for this task, not an O(1) theorem.

## Why the boundary ranking matters

The sparse bins were not selected by inspecting each internal bond. The successful ranking uses quantities available at the ports:

```text
external/source spectrum
        x
soma-return / objective-derivative source spectrum.
```

So the global controller can choose a common small frequency set, while every tunable element accumulates only those same phasors locally.

That is much closer to a fabricable compiler architecture than an oracle-selected per-bond frequency set.

## Relation to TRIM / in-situ backprop

This result does **not** claim that an intensity-only two-pass TRIM circuit for the broadband task has already been built here.

What it does is bridge the time-domain problem toward the steady-state photonic picture:

```text
long dynamic correlation
    -> Parseval / DFT decomposition
    -> small set of coherent phasor correlations
    -> candidate lock-in / interferometric local measurement
```

A physical implementation can now target roughly O(10) coherent channels rather than O(T) raw local storage for this registered task.

For real physical fields, conjugate-frequency pairing and the exact intensity/phase-stepping implementation must be derived explicitly; the analytic complex field used here should not be mapped one-for-one onto one-sided hardware frequencies without that step.

## Next test

Map accuracy is not enough. The next decisive check is **closed-loop learning with the compressed physical gradient**.

At every relinearization step:

1. run the normal forward task;
2. construct the soma derivative and reciprocal retro source;
3. select K bins from boundary spectra only;
4. build the approximate local gradient from only those bins;
5. update the same graded frontier conductances;
6. repeat.

Compare K=8 and K=16 with the exact adjoint under the identical candidate set, step normalization, and number of iterations.

If K=8/16 preserve most of the exact learner's objective gain, the compression has crossed from **map reconstruction** into **useful hardware training**.

## Wall sentence

> **The time-domain adjoint does not require every tuner to remember the whole waveform: for this task, the required forward×backward correlation is spectrally sparse, and a common frequency set chosen from boundary signals reconstructs the internal gradient with high fidelity.**
