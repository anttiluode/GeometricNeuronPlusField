# Soma mode-pair mechanism v0.1 — sparse cross-mode mixing and a blind common-mode reference

This document records three linked experiments:

1. `MODE_PAIR_DISCOVERY_PREREG_V01.md` / `mode_pair_discovery_probe.py`, discovery seeds 72-83;
2. `MODE_PAIR_CONFIRM_PREREG_V01.md` / `mode_pair_confirm_probe.py`, held-out seeds 84-95;
3. `COMMON_MODE_MIXER_CONFIRM_PREREG_V01.md` / `common_mode_mixer_confirm_probe.py`, held-out seeds 96-107 after an exploratory development run on reused seeds 72-95.

The question was:

> **Which geometry-defined mode pairs create the soma's order-sensitive coherent cross term?**

## 1. Exact pair decomposition

For each body the full graph-Laplacian basis was used. Single-source fields were projected into modal coefficients

```text
q_A,n(t), q_B,n(t)
```

and reconstructed at the soma through

```text
u_A,n(t) = q_A,n(t) phi_n(s)
u_B,n(t) = q_B,n(t) phi_n(s).
```

At the target and distractor soma-power peak times the order-sensitive cross-source pair matrix was

```text
M_nm = 2 Re[u_A,n(t_T) conj(u_B,m(t_T-tau))]
       - 2 Re[u_B,n(t_D) conj(u_A,m(t_D-tau))].
```

Modal reconstruction errors were at machine precision in both discovery and confirmation (`~1e-15` relative for single-source reconstruction; `~1e-14` or better for the pair sum).

So the soma cross term can be decomposed exactly into graph-mode pair interactions.

## 2. The pair support is sparse

Discovery seeds 72-83:

```text
mean unordered-pair fraction for 50% absolute mass     0.01066
mean unordered-pair fraction for 80% absolute mass     0.04638
mean participation-ratio effective pair fraction       0.02637
```

Held-out seeds 84-95:

```text
mean unordered-pair fraction for 50% absolute mass     0.01154
mean unordered-pair fraction for 80% absolute mass     0.04946
mean participation-ratio effective pair fraction       0.02738
```

All registered sparsity thresholds passed.

With 70 graph modes there are 2485 unordered pairs. Roughly 1.1% of those pairs carry half of the absolute order-sensitive cross-term mass and about 5% carry 80%.

This is not a broad all-pairs soup.

## 3. The computation is overwhelmingly cross-mode

Discovery mean diagonal `n=n` absolute-pair fraction:

```text
0.02948
```

Held-out confirmation:

```text
0.02167
```

The preregistered `<0.06` prediction passed.

So about 97-98% of absolute pair mass comes from interactions between *different* graph modes, not A and B occupying the same mode and simply interfering there.

This is the central structural result:

> **The local soma power readout is acting as a mode mixer.**

A point readout at soma `s` is

```text
y_s = |sum_n q_n phi_n(s)|^2.
```

In modal coordinates this is a quadratic form

```text
y_s = q^* K_s q
K_s = v_s v_s^T
v_s = [phi_0(s), phi_1(s), ...].
```

The off-diagonal entries of the rank-1 kernel `K_s` are exactly the mode-pair mixing terms. The soma is therefore not a diagonal mode selector. A local point measurement is intrinsically a dense modal mixer once it is squared.

## 4. The surprising hub is mode 0

The unique constant graph mode was the highest mean-involvement mode in both discovery and held-out confirmation.

Discovery:

```text
mode-0 involvement fraction                  0.1891
approx. absolute pair mass involving mode 0  ~0.36
```

Held-out confirmation:

```text
mode-0 involvement fraction                  0.1700
absolute pair mass involving mode 0           0.3318
highest mean-involvement mode index           0
```

All registered common-mode-hub predictions passed.

This is striking because mode 0 was already known to be almost perfectly blind to temporal order when its own coherent power was read directly.

The dominant discovery pairs were mostly of the form

```text
(0,3), (0,2), (0,7), (0,8), (0,6), (0,11), (0,10), (0,13), ...
```

rather than pairs within the old directly informative band.

## 5. The old 18-20 band is not the dominant soma mixer

Modes 18-20 were the held-out positive band in the earlier *direct coherent mode readout* experiment.

But their soma pair involvement was not enriched.

Discovery mean enrichment relative to dimensional share:

```text
0.727
```

Held-out confirmation:

```text
0.464
```

The preregistered `<1.0` confirmation prediction passed.

Meanwhile modes 0-17 carried:

```text
discovery      0.753 of total mode involvement
confirmation   0.763
```

This resolves an apparent contradiction. `|q_n|^2` asks whether a mode is informative when observed *alone*. Soma point power asks what happens after a spatial point projection mixes modes and squares the result. These are different quadratic readouts, so the directly informative mode band need not dominate the local soma interaction.

## 6. Absolute interaction mass is not the same as causal necessity

The held-out confirmation also made one important correction to the tempting story.

Removing only the cross-source interactions between mode 0 and nonzero modes did **not** reduce mean absolute soma contrast:

```text
mean |C| full                         0.16597
mean |C| without mode0/nonzero mix    0.16750
```

The signed result did change (`mean signed error ~0.0194`), but the magnitude did not collapse.

So mode 0 is a large interaction hub, but it is not legitimate to say that the task disappears without those particular cross-source pairs. The pair system contains substantial signed cancellation and redundancy.

Likewise:

```text
mean |C| diagonal cross only          0.10249
mean |C| off-diagonal cross only      0.18039
```

The off-diagonal structure is much closer to the full computation, but different subsets can compensate for one another.

## 7. Common mode as a reference component

The confirmed pair hub suggested a more precise test. Decompose the total soma field into

```text
psi_soma = c + r
```

where `c` is graph mode 0 and `r` is the sum of all nonconstant modes. Then

```text
|psi_soma|^2 = |c|^2 + 2 Re[c conj(r)] + |r|^2.
```

The exploratory 24-body run found:

```text
mean |C| common mode alone        0.000068
corr(first-order mixer, full C)   0.9813
mean |C_ref - C_full|             0.0576
sign match                        23/24
```

A fresh preregistered confirmation on seeds 96-107 then gave:

```text
mean |C| full                     0.16173
mean |C| common only              0.0000748
mean |C| first-order reference    0.11578
mean |C| residual only            0.23342

corr(reference, full)             0.95771
corr(residual, full)              0.97087

mean absolute error
reference -> full                 0.06768
residual  -> full                 0.09283

reference lower error            10/12 bodies
reference sign match             10/12 bodies
```

All three preregistered common-mode-mixer predictions passed.

The nuance matters: residual-only contrast is itself highly informative and its cross-body correlation is slightly higher in this confirmation. Mode 0 is therefore **not the sole carrier of the code** and not a magical oscillator that creates information from nothing.

The cleaner interpretation is:

> **Mode 0 is almost order-blind alone, while its first-order cross term with nonconstant modes acts as a reference/calibration component that substantially improves reconstruction of the full local quadratic readout.**

This is analogous to reference-assisted square-law detection, but the analogy should not be promoted into a biological claim.

## Current mechanism

The strongest current reduction is now:

```text
anatomy G
   -> graph spectrum {lambda_n, phi_n}
   -> source-specific modal histories q_A,n(t), q_B,n(t)
   -> local soma vector v_s = phi_n(s)
   -> rank-1 quadratic kernel K_s = v_s v_s^T
   -> sparse, mostly off-diagonal mode-pair mixing
   -> coherent source-source cross term
   -> temporal-order contrast
```

Amplitude balance remains an opportunity variable: it controls where both source histories have enough local leverage to interact. The actual order-sensitive value depends on the geometry-shaped temporal/complex histories and the quadratic modal mixing at the readout.

## Wall sentence

> **The soma does not read a winning graph mode. It performs a local point projection whose squared magnitude mixes graph modes pairwise. Temporal-order computation is concentrated in a sparse set of mostly off-diagonal interactions; the constant mode is blind alone but acts as a large reference/calibration hub for nonconstant geometry-shaped modes.**

## Next clean question

Every point readout has a geometry-only modal vector `v_x` and therefore a rank-1 quadratic kernel

```text
K_x = v_x v_x^T.
```

The next question is therefore no longer "which soma taps?" It is:

> **Why is the soma's quadratic kernel aligned with the temporal-order interaction produced by the two sources better than the kernels of most other cells?**

That can be tested by constructing a task interaction operator in modal coordinates from the two single-source histories, then measuring its alignment with `K_x` across every body cell. It would join the amplitude-balance result and the mode-pair result in one object, without adding another biological compartment.
