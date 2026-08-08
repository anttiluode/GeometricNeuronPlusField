# AIS selection-information test v0.1 — preregistration

## Why this test exists

The active-boundary lineage now has a consistent negative timing result: the
full HH-like eventizer is less phase-precise than its own rate-matched
linearization, including when the phase-bearing `Re(psi_soma)` signal is supplied.

What remains is a narrower positive hypothesis suggested by Claude's audit:

> **the nonlinear boundary may buy event selection at a cost in timing precision.**

The historical power-interface runs visibly allocated their finite spike budget
very unevenly across modulation frequencies.  This test asks whether that
allocation carries more information about the input frequency than a
memoryless detector or the same membrane's linearization when all three have
**exactly the same total event budget over the full frequency battery**.

This is a mechanism/utility test on the existing model family, not an
independent replication: the frequency battery and qualitative power-allocation
phenomenon have already been seen.

## Frozen pieces

- same 24 frozen FunctionalArbor seeds;
- same carrier physics and modulation frequencies;
- same 640-frame traces and burn-in 160;
- same HH parameters, integration, gain and clipping magnitude;
- same small-signal linearization of that HH membrane;
- same `power`, `magnitude`, and `real` interface definitions and per-body/feed
  q95 absolute normalization used by AIS_FINAL_PHASE_V01.

No HH or upstream parameter is tuned.

## Event-budget matching

For each body and feed:

1. run the full active boundary naturally on all six frequency conditions;
2. let `K` be its total number of post-burn spikes across the six equal-duration
   conditions;
3. pool all allowed post-burn frames from all six conditions;
4. select exactly the top `K` frames for the memoryless current score;
5. select exactly the top `K` frames for the linearized-membrane score.

Thus total event count is identical while **allocation across input frequencies
is free to differ**.

## Primary feed

`Re(psi_soma)` is primary because it is the only tested scalar that preserves
the signed carrier and is the most voltage-like interface in this toy.

`|psi|^2` is retained as the historical reference and `|psi|` as an envelope
control.

## Primary statistic: frequency information per spike

For an encoder with counts `n_f` over the six equally long/equally weighted
frequency conditions, define

```text
p(f | spike) = n_f / sum_f n_f

I_spike = KL[p(f | spike) || Uniform(6)]
        = log2(6) - H[p(f | spike)]
```

`I_spike` is in **bits per spike occurrence**.  It is zero if spikes are evenly
allocated across frequencies and reaches `log2(6)` if all spikes occur under
one frequency.

A body is valid if the active arm emits at least 12 total post-burn spikes over
the battery.  The control has the same `K` by construction.

Primary paired quantity per body:

```text
Delta I = I_spike(active) - I_spike(linearized)
```

Registered success rule for saying the active boundary earns a positive
frequency-selection role:

```text
>= 12 valid Re(psi) bodies
median Delta I > 0
one-sided Wilcoxon(active > linearized) p < 0.05
```

Also report the exact sign test and active-vs-memoryless comparison.

## Secondary statistic: full binary event information

Treat a uniformly sampled post-burn frame as `Y = spike/no spike` and frequency
condition as `F`.  Compute `I(F;Y)` in bits/frame.  Because all encoders are
exactly total-rate matched, differences reflect frequency-dependent allocation,
including information carried by both spikes and silences.

This is secondary; it cannot replace the registered bits-per-spike test.

## Interpretation boundary

A positive result means only that this nonlinear boundary makes **event
occurrence more frequency-selective per emitted spike** under the frozen toy
battery.  It does not establish biological optimality, task usefulness, or AIS
geometry.

A null/negative result means the surviving "selection at the cost of precision"
story has not earned a positive computational role.  Do not proceed to AIS
position/length co-adaptation from a null.
