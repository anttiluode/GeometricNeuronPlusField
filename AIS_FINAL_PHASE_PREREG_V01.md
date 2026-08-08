# AIS final phase/interface test v0.1 — preregistration

## Why this test exists

The active-boundary v0.2 result used the historical FunctionalArbor soma readout
`|psi_soma|^2`. That is a power/envelope observable. Claude's audit pointed out
that both `|psi|^2` and `|psi|` discard the complex carrier phase, while
`Re(psi)` preserves a signed carrier waveform.

This test changes **only that interface variable**. It is the last planned phase
test before any AIS geometry/position/extent co-adaptation.

## Frozen pieces

No retuning of:

- frozen FunctionalArbor body generation or carrier physics;
- HH conductances, reversal potentials, capacitance or integration step;
- AIS input gain (`15`) or clipping magnitude (`2.5`);
- frequency battery (`0.00625, 0.0125, 0.025, 0.05, 0.0833333, 0.125`);
- burn-in (`160`) or trace length (`640`);
- v0.2 exact per-frequency event-count matching;
- minimum timing exposure (`>=4` active events);
- registered upper band (`f >= 0.05` cycles/frame).

## Experimental factor

The complex soma value is recorded after each frozen field update and converted
to one of three scalar feeds:

```text
power       |psi_soma|^2
magnitude   |psi_soma|
real        Re(psi_soma)
```

Each feed is normalized **within body and feed** by the 95th percentile of its
absolute value over the whole registered frequency battery.

Power and magnitude are clipped to `[0, 2.5]`. `Re(psi)` is signed and clipped
to `[-2.5, 2.5]`. The same fixed input gain then converts the normalized feed to
injected current.

The linearized control receives that exact same clipped current and uses the
measured small-signal impulse response of the unchanged HH compartment.
The memoryless control ranks the exact same current samples.

## Primary test

The primary arm is `Re(psi)`.

For every body and upper-band frequency with at least four active spikes:

1. exactly match the linearized control to the active event count;
2. compute envelope-phase PPC for active and linearized events;
3. within each body, average `PPC_active - PPC_linearized` over its valid
   upper-band frequencies.

This yields **one primary number per body**, avoiding the old pseudo-replication
problem of treating multiple frequencies from one body as independent
organisms.

Registered success rule for saying that phase earns a special role in active
eventization:

```text
>= 8 valid bodies
median body-level delta PPC > 0
one-sided Wilcoxon(active > linearized) p < 0.05
```

If this rule fails, do not tune HH parameters or move the frequency window.

## Secondary receipts

- body-level sign test;
- per-frequency PPC for all three feeds;
- carrier-phase PPC using `carrier_omega / 2pi`;
- active event counts by feed/frequency;
- memoryless comparison.

The carrier-phase metric is descriptive. The primary claim remains whether the
nonlinear active boundary extracts more envelope timing precision from a
phase-preserving input than its own linearization.

## Interpretation boundary

A failure does **not** prove that physical neuronal phase can never matter. It
means that in this lineage, after frozen growth and the soma bottleneck, giving
the present active boundary the remaining phase-bearing scalar does not earn a
timing advantage over the corresponding linear filter.

No AIS position/length/extent co-adaptation happens before this receipt is
read.
