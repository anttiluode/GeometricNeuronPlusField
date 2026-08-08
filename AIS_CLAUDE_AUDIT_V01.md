# Claude HH-ablation audit — correction note

Claude independently checked the instantaneous-gate experiment using standard
Hodgkin-Huxley dynamics over a broad constant-current sweep. The important
correction is that the old `h_instant -> 0/24` and `n_instant -> 0/24` result is
a **dead-arm result**, not evidence that those histories specifically create the
observed frequency allocation.

The independent check reproduces the qualitative asymmetry:

```text
m instantaneous   remains strongly excitable
h instantaneous   essentially/silently loses repetitive firing
n instantaneous   loses repetitive firing
```

Thus the old 0/24 ablation cannot identify a frequency-selection mechanism
without first restoring a comparable operating/firing regime. The historical
AIS_GATE_ABLATION_V01 receipt is retained, but its boxed interpretation is
superseded by this note.

The audit also caught a statistical wording issue in AIS_ACTIVE_V02. The pooled
21 body/frequency pairs gave a sign-test p of 0.078, so the win/loss count alone
was directional rather than formally significant. Re-analysis of the stored
paired magnitudes gives:

```text
21 pair Wilcoxon, two-sided p ~= 0.0175
```

but those 21 observations contain repeated frequencies from the same bodies.
The cleaner organism-level reduction averages the valid upper-band deltas
within each of the 11 bodies:

```text
mean body delta PPC (active-linearized)  ~= -0.1348
median body delta                        ~= -0.0849
active better / worse bodies             2 / 9
body-level Wilcoxon, two-sided p          ~= 0.0322
```

So the defensible statement is not simply "the active compartment failed to
beat the linearization." In this registered upper-band sample its nonlinear
eventization **degraded timing precision relative to its own linearized
response**, with the body-level paired magnitudes supporting that direction.

Claude's proposed working interpretation is therefore:

> the active compartment may buy event selection at a cost in timing precision.

That tradeoff still needs a direct information-coding test. Before doing that,
the more fundamental interface issue is tested first: `|psi|^2` and `|psi|`
remove carrier phase, whereas `Re(psi)` preserves it. The next registered
experiment is therefore AIS_FINAL_PHASE_PREREG_V01.
