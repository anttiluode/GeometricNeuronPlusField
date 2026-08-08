# Benchmark scope v0.1 — what a graph wave can and cannot beat

The exact adjoint makes a matched benchmark practical, but one conceptual point has to be fixed before writing it.

## A fully free pole bank is a superset of the linear graph wave

For frozen reciprocal linear coupling,

```text
psi'' + gamma psi' + (rho I - c L) psi = B u
```

with symmetric `L`, diagonalize

```text
L = Phi Lambda Phi^T.
```

In modal coordinates `q=Phi^T psi`, every mode is an independent damped second-order oscillator:

```text
q_n'' + gamma q_n' + (rho + c lambda_n) q_n = b_n^T u.
```

A local linear observation is a weighted sum of those modal states, followed here by the square-law readout.

Therefore the complete frozen graph transfer function is already a pole/residue bank. A sufficiently unconstrained free-pole model with the same number of dynamical modes and free residues contains the graph model as a special case.

So this question is **not scientifically meaningful**:

> “Can the graph beat a fully free pole bank at the same number of modes under globally optimal training?”

The free bank cannot be less expressive if it truly includes all graph poles and residues. At best the graph ties it.

## The meaningful claim is parameterization efficiency

The graph can still win a different and much more hardware-relevant contest:

> **How much useful temporal computation is obtained per trainable physical scalar / tuner?**

One local conductance changes many eigenvalues and residues coherently because the poles are tied together by geometry. That structured parameter sharing is exactly what a physical scattering mesh buys.

The benchmark must therefore distinguish:

```text
state count / number of poles
trainable scalar count / number of tuners
fixed inherited structure
training evaluations / passes
hardware operations
```

## Benchmark ladder

### B0 — containment sanity check

Construct the exact pole/residue decomposition of the same reciprocal operator and verify that the full free-pole reconstruction reproduces the graph response. This is expected to pass and prevents us from making an impossible expressivity claim.

### B1 — matched trainable-scalar benchmark

Use the same frozen base body and the same task distribution. Give each arm exactly `P=8` trainable real scalars and the same number of relinearized update iterations / objective evaluations.

Candidate arms:

```text
G8   eight local frontier conductance tuners
F8   eight directly tunable modal pole-frequency coordinates
R8   eight directly tunable modal residue/gain coordinates
FR8  four pole-frequency + four residue coordinates
```

The free coordinates should be chosen by a fixed discovery rule before held-out bodies (for example, largest base sensitivity), so the baseline is not deliberately weak.

The graph wins only if it beats the strongest free-coordinate arm on held-out task performance at the same tuner count and evaluation budget.

### B2 — generalization, not one-lag memorization

A single lag is too easy to overfit. Train on a fixed set of temporal separations and evaluate on interleaved held-out separations and amplitude perturbations. Both systems receive identical labelled source histories and the same square-law output convention.

### B3 — hardware-cost benchmark

Only after B1/B2 should we compare physical cost:

```text
number of tuners
number of local monitors
forward/backward passes per update
loss / attenuation budget
phase-setting precision
thermal drift sensitivity
latency and stored state
```

This is where a local reciprocal mesh can plausibly beat a digitally parameterized free-pole model even when the latter is more expressive in abstract function space.

## What would count as the architectural result

Not:

> “gradient ascent improves the graph.”

That is now a sanity check.

The interesting result would be:

> **At the same number of physical tuners and the same training budget, a locally coupled reciprocal geometry reaches better held-out temporal performance or robustness than direct tuning of an equally sized set of unconstrained modal coordinates.**

That would establish a useful inductive bias / physical parameter-sharing advantage, not magical expressivity beyond a general pole bank.

## Immediate next step

Implement B0 first. It should be an identity test. Then use that exact modal representation as the common substrate for B1 so that graph-vs-pole differences cannot be blamed on different simulators.
