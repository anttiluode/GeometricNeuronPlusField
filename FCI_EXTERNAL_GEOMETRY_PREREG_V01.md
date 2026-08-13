# External FCI geometry gate v0.1 — preregistration

## Status

**Frozen before retrieving or transcribing the 24 cell-level FCI labels.**

External source:

> Aizenbud et al. (2026), *Dendritic morphology and synaptic nonlinearities enhance
> functional complexity in human cortical neurons*, PNAS 123(28):e2533168123.

The authors released the reconstructed morphologies and neuron models in
`ido4848/FCI`.

This is an external test of the surviving upstream GeometricNeuronPlusField claim:

> frozen morphology can be usefully treated as an operator that constrains the modes,
> time scales, and local mixtures available to a moving field.

It is **not** a test of the failed generic-HH/AIS claims and must not be used to rescue
them.

---

## 1. Why this dataset is unusually hostile/useful

Aizenbud et al. already report strong simple morphology baselines for their Functional
Complexity Index (FCI):

```text
total dendritic area                         R^2 ~ 0.74
longest bifurcation branch                   R^2 ~ 0.44
number of bifurcation branches               R^2 ~ 0.29
area + longest bifurcation branch            R^2 ~ 0.81
best 3 reported morphology features          R^2 ~ 0.85
best 4 reported morphology features          R^2 ~ 0.88
```

These are in-sample explanatory values in the paper, not our predictive scores.

Therefore the question is **not**:

> does geometry correlate with complexity?

That is already established in the paper.

The question is:

> **does a frozen graph/operator description add held-out information about neuronal
> functional complexity beyond ordinary size/path morphometrics?**

If not, this branch loses cleanly.

---

## 2. Data firewall

Before any cell-level FCI values are loaded, freeze:

1. morphology parser and graph reduction;
2. baseline features;
3. Geometric Neuron features;
4. prediction model;
5. cross-validation scheme;
6. verdict rule.

Allowed before labels:

- inspect morphology file formats;
- confirm the number and identity of available morphologies;
- debug topology reconstruction;
- compare our basic morphometrics against NeuroM outputs without FCI;
- inspect species/layer metadata.

Not allowed before freezing:

- rank cells by FCI;
- choose graph metrics because they correlate with FCI;
- choose mode count by FCI performance;
- discard cells because they hurt the result.

The paper's published aggregate claims and two exemplar FCI values are already known;
those do not determine the 24-cell mapping and are not to be used for feature tuning.

---

## 3. Frozen graph construction

For each reconstructed dendritic morphology:

1. identify the soma/root;
2. retain dendritic cable geometry only for the primary morphology graph;
3. collapse consecutive reconstruction samples between structural events into one
   branch edge;
4. structural nodes are:
   - soma/root,
   - bifurcation points,
   - terminal tips;
5. for each collapsed edge store:
   - cable length `ell` = sum of Euclidean segment lengths;
   - approximate membrane area `A` = sum of truncated-cylinder lateral areas using
     local diameters;
6. exclude axon from the primary morphology graph if present, matching the paper's
   dendritic morphology analysis.

This collapse is important: graph spectra must not depend mainly on how densely a
reconstructor sampled points along the same cable.

### Two operators

Primary cable operator:

```text
w_ij = 1 / ell_ij
L = D - W
```

with a node mass equal to one half of each incident edge membrane area assigned to the
node.  Solve the symmetric mass-normalized operator

```text
L_M = M^(-1/2) L M^(-1/2)
```

where zero/degenerate masses are rejected as parser failures, not silently repaired.

Negative-control operator:

```text
unweighted topology-only graph
```

The topology-only arm is descriptive and cannot promote the primary gate.

---

## 4. Frozen ordinary morphology baselines

Compute independently from the same collapsed tree:

```text
log_total_dendritic_area
log_total_dendritic_length
log_longest_root_to_tip_path
number_of_bifurcations
```

The **primary baseline** is intentionally small and paper-motivated:

```text
B2 = [log_total_dendritic_area,
      log_longest_root_to_tip_path]
```

`longest_root_to_tip_path` is our reproducible cable-graph analogue of the paper's
strong `longest bifurcation branch / maximal path-distance` family.  We will also report
single-area and B4 baselines, but B2 is the comparison frozen for the gate.

If later NeuroM comparison shows our path quantity is not the same quantity as the
paper's reported longest bifurcation branch, rename it honestly; do not tune its
definition against FCI.

---

## 5. Frozen Geometric Neuron features

Use the first `K = 16` **nonzero** eigenpairs of `L_M`, or all available nonzero modes
when a morphology has fewer than 16.  `K=16` is frozen from the existing small-bank
Geometric Neuron style; it is not selected against this dataset.

Exactly three primary operator features are permitted:

### G1 — low-spectrum entropy

Normalize the first-K nonzero eigenvalues by their sum:

```text
p_k = lambda_k / sum(lambda_1 ... lambda_K)
G1 = -sum p_k log(p_k) / log(K)
```

Interpretation: how evenly the low operator spectrum allocates characteristic spatial
scales.

### G2 — root modal participation entropy

At the soma/root node, use squared mass-normalized eigenvector loading:

```text
q_k = phi_k(root)^2 / sum_j phi_j(root)^2
G2 = -sum q_k log(q_k) / log(K)
```

Interpretation: whether the root is dominated by a few low modes or mixes many of
them.  This is the external analogue of the repo's local-observability/root-mixture
question.

### G3 — low-mode spacing irregularity

For adjacent first-K nonzero eigenvalues:

```text
d_k = log(lambda_(k+1)) - log(lambda_k)
G3 = std(d_k) / (abs(mean(d_k)) + eps)
```

Interpretation: irregularity of the low resonator spacing after using log spacing to
reduce absolute size scaling.

No additional spectral/topological feature may be added to the confirmatory gate after
FCI values are opened.  Additional metrics belong to a later explicitly exploratory
analysis.

---

## 6. Target hierarchy

### Preferred primary target: morphology-isolated FCI

If the paper/SI provides **per-cell FCI for the experiment in which the same rat-type
synapses are assigned to both human and rat morphologies**, use that as the primary
target.  This isolates morphology from species-specific synaptic parameters and is the
cleanest question for this repo.

### Fallback target: main FCI

If individual morphology-isolated FCI values are not publicly recoverable without
rerunning the full biophysical/DNN pipeline, use the main per-cell FCI values only as a
secondary/fallback target and include species as an explicit nuisance/baseline feature.
Do not pretend a morphology-only metric should explain species-specific NMDA effects.

The choice is determined by data availability, not performance.

---

## 7. Frozen predictive models

Because `n = 24` is tiny, no nonlinear learner and no hyperparameter search.

All features are standardized inside each training fold.

Use ridge regression with fixed `alpha = 1.0` for:

```text
M0: intercept / nuisance only
B2: two-feature morphology baseline
B2+G: B2 + G1 + G2 + G3
B4: four ordinary morphology features
```

For fallback main-FCI analysis, species is included in every model including M0.
Layer is reported descriptively but not one-hot-expanded into the primary 24-cell fit.

---

## 8. Frozen validation

Primary: leave-one-cell-out cross-validation (LOOCV), producing one out-of-sample
prediction for every morphology.

Report:

```text
CV R^2
MAE
Spearman rho(prediction, target)
```

The primary comparison is paired absolute error per held-out cell:

```text
Delta_i = |y_i - yhat_B2_i| - |y_i - yhat_B2+G_i|
```

Positive `Delta` means graph features improved the held-out prediction.

Inference:

- mean and median Delta;
- bootstrap 95% CI of mean Delta;
- two-sided exact/sign-flip test across the 24 held-out errors;
- number of cells with Delta > 0.

Secondary robustness checks, reported but not used to rescue the primary gate:

- leave-one-layer-out;
- human-only and rat-only descriptive fits when sample size permits;
- topology-only operator features;
- B4 versus B4+G.

---

## 9. Gate

Call only

```text
OPERATOR_ADDS_EXTERNAL_COMPLEXITY_SIGNAL
```

if all are true for the primary target:

1. `CV_R2(B2+G) > CV_R2(B2)`;
2. `MAE(B2+G) <= 0.90 * MAE(B2)` (at least 10% relative improvement);
3. mean paired error improvement `Delta > 0`;
4. bootstrap 95% CI of mean Delta is entirely above 0;
5. two-sided sign-flip `p < 0.05`.

Otherwise:

```text
NO_EXTERNAL_OPERATOR_ADVANTAGE
```

A correlation of any single G feature with FCI does not pass the gate.
A human-vs-rat separation does not pass the gate.
An in-sample R^2 does not pass the gate.

---

## 10. What each outcome means

### If the gate passes

Earned statement:

> a small, frozen operator description of real dendritic trees carries predictive
> information about single-neuron I/O complexity beyond strong ordinary morphology
> baselines in this external 24-cell set.

Still **not** earned:

- biological neurons calculate graph eigenvectors;
- our synthetic wave equations are the biological mechanism;
- the original checkerboard story was biologically correct;
- AIS/eventization claims;
- a general learning advantage.

### If the gate fails

Earned statement:

> on this external set, the Geometric Neuron modal language does not add useful
> held-out information beyond ordinary dendritic size/path descriptors.

That is a valuable result.  It would tell us that the real external morphology result
is already captured at a simpler level than our operator description.

---

## 11. Why this is a good next test

This dataset was not designed by us, its neuron morphologies were not generated by our
code, its target FCI was not invented for Geometric Neuron, and its simple baselines are
already strong.

That is exactly the kind of external pressure this repo needs.
