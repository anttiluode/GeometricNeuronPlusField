# Random frontier control v0.1 — graded coupling survives the selection-effect test

`BOND_RESPONSE_V01.md` reported that 80.5% of held-out additions selected because the base adjoint derivative was positive had an interior optimum on the bath-to-arbor conductance path. Claude correctly identified the missing control: that conditioning itself forces the curve to rise initially and could manufacture an apparent interior-optimum enrichment.

`RANDOM_FRONTIER_CONTROL_PREREG_V01.md` therefore fixed fresh seeds 276-287 and selected frontier additions **before looking at the gradient sign**.

## Fresh preregistered result

Across 12 bodies, 72 uniformly sampled legal tip-like additions were swept on the same fixed alpha grid.

```text
random additions                           72
interior optimum                           53 / 72 = 0.7361
alpha_best = 0                             12 / 72 = 0.1667
alpha_best = 1                              7 / 72 = 0.0972
median alpha_best among interior events                  0.15
slope-sign reversal                        63 / 72 = 0.8750
mean binary regret                                      0.01459
positive binary regret                    64 / 72 = 0.8889
```

The preregistered `BROAD_GRADED` threshold was random interior fraction `>=0.70`.

**Verdict: BROAD_GRADED.**

The missing control does **not** collapse the interior-optimum effect.

## There is a selection effect, but it is not the explanation

After random selection, the events happened to split exactly 36/36 by the sign of the base directional derivative.

Positive-derivative subset:

```text
n                                           36
interior optimum                            31 / 36 = 0.8611
alpha_best = 0                               0 / 36 = 0
alpha_best = 1                               5 / 36 = 0.1389
median interior alpha_best                                0.03
slope reversal                              32 / 36 = 0.8889
mean binary regret                                        0.01556
```

Negative-derivative subset:

```text
n                                           36
interior optimum                            22 / 36 = 0.6111
alpha_best = 0                              12 / 36 = 0.3333
alpha_best = 1                               2 / 36 = 0.0556
median interior alpha_best                                0.20
slope reversal                              31 / 36 = 0.8611
mean binary regret                                        0.01363
```

So gradient conditioning raises the interior fraction on these fresh bodies by about `0.125` (`0.8611 - 0.7361`). That is a real selection effect.

But the unconditioned random fraction remains 73.6%, above the preregistered broad-graded threshold, and even the initially negative subset has 61.1% interior optima.

The earlier result was therefore **enriched by selection, not created by selection**.

## A stronger nonmonotonicity appears

Descriptively, across all 72 random additions:

```text
some alpha has positive gain               60 / 72 = 0.8333
full binary endpoint has positive gain      35 / 72 = 0.4861
mean best gain                                          +0.01248
mean full-binary gain                                   -0.00212
```

The initially negative-gradient subset is especially informative. Two thirds of those events eventually attain a positive best gain somewhere on the path, despite beginning by moving in the wrong direction. Their median interior optimum is farther out (`alpha=0.20`) than the positive-gradient subset (`alpha=0.03`).

That means the response landscape is not merely “positive slope followed by saturation.” It commonly contains genuine slope reversals and, for some bonds, valleys before a better finite-coupling state.

This does **not** mean a local gradient learner can discover those negative-first valleys without exploration; it means the structural response surface itself is strongly nonmonotone.

## What is now earned

The graded-medium interpretation survives its most obvious selection confound:

> **Intermediate coupling is commonly preferred even for frontier bonds chosen without gradient conditioning. The binary arbor is therefore not merely a convenient endpoint representation; in this model, conductance magnitude itself is often a computational degree of freedom.**

The stronger formulation should still stay local to this model. We have not shown that arbitrary material pixels in arbitrary wave devices prefer intermediate values, nor that a final hardware implementation should use exactly this alpha scale.

## Consequence for the benchmark

This makes thresholding a continuously trained structure back to binary particularly inappropriate as the default benchmark. The spatial model should be compared to a conventional temporal/filter model in its **native continuous-coupling state space**, with the same number of trainable degrees of freedom and the same training/evaluation budget.

That is now the next architectural test.
