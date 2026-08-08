# Modal locality audit v0.1

Canonical GitHub Actions run: `31241224985`

This is the Claude-branch question produced by the modal reduction:

> If the mature arbor computes in global graph modes, is a one-cell structural change really a *local* event in those computational coordinates?

Twelve frozen 70-cell bodies were each given 25 legal one-cell additions.  For each perturbation the body Laplacian was recomputed and its spectrum matched back to the unperturbed modes.  No wave simulation, learning, or fitted parameter appears in this audit.

## Receipt

```text
bodies                                      12
single-cell additions                       300
modes/body                                  70

modes with best overlap < 0.90
    mean                                    16.86 / 70
    median                                  17
    range                                   5 - 43
    fraction of spectrum                    0.2409

mean best mode overlap                      0.9308
confirmed band 18-20 identity loss          0.0804
||delta lambda|| / median eigenvalue gap    11.09
near-degenerate neighboring pairs/body      13.25

corr(distance from soma, modes scrambled)  -0.449
corr(distance from soma, band loss)         -0.191
```

The spatial distribution of the eigenvector change is only weakly localized around the added cell:

```text
radius    delta-phi mass    cell share    enrichment
1           0.0486           0.0327         1.486x
2           0.0966           0.0617         1.567x
3           0.1403           0.0911         1.541x
5           0.2287           0.1562         1.465x
```

The roughly flat ~1.5x enrichment is the important shape: there is a local excess, but most of the modal perturbation is distributed across the body rather than decaying rapidly away from the structural event.

Near-degeneracies can rotate individual eigenvectors and therefore inflate a raw mode-scramble count.  They do **not** explain away the 11x eigenvalue displacement relative to the typical gap, the task-band identity loss, or the broad spatial distribution of `|delta phi|^2`.

## Interpretation

The strong claim should be phrased carefully.

This audit does not mathematically prove that *no conceivable local biochemical eligibility signal* could ever support useful structural learning.  It does show that in this particular fixed-K graph-wave toy, the computational coordinates identified by the modal reduction are globally entangled with local anatomy.

A local event such as "this cell was born" is therefore not a local event in the coordinates that determine the field transfer function.

That makes the v0.8/v0.9 credit-assignment null much less mysterious:

```text
local anatomical event
       |
       v
many simultaneous modal changes
       |
       v
global task consequence
```

v0.5's atomic accept/reject could succeed because it evaluated the **whole structural counterfactual** and either kept or reverted it.  v0.7-v0.9 tried to distribute delayed scalar consequence back onto specific locations.  The modal audit shows why that attribution is intrinsically difficult in the current representation.

The location result is also interesting: perturbations closer to the soma tend to scramble more of the spectrum (`r = -0.449`).  That does not make the effect local, but it says the soma-near region has greater **spectral leverage** over the global operator.

This is one reason the AIS bridge is worth testing: biology places a small active and plastic output compartment near the convergence point rather than trying to make every dendritic compartment an equally powerful global controller.

## Wall sentence

> **In the current toy, one local anatomical edit changes a global filter bank: about a quarter of the modal identities are substantially disturbed, with only weak spatial localization of the eigenvector change.  The credit problem is therefore partly a coordinate problem — local anatomy and global computational modes are not aligned.**
