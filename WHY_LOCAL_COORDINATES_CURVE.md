# Why local bond coordinates curve — exact spectral perturbation geometry

`COORDINATE_CURVATURE_CONFIRM_V01.md` found an empirical fact that needs an analytical explanation:

> At matched task-space trust radius, a local bond coordinate departs from its own first-order prediction about twenty times faster than a direct modal coordinate.

For a reciprocal symmetric wave operator, the reason is visible directly from eigenvalue perturbation theory.

## 1. One local bond is rank one in physical coordinates

Let the positive stiffness/Laplacian operator be

```text
K = sum_e k_e b_e b_e^T
```

where for bond `e=(a,b)`

```text
b_e = e_a - e_b.
```

Changing one local coupling by `delta k` gives

```text
Delta K_e = delta k b_e b_e^T.
```

This is extremely local and rank one in the physical node basis.

## 2. The same local edit is dense in the modal basis

Diagonalize the reciprocal operator:

```text
K Phi = Phi Lambda
Phi^T Phi = I.
```

Define the edge-difference participation of mode `n`

```text
z_n = phi_n(a) - phi_n(b).
```

Then in modal coordinates

```text
Phi^T Delta K_e Phi
    = delta k z z^T.
```

Every modal pair `(n,m)` with nonzero edge differences receives coupling

```text
Delta K_nm = delta k z_n z_m.
```

So the physical/local and modal/direct coordinate systems have exactly opposite sparsity:

```text
local bond coordinate
    sparse rank-1 in node basis
    dense rank-1 in modal basis

direct modal pole coordinate
    sparse/diagonal in modal basis
    dense phi_n phi_n^T in node basis.
```

This is the algebraic form of the compiler trade.

## 3. First-order pole motion is simple

For a nondegenerate mode,

```text
d lambda_n / d k_e = z_n^2.
```

Thus one physical bond moves **many modal frequencies at once**, with nonnegative first-order shifts weighted by their squared strain across the bond.

A direct spectral `F_n` coordinate instead moves one chosen pole while holding the others fixed by construction.

## 4. The expensive part is eigenvector rotation

The first-order change of eigenvector `phi_n` is

```text
d phi_n / d k_e
  = sum_{m != n}
      [ z_m z_n / (lambda_n-lambda_m) ] phi_m.
```

This is the crucial term.

A local bond does not merely move eigenvalues. It **rotates the modal basis itself**. The rotation is strongest when:

```text
|z_n z_m| is large
and
|lambda_n-lambda_m| is small.
```

Near-degenerate modes therefore produce large coordinate curvature.

Once the eigenvectors rotate, all source and readout residues also change:

```text
b_A,n = phi_n(A)
b_B,n = phi_n(B)
c_n   = phi_n(soma).
```

So one local bond simultaneously changes

```text
modal frequencies
+ source-A residues
+ source-B residues
+ soma/output residues
+ their pairwise interference relationships.
```

That is precisely the bundle of quantities the free `F/A/B/C` parameterization was allowed to address independently.

## 5. Why the local tangent can still be excellent

At one operating point the first derivative of that bundled local perturbation can align almost perfectly with the desired task direction. `SPECTRAL_TO_LOCAL_COMPILER_V01.md` found that one feasible local bond reproduces the 25-lag **shape** of 90.6% of winning free-coordinate tangents to <10% relative error.

There is no contradiction:

```text
at k0:
    bundled local derivative
    can align with
    desired spectral derivative

but after Delta k:
    the eigenvectors and residues have moved
    so the bundle itself changes direction/scale.
```

The local coordinate vector field is state-dependent.

## 6. Why direct spectral coordinates are straighter

The benchmark's direct coordinates deliberately factor the dynamics into near-normal coordinates:

```text
F_n   change one modal stiffness/pole
A_n   change one source-A residue
B_n   change one source-B residue
C_n   change one soma/output residue.
```

They are not perfectly linear coordinates because the output objective is nonlinear and the pole dynamics depend nonlinearly on stiffness. But they suppress the strongest source of geometric curvature: **eigenbasis rotation induced by local operator perturbation**.

That matches the fresh curvature measurement:

```text
median relative finite-step error at delta=.01
local bond       0.7997
free spectral    0.0363.
```

## 7. Connection to the earlier modal-locality result

`MODAL_LOCALITY_V01.md` found that a one-cell structural edit reorganizes a broad fraction of the graph spectrum. That was the finite/discrete version of the same phenomenon.

The formula above supplies the infinitesimal version:

```text
local physical change
   -> dense off-diagonal modal perturbation z z^T
   -> gap-amplified eigenvector rotation
   -> globally changed residues/interference.
```

So the old “one cell changes global modes” observation and the new “local coordinates are ~20x more curved” result are the same spectral geometry viewed at two scales.

## 8. The compiler problem is now mathematically explicit

A physical compiler has two opposing coordinate systems:

### Local hardware chart

```text
{k_e}
```

Advantages:
- one coordinate = one local physical element;
- exact local in-situ gradient is available;
- fabrication/control can be spatially local.

Cost:
- one coordinate perturbs many poles/residues together;
- eigenbasis rotates;
- response is strongly curved/nonmonotone;
- gradients must be refreshed frequently.

### Spectral chart

```text
{lambda_n, b_A,n, b_B,n, c_n}
```

Advantages:
- task-aligned normal coordinates;
- substantially straighter finite-step response;
- larger useful optimization steps.

Cost:
- individual coordinates are globally nonlocal in the physical basis;
- synthesizing them requires a compiler / many local controls / generic programmable mesh.

## Wall sentence

> **A local edge is rank one in space but dense in mode space. Its off-diagonal modal couplings rotate the eigenbasis with strength inversely proportional to spectral gaps. That is why local physical tuning can have an excellent instantaneous task gradient yet become strongly curved after finite movement, while direct modal tuning remains comparatively straight.**
