# Plan: A minimal TinyAD reproduction of the Bending-Active Cosserat (BAC) shell model

Reference paper: Chen, Vouga, Kaufman, *"Better Bending: Analysis, Construction and
Verification of Discrete Bending Models for Kirchhoff–Love Shells,"* SIGGRAPH 2026.

BAC is the midedge-director shell energy whose discrete second fundamental form uses a
`tan(alpha)` hinge measure (per libshell it is `MidedgeAngleTanFormulation`). The
`tan` gives an energy barrier as a hinge approaches a 180° fold — the key property the
paper highlights for under-resolved high-curvature regions.

This repo re-implements BAC from scratch on top of **TinyAD** (autodiff gradients +
sparse Hessians). It does **not** link libshell, libshell's energies, or better-bending.
The upstream repos are cloned into `reference/` (git-ignored) only as a *development
oracle* for validating our energy numerically.

## The model (what we implement)

Degrees of freedom: vertex positions `q ∈ R^{3V}` plus one director angle `θ_e` per edge.

Per face `f` with vertices `q0,q1,q2` and local edges `i=0,1,2` (edge `i` opposite vertex `i`):
- First fundamental form `a` (edge Gram matrix).
- Altitude `h_i` = distance from vertex `i` to opposite edge.
- Dihedral angle `θ_i` at edge `i` (signed angle between the two incident face normals).
- Director angle `α_i = θ_i/2 + orient_i · θ_{edge(i)}`.
- Hinge measure `II_i = 2 h_i tan(α_i)`  ← **the BAC `tan` barrier**.
- Second fundamental form `b = [[II0+II1, II0],[II0, II0+II2]]`.

StVK energy per face with rest forms `abar, bbar`, thickness `t`, Lamé `α,β`:
- `M_s = abar⁻¹(a - abar)`, membrane `= (t/4)·dA·(½α tr(M_s)² + β tr(M_s²))`.
- `M_b = abar⁻¹(b - bbar)`, bending  `= (t³/12)·dA·(½α tr(M_b)² + β tr(M_b²))`.
- `dA = ½√det(abar)`.

2D Lamé from Young `Y`, Poisson `ν`: `α = Yν/(1-ν²)`, `β = Y/(2(1+ν))`.

## Milestones (each has a concrete test)

- **M0 Scaffold + deps.** CMake fetches Eigen + TinyAD + libigl + polyscope.
  *Test:* trivial TinyAD energy builds & runs; polyscope EGL headless writes a PNG.
- **M1 Connectivity.** `MeshConnectivity` (EV/EF/EOpp/FE/FEorient) matching libshell
  conventions so we can cross-check.
  *Test:* adjacency on a 2-triangle mesh matches a libshell dump.
- **M2 Energies (TinyAD).** StVK membrane (9-dof stencil) + BAC bending (21-dof stencil).
  *Tests:* (a) flat config, flat rest → bending E=0, ∇E=0; (b) rigid-motion invariance;
  (c) TinyAD finite-difference gradient/Hessian check; (d) **golden**: match a libshell
  energy+gradient dump on a random crumpled config to ~1e-8.
- **M3 Projected Newton.** Fixed (pinned) DOFs, gravity potential, per-element PSD Hessian
  projection, backtracking line search, edge-DOF validity guard.
  *Test:* a single creased hinge and a cantilever relax; energy decreases monotonically;
  ‖∇E‖ → 0.
- **M4 Rendering.** polyscope headless (EGL) screenshots → PNG frames → GIF.
  *Test:* a PNG of a deformed mesh is produced offscreen.
- **M5 Flagship — lateral buckling (paper §10.3, Romero et al.).** Cantilevered strip
  (L=1, width W, thickness 1e-3), transverse gravity load γ* swept high→low with
  Newton continuation; measure max out-of-plane deflection → buckling bifurcation curve.
  *Test:* flat below a critical γ*, buckles above; curve + rendered shapes.
- **M6 Bonus — dynamic drape.** Inertia (implicit Euler) + projected Newton per step;
  a hanging sheet / sharp fold GIF showcasing the `tan` barrier.
- **M7 Polish.** README with results, figures, build/run instructions, verification notes.

## Layout
```
src/    MeshConnectivity, BAC energy (TinyAD), projected-Newton solver, sim utilities
apps/   lateral_buckling, drape, validate (FD + golden), render helper
tests/  unit tests
docs/   committed result figures/GIFs
reference/  (git-ignored) upstream repos used only as a validation oracle
```
