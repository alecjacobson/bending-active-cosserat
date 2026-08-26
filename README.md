# Bending-Active Cosserat (BAC) — a minimal TinyAD reproduction

A from-scratch, minimal reproduction of the **Bending-Active Cosserat (BAC)** thin-shell
bending model from:

> Zhen Chen, Etienne Vouga, Danny M. Kaufman.
> *Better Bending: Analysis, Construction and Verification of Discrete Bending Models for
> Kirchhoff–Love Shells.* SIGGRAPH 2026.

BAC is the midedge-director discretization whose second fundamental form uses a
`tan(α)` hinge measure, giving an energy barrier as a hinge folds toward 180°.

This implementation uses **[TinyAD](https://github.com/patr-schm/TinyAD)** for gradients
and (projected) sparse Hessians, **libigl** for mesh IO, and **polyscope** for headless
rendering. It does **not** depend on libshell / better-bending / any physics-energy
library — those are used only as a numerical validation oracle during development.

See [`PLAN.md`](PLAN.md) for the model description and milestone/testing plan.

## Status

Work in progress — see PLAN.md milestones.
