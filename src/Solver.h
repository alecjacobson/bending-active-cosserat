// Projected-Newton static equilibrium solver for the BAC shell.
//
// Minimizes  E_total(x) = E_elastic(x) - extForce . x  over the free DOFs
// (fixed DOFs are held via a selection matrix), using TinyAD's per-element
// PSD-projected Hessian, with backtracking line search and a validity guard
// (rejects steps that push a hinge past the tan barrier).
#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <cstdio>
#include <vector>
#include "BACModel.h"

namespace bac {

struct NewtonOptions {
    int maxIters = 2000;
    double tol = 1e-6;       // stop when ||reduced gradient|| < tol
    double reg0 = 1e-9;      // initial Tikhonov regularization
    bool verbose = false;
};

struct NewtonResult {
    int iters = 0;
    double gradNorm = 0;
    double energy = 0;
    bool converged = false;
};

// fixed: length nDofs, nonzero => DOF held fixed. extForce: length nDofs (constant).
template <class Func>
NewtonResult projectedNewton(const BACModel& model, Func& func,
                             Eigen::VectorXd& x,
                             const std::vector<char>& fixed,
                             const Eigen::VectorXd& extForce,
                             const NewtonOptions& opt = {}) {
    const int n = (int)x.size();

    // Selection matrix P (nFree x n) mapping free DOFs -> full.
    std::vector<int> freeIdx;
    freeIdx.reserve(n);
    for (int i = 0; i < n; i++) if (!fixed[i]) freeIdx.push_back(i);
    const int nf = (int)freeIdx.size();
    std::vector<Eigen::Triplet<double>> Pc;
    Pc.reserve(nf);
    for (int i = 0; i < nf; i++) Pc.emplace_back(i, freeIdx[i], 1.0);
    Eigen::SparseMatrix<double> P(nf, n);
    P.setFromTriplets(Pc.begin(), Pc.end());
    Eigen::SparseMatrix<double> PT = P.transpose();

    auto objective = [&](const Eigen::VectorXd& xx, bool& valid) -> double {
        Eigen::MatrixXd V; Eigen::VectorXd th;
        model.unpack(xx, V, th);
        valid = model.dofsValid(V, th);
        if (!valid) return std::numeric_limits<double>::infinity();
        return func.eval(xx) - extForce.dot(xx);
    };

    NewtonResult res;
    double reg = opt.reg0;
    for (int it = 0; it < opt.maxIters; it++) {
        auto [f, g, H] = func.eval_with_hessian_proj(x);
        g -= extForce;                 // gradient of the linear gravity potential
        double E = f - extForce.dot(x);

        Eigen::VectorXd gFree = P * g;
        res.iters = it; res.gradNorm = gFree.norm(); res.energy = E;
        if (opt.verbose && (it % 20 == 0))
            std::printf("  it %4d  E=%.8e  |g|=%.3e  reg=%.1e\n", it, E, gFree.norm(), reg);
        if (gFree.norm() < opt.tol) { res.converged = true; break; }

        Eigen::SparseMatrix<double> Hff = (P * H * PT).pruned();
        Eigen::SparseMatrix<double> I(nf, nf); I.setIdentity();

        // Solve (Hff + reg I) d = -gFree, growing reg until SPD factorization succeeds
        // and yields a descent direction.
        Eigen::VectorXd dFree;
        bool ok = false;
        for (int tries = 0; tries < 40; tries++) {
            Eigen::SparseMatrix<double> M = Hff + reg * I;
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(M);
            if (solver.info() == Eigen::Success) {
                dFree = solver.solve(-gFree);
                if (solver.info() == Eigen::Success && dFree.dot(gFree) < 0) { ok = true; break; }
            }
            reg *= 4.0;
        }
        if (!ok) { break; }

        Eigen::VectorXd d = PT * dFree;

        // Backtracking Armijo line search with validity guard.
        double gd = g.dot(d);          // directional derivative (< 0)
        double t = 1.0;
        bool stepped = false;
        for (int ls = 0; ls < 60; ls++) {
            Eigen::VectorXd xn = x + t * d;
            bool valid;
            double En = objective(xn, valid);
            if (valid && En <= E + 1e-4 * t * gd) { x = xn; stepped = true; break; }
            t *= 0.5;
        }
        if (!stepped) { break; }        // stalled
        reg = std::max(opt.reg0, reg * 0.5);
    }

    // final gradient
    {
        auto [f, g] = func.eval_with_gradient(x);
        g -= extForce;
        res.energy = f - extForce.dot(x);
        res.gradNorm = (P * g).norm();
        res.converged = res.gradNorm < opt.tol;
    }
    return res;
}

} // namespace bac
