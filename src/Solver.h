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
    double tol = 1e-6;       // stop when ||reduced gradient|| < tol (absolute)
    double relTol = 1e-8;    // ... or < relTol * ||reduced external force||
    double reg0 = 1e-9;      // initial Tikhonov regularization
    bool projected = true;   // per-element PSD-projected Hessian (else true Hessian)
    bool verbose = false;

    // Optional implicit-Euler inertia term:  + 0.5*coeff*(x-target)^T diag(mass) (x-target)
    const Eigen::VectorXd* inertiaMass = nullptr;    // per-DOF lumped mass (0 for massless)
    const Eigen::VectorXd* inertiaTarget = nullptr;  // x* = x_n + dt*v_n
    double inertiaCoeff = 0.0;                        // 1/dt^2
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

    const bool haveInertia = opt.inertiaCoeff > 0 && opt.inertiaMass && opt.inertiaTarget;
    auto inertiaEnergy = [&](const Eigen::VectorXd& xx) -> double {
        if (!haveInertia) return 0.0;
        Eigen::VectorXd d = xx - *opt.inertiaTarget;
        return 0.5 * opt.inertiaCoeff * (opt.inertiaMass->array() * d.array() * d.array()).sum();
    };
    auto objective = [&](const Eigen::VectorXd& xx, bool& valid) -> double {
        Eigen::MatrixXd V; Eigen::VectorXd th;
        model.unpack(xx, V, th);
        valid = model.dofsValid(V, th);
        if (!valid) return std::numeric_limits<double>::infinity();
        return func.eval(xx) - extForce.dot(xx) + inertiaEnergy(xx);
    };

    const double stopTol = std::max(opt.tol, opt.relTol * (P * extForce).norm());

    NewtonResult res;
    double reg = opt.reg0;
    for (int it = 0; it < opt.maxIters; it++) {
        double f; Eigen::VectorXd g; Eigen::SparseMatrix<double> H;
        if (opt.projected) {
            auto r = func.eval_with_hessian_proj(x);
            f = std::get<0>(r); g = std::get<1>(r); H = std::get<2>(r);
        } else {
            auto r = func.eval_with_derivatives(x);
            f = std::get<0>(r); g = std::get<1>(r); H = std::get<2>(r);
        }
        g -= extForce;                 // gradient of the linear gravity potential
        double E = f - extForce.dot(x);
        if (haveInertia) {
            Eigen::VectorXd d = x - *opt.inertiaTarget;
            g += opt.inertiaCoeff * (opt.inertiaMass->array() * d.array()).matrix();
            E += inertiaEnergy(x);
            // add coeff*mass to the Hessian diagonal
            for (int i = 0; i < n; i++) {
                double m = opt.inertiaCoeff * (*opt.inertiaMass)[i];
                if (m != 0.0) H.coeffRef(i, i) += m;
            }
        }

        Eigen::VectorXd gFree = P * g;
        res.iters = it; res.gradNorm = gFree.norm(); res.energy = E;
        if (opt.verbose && (it % 20 == 0))
            std::printf("  it %4d  E=%.8e  |g|=%.3e  reg=%.1e\n", it, E, gFree.norm(), reg);
        if (gFree.norm() < stopTol) { res.converged = true; break; }

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
                // Iterative refinement: recover precision lost to the huge
                // membrane/bending conditioning (cond ~ (L/h)^2).
                for (int ref = 0; ref < 3; ref++) {
                    Eigen::VectorXd resid = -gFree - M * dFree;
                    dFree += solver.solve(resid);
                }
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
        if (opt.verbose)
            std::printf("  it %4d  E=%.8e |g|=%.3e reg=%.1e t=%.3e |d|=%.2e\n",
                        it, E, gFree.norm(), reg, t, d.norm());
        if (!stepped) { break; }        // stalled
        reg = std::max(opt.reg0, reg * 0.5);
    }

    // final gradient
    {
        auto [f, g] = func.eval_with_gradient(x);
        g -= extForce;
        res.energy = f - extForce.dot(x);
        if (haveInertia) {
            Eigen::VectorXd d = x - *opt.inertiaTarget;
            g += opt.inertiaCoeff * (opt.inertiaMass->array() * d.array()).matrix();
            res.energy += inertiaEnergy(x);
        }
        res.gradNorm = (P * g).norm();
        res.converged = res.gradNorm < stopTol;
    }
    return res;
}

} // namespace bac
