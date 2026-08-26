// M3 test: a cantilevered strip relaxes under gravity via projected Newton.
// Checks: solver converges (||grad|| small), total energy decreases, tip deflects.
#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include "BACModel.h"
#include "Solver.h"
#include "Meshes.h"

using namespace bac;

int main() {
    double W = 0.2, L = 1.0, thickness = 2e-3;
    int nx = 6, ny = 30;
    Eigen::MatrixXd V; Eigen::MatrixXi F;
    makeRectangleGrid(W, L, nx, ny, V, F);

    double Y = 1e6, nu = 0.3;
    double lameAlpha = Y * nu / (1 - nu * nu), lameBeta = Y / 2 / (1 + nu);
    BACModel model(V, F, thickness, lameAlpha, lameBeta);

    int nverts = V.rows();
    // Pin the y=0 edge (clamped cantilever root).
    std::vector<char> fixed(model.nDofs(), 0);
    for (int v = 0; v < nverts; v++)
        if (V(v, 1) < 1e-9)
            for (int c = 0; c < 3; c++) fixed[3 * v + c] = 1;

    // Gravity in -z: per-vertex force = -rho*t*g * (vertex area).
    Eigen::VectorXd vArea = Eigen::VectorXd::Zero(nverts);
    for (int f = 0; f < F.rows(); f++) {
        Eigen::Vector3d a = V.row(F(f, 0)), b = V.row(F(f, 1)), c = V.row(F(f, 2));
        double ar = 0.5 * (b - a).cross(c - a).norm();
        for (int k = 0; k < 3; k++) vArea[F(f, k)] += ar / 3.0;
    }
    double rho = 1000.0, g = 9.8;
    Eigen::VectorXd extForce = Eigen::VectorXd::Zero(model.nDofs());
    for (int v = 0; v < nverts; v++)
        extForce[3 * v + 2] = -rho * thickness * g * vArea[v];

    auto func = model.makeEnergyFunction();
    Eigen::VectorXd x = model.pack(V, model.restEdgeThetas());

    bool v0; double E0;
    {
        Eigen::MatrixXd Vv; Eigen::VectorXd th; model.unpack(x, Vv, th);
        E0 = func.eval(x) - extForce.dot(x);
    }

    NewtonOptions opt; opt.tol = 1e-6; opt.maxIters = 2000; opt.verbose = true;
    NewtonResult r = projectedNewton(model, func, x, fixed, extForce, opt);

    Eigen::MatrixXd Vf; Eigen::VectorXd thf; model.unpack(x, Vf, thf);
    double maxDrop = 0, tipZ = 0, maxY = 0;
    for (int v = 0; v < nverts; v++) {
        maxDrop = std::max(maxDrop, -Vf(v, 2));
        if (V(v, 1) > maxY) { maxY = V(v, 1); tipZ = Vf(v, 2); }
    }

    std::printf("\nInitial total energy: %.6e\n", E0);
    std::printf("Final   total energy: %.6e\n", r.energy);
    std::printf("Iterations: %d  converged=%d  |grad|=%.3e\n", r.iters, r.converged, r.gradNorm);
    std::printf("Max downward deflection: %.4f m   tip z: %.4f m\n", maxDrop, tipZ);

    writeOBJ("cantilever_final.obj", Vf, F);

    int fail = 0;
    auto ck = [&](bool ok, const char* n) { std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", n); if (!ok) fail++; };
    ck(r.converged, "converged to ||grad|| < tol");
    ck(r.energy < E0, "total energy decreased");
    ck(maxDrop > 1e-3, "cantilever deflected downward");
    ck(tipZ < -1e-3, "tip dropped below rest plane");
    std::printf("\n%s\n", fail ? "FAILURES" : "ALL PASS");
    return fail ? 1 : 0;
}
