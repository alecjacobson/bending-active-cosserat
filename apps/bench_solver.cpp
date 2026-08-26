// Quick benchmark: projected vs true-Hessian Newton on one buckling solve.
#include <cstdio>
#include <chrono>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "BACModel.h"
#include "Solver.h"
#include "Meshes.h"
using namespace bac;

int main(int argc, char** argv) {
    int ny = argc > 1 ? atoi(argv[1]) : 24;
    double W = argc > 2 ? atof(argv[2]) : 0.2;
    double gamma = argc > 3 ? atof(argv[3]) : 40.0;
    int projected = argc > 4 ? atoi(argv[4]) : 1;

    double nu = 0.35, L = 1, D = 1, h = 1e-3, density = 1;
    double Y = D * 12 * (1 - nu * nu) / (h * h * h);
    double la = Y * nu / (1 - nu * nu), lb = Y / 2 / (1 + nu);
    int nx = std::max(2, (int)llround(W * ny / L));
    Eigen::MatrixXd V; Eigen::MatrixXi F;
    makeRectangleGrid(W, L, nx, ny, V, F);
    int nverts = V.rows();
    BACModel model(V, F, h, la, lb); model.setThreadCount(8);
    auto func = model.makeEnergyFunction();
    double tol = 1.5 * L / ny;
    std::vector<char> fixed(model.nDofs(), 0);
    for (int v = 0; v < nverts; v++) if (V(v, 1) < tol) for (int c = 0; c < 3; c++) fixed[3 * v + c] = 1;
    Eigen::VectorXd vArea = Eigen::VectorXd::Zero(nverts);
    for (int f = 0; f < F.rows(); f++) {
        Eigen::Vector3d a = V.row(F(f, 0)), b = V.row(F(f, 1)), c = V.row(F(f, 2));
        double ar = 0.5 * (b - a).cross(c - a).norm();
        for (int k = 0; k < 3; k++) vArea[F(f, k)] += ar / 3.0;
    }
    for (int v = 0; v < nverts; v++) if (fixed[3 * v]) vArea[v] = 0;
    std::mt19937 rng(7); std::uniform_real_distribution<double> U(0, 1);
    Eigen::MatrixXd curV = V;
    for (int v = 0; v < nverts; v++) if (!fixed[3 * v]) curV(v, 2) += 1e-4 * U(rng);
    Eigen::VectorXd x = model.pack(curV, model.restEdgeThetas());
    double g = D * gamma / (density * L * L * L * h);
    Eigen::VectorXd extForce = Eigen::VectorXd::Zero(model.nDofs());
    for (int v = 0; v < nverts; v++) extForce[3 * v + 0] = g * h * vArea[v] * density;

    NewtonOptions opt; opt.tol = 1e-6; opt.maxIters = 2000; opt.projected = projected;
    auto t0 = std::chrono::high_resolution_clock::now();
    NewtonResult r = projectedNewton(model, func, x, fixed, extForce, opt);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    Eigen::MatrixXd Vf; Eigen::VectorXd th; model.unpack(x, Vf, th);
    double maxLat = 0; for (int v = 0; v < nverts; v++) maxLat = std::max(maxLat, std::abs(Vf(v, 2)));
    printf("projected=%d  iters=%d conv=%d  |g|=%.2e  E=%.4e  max|z|=%.4e  time=%.0fms\n",
           projected, r.iters, r.converged, r.gradNorm, r.energy, maxLat, ms);
    return 0;
}
