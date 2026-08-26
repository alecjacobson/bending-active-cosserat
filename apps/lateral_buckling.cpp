// Flagship reproduction: lateral (lateral-torsional) buckling of a cantilevered
// strip under gravity, following "Better Bending" Section 10.3 (Romero et al.).
//
// A thin strip of length L (along y, clamped at y=0), depth W (along x), thickness
// h (along z) is loaded by a body force of dimensionless magnitude gamma* along +x
// (its stiff, in-plane direction). Below a critical gamma* it stays planar; above it,
// it buckles out of plane (z) via twisting. We sweep gamma* from high to low with
// Newton continuation and record the maximum out-of-plane deflection max|z|, tracing
// the buckling bifurcation branch.
//
// Usage: lateral_buckling [ny] [gammaMax] [gammaMin] [dGamma] [W1 W2 ...]
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "BACModel.h"
#include "Solver.h"
#include "Meshes.h"

using namespace bac;

int main(int argc, char** argv) {
    // ---- physical parameters (match the paper) ----
    const double nu = 0.35;
    const double L = 1.0;
    const double D = 1.0;         // bending stiffness
    const double h = 1e-3;        // thickness
    const double density = 1.0;   // divides out of gamma*
    const double Y = D * 12.0 * (1.0 - nu * nu) / (h * h * h);
    const double lameAlpha = Y * nu / (1.0 - nu * nu);
    const double lameBeta = Y / 2.0 / (1.0 + nu);

    // ---- sweep configuration (overridable for quick runs) ----
    int ny = (argc > 1) ? std::atoi(argv[1]) : 48;
    double gammaMax = (argc > 2) ? std::atof(argv[2]) : 40.0;
    double gammaMin = (argc > 3) ? std::atof(argv[3]) : 8.0;
    double dGamma = (argc > 4) ? std::atof(argv[4]) : 1.0;
    std::vector<double> widths;
    for (int i = 5; i < argc; i++) widths.push_back(std::atof(argv[i]));
    if (widths.empty()) widths = {0.1, 0.2, 0.3};

    system("mkdir -p results");
    FILE* csv = fopen("results/buckling.csv", "w");
    fprintf(csv, "width,gamma,max_lateral,energy,iters,converged\n");

    for (double W : widths) {
        int nx = std::max(2, (int)std::llround(W * ny / L));
        Eigen::MatrixXd V; Eigen::MatrixXi F;
        makeRectangleGrid(W, L, nx, ny, V, F);
        int nverts = V.rows();
        printf("\n=== W=%.3f  mesh %dx%d  (%d verts, %ld faces) ===\n", W, nx, ny, nverts,
               (long)F.rows());

        BACModel model(V, F, h, lameAlpha, lameBeta);
        model.setThreadCount(8);
        auto func = model.makeEnergyFunction();

        double tol = 1.5 * L / ny;
        std::vector<char> fixed(model.nDofs(), 0);
        for (int v = 0; v < nverts; v++)
            if (V(v, 1) < tol)
                for (int c = 0; c < 3; c++) fixed[3 * v + c] = 1;

        Eigen::VectorXd vArea = Eigen::VectorXd::Zero(nverts);
        for (int f = 0; f < F.rows(); f++) {
            Eigen::Vector3d a = V.row(F(f, 0)), b = V.row(F(f, 1)), c = V.row(F(f, 2));
            double ar = 0.5 * (b - a).cross(c - a).norm();
            for (int k = 0; k < 3; k++) vArea[F(f, k)] += ar / 3.0;
        }
        for (int v = 0; v < nverts; v++) if (fixed[3 * v]) vArea[v] = 0;

        // seed the buckling mode with a tiny out-of-plane perturbation
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> U(0, 1);
        Eigen::MatrixXd curV = V;
        for (int v = 0; v < nverts; v++)
            if (!fixed[3 * v]) curV(v, 2) += 1e-4 * U(rng);

        Eigen::VectorXd x = model.pack(curV, model.restEdgeThetas());

        char dir[128]; snprintf(dir, sizeof(dir), "results/W%03d", (int)std::llround(W * 100));
        { std::string cmd = std::string("mkdir -p ") + dir; system(cmd.c_str()); }

        for (double gamma = gammaMax; gamma >= gammaMin - 1e-9; gamma -= dGamma) {
            double g = D * gamma / (density * L * L * L * h);
            Eigen::VectorXd extForce = Eigen::VectorXd::Zero(model.nDofs());
            for (int v = 0; v < nverts; v++)
                extForce[3 * v + 0] = g * h * vArea[v] * density;  // +x load

            NewtonOptions opt;
            opt.tol = 1e-6;
            opt.maxIters = 1500;
            // True-Hessian Newton with adaptive PD regularization: efficient on the
            // stiff, near-inextensible membrane. (The per-element PSD-projected mode,
            // opt.projected=true, reaches the same equilibria but needs more iters.)
            opt.projected = false;
            opt.verbose = false;
            NewtonResult r = projectedNewton(model, func, x, fixed, extForce, opt);

            Eigen::MatrixXd Vf; Eigen::VectorXd thf; model.unpack(x, Vf, thf);
            double maxLat = 0;
            for (int v = 0; v < nverts; v++) maxLat = std::max(maxLat, std::abs(Vf(v, 2)));

            printf("  gamma=%5.1f  max|z|=%.5e  E=%.4e  it=%d conv=%d\n",
                   gamma, maxLat, r.energy, r.iters, r.converged);
            fprintf(csv, "%.4f,%.4f,%.8e,%.8e,%d,%d\n", W, gamma, maxLat, r.energy,
                    r.iters, r.converged);
            fflush(csv);

            char fn[256]; snprintf(fn, sizeof(fn), "%s/g%05.1f.obj", dir, gamma);
            writeOBJ(fn, Vf, F);
        }
    }
    fclose(csv);
    printf("\nWrote results/buckling.csv and per-width OBJ frames.\n");
    return 0;
}
