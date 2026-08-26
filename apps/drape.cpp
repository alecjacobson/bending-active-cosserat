// Bonus dynamic demo: a square sheet pinned at two top corners is released from
// horizontal and drapes under gravity. Backward-Euler time stepping, each step a
// projected-Newton solve of the incremental potential
//     E_elastic(x) - g.x + (1/2 dt^2) (x - x*)^T M (x - x*),   x* = x_n + dt v_n.
// Uses StVK membrane + BAC bending + lumped inertia. Renders frames via headless EGL.
//
// Usage: drape [nx] [steps] [dt]
#include <cstdio>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "BACModel.h"
#include "Solver.h"
#include "Meshes.h"
#include "Render.h"

using namespace bac;

int main(int argc, char** argv) {
    int n = (argc > 1) ? std::atoi(argv[1]) : 24;
    int steps = (argc > 2) ? std::atoi(argv[2]) : 140;
    double dt = (argc > 3) ? std::atof(argv[3]) : 0.02;

    double W = 1.0, L = 1.0, thickness = 5e-3;
    double Y = 2e5, nu = 0.3;
    double lameAlpha = Y * nu / (1 - nu * nu), lameBeta = Y / 2 / (1 + nu);
    double rho = 200.0, g = 9.8;   // areal mass rho*thickness = 1 kg/m^2

    Eigen::MatrixXd V; Eigen::MatrixXi F;
    makeRectangleGrid(W, L, n, n, V, F);
    int nverts = V.rows();
    BACModel model(V, F, thickness, lameAlpha, lameBeta);
    model.setThreadCount(8);
    auto func = model.makeEnergyFunction();

    // Pin the two top corners (y=L, x=0 and x=W).
    std::vector<char> fixed(model.nDofs(), 0);
    auto pinCorner = [&](double tx, double ty) {
        int best = 0; double bd = 1e18;
        for (int v = 0; v < nverts; v++) {
            double d = std::hypot(V(v, 0) - tx, V(v, 1) - ty);
            if (d < bd) { bd = d; best = v; }
        }
        for (int c = 0; c < 3; c++) fixed[3 * best + c] = 1;
    };
    pinCorner(0, L); pinCorner(W, L);

    // Lumped vertex mass -> per-DOF mass vector (0 on director DOFs).
    Eigen::VectorXd vArea = Eigen::VectorXd::Zero(nverts);
    for (int f = 0; f < F.rows(); f++) {
        Eigen::Vector3d a = V.row(F(f, 0)), b = V.row(F(f, 1)), c = V.row(F(f, 2));
        double ar = 0.5 * (b - a).cross(c - a).norm();
        for (int k = 0; k < 3; k++) vArea[F(f, k)] += ar / 3.0;
    }
    Eigen::VectorXd massDiag = Eigen::VectorXd::Zero(model.nDofs());
    Eigen::VectorXd extForce = Eigen::VectorXd::Zero(model.nDofs());
    for (int v = 0; v < nverts; v++) {
        double m = rho * thickness * vArea[v];
        for (int c = 0; c < 3; c++) massDiag[3 * v + c] = m;
        extForce[3 * v + 2] = -g * m;   // gravity in -z
    }

    Eigen::VectorXd x = model.pack(V, model.restEdgeThetas());
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nDofs());

    renderInit(1280, 960);
    Camera cam;
    cam.eye = {2.4, -1.9, 0.15}; cam.target = {0.5, 0.75, -0.45}; cam.up = {0, 0, 1};
    system("mkdir -p results/drape");

    int frame = 0;
    for (int s = 0; s < steps; s++) {
        Eigen::VectorXd xPrev = x;
        Eigen::VectorXd xStar = x + dt * v;

        NewtonOptions opt;
        opt.projected = true;              // projected Newton (well-conditioned with inertia)
        opt.tol = 1e-7; opt.maxIters = 200;
        opt.inertiaMass = &massDiag; opt.inertiaTarget = &xStar; opt.inertiaCoeff = 1.0 / (dt * dt);
        NewtonResult r = projectedNewton(model, func, x, fixed, extForce, opt);

        v = (x - xPrev) / dt;

        if (s % 2 == 0) {
            Eigen::MatrixXd Vc; Eigen::VectorXd th; model.unpack(x, Vc, th);
            char fn[128]; snprintf(fn, sizeof(fn), "results/drape/f%04d.png", frame++);
            renderMeshPNG(Vc, F, fn, cam, 2);
            printf("step %3d t=%.2f  it=%d conv=%d  Etot=%.4e\n", s, s * dt, r.iters,
                   r.converged, r.energy);
        }
    }
    printf("done: %d frames in results/drape/\n", frame);
    return 0;
}
