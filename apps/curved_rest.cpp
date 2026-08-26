// Non-flat rest-state example: a curved cylindrical arch (barrel vault).
//
// The rest configuration is a half-cylinder (rest second fundamental form b̄ != 0),
// set via BACModel::setCurvedRest. We first verify that the elastic energy at the
// exact rest configuration is ~0 (the key correctness check for a curved rest), then
// clamp the two feet and release the arch under gravity: it sags to a new equilibrium
// while "remembering" its curved rest shape. Backward-Euler dynamics; writes OBJ frames.
//
// This mirrors the paper's curved test geometries (e.g. the HalfCylinder benchmark).
//
// Usage: curved_rest [nu] [nv] [steps] [dt]
#include <cstdio>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "BACModel.h"
#include "Solver.h"
#include "Meshes.h"

using namespace bac;

int main(int argc, char** argv) {
    int nu = (argc > 1) ? std::atoi(argv[1]) : 24;
    int nv = (argc > 2) ? std::atoi(argv[2]) : 36;
    int steps = (argc > 3) ? std::atoi(argv[3]) : 150;
    double dt = (argc > 4) ? std::atof(argv[4]) : 0.02;

    const double R = 0.5, Lv = 1.5, phiSpan = M_PI;   // half-cylinder arch
    const double thickness = 8e-3, Y = 2e6, nu_p = 0.3;
    double lameAlpha = Y * nu_p / (1 - nu_p * nu_p), lameBeta = Y / 2 / (1 + nu_p);
    double rho = 45.0, g = 9.8;                         // gentle load: sags, no snap-through

    Eigen::MatrixXd V; Eigen::MatrixXi F;
    makeCylinderArch(R, Lv, nu, nv, phiSpan, V, F);
    int nverts = V.rows();

    BACModel model(V, F, thickness, lameAlpha, lameBeta);
    model.setThreadCount(8);
    // Define the CURVED rest state (rest ā, b̄ from the arch; rest directors = 0).
    Eigen::VectorXd restThetas = Eigen::VectorXd::Zero(model.nEdges());
    model.setCurvedRest(V, restThetas);

    auto func = model.makeEnergyFunction();

    // --- key check: energy at the exact curved rest configuration must be ~0 ---
    Eigen::VectorXd xRest = model.pack(V, restThetas);
    double Erest = func.eval(xRest);
    auto [fr, gr] = func.eval_with_gradient(xRest);
    std::printf("Curved-rest sanity: E(rest) = %.3e   |grad(rest)| = %.3e  (both ~0 expected)\n",
                Erest, gr.norm());

    // Clamp the two feet (the u=0 and u=nu columns, at z=0).
    std::vector<char> fixed(model.nDofs(), 0);
    int Nu = nu + 1;
    for (int j = 0; j <= nv; j++) {
        for (int i : {0, nu}) {
            int v = j * Nu + i;
            for (int c = 0; c < 3; c++) fixed[3 * v + c] = 1;
        }
    }

    // Lumped mass + gravity (-z).
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
        extForce[3 * v + 2] = -g * m;
    }

    Eigen::VectorXd x = xRest;
    Eigen::VectorXd vel = Eigen::VectorXd::Zero(model.nDofs());
    system("mkdir -p results/arch");
    int frame = 0;
    for (int s = 0; s < steps; s++) {
        Eigen::VectorXd xPrev = x, xStar = x + dt * vel;
        NewtonOptions opt;
        opt.projected = true; opt.tol = 1e-7; opt.maxIters = 200;
        opt.inertiaMass = &massDiag; opt.inertiaTarget = &xStar; opt.inertiaCoeff = 1.0 / (dt * dt);
        NewtonResult r = projectedNewton(model, func, x, fixed, extForce, opt);
        vel = (x - xPrev) / dt;
        if (s % 2 == 0) {
            Eigen::MatrixXd Vc; Eigen::VectorXd th; model.unpack(x, Vc, th);
            char fn[128]; snprintf(fn, sizeof(fn), "results/arch/f%04d.obj", frame++);
            writeOBJ(fn, Vc, F);
            printf("step %3d t=%.2f  it=%d conv=%d  Etot=%.4e\n", s, s * dt, r.iters,
                   r.converged, r.energy);
        }
    }
    printf("done: %d OBJ frames in results/arch/\n", frame);
    return 0;
}
