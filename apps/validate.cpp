// Validation harness for the BAC energy.
//   1. Connectivity sanity on a tiny mesh.
//   2. Flat rest + flat config -> bending energy ~ 0, gradient ~ 0.
//   3. Rigid-motion invariance of the energy.
//   4. Finite-difference check of TinyAD gradient (and Hessian) on a random config.
//   5. Golden check vs a libshell dump (if reference/golden.txt exists).
#include <cstdio>
#include <random>
#include <Eigen/Dense>
#include "MeshConnectivity.h"
#include "BACModel.h"
#include "Meshes.h"

using namespace bac;

static int g_fail = 0;
static void check(bool ok, const char* name, double val = 0) {
    std::printf("  [%s] %s (%.3e)\n", ok ? "PASS" : "FAIL", name, val);
    if (!ok) g_fail++;
}

int main() {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> U(-1, 1);

    // ---------------- M1: connectivity on 2 triangles ----------------
    {
        std::printf("== M1 connectivity ==\n");
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0;
        Eigen::MatrixXi F(2, 3);
        F << 0, 1, 2,  0, 2, 3;
        MeshConnectivity mesh(F);
        check(mesh.nEdges() == 5, "nEdges == 5", mesh.nEdges());
        // the shared diagonal edge (0,2) must have two incident faces and two apexes
        int shared = -1;
        for (int e = 0; e < mesh.nEdges(); e++) {
            int a = mesh.edgeVertex(e, 0), b = mesh.edgeVertex(e, 1);
            if ((a == 0 && b == 2)) shared = e;
        }
        check(shared != -1, "found shared diagonal edge");
        if (shared != -1) {
            bool interior = mesh.edgeFace(shared, 0) != -1 && mesh.edgeFace(shared, 1) != -1 &&
                            mesh.edgeOppositeVertex(shared, 0) != -1 &&
                            mesh.edgeOppositeVertex(shared, 1) != -1;
            check(interior, "shared edge is interior with 2 apexes");
            // apexes should be vertices 1 and 3
            int o0 = mesh.edgeOppositeVertex(shared, 0), o1 = mesh.edgeOppositeVertex(shared, 1);
            check((o0 == 1 && o1 == 3) || (o0 == 3 && o1 == 1), "apexes are {1,3}");
        }
    }

    // ---------------- Build a test mesh + model ----------------
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    makeRectangleGrid(0.3, 1.0, 3, 6, V, F);  // small strip
    double Y = 1e5, nu = 0.35, thickness = 1e-2;
    double lameAlpha = Y * nu / (1.0 - nu * nu);
    double lameBeta = Y / 2.0 / (1.0 + nu);
    BACModel model(V, F, thickness, lameAlpha, lameBeta);
    auto func = model.makeEnergyFunction();

    // ---------------- M2a: flat config -> ~0 energy & gradient ----------------
    {
        std::printf("== M2a flat rest, flat config ==\n");
        Eigen::VectorXd x = model.pack(V, model.restEdgeThetas());
        auto [fval, grad] = func.eval_with_gradient(x);
        check(std::abs(fval) < 1e-18, "energy ~ 0", fval);
        check(grad.norm() < 1e-9, "gradient ~ 0", grad.norm());
    }

    // ---------------- M2b: rigid-motion invariance ----------------
    {
        std::printf("== M2b rigid invariance ==\n");
        // deform to a nontrivial shape first
        Eigen::MatrixXd Vd = V;
        for (int i = 0; i < Vd.rows(); i++) Vd(i, 2) += 0.1 * std::sin(3 * Vd(i, 1)) * Vd(i, 0);
        Eigen::VectorXd th = model.restEdgeThetas();
        for (int e = 0; e < th.size(); e++) th[e] = 0.05 * U(rng);
        double E0 = func.eval(model.pack(Vd, th));

        // apply random rotation + translation
        Eigen::Vector3d axis(U(rng), U(rng), U(rng)); axis.normalize();
        double ang = 0.7;
        Eigen::Matrix3d R = Eigen::AngleAxisd(ang, axis).toRotationMatrix();
        Eigen::Vector3d t(0.3, -0.2, 0.5);
        Eigen::MatrixXd Vr = (Vd * R.transpose()).rowwise() + t.transpose();
        double E1 = func.eval(model.pack(Vr, th));
        check(std::abs(E1 - E0) < 1e-6 * (1 + std::abs(E0)), "energy invariant under rigid motion",
              std::abs(E1 - E0));
    }

    // ---------------- M2c: finite-difference gradient & Hessian ----------------
    {
        std::printf("== M2c finite-difference derivative check ==\n");
        Eigen::MatrixXd Vd = V;
        for (int i = 0; i < Vd.rows(); i++) {
            Vd(i, 0) += 0.02 * U(rng);
            Vd(i, 1) += 0.02 * U(rng);
            Vd(i, 2) += 0.05 * U(rng);
        }
        Eigen::VectorXd th = model.restEdgeThetas();
        for (int e = 0; e < th.size(); e++) th[e] = 0.1 * U(rng);
        Eigen::VectorXd x = model.pack(Vd, th);

        auto [fval, grad, Hess] = func.eval_with_derivatives(x);

        // central finite differences of the energy -> gradient
        double h = 1e-6;
        double maxGErr = 0, maxGMag = 0;
        for (int i = 0; i < x.size(); i++) {
            Eigen::VectorXd xp = x, xm = x;
            xp[i] += h; xm[i] -= h;
            double gnum = (func.eval(xp) - func.eval(xm)) / (2 * h);
            maxGErr = std::max(maxGErr, std::abs(gnum - grad[i]));
            maxGMag = std::max(maxGMag, std::abs(grad[i]));
        }
        check(maxGErr < 1e-4 * (1 + maxGMag), "grad matches finite differences", maxGErr);

        // FD of gradient -> Hessian (spot check a few columns)
        double maxHErr = 0, maxHMag = 0;
        Eigen::MatrixXd Hd(Hess);
        for (int i = 0; i < x.size(); i += std::max(1, (int)x.size() / 12)) {
            Eigen::VectorXd xp = x, xm = x;
            xp[i] += h; xm[i] -= h;
            auto [fp, gp] = func.eval_with_gradient(xp);
            auto [fm, gm] = func.eval_with_gradient(xm);
            Eigen::VectorXd col = (gp - gm) / (2 * h);
            for (int j = 0; j < x.size(); j++) {
                maxHErr = std::max(maxHErr, std::abs(col[j] - Hd(j, i)));
                maxHMag = std::max(maxHMag, std::abs(Hd(j, i)));
            }
        }
        check(maxHErr < 1e-3 * (1 + maxHMag), "Hessian matches finite differences", maxHErr);
    }

    // ---------------- M2d: golden check vs libshell ----------------
    {
        std::printf("== M2d golden check (reference/golden.txt) ==\n");
        std::ifstream gf("reference/golden.txt");
        if (!gf.good()) {
            std::printf("  [SKIP] no reference/golden.txt (run reference dumper first)\n");
        } else {
            int gnv, gne, gnf;
            gf >> gnv >> gne >> gnf;
            Eigen::MatrixXd Vrest(gnv, 3);
            for (int i = 0; i < gnv; i++) gf >> Vrest(i, 0) >> Vrest(i, 1) >> Vrest(i, 2);
            Eigen::MatrixXi Fg(gnf, 3);
            for (int i = 0; i < gnf; i++) gf >> Fg(i, 0) >> Fg(i, 1) >> Fg(i, 2);
            Eigen::MatrixXd Vcur(gnv, 3);
            for (int i = 0; i < gnv; i++) gf >> Vcur(i, 0) >> Vcur(i, 1) >> Vcur(i, 2);
            Eigen::VectorXd thg(gne);
            for (int i = 0; i < gne; i++) gf >> thg[i];
            double gy, gnu, gth, genergy;
            gf >> gy >> gnu >> gth >> genergy;
            Eigen::VectorXd ggrad(3 * gnv + gne);
            for (int i = 0; i < ggrad.size(); i++) gf >> ggrad[i];

            double la = gy * gnu / (1 - gnu * gnu), lb = gy / 2 / (1 + gnu);
            BACModel gm(Vrest, Fg, gth, la, lb);       // rest state from Vrest
            auto gfunc = gm.makeEnergyFunction();
            Eigen::VectorXd xg = gm.pack(Vcur, thg);   // evaluate at deformed Vcur
            auto [fv, gr] = gfunc.eval_with_gradient(xg);
            check(std::abs(fv - genergy) < 1e-7 * (1 + std::abs(genergy)), "energy matches libshell",
                  std::abs(fv - genergy));
            // libshell gradient ordering: positions (row-major 3*nv) then edge thetas -> same as ours
            check((gr - ggrad).norm() < 1e-6 * (1 + ggrad.norm()), "gradient matches libshell",
                  (gr - ggrad).norm());
        }
    }

    std::printf("\n%s (%d failures)\n", g_fail == 0 ? "ALL PASS" : "FAILURES", g_fail);
    return g_fail == 0 ? 0 : 1;
}
