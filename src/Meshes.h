// Simple structured mesh generators + OBJ writer (no external mesh IO needed).
#pragma once
#include <Eigen/Core>
#include <fstream>
#include <string>
#include <vector>

namespace bac {

// A W x L rectangle in the xy-plane, (nx+1)x(ny+1) vertices, 2*nx*ny triangles.
// x in [0,W] (width), y in [0,L] (length). z = 0.
inline void makeRectangleGrid(double W, double L, int nx, int ny,
                              Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    int Nx = nx + 1, Ny = ny + 1;
    V.resize(Nx * Ny, 3);
    auto vid = [&](int i, int j) { return j * Nx + i; };
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            V.row(vid(i, j)) << W * double(i) / nx, L * double(j) / ny, 0.0;
        }
    F.resize(2 * nx * ny, 3);
    int f = 0;
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            int a = vid(i, j), b = vid(i + 1, j), c = vid(i + 1, j + 1), d = vid(i, j + 1);
            // consistent CCW winding (normals +z)
            F.row(f++) << a, b, c;
            F.row(f++) << a, c, d;
        }
}

inline void writeOBJ(const std::string& path, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    std::ofstream o(path);
    o.setf(std::ios::fixed);
    o.precision(8);
    for (int i = 0; i < V.rows(); i++)
        o << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    for (int i = 0; i < F.rows(); i++)
        o << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << "\n";
}

} // namespace bac
