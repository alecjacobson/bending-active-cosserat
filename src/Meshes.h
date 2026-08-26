// Simple structured mesh generators + OBJ writer (no external mesh IO needed).
#pragma once
#include <Eigen/Core>
#include <array>
#include <cstdio>
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

// A curved cylindrical arch: an arc of radius R spanning angle phiSpan in the xz-plane
// (crown up in +z, feet at z=0), extruded along y for length Lv. Developable, so its
// rest metric matches a flat sheet -- a clean non-flat rest (bbar != 0) test shape.
inline void makeCylinderArch(double R, double Lv, int nu, int nv, double phiSpan,
                             Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    int Nu = nu + 1, Nv = nv + 1;
    V.resize(Nu * Nv, 3);
    auto vid = [&](int i, int j) { return j * Nu + i; };
    for (int j = 0; j < Nv; j++)
        for (int i = 0; i < Nu; i++) {
            double u = double(i) / nu;                 // 0..1 across the arc
            double phi = (u - 0.5) * phiSpan;          // centered angle
            double x = R * std::sin(phi), z = R * std::cos(phi);
            double y = Lv * double(j) / nv;
            V.row(vid(i, j)) << x, y, z;
        }
    F.resize(2 * nu * nv, 3);
    int f = 0;
    for (int j = 0; j < nv; j++)
        for (int i = 0; i < nu; i++) {
            int a = vid(i, j), b = vid(i + 1, j), c = vid(i + 1, j + 1), d = vid(i, j + 1);
            F.row(f++) << a, b, c;
            F.row(f++) << a, c, d;
        }
}

inline bool readOBJ(const std::string& path, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    std::ifstream in(path);
    if (!in.good()) return false;
    std::vector<std::array<double, 3>> verts;
    std::vector<std::array<int, 3>> faces;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2) continue;
        if (line[0] == 'v' && line[1] == ' ') {
            double x, y, z; sscanf(line.c_str() + 2, "%lf %lf %lf", &x, &y, &z);
            verts.push_back({x, y, z});
        } else if (line[0] == 'f' && line[1] == ' ') {
            int a, b, c; sscanf(line.c_str() + 2, "%d %d %d", &a, &b, &c);
            faces.push_back({a - 1, b - 1, c - 1});
        }
    }
    V.resize(verts.size(), 3);
    for (size_t i = 0; i < verts.size(); i++) V.row(i) << verts[i][0], verts[i][1], verts[i][2];
    F.resize(faces.size(), 3);
    for (size_t i = 0; i < faces.size(); i++) F.row(i) << faces[i][0], faces[i][1], faces[i][2];
    return true;
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
