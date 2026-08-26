// Bending-Active Cosserat (BAC) shell energy, implemented on TinyAD.
//
// DOF layout (one flat vector x of length 3*nVerts + nEdges):
//   x[3*v + c]              : component c of vertex v's position
//   x[3*nVerts + e]         : director angle theta_e on edge e
//
// Energy = StVK membrane (stretching) + BAC bending, summed over faces.
// The BAC bending uses the tan(alpha) hinge measure of the paper's
// MidedgeAngleTanFormulation: II_i = 2 * altitude_i * tan(alpha_i),
// alpha_i = theta_i/2 + orient_i * theta_edge(i), giving an energy barrier as a
// hinge folds toward 180 degrees.
#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <TinyAD/Scalar.hh>
#include <TinyAD/ScalarFunction.hh>
#include "MeshConnectivity.h"

namespace bac {

struct RestState {
    std::vector<double> thickness;    // per face
    std::vector<double> lameAlpha;    // per face
    std::vector<double> lameBeta;     // per face
    std::vector<Eigen::Matrix2d> abar; // per face rest first fundamental form
    std::vector<Eigen::Matrix2d> bbar; // per face rest second fundamental form
};

// ---- small templated geometry helpers (work for double and TinyAD::Scalar) ----

template <class T>
inline Eigen::Matrix<T, 3, 1> cross3(const Eigen::Matrix<T, 3, 1>& a,
                                     const Eigen::Matrix<T, 3, 1>& b) {
    return a.cross(b);
}

// StVK quadratic form: 0.5*alpha*tr(M)^2 + beta*tr(M^2)
template <class T>
inline T stvk(const Eigen::Matrix<T, 2, 2>& M, double lameAlpha, double lameBeta) {
    T tr = M(0, 0) + M(1, 1);
    Eigen::Matrix<T, 2, 2> M2 = M * M;
    T trM2 = M2(0, 0) + M2(1, 1);
    return T(0.5) * lameAlpha * tr * tr + lameBeta * trM2;
}

template <class T>
inline Eigen::Matrix<T, 2, 2> castT(const Eigen::Matrix2d& A) {
    Eigen::Matrix<T, 2, 2> R;
    R << T(A(0, 0)), T(A(0, 1)), T(A(1, 0)), T(A(1, 1));
    return R;
}

class BACModel {
public:
    // Flat rest state: rest first fundamental forms from V, rest bbar = 0.
    BACModel(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
             double thickness, double lameAlpha, double lameBeta)
        : mesh_(F), V0_(V), F_(F) {
        nVerts_ = (int)V.rows();
        nEdges_ = mesh_.nEdges();
        int nf = (int)F.rows();
        rest_.thickness.assign(nf, thickness);
        rest_.lameAlpha.assign(nf, lameAlpha);
        rest_.lameBeta.assign(nf, lameBeta);
        rest_.abar.resize(nf);
        rest_.bbar.assign(nf, Eigen::Matrix2d::Zero());
        for (int f = 0; f < nf; f++)
            rest_.abar[f] = faceFirstFundamentalForm(V, F, f);
        restEdgeThetas_ = Eigen::VectorXd::Zero(nEdges_);
    }

    const MeshConnectivity& mesh() const { return mesh_; }
    RestState& rest() { return rest_; }
    const RestState& rest() const { return rest_; }
    int nVerts() const { return nVerts_; }
    int nEdges() const { return nEdges_; }
    int nDofs() const { return 3 * nVerts_ + nEdges_; }
    const Eigen::VectorXd& restEdgeThetas() const { return restEdgeThetas_; }

    // Pack/unpack the flat DOF vector.
    Eigen::VectorXd pack(const Eigen::MatrixXd& V, const Eigen::VectorXd& thetas) const {
        Eigen::VectorXd x(nDofs());
        for (int v = 0; v < nVerts_; v++)
            for (int c = 0; c < 3; c++) x[3 * v + c] = V(v, c);
        x.segment(3 * nVerts_, nEdges_) = thetas;
        return x;
    }
    void unpack(const Eigen::VectorXd& x, Eigen::MatrixXd& V, Eigen::VectorXd& thetas) const {
        V.resize(nVerts_, 3);
        for (int v = 0; v < nVerts_; v++)
            for (int c = 0; c < 3; c++) V(v, c) = x[3 * v + c];
        thetas = x.segment(3 * nVerts_, nEdges_);
    }

    // First fundamental form (edge Gram matrix) of face f from positions V.
    static Eigen::Matrix2d faceFirstFundamentalForm(const Eigen::MatrixXd& V,
                                                    const Eigen::MatrixXi& F, int f) {
        Eigen::Vector3d q0 = V.row(F(f, 0)), q1 = V.row(F(f, 1)), q2 = V.row(F(f, 2));
        Eigen::Vector3d e1 = q1 - q0, e2 = q2 - q0;
        Eigen::Matrix2d a;
        a << e1.dot(e1), e1.dot(e2), e2.dot(e1), e2.dot(e2);
        return a;
    }

    // Per-face BAC second-fundamental-form entries (II0,II1,II2), double precision.
    // Mirrors the TinyAD bending element exactly (validated by the "rest energy = 0"
    // check when used to define a curved rest state).
    Eigen::Vector3d faceIIEntries(const Eigen::MatrixXd& V, const Eigen::VectorXd& thetas,
                                  int f) const {
        Eigen::Vector3d qf[3];
        for (int k = 0; k < 3; k++) qf[k] = V.row(mesh_.faceVertex(f, k));
        Eigen::Vector3d II;
        for (int i = 0; i < 3; i++) {
            const Eigen::Vector3d& a0 = qf[i];
            const Eigen::Vector3d& a1 = qf[(i + 1) % 3];
            const Eigen::Vector3d& a2 = qf[(i + 2) % 3];
            double altitude = (a1 - a0).cross(a2 - a0).norm() / (a2 - a1).norm();
            int edge = mesh_.faceEdge(f, i);
            double orient = (mesh_.faceEdgeOrientation(f, i) == 0) ? 1.0 : -1.0;
            double halfTheta = 0.0;
            int v2 = mesh_.edgeOppositeVertex(edge, 0);
            int v3 = mesh_.edgeOppositeVertex(edge, 1);
            if (v2 != -1 && v3 != -1) {
                Eigen::Vector3d Q0 = V.row(mesh_.edgeVertex(edge, 0));
                Eigen::Vector3d Q1 = V.row(mesh_.edgeVertex(edge, 1));
                Eigen::Vector3d Q2 = V.row(v2), Q3 = V.row(v3);
                Eigen::Vector3d n0 = (Q0 - Q2).cross(Q1 - Q2);
                Eigen::Vector3d n1 = (Q1 - Q3).cross(Q0 - Q3);
                Eigen::Vector3d axis = Q1 - Q0;
                double s = n0.cross(n1).dot(axis) / axis.norm();
                double c = n0.dot(n1) + n0.norm() * n1.norm();
                halfTheta = std::atan2(s, c);
            }
            double alpha = halfTheta + orient * thetas[edge];
            II[i] = 2.0 * altitude * std::tan(alpha);
        }
        return II;
    }

    Eigen::Matrix2d faceSecondFundamentalForm(const Eigen::MatrixXd& V,
                                              const Eigen::VectorXd& thetas, int f) const {
        Eigen::Vector3d II = faceIIEntries(V, thetas, f);
        Eigen::Matrix2d b;
        b << II[0] + II[1], II[0], II[0], II[0] + II[2];
        return b;
    }

    // Define a curved (non-flat) rest state from a rest configuration: rest first and
    // second fundamental forms are taken from (Vrest, thetasRest). Evaluating the energy
    // at exactly (Vrest, thetasRest) then gives ~0.
    void setCurvedRest(const Eigen::MatrixXd& Vrest, const Eigen::VectorXd& thetasRest) {
        for (int f = 0; f < mesh_.nFaces(); f++) {
            rest_.abar[f] = faceFirstFundamentalForm(Vrest, F_, f);
            rest_.bbar[f] = faceSecondFundamentalForm(Vrest, thetasRest, f);
        }
        restEdgeThetas_ = thetasRest;
    }

    // Validity of director DOFs: every hinge angle alpha must stay strictly inside
    // (-pi/2, pi/2) so the tan() hinge measure stays finite (the BAC barrier).
    bool dofsValid(const Eigen::MatrixXd& V, const Eigen::VectorXd& thetas,
                   double margin = 1e-6) const {
        const double lim = M_PI / 2.0 - margin;
        for (int f = 0; f < mesh_.nFaces(); f++) {
            Eigen::Vector3d qf[3];
            for (int k = 0; k < 3; k++) qf[k] = V.row(mesh_.faceVertex(f, k));
            for (int i = 0; i < 3; i++) {
                int edge = mesh_.faceEdge(f, i);
                double orient = (mesh_.faceEdgeOrientation(f, i) == 0) ? 1.0 : -1.0;
                double halfTheta = 0.0;
                int v2 = mesh_.edgeOppositeVertex(edge, 0);
                int v3 = mesh_.edgeOppositeVertex(edge, 1);
                if (v2 != -1 && v3 != -1) {
                    Eigen::Vector3d Q0 = V.row(mesh_.edgeVertex(edge, 0));
                    Eigen::Vector3d Q1 = V.row(mesh_.edgeVertex(edge, 1));
                    Eigen::Vector3d Q2 = V.row(v2), Q3 = V.row(v3);
                    Eigen::Vector3d n0 = (Q0 - Q2).cross(Q1 - Q2);
                    Eigen::Vector3d n1 = (Q1 - Q3).cross(Q0 - Q3);
                    Eigen::Vector3d axis = Q1 - Q0;
                    double s = n0.cross(n1).dot(axis) / axis.norm();
                    double c = n0.dot(n1) + n0.norm() * n1.norm();
                    halfTheta = std::atan2(s, c);
                }
                double alpha = halfTheta + orient * thetas[edge];
                if (!(std::abs(alpha) < lim)) return false;
            }
        }
        return true;
    }

    // Build the TinyAD scalar function for the total elastic energy.
    // The returned function references this model; the model must outlive it.
    auto makeEnergyFunction() const {
        auto func = TinyAD::scalar_function<1>(TinyAD::range(nDofs()));
        // Cap OpenMP threads: TinyAD's default (cores-1) oversubscribes on many-core
        // machines and can hang. A modest bounded count is plenty for these problems.
        func.settings.n_threads = threadCount_;

        const MeshConnectivity& mesh = mesh_;
        const RestState& rest = rest_;
        const Eigen::MatrixXi& F = F_;
        const int nVerts = nVerts_;

        // ---- Membrane (stretching): one element per face, 9 scalar variables ----
        func.template add_elements<9>(TinyAD::range(F.rows()),
            [&mesh, &rest, &F](auto& element) -> TINYAD_SCALAR_TYPE(element) {
                using T = TINYAD_SCALAR_TYPE(element);
                int f = element.handle;
                Eigen::Matrix<T, 3, 1> q[3];
                for (int k = 0; k < 3; k++) {
                    int gv = F(f, k);
                    q[k] = Eigen::Matrix<T, 3, 1>(element.variable(3 * gv + 0),
                                                  element.variable(3 * gv + 1),
                                                  element.variable(3 * gv + 2));
                }
                Eigen::Matrix<T, 2, 2> a;
                Eigen::Matrix<T, 3, 1> e1 = q[1] - q[0];
                Eigen::Matrix<T, 3, 1> e2 = q[2] - q[0];
                a << e1.dot(e1), e1.dot(e2), e2.dot(e1), e2.dot(e2);
                Eigen::Matrix<T, 2, 2> abar = castT<T>(rest.abar[f]);
                Eigen::Matrix<T, 2, 2> abarinv = castT<T>(rest.abar[f].inverse());
                Eigen::Matrix<T, 2, 2> M = abarinv * (a - abar);
                double dA = 0.5 * std::sqrt(rest.abar[f].determinant());
                double coeff = rest.thickness[f] / 4.0;
                return coeff * dA * stvk<T>(M, rest.lameAlpha[f], rest.lameBeta[f]);
            });

        // ---- BAC bending: one element per face, up to 21 scalar variables ----
        func.template add_elements<21>(TinyAD::range(F.rows()),
            [&mesh, &rest, &F, nVerts](auto& element) -> TINYAD_SCALAR_TYPE(element) {
                using T = TINYAD_SCALAR_TYPE(element);
                using Vec3 = Eigen::Matrix<T, 3, 1>;
                int f = element.handle;

                // Accessor for a vertex position by global index (dedup handled by TinyAD).
                auto P = [&](int gv) -> Vec3 {
                    return Vec3(element.variable(3 * gv + 0),
                                element.variable(3 * gv + 1),
                                element.variable(3 * gv + 2));
                };

                Vec3 qf[3];
                for (int k = 0; k < 3; k++) qf[k] = P(F(f, k));

                Eigen::Matrix<T, 3, 1> II;
                for (int i = 0; i < 3; i++) {
                    // altitude from local vertex i to opposite edge
                    const Vec3& a0 = qf[i];
                    const Vec3& a1 = qf[(i + 1) % 3];
                    const Vec3& a2 = qf[(i + 2) % 3];
                    Vec3 n = (a1 - a0).cross(a2 - a0);
                    Vec3 ee = a2 - a1;
                    T altitude = n.norm() / ee.norm();

                    int edge = mesh.faceEdge(f, i);
                    double orient = (mesh.faceEdgeOrientation(f, i) == 0) ? 1.0 : -1.0;
                    T thetaEdge = element.variable(3 * nVerts + edge);

                    // dihedral half-angle theta_i/2 = atan2(s, c); 0 on the boundary.
                    T halfTheta = T(0.0);
                    int v2 = mesh.edgeOppositeVertex(edge, 0);
                    int v3 = mesh.edgeOppositeVertex(edge, 1);
                    if (v2 != -1 && v3 != -1) {
                        int v0 = mesh.edgeVertex(edge, 0);
                        int v1 = mesh.edgeVertex(edge, 1);
                        Vec3 Q0 = P(v0), Q1 = P(v1), Q2 = P(v2), Q3 = P(v3);
                        Vec3 n0 = (Q0 - Q2).cross(Q1 - Q2);
                        Vec3 n1 = (Q1 - Q3).cross(Q0 - Q3);
                        Vec3 axis = Q1 - Q0;
                        T s = n0.cross(n1).dot(axis) / axis.norm();
                        T c = n0.dot(n1) + n0.norm() * n1.norm();
                        halfTheta = atan2(s, c);
                    }
                    T alpha = halfTheta + orient * thetaEdge;
                    II[i] = T(2.0) * altitude * tan(alpha);
                }

                Eigen::Matrix<T, 2, 2> b;
                b << II[0] + II[1], II[0], II[0], II[0] + II[2];
                Eigen::Matrix<T, 2, 2> bbar = castT<T>(rest.bbar[f]);
                Eigen::Matrix<T, 2, 2> abarinv = castT<T>(rest.abar[f].inverse());
                Eigen::Matrix<T, 2, 2> M = abarinv * (b - bbar);
                double dA = 0.5 * std::sqrt(rest.abar[f].determinant());
                double coeff = std::pow(rest.thickness[f], 3) / 12.0;
                return coeff * dA * stvk<T>(M, rest.lameAlpha[f], rest.lameBeta[f]);
            });

        return func;
    }

private:
    MeshConnectivity mesh_;
    Eigen::MatrixXd V0_;
    Eigen::MatrixXi F_;
    RestState rest_;
    Eigen::VectorXd restEdgeThetas_;
    int nVerts_ = 0, nEdges_ = 0;
    int threadCount_ = 8;

public:
    void setThreadCount(int n) { threadCount_ = n; }
};

} // namespace bac
