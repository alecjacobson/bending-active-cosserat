// Mesh connectivity for shell simulation.
//
// Reproduces the edge/face adjacency conventions used by libshell so that our
// from-scratch BAC energy can be validated numerically against it. Conventions:
//   - faceEdge(f,i)         : the edge opposite local vertex i of face f
//   - faceEdgeOrientation   : 0 if f == edgeFace(edge,0), else 1
//   - edgeVertex(e,0/1)     : the two endpoints of edge e (sorted ascending)
//   - edgeFace(e,0/1)       : the two incident faces (or -1 on the boundary)
//   - edgeOppositeVertex    : the apex of edgeFace(e,0/1) not on the edge (or -1)
#pragma once
#include <Eigen/Core>

namespace bac {

class MeshConnectivity {
public:
    MeshConnectivity() = default;
    explicit MeshConnectivity(const Eigen::MatrixXi& F);

    int nFaces() const { return (int)F_.rows(); }
    int nEdges() const { return (int)EV_.rows(); }

    int faceVertex(int face, int j) const { return F_(face, j); }
    int faceEdge(int face, int j) const { return FE_(face, j); }
    int faceEdgeOrientation(int face, int j) const { return FEorient_(face, j); }

    int edgeVertex(int edge, int j) const { return EV_(edge, j); }
    int edgeFace(int edge, int j) const { return EF_(edge, j); }
    int edgeOppositeVertex(int edge, int j) const { return EOpp_(edge, j); }

    // Vertex opposite face f across its edge j (the apex of the neighbor triangle),
    // or -1 if that edge is on the boundary.
    int vertexOppositeFaceEdge(int face, int j) const;

    const Eigen::MatrixXi& faces() const { return F_; }

private:
    Eigen::MatrixXi F_, FE_, FEorient_, EV_, EF_, EOpp_;
};

} // namespace bac
