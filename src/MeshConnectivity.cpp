#include "MeshConnectivity.h"
#include <map>
#include <utility>

namespace bac {

MeshConnectivity::MeshConnectivity(const Eigen::MatrixXi& F) : F_(F) {
    const int nfaces = (int)F.rows();

    // Collect edges. For each face, local edge j is opposite vertex j and spans
    // (v_{j+1}, v_{j+2}). Store which side (idx) the incident face sits on so that
    // edgeFace(e,0)/edgeFace(e,1) are stable.
    std::map<std::pair<int, int>, Eigen::Vector2i> edgeFaces;
    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int v0 = F(i, (j + 1) % 3);
            int v1 = F(i, (j + 2) % 3);
            int idx = 0;
            if (v0 > v1) { std::swap(v0, v1); idx = 1; }
            std::pair<int, int> p(v0, v1);
            auto it = edgeFaces.find(p);
            if (it == edgeFaces.end()) {
                edgeFaces[p][idx] = i;
                edgeFaces[p][1 - idx] = -1;
            } else {
                edgeFaces[p][idx] = i;
            }
        }
    }

    const int nedges = (int)edgeFaces.size();
    FE_.resize(nfaces, 3);
    FEorient_.resize(nfaces, 3);
    EV_.resize(nedges, 2);
    EF_.resize(nedges, 2);
    EOpp_.resize(nedges, 2);

    std::map<std::pair<int, int>, int> edgeIndices;
    int idx = 0;
    for (auto& it : edgeFaces) {
        edgeIndices[it.first] = idx;
        EV_(idx, 0) = it.first.first;
        EV_(idx, 1) = it.first.second;
        EF_(idx, 0) = it.second[0];
        EF_(idx, 1) = it.second[1];
        idx++;
    }

    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int v0 = F(i, (j + 1) % 3);
            int v1 = F(i, (j + 2) % 3);
            if (v0 > v1) std::swap(v0, v1);
            FE_(i, j) = edgeIndices[std::make_pair(v0, v1)];
        }
    }

    // Opposite vertex (apex) of each incident face.
    for (int e = 0; e < nedges; e++) {
        for (int s = 0; s < 2; s++) {
            int face = EF_(e, s);
            EOpp_(e, s) = -1;
            if (face == -1) continue;
            for (int j = 0; j < 3; j++) {
                int v = F(face, j);
                if (v != EV_(e, 0) && v != EV_(e, 1)) { EOpp_(e, s) = v; break; }
            }
        }
    }

    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int edge = FE_(i, j);
            FEorient_(i, j) = (EF_(edge, 0) == i) ? 0 : 1;
        }
    }
}

int MeshConnectivity::vertexOppositeFaceEdge(int face, int j) const {
    int edge = FE_(face, j);
    int orient = FEorient_(face, j);
    return EOpp_(edge, 1 - orient);
}

} // namespace bac
