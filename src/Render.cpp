#include "Render.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include <glm/glm.hpp>
#include <vector>

namespace bac {

static bool g_inited = false;

void renderInit(int width, int height) {
    if (g_inited) return;
    polyscope::options::allowHeadlessBackends = true;
    polyscope::options::programName = "bac-render";
    polyscope::options::verbosity = 0;
    polyscope::options::ssaaFactor = 2;                 // anti-aliasing
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::init("openGL3_egl");                     // headless EGL
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setWindowSize(width, height);
    g_inited = true;
}

void renderMeshPNG(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                   const std::string& pngPath, const Camera& cam, int heightAxis) {
    renderInit();
    auto* m = polyscope::registerSurfaceMesh("shell", V, F);
    m->setSurfaceColor({0.85, 0.85, 0.90});
    m->setEdgeWidth(1.0);
    m->setEdgeColor({0.2, 0.2, 0.25});
    m->setSmoothShade(true);
    if (heightAxis >= 0 && heightAxis < 3) {
        std::vector<double> h(V.rows());
        for (int i = 0; i < V.rows(); i++) h[i] = V(i, heightAxis);
        auto* q = m->addVertexScalarQuantity("coord", h);
        q->setColorMap("coolwarm");
        q->setEnabled(true);
    }
    polyscope::view::lookAt(
        glm::vec3(cam.eye[0], cam.eye[1], cam.eye[2]),
        glm::vec3(cam.target[0], cam.target[1], cam.target[2]),
        glm::vec3(cam.up[0], cam.up[1], cam.up[2]), false);
    polyscope::screenshot(pngPath, false);
    polyscope::removeStructure("shell");
}

} // namespace bac
