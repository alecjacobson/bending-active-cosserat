// Headless polyscope rendering helpers (EGL offscreen). Compiled only when
// BAC_WITH_RENDER is defined.
#pragma once
#include <Eigen/Core>
#include <array>
#include <string>

namespace bac {

// Initialize polyscope with the headless EGL backend. Call once.
void renderInit(int width = 1280, int height = 960);

// Camera specification in world coordinates.
struct Camera {
    std::array<double, 3> eye{2.0, -1.5, 1.2};
    std::array<double, 3> target{0.0, 0.5, 0.0};
    std::array<double, 3> up{0.0, 0.0, 1.0};
};

// Render one mesh to a PNG. `heightAxis` (0/1/2) colors the surface by that
// coordinate for depth cueing; pass -1 for a flat color.
void renderMeshPNG(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                   const std::string& pngPath, const Camera& cam, int heightAxis = 2);

} // namespace bac
