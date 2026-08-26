// Headless polyscope rendering helpers (EGL offscreen). Compiled only when
// BAC_WITH_RENDER is defined.
#pragma once
#include <Eigen/Core>
#include <array>
#include <string>

namespace bac {

// Initialize polyscope with the headless EGL backend. Call once.
void renderInit(int width = 1280, int height = 960);

// Pin the shadow/ground plane to a fixed world-up (Z) height, so it does not move
// between frames of an animation. Call once before rendering a sequence.
void renderSetFixedGround(double heightZ);

// Camera specification in world coordinates.
struct Camera {
    std::array<double, 3> eye{2.0, -1.5, 1.2};
    std::array<double, 3> target{0.0, 0.5, 0.0};
    std::array<double, 3> up{0.0, 0.0, 1.0};
};

// Render one mesh to a PNG. Positions are first rotated by R (world = R * v) for
// display, e.g. to point the gravity/load direction downward. `heightAxis` (0/1/2)
// colors the surface by that ORIGINAL (pre-rotation) coordinate; pass -1 for flat.
void renderMeshPNG(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                   const std::string& pngPath, const Camera& cam, int heightAxis = 2,
                   const Eigen::Matrix3d& R = Eigen::Matrix3d::Identity());

} // namespace bac
