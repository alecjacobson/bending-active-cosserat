// Render a list of OBJ files to PNGs via headless polyscope (EGL), with a fixed
// shadow floor across the whole sequence and an optional world-up remap so the
// gravity/load direction points downward.
//
// Usage:
//   render_frames <outdir> <eye_x,y,z> <target_x,y,z> <up_axis> <heightAxis> <obj...>
// where <up_axis> is a signed sim axis that should point up, e.g. z+ (default),
// x-, y+ ... and <heightAxis> in {0,1,2,-1} selects the coloring coordinate.
#include <cstdio>
#include <string>
#include <vector>
#include <limits>
#include <Eigen/Core>
#include "Render.h"
#include "Meshes.h"

using namespace bac;

static bool parse3(const std::string& s, std::array<double, 3>& out) {
    return sscanf(s.c_str(), "%lf,%lf,%lf", &out[0], &out[1], &out[2]) == 3;
}

// Rotation mapping the given signed sim axis (e.g. "x-") to world +Z (up).
static Eigen::Matrix3d upRotation(const std::string& axis) {
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    if (axis == "z+" || axis.empty()) R << 1,0,0, 0,1,0, 0,0,1;
    else if (axis == "z-") R << 1,0,0, 0,-1,0, 0,0,-1;
    else if (axis == "x+") R << 0,0,-1, 0,1,0, 1,0,0;
    else if (axis == "x-") R << 0,0,1, 0,1,0, -1,0,0;
    else if (axis == "y+") R << 1,0,0, 0,0,-1, 0,1,0;
    else if (axis == "y-") R << 1,0,0, 0,0,1, 0,-1,0;
    return R;
}

int main(int argc, char** argv) {
    if (argc < 7) {
        std::printf("usage: %s <outdir> <eye> <target> <up_axis> <heightAxis> <obj...>\n", argv[0]);
        return 1;
    }
    std::string outdir = argv[1];
    Camera cam;
    parse3(argv[2], cam.eye);
    parse3(argv[3], cam.target);
    Eigen::Matrix3d R = upRotation(argv[4]);
    int heightAxis = std::atoi(argv[5]);

    renderInit(1280, 960);

    // First pass: find the global lowest world-Z across all frames for a fixed floor.
    double minZ = std::numeric_limits<double>::infinity();
    std::vector<std::string> objs;
    for (int i = 6; i < argc; i++) {
        Eigen::MatrixXd V; Eigen::MatrixXi F;
        if (!readOBJ(argv[i], V, F)) { std::printf("skip %s\n", argv[i]); continue; }
        Eigen::MatrixXd Vr = V * R.transpose();
        minZ = std::min(minZ, Vr.col(2).minCoeff());
        objs.push_back(argv[i]);
    }
    double span = 1.0;
    renderSetFixedGround(minZ - 0.02 * span);

    for (const auto& path : objs) {
        Eigen::MatrixXd V; Eigen::MatrixXi F;
        if (!readOBJ(path, V, F)) continue;
        std::string base = path;
        auto slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        auto dot = base.find_last_of('.');
        if (dot != std::string::npos) base = base.substr(0, dot);
        std::string png = outdir + "/" + base + ".png";
        renderMeshPNG(V, F, png, cam, heightAxis, R);
        std::printf("rendered %s\n", png.c_str());
    }
    return 0;
}
