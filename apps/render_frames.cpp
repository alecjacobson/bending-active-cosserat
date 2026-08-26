// Render a list of OBJ files to PNGs via headless polyscope (EGL).
// Usage: render_frames <out_dir> <cam_eye_x,y,z> <cam_target_x,y,z> <obj1> <obj2> ...
// If cam args are "-" a default camera is used.
#include <cstdio>
#include <string>
#include <vector>
#include <Eigen/Core>
#include "Render.h"
#include "Meshes.h"

using namespace bac;

static bool parse3(const std::string& s, std::array<double, 3>& out) {
    return sscanf(s.c_str(), "%lf,%lf,%lf", &out[0], &out[1], &out[2]) == 3;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::printf("usage: %s <outdir> <eye|-> <target|-> <obj...>\n", argv[0]);
        return 1;
    }
    std::string outdir = argv[1];
    Camera cam;
    parse3(argv[2], cam.eye);
    parse3(argv[3], cam.target);

    renderInit(1280, 960);
    for (int i = 4; i < argc; i++) {
        Eigen::MatrixXd V; Eigen::MatrixXi F;
        if (!readOBJ(argv[i], V, F)) { std::printf("skip %s\n", argv[i]); continue; }
        std::string base = argv[i];
        auto slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        auto dot = base.find_last_of('.');
        if (dot != std::string::npos) base = base.substr(0, dot);
        std::string png = outdir + "/" + base + ".png";
        renderMeshPNG(V, F, png, cam, 2);
        std::printf("rendered %s\n", png.c_str());
    }
    return 0;
}
