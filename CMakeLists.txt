cmake_minimum_required(VERSION 3.16)
project(bending_active_cosserat LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(BAC_WITH_RENDER "Build polyscope headless rendering apps" OFF)

include(FetchContent)

# ---- Eigen (header only) ----
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)

# ---- TinyAD (header only; needs Eigen3::Eigen) ----
FetchContent_Declare(tinyad
  GIT_REPOSITORY https://github.com/patr-schm/TinyAD.git
  GIT_TAG main
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(tinyad)

# ---- core library ----
add_library(bac_core STATIC src/MeshConnectivity.cpp)
target_include_directories(bac_core PUBLIC src)
target_link_libraries(bac_core PUBLIC Eigen3::Eigen TinyAD)

# ---- validation / test app (finite-diff + golden checks) ----
add_executable(validate apps/validate.cpp)
target_link_libraries(validate PRIVATE bac_core)

# ---- solver test ----
add_executable(test_solver apps/test_solver.cpp)
target_link_libraries(test_solver PRIVATE bac_core)

# ---- dynamic drape (pure physics; writes OBJ frames, render with render_frames) ----
add_executable(drape apps/drape.cpp)
target_link_libraries(drape PRIVATE bac_core)

# ---- curved (non-flat) rest-state arch ----
add_executable(curved_rest apps/curved_rest.cpp)
target_link_libraries(curved_rest PRIVATE bac_core)

# ---- lateral buckling flagship ----
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/apps/lateral_buckling.cpp)
  add_executable(lateral_buckling apps/lateral_buckling.cpp)
  target_link_libraries(lateral_buckling PRIVATE bac_core)
endif()

if(BAC_WITH_RENDER)
  # ---- polyscope (headless via EGL on Linux) ----
  set(POLYSCOPE_BACKEND_OPENGL3_GLFW OFF CACHE BOOL "" FORCE)
  set(POLYSCOPE_BACKEND_OPENGL3_EGL  ON  CACHE BOOL "" FORCE)
  # polyscope v2.4.0 only adds its bundled `glad` under the GLFW-backend branch, but
  # the EGL (headless) backend needs it too. Patch deps/CMakeLists.txt to add glad
  # unconditionally when we build EGL-only. (Idempotent.)
  FetchContent_Declare(polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG v2.4.0
    GIT_SHALLOW TRUE
    PATCH_COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_polyscope_glad.sh)
  FetchContent_MakeAvailable(polyscope)

  add_library(bac_render src/Render.cpp)
  target_include_directories(bac_render PUBLIC src)
  target_link_libraries(bac_render PUBLIC bac_core polyscope)
  target_compile_definitions(bac_render PUBLIC BAC_WITH_RENDER)

  add_executable(render_frames apps/render_frames.cpp)
  target_link_libraries(render_frames PRIVATE bac_render)
endif()

add_executable(bench_solver apps/bench_solver.cpp)
target_link_libraries(bench_solver PRIVATE bac_core)
