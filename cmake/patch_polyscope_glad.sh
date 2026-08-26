#!/usr/bin/env bash
# polyscope v2.4.0 only adds its bundled `glad` under the GLFW-backend branch, but the
# EGL (headless) backend needs it too. Prepend an unconditional glad add to the deps
# CMake so EGL-only builds link. Idempotent. Run from the polyscope source root.
set -e
f="deps/CMakeLists.txt"
if grep -q "BAC_GLAD_PATCH" "$f"; then
  exit 0
fi
{
  echo "# BAC_GLAD_PATCH"
  echo "if(NOT TARGET glad)"
  echo "  add_subdirectory(glad)"
  echo "endif()"
  cat "$f"
} > "$f.patched"
mv "$f.patched" "$f"
