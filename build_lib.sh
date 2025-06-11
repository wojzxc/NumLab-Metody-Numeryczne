#!/usr/bin/env bash
set -e                     # przerwij przy pierwszym b³êdzie

# 1. katalog budowania
BUILD_DIR="build"
rm -rf "${BUILD_DIR}"
mkdir  "${BUILD_DIR}"

# 2. konfiguracja – wy³¹czamy examples i tests
cmake -S . -B "${BUILD_DIR}" \
      -DCMAKE_BUILD_TYPE=Debug \
      -DNUMLAB_BUILD_EXAMPLES=OFF \
      -DNUMLAB_BUILD_TESTS=OFF

# 3. kompilacja wy³¹cznie celu NumLab
cmake --build "${BUILD_DIR}" --target NumLab --config Debug

echo "========================================================"
echo "Biblioteka zbudowana:"
if [ -f "${BUILD_DIR}/libNumLab.a" ];        then echo "  ${BUILD_DIR}/libNumLab.a"; fi
if [ -f "${BUILD_DIR}/Debug/NumLab.lib" ]; then echo "  ${BUILD_DIR}/Debug/NumLab.lib"; fi
echo "Nag³ówki:  include/  (kopiuj lub do³¹cz do -I)"
