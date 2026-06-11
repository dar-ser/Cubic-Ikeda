#!/usr/bin/env bash
set -e

echo "=== Installing system dependencies ==="
sudo apt-get update
sudo apt-get install -y build-essential cmake g++ libeigen3-dev

echo "=== Initializing CAPD submodule ==="
git submodule update --init --recursive

echo "=== Building CAPD ==="
CAPD_DIR="external/capd"
mkdir -p ${CAPD_DIR}/build
cd ${CAPD_DIR}/build

cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

if [ ! -f "./bin/capd-config" ]; then
    echo "CAPD build failed: capd-config not found"
    exit 1
fi

cd ../../..

echo "=== Exporting environment ==="
export CAPDBINDIR="$(pwd)/external/capd/build/bin/"
export PATH="$CAPDBINDIR:$PATH"

echo "CAPDBINDIR set to: $CAPDBINDIR"

echo "=== Building project ==="
make clean
make -j$(nproc)

echo "=== DONE ==="
echo "Run:"
echo "  make run: executes the code that proofs existence of a periodic orbit"