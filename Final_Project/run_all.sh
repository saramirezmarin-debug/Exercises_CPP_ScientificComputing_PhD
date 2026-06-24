#!/bin/bash

set -e

echo "=============================="
echo "Cleaning previous build..."
echo "=============================="

rm -rf build

echo ""
echo "=============================="
echo "Configuring project with CMake + Ninja..."
echo "=============================="

cmake -S . -B build -G Ninja

echo ""
echo "=============================="
echo "Building project..."
echo "=============================="

cmake --build build

echo ""
echo "=============================="
echo "Generating Doxygen docs..."
echo "=============================="

cmake --build build --target doc

echo ""
echo "=============================="
echo "Opening documentation..."
echo "=============================="

DOC_INDEX="$PWD/build/docs/html/index.html"

if [ -f "$DOC_INDEX" ]; then
    WIN_DOC_INDEX=$(cygpath -w "$DOC_INDEX")
    powershell.exe -NoProfile -ExecutionPolicy Bypass -Command "Start-Process -FilePath '$WIN_DOC_INDEX'"
else
    echo "Documentation index not found at:"
    echo "$DOC_INDEX"
    exit 1
fi

echo ""
echo "=============================="
echo "Running executable..."
echo "=============================="

./exe/CSC.exe