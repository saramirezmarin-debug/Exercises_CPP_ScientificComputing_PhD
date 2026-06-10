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

cmake --build build --target docs

echo ""
echo "=============================="
echo "Opening documentation..."
echo "=============================="

start docs/html/index.html

echo ""
echo "=============================="
echo "Running executable..."
echo "=============================="

./exe/Exercise4.exe