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
echo "Running executable..."
echo "=============================="

./build/Exercise2.exe