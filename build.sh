#!/bin/bash

# Create a build directory
mkdir build

# Change to the build directory
cd build

# Run CMake to generate the build files
cmake ..

# Build the project
make

# Package the executable and documentation into a zip file
make package
