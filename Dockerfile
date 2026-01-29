# Use Ubuntu 22.04 LTS as base image for stability and compatibility
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install dependencies
# - build-essential: gcc, g++, make
# - mpi-default-dev: MPI compiler headers and libraries
# - libpetsc-real-dev: PETSc library (real numbers)
# - libslepc-real-dev: SLEPc library (real numbers)
# - python3, python3-pip: For visualization scripts
# - git: For version control (optional but good practice)
RUN apt-get update && apt-get install -y \
    build-essential \
    mpi-default-dev \
    libpetsc-real-dev \
    libslepc-real-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python Python 3 dependencies for visualization
# - matplotlib: For plotting
# - numpy: For numerical operations
RUN pip3 install --no-cache-dir matplotlib numpy

# Set Environment Variables for PETSc and SLEPc
# Default Installation paths on Ubuntu 22.04
ENV PETSC_DIR=/usr/lib/petsc
ENV SLEPC_DIR=/usr/lib/slepc

# Set working directory
WORKDIR /app

# Create output directory structure
RUN mkdir -p /app/output/figures

# Copy source code
COPY src /app/src

# Set working directory to src for building
WORKDIR /app/src

# Build the application
# We override variables to ensure they use the container's paths if needed, 
# though Env vars should handle it.
RUN make clean && make all

# Command to run by default. 
# We'll default to running the simulation and then the visualizer.
# Using a shell script to chain commands is often cleaner.
CMD ["/bin/bash", "-c", "./Laplace && python3 visualize_mesh.py ../mesh"]
