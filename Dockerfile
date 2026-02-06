# Use Ubuntu 20.04 LTS with older PETSc 3.12 (closer to historical dev environment)
FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install dependencies
# - build-essential: gcc, g++, make
# - mpi-default-dev: MPI compiler headers and libraries
# - libpetsc-real3.12-dev: PETSc 3.12 library (real numbers)
# - libslepc-real3.12-dev: SLEPc 3.12 library (real numbers)
# - python3, python3-pip: For visualization scripts
RUN apt-get update && apt-get install -y \
    build-essential \
    mpi-default-dev \
    libpetsc-real3.12-dev \
    libslepc-real3.12-dev \
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

# Copy source code and configuration
COPY src /app/src
COPY config /app/config
COPY mesh /app/mesh
COPY mesh_4 /app/mesh_4

# Set working directory to src for building
WORKDIR /app/src

# Build the application
# We override variables to ensure they use the container's paths if needed
# Use -O1 instead of -O2 to avoid optimizer-triggered memory bugs
RUN make clean && make all ASAN_FLAGS="-g -O1"

# Ensure scripts are executable
RUN chmod +x run_emr_sweep.sh

# Command to run by default.
# Execute the full EMR simulation sweep and plotting pipeline
CMD ["/bin/bash", "-c", "./run_emr_sweep.sh && python3 plot_emr_vs_h.py"]
