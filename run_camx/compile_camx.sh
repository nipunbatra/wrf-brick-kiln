#!/bin/bash
#
# CAMx v7.32 Compilation Script
# For Brick Kiln Air Quality Impact Analysis
#

set -e

# =============================================================================
# Configuration
# =============================================================================
CAMX_SRC="/home/nipun.batra/git/wrf-brick-kiln/camx_src"
NCF_INST="/home/deenalad/WRF_env/BUILD_WRF/LIBRARIES/netcdf"
MPI_INST="/home/deenalad/WRF_env/BUILD_WRF/LIBRARIES/mpich"

# Output executable location
OUTPUT_DIR="/home/nipun.batra/git/wrf-brick-kiln/run_camx"

echo "=============================================="
echo "CAMx v7.32 Compilation Script"
echo "=============================================="

# Check if source exists
if [ ! -d "$CAMX_SRC" ]; then
    echo "ERROR: CAMx source directory not found at $CAMX_SRC"
    exit 1
fi

cd "$CAMX_SRC"

# Clean previous build
echo "Cleaning previous build..."
make clean 2>/dev/null || true

# Set environment variables
export LD_LIBRARY_PATH="${NCF_INST}/lib:${MPI_INST}/lib:$LD_LIBRARY_PATH"

# =============================================================================
# Compile CAMx with gfortran, no MPI, with NetCDF4 support
# =============================================================================
echo ""
echo "Compiling CAMx with gfortran (no MPI, NetCDF4 with compression)..."
echo "This may take 10-15 minutes..."
echo ""

# Build command - using gfortran without MPI (as per Makefile note about OMP issues)
make COMPILER=gfort NCF=NCF4_C MPI=false 2>&1 | tee compile.log

# Check if compilation succeeded
EXEC_NAME="CAMx.v7.32.noMPI.NCF4.gfort"
if [ -f "$EXEC_NAME" ]; then
    echo ""
    echo "=============================================="
    echo "Compilation SUCCESSFUL!"
    echo "Executable: $CAMX_SRC/$EXEC_NAME"
    echo "=============================================="

    # Copy to output directory
    mkdir -p "$OUTPUT_DIR"
    cp "$EXEC_NAME" "$OUTPUT_DIR/"
    chmod +x "$OUTPUT_DIR/$EXEC_NAME"

    echo "Copied to: $OUTPUT_DIR/$EXEC_NAME"
else
    echo ""
    echo "=============================================="
    echo "ERROR: Compilation FAILED!"
    echo "Check compile.log for details"
    echo "=============================================="
    exit 1
fi
