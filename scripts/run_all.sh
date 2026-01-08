#!/bin/bash
#
# Master script to run the complete brick kiln air quality impact analysis
#
# Usage: ./run_all.sh [--skip-met] [--skip-emissions] [--scenario SCENARIO]
#
# Steps:
#   1. Generate meteorology (WRFCAMx)
#   2. Process emissions (EDGAR + brick kilns)
#   3. Run CAMx air quality model
#   4. Analyze results
#

set -e  # Exit on error

# ============================================================================
# Configuration - EDIT THESE PATHS FOR YOUR SYSTEM
# ============================================================================
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Path to WRF output files (wrfout_d01_*)
# Default: use included sample data (1 day)
# Override: WRF_INPUT=/path/to/full/wrf/output ./run_all.sh
WRF_INPUT="${WRF_INPUT:-$PROJECT_ROOT/wrf_sample}"

# Path to CAMx v7.32 source
# Default: use included source
CAMX_ROOT="${CAMX_ROOT:-$PROJECT_ROOT/camx_src}"

# Default options
SKIP_MET=false
SKIP_EMISSIONS=false
SCENARIO="all"  # base, zigzag_50, or all

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-met)
            SKIP_MET=true
            shift
            ;;
        --skip-emissions)
            SKIP_EMISSIONS=true
            shift
            ;;
        --scenario)
            SCENARIO="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--skip-met] [--skip-emissions] [--scenario SCENARIO]"
            echo ""
            echo "Options:"
            echo "  --skip-met       Skip meteorology generation (use existing)"
            echo "  --skip-emissions Skip emissions processing"
            echo "  --scenario       Run specific scenario (base, zigzag_50, or all)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=============================================================="
echo "  BRICK KILN AIR QUALITY IMPACT ANALYSIS"
echo "  WRF-CAMx Modeling Pipeline"
echo "=============================================================="
echo ""
echo "Project root: $PROJECT_ROOT"
echo "Scenario: $SCENARIO"
echo ""

# ============================================================================
# Step 1: Generate Meteorology with WRFCAMx
# ============================================================================
if [ "$SKIP_MET" = false ]; then
    echo "=============================================================="
    echo "Step 1: Generating CAMx meteorology from WRF output"
    echo "=============================================================="

    cd "$PROJECT_ROOT/wrfcamx_v5.2"

    # Check if WRF output exists
    if [ ! -d "$WRF_INPUT" ]; then
        echo "ERROR: WRF output directory not found: $WRF_INPUT"
        exit 1
    fi

    # Generate snow files if needed
    if [ ! -f snow.file.240131 ]; then
        echo "Generating snow age files..."
        gfortran -o make_snow make_snow.f90
        ./make_snow
    fi

    # Run WRFCAMx
    echo "Running WRFCAMx..."
    ./test_camx.job

    echo "Meteorology files generated in: $PROJECT_ROOT/wrfcamx_v5.2/camx_input/"
else
    echo "Step 1: Skipping meteorology generation (--skip-met)"
fi

# ============================================================================
# Step 2: Process Emissions
# ============================================================================
if [ "$SKIP_EMISSIONS" = false ]; then
    echo ""
    echo "=============================================================="
    echo "Step 2: Processing emissions (EDGAR + brick kilns)"
    echo "=============================================================="

    cd "$PROJECT_ROOT/edgar_emissions"

    # Check for EDGAR data
    if [ ! -d "raw_data" ] || [ -z "$(ls -A raw_data 2>/dev/null)" ]; then
        echo "NOTE: EDGAR data not found in raw_data/"
        echo "      Download from: https://edgar.jrc.ec.europa.eu/"
        echo "      Running with synthetic emissions for demonstration..."
    fi

    # Process emissions
    python3 process_edgar_to_camx.py

    echo "Emissions files generated in: $PROJECT_ROOT/edgar_emissions/camx_emissions/"
else
    echo "Step 2: Skipping emissions processing (--skip-emissions)"
fi

# ============================================================================
# Step 3: Run CAMx Air Quality Model
# ============================================================================
echo ""
echo "=============================================================="
echo "Step 3: Running CAMx air quality simulations"
echo "=============================================================="

cd "$PROJECT_ROOT/run_camx"

# Check if CAMx executable exists
CAMX_EXEC="$CAMX_ROOT/CAMx.v7.32.noMPI.gfort"
if [ ! -f "$CAMX_EXEC" ]; then
    echo "WARNING: CAMx executable not found at $CAMX_EXEC"
    echo "         Skipping CAMx simulation. To run:"
    echo "         1. Compile CAMx: cd $CAMX_ROOT && make"
    echo "         2. Re-run this script"
    SKIP_CAMX=true
else
    SKIP_CAMX=false
fi

if [ "$SKIP_CAMX" = false ]; then
    if [ "$SCENARIO" = "all" ] || [ "$SCENARIO" = "base" ]; then
        echo "Running base case (100% FCBTK)..."
        sed -i 's/set SCENARIO = .*/set SCENARIO = "base"/' CAMx.india.job
        ./CAMx.india.job
    fi

    if [ "$SCENARIO" = "all" ] || [ "$SCENARIO" = "zigzag_50" ]; then
        echo "Running scenario: 50% Zigzag conversion..."
        sed -i 's/set SCENARIO = .*/set SCENARIO = "zigzag_50"/' CAMx.india.job
        ./CAMx.india.job
    fi
fi

# ============================================================================
# Step 4: Run Brick Kiln Impact Analysis
# ============================================================================
echo ""
echo "=============================================================="
echo "Step 4: Analyzing brick kiln impact on Delhi air quality"
echo "=============================================================="

cd "$PROJECT_ROOT/brick_kiln_analysis"
python3 brick_kiln_simulation.py

# Move figures to figures directory
mv -f *.png "$PROJECT_ROOT/figures/" 2>/dev/null || true

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=============================================================="
echo "  ANALYSIS COMPLETE"
echo "=============================================================="
echo ""
echo "Output locations:"
echo "  Meteorology:  $PROJECT_ROOT/wrfcamx_v5.2/camx_input/"
echo "  Emissions:    $PROJECT_ROOT/edgar_emissions/camx_emissions/"
echo "  CAMx output:  $PROJECT_ROOT/output/camx_results/"
echo "  Figures:      $PROJECT_ROOT/figures/"
echo ""
echo "Key results:"
echo "  - brick_kiln_analysis.png    : Spatial emissions and concentration maps"
echo "  - summary_infographic.png    : Summary of Delhi PM2.5 impact"
echo ""
