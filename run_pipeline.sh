#!/bin/bash
#
# =============================================================================
# WRF-CAMx Brick Kiln Impact Analysis Pipeline
# =============================================================================
#
# This master script runs the complete pipeline:
#   1. Generate brick kiln locations and emissions (kiln_camx_analysis.py)
#   2. Create CAMx input files (create_camx_inputs.py)
#   3. Compile CAMx if needed (compile_camx.sh)
#   4. Run CAMx for both 2020 and 2025 scenarios
#   5. Analyze and visualize results (analyze_camx_output.py)
#
# Usage: ./run_pipeline.sh [step]
#        step = all, emissions, inputs, compile, run, analyze
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

# =============================================================================
# STEP 1: Generate Emissions
# =============================================================================
run_emissions() {
    print_header "STEP 1: Generating Brick Kiln Emissions"

    cd "$SCRIPT_DIR/brick_kiln_analysis"

    if [ -f kiln_camx_analysis.py ]; then
        python3 kiln_camx_analysis.py
        print_success "Emissions generated"
    else
        print_error "kiln_camx_analysis.py not found"
        exit 1
    fi

    cd "$SCRIPT_DIR"
}

# =============================================================================
# STEP 2: Create CAMx Inputs
# =============================================================================
run_inputs() {
    print_header "STEP 2: Creating CAMx Input Files"

    cd "$SCRIPT_DIR/run_camx"

    if [ -f create_camx_inputs.py ]; then
        python3 create_camx_inputs.py
        print_success "CAMx inputs created"
    else
        print_error "create_camx_inputs.py not found"
        exit 1
    fi

    cd "$SCRIPT_DIR"
}

# =============================================================================
# STEP 3: Compile CAMx
# =============================================================================
run_compile() {
    print_header "STEP 3: Compiling CAMx v7.32"

    cd "$SCRIPT_DIR/run_camx"

    # Check if executable already exists
    if [ -f "CAMx.v7.32.noMPI.NCF4.gfort" ]; then
        print_warning "CAMx executable already exists, skipping compilation"
        print_success "Using existing executable"
    else
        if [ -f compile_camx.sh ]; then
            ./compile_camx.sh
            print_success "CAMx compiled"
        else
            print_error "compile_camx.sh not found"
            exit 1
        fi
    fi

    cd "$SCRIPT_DIR"
}

# =============================================================================
# STEP 4: Run CAMx
# =============================================================================
run_camx() {
    print_header "STEP 4: Running CAMx for Both Scenarios"

    cd "$SCRIPT_DIR/run_camx"

    if [ ! -f "CAMx.brick_kiln.job" ]; then
        print_error "CAMx.brick_kiln.job not found"
        exit 1
    fi

    # Run 2020 scenario (base case - 100% FCBK)
    echo ""
    echo -e "${YELLOW}Running 2020 scenario (100% FCBK)...${NC}"
    ./CAMx.brick_kiln.job 2020
    print_success "2020 scenario complete"

    # Run 2025 scenario (converted - 50% Zigzag)
    echo ""
    echo -e "${YELLOW}Running 2025 scenario (50% Zigzag)...${NC}"
    ./CAMx.brick_kiln.job 2025
    print_success "2025 scenario complete"

    cd "$SCRIPT_DIR"
}

# =============================================================================
# STEP 5: Analyze Results
# =============================================================================
run_analyze() {
    print_header "STEP 5: Analyzing and Visualizing Results"

    cd "$SCRIPT_DIR/run_camx"

    if [ -f analyze_camx_output.py ]; then
        python3 analyze_camx_output.py
        print_success "Analysis complete"
    else
        print_error "analyze_camx_output.py not found"
        exit 1
    fi

    cd "$SCRIPT_DIR"
}

# =============================================================================
# MAIN
# =============================================================================

print_header "WRF-CAMx BRICK KILN IMPACT ANALYSIS PIPELINE"

STEP=${1:-all}

case $STEP in
    emissions)
        run_emissions
        ;;
    inputs)
        run_inputs
        ;;
    compile)
        run_compile
        ;;
    run)
        run_camx
        ;;
    analyze)
        run_analyze
        ;;
    all)
        run_emissions
        run_inputs
        run_compile
        run_camx
        run_analyze
        ;;
    *)
        echo "Usage: $0 [step]"
        echo "  step = all, emissions, inputs, compile, run, analyze"
        echo ""
        echo "Steps:"
        echo "  emissions - Generate brick kiln locations and emission files"
        echo "  inputs    - Create CAMx input files (IC/BC/photolysis)"
        echo "  compile   - Compile CAMx v7.32"
        echo "  run       - Run CAMx for both 2020 and 2025 scenarios"
        echo "  analyze   - Analyze and visualize CAMx output"
        echo "  all       - Run complete pipeline"
        exit 1
        ;;
esac

print_header "PIPELINE COMPLETE!"

echo ""
echo "Output locations:"
echo "  Figures:     $SCRIPT_DIR/figures/"
echo "  Emissions:   $SCRIPT_DIR/camx_emissions/"
echo "  CAMx Output: $SCRIPT_DIR/camx_output/"
echo ""
echo "Key results:"
echo "  - figures/brick_kiln_camx_analysis.png"
echo "  - figures/camx_pm25_comparison.png"
echo "  - figures/camx_delhi_comparison.png"
echo ""
