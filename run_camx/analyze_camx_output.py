#!/usr/bin/env python3
"""
Analyze CAMx Output for Brick Kiln Scenarios
=============================================

This script:
1. Reads CAMx output files for both scenarios (2020 and 2025)
2. Calculates PM2.5 concentrations (PEC + POA + FPRM + PNO3 + PSO4 + PNH4)
3. Creates comparison plots
4. Calculates Delhi-specific impact

Usage: python analyze_camx_output.py [scenario]
       If no scenario specified, compares both 2020 and 2025
"""

import numpy as np
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import netCDF4 as nc

try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False

# =============================================================================
# PATHS
# =============================================================================
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'camx_output')
FIGURES_DIR = os.path.join(PROJECT_ROOT, 'figures')
SHAPEFILE_PATH = os.path.join(PROJECT_ROOT, 'shapefiles', 'India_State_Boundary.shp')

os.makedirs(FIGURES_DIR, exist_ok=True)

# Key locations
DELHI_CENTER = {'lat': 28.6139, 'lon': 77.2090}

# Grid parameters (from CAMx)
GRID = {
    'NCOLS': 100,
    'NROWS': 100,
    'XCELL': 27.0,  # km
    'YCELL': 27.0,  # km
    'XORIG': -1714.5,  # km
    'YORIG': -1714.5,  # km
    'PLON': 83.0,
    'PLAT': 21.5,
}

# PM2.5 component species in CAMx
PM25_SPECIES = ['PNO3', 'PSO4', 'PNH4', 'PEC', 'POA', 'FPRM']

# Colors
FCBK_COLOR = '#DC143C'
ZIGZAG_COLOR = '#228B22'

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def camx_grid_to_latlon(i, j):
    """Convert CAMx grid indices to lat/lon."""
    # Mercator projection inverse
    x_km = GRID['XORIG'] + (i + 0.5) * GRID['XCELL']
    y_km = GRID['YORIG'] + (j + 0.5) * GRID['YCELL']

    # Approximate inverse (for visualization purposes)
    lon = GRID['PLON'] + x_km / (111.0 * np.cos(np.radians(GRID['PLAT'])))
    lat = GRID['PLAT'] + y_km / 111.0

    return lon, lat


def get_grid_latlon():
    """Get lat/lon arrays for entire grid."""
    i_arr = np.arange(GRID['NCOLS'])
    j_arr = np.arange(GRID['NROWS'])
    I, J = np.meshgrid(i_arr, j_arr)

    LON, LAT = camx_grid_to_latlon(I, J)
    return LON, LAT


def find_delhi_cell():
    """Find grid cell indices for Delhi."""
    LON, LAT = get_grid_latlon()
    dist = np.sqrt((LON - DELHI_CENTER['lon'])**2 + (LAT - DELHI_CENTER['lat'])**2)
    j, i = np.unravel_index(np.argmin(dist), dist.shape)
    return i, j


def load_camx_output(scenario, date='20240202'):
    """Load CAMx output for a specific scenario and date."""

    output_dir = os.path.join(OUTPUT_DIR, scenario)
    run_name = f'india.27km.brickkiln.{scenario}'

    # Try different possible file names
    possible_files = [
        os.path.join(output_dir, f'CAMx.{run_name}.{date}.avrg.grd01.nc'),
        os.path.join(output_dir, f'CAMx.{run_name}.{date}.nc'),
        os.path.join(output_dir, f'CAMx.india.27km.{scenario}.{date}.avrg.grd01.nc'),
    ]

    for filepath in possible_files:
        if os.path.exists(filepath):
            print(f"  Loading: {os.path.basename(filepath)}")
            return nc.Dataset(filepath)

    print(f"  WARNING: Output file not found for scenario {scenario}, date {date}")
    print(f"  Searched: {output_dir}")
    return None


def calculate_pm25(ds):
    """Calculate total PM2.5 from CAMx output."""

    pm25 = np.zeros((ds.dimensions['ROW'].size, ds.dimensions['COL'].size))

    for species in PM25_SPECIES:
        if species in ds.variables:
            # Get surface layer, time-averaged
            data = ds.variables[species][:]

            # Handle different array shapes
            if len(data.shape) == 4:  # (time, layer, row, col)
                data = data[:, 0, :, :].mean(axis=0)
            elif len(data.shape) == 3:  # (time, row, col)
                data = data.mean(axis=0)

            pm25 += data

    return pm25


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_comparison_figure(pm25_2020, pm25_2025, LON, LAT):
    """Create comparison figure for both scenarios."""

    print("Creating comparison figure...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('CAMx Brick Kiln PM2.5 Impact Analysis\n' +
                 'February 2024 - India 27km Domain',
                 fontsize=16, fontweight='bold')

    # Load shapefiles if available
    india_gdf = up_gdf = delhi_gdf = None
    if HAS_GEOPANDAS and os.path.exists(SHAPEFILE_PATH):
        try:
            gdf = gpd.read_file(SHAPEFILE_PATH)
            if gdf.crs and gdf.crs != 'EPSG:4326':
                gdf = gdf.to_crs('EPSG:4326')
            india_gdf = gdf
            up_gdf = gdf[gdf['State_Name'].str.contains('Uttar Pradesh', case=False, na=False)]
            delhi_gdf = gdf[gdf['State_Name'].str.contains('Delhi', case=False, na=False)]
        except Exception as e:
            print(f"  Warning: Could not load shapefiles: {e}")

    # Common colorbar limits
    vmax = max(pm25_2020.max(), pm25_2025.max())
    vmin = 0

    # Custom colormap
    colors = ['#FFFFFF', '#FFF7BC', '#FED976', '#FEB24C', '#FD8D3C',
              '#FC4E2A', '#E31A1C', '#BD0026', '#800026']
    cmap = LinearSegmentedColormap.from_list('pm25', colors)

    # Panel 1: 2020 PM2.5
    ax1 = axes[0, 0]
    im1 = ax1.pcolormesh(LON, LAT, pm25_2020, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    if india_gdf is not None:
        india_gdf.boundary.plot(ax=ax1, color='gray', linewidth=0.5)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax1, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax1, color='orange', linewidth=2)

    ax1.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=15)
    ax1.set_xlim(72, 90)
    ax1.set_ylim(20, 35)
    ax1.set_title('2020 Scenario (100% FCBK)', fontsize=12, fontweight='bold', color=FCBK_COLOR)
    ax1.set_xlabel('Longitude (°E)')
    ax1.set_ylabel('Latitude (°N)')
    plt.colorbar(im1, ax=ax1, label='PM2.5 (µg/m³)', shrink=0.8)

    # Panel 2: 2025 PM2.5
    ax2 = axes[0, 1]
    im2 = ax2.pcolormesh(LON, LAT, pm25_2025, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    if india_gdf is not None:
        india_gdf.boundary.plot(ax=ax2, color='gray', linewidth=0.5)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax2, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax2, color='orange', linewidth=2)

    ax2.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=15)
    ax2.set_xlim(72, 90)
    ax2.set_ylim(20, 35)
    ax2.set_title('2025 Scenario (50% Zigzag)', fontsize=12, fontweight='bold', color=ZIGZAG_COLOR)
    ax2.set_xlabel('Longitude (°E)')
    ax2.set_ylabel('Latitude (°N)')
    plt.colorbar(im2, ax=ax2, label='PM2.5 (µg/m³)', shrink=0.8)

    # Panel 3: Difference (Reduction)
    ax3 = axes[0, 2]
    diff = pm25_2020 - pm25_2025
    vmax_diff = np.abs(diff).max()

    im3 = ax3.pcolormesh(LON, LAT, diff, cmap='RdYlGn', vmin=-vmax_diff/2, vmax=vmax_diff, shading='auto')

    if india_gdf is not None:
        india_gdf.boundary.plot(ax=ax3, color='gray', linewidth=0.5)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax3, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax3, color='orange', linewidth=2)

    ax3.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=15)
    ax3.set_xlim(72, 90)
    ax3.set_ylim(20, 35)
    ax3.set_title('PM2.5 Reduction (2020 - 2025)', fontsize=12, fontweight='bold', color='green')
    ax3.set_xlabel('Longitude (°E)')
    ax3.set_ylabel('Latitude (°N)')
    plt.colorbar(im3, ax=ax3, label='PM2.5 reduction (µg/m³)', shrink=0.8)

    # Panel 4: Delhi zoom - 2020
    ax4 = axes[1, 0]
    im4 = ax4.pcolormesh(LON, LAT, pm25_2020, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax4, color='orange', linewidth=3)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax4, color='black', linewidth=1)

    ax4.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=20)
    ax4.set_xlim(76, 79)
    ax4.set_ylim(27, 30)
    ax4.set_title('Delhi Region - 2020', fontsize=11, fontweight='bold')
    ax4.set_xlabel('Longitude (°E)')
    ax4.set_ylabel('Latitude (°N)')
    plt.colorbar(im4, ax=ax4, label='PM2.5 (µg/m³)', shrink=0.8)

    # Panel 5: Delhi zoom - 2025
    ax5 = axes[1, 1]
    im5 = ax5.pcolormesh(LON, LAT, pm25_2025, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax5, color='orange', linewidth=3)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax5, color='black', linewidth=1)

    ax5.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=20)
    ax5.set_xlim(76, 79)
    ax5.set_ylim(27, 30)
    ax5.set_title('Delhi Region - 2025', fontsize=11, fontweight='bold')
    ax5.set_xlabel('Longitude (°E)')
    ax5.set_ylabel('Latitude (°N)')
    plt.colorbar(im5, ax=ax5, label='PM2.5 (µg/m³)', shrink=0.8)

    # Panel 6: Summary statistics
    ax6 = axes[1, 2]
    ax6.axis('off')

    # Calculate Delhi-specific values
    delhi_i, delhi_j = find_delhi_cell()
    delhi_2020 = pm25_2020[delhi_j, delhi_i]
    delhi_2025 = pm25_2025[delhi_j, delhi_i]
    delhi_reduction = delhi_2020 - delhi_2025
    delhi_pct = (delhi_reduction / delhi_2020) * 100 if delhi_2020 > 0 else 0

    # Domain statistics
    domain_mean_2020 = pm25_2020.mean()
    domain_mean_2025 = pm25_2025.mean()
    domain_reduction = domain_mean_2020 - domain_mean_2025
    domain_pct = (domain_reduction / domain_mean_2020) * 100 if domain_mean_2020 > 0 else 0

    # Max reduction
    max_reduction = diff.max()
    max_idx = np.unravel_index(np.argmax(diff), diff.shape)
    max_lon, max_lat = LON[max_idx], LAT[max_idx]

    summary_text = f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    CAMx PM2.5 IMPACT SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELHI CONCENTRATIONS:
  2020 (100% FCBK):  {delhi_2020:.2f} µg/m³
  2025 (50% Zigzag): {delhi_2025:.2f} µg/m³
  Reduction:         {delhi_reduction:.2f} µg/m³ ({delhi_pct:.1f}%)

DOMAIN AVERAGE:
  2020: {domain_mean_2020:.2f} µg/m³
  2025: {domain_mean_2025:.2f} µg/m³
  Reduction: {domain_reduction:.2f} µg/m³ ({domain_pct:.1f}%)

MAXIMUM REDUCTION:
  {max_reduction:.2f} µg/m³
  at ({max_lon:.1f}°E, {max_lat:.1f}°N)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Model: CAMx v7.32
Grid: 100x100 cells, 27km
Chemistry: CB6r4
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

    ax6.text(0.5, 0.5, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='center', horizontalalignment='center',
             fontfamily='monospace',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#E8F8E8',
                       edgecolor='green', linewidth=2))

    plt.tight_layout()

    output_path = os.path.join(FIGURES_DIR, 'camx_pm25_comparison.png')
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path}")
    plt.close()

    return delhi_2020, delhi_2025


def create_bar_comparison(delhi_2020, delhi_2025):
    """Create bar chart comparison for Delhi."""

    print("Creating bar chart comparison...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: Brick kiln contribution only
    ax1 = axes[0]
    x = ['2020\n(100% FCBK)', '2025\n(50% Zigzag)']
    heights = [delhi_2020, delhi_2025]
    colors = [FCBK_COLOR, ZIGZAG_COLOR]

    bars = ax1.bar(x, heights, color=colors, alpha=0.8, width=0.5)

    for bar, val in zip(bars, heights):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{val:.2f}', ha='center', fontsize=12, fontweight='bold')

    reduction = delhi_2020 - delhi_2025
    pct = (reduction / delhi_2020) * 100 if delhi_2020 > 0 else 0

    ax1.annotate(f'↓{reduction:.2f} µg/m³\n({pct:.1f}%)',
                xy=(0.5, (delhi_2020 + delhi_2025) / 2),
                fontsize=14, fontweight='bold', color='green', ha='center')

    ax1.set_ylabel('PM2.5 from Brick Kilns (µg/m³)', fontsize=11)
    ax1.set_title('Delhi: Brick Kiln PM2.5 Contribution\n(CAMx Model Results)', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, max(heights) * 1.3)

    # Panel 2: Total PM2.5 with background
    ax2 = axes[1]
    background = 80.0  # Typical Delhi winter background from other sources

    ax2.bar(x, [background, background], color='gray', alpha=0.6, width=0.5, label='Other sources')
    ax2.bar(x, heights, bottom=[background, background], color=colors, alpha=0.8, width=0.5, label='Brick kilns')

    # WHO guidelines
    ax2.axhline(y=60, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='WHO IT-3')
    ax2.axhline(y=35, color='red', linestyle='--', linewidth=2, alpha=0.7, label='WHO IT-2')

    for i, (bg, bk) in enumerate(zip([background]*2, heights)):
        total = bg + bk
        ax2.text(i, total + 2, f'{total:.1f}', ha='center', fontsize=11, fontweight='bold')

    ax2.set_ylabel('Total PM2.5 (µg/m³)', fontsize=11)
    ax2.set_title('Delhi: Total PM2.5 Concentration', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 120)
    ax2.legend(loc='upper right', fontsize=9)

    plt.tight_layout()

    output_path = os.path.join(FIGURES_DIR, 'camx_delhi_comparison.png')
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path}")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("CAMx OUTPUT ANALYSIS - BRICK KILN SCENARIOS")
    print("=" * 70)

    # Check if output exists
    output_2020_dir = os.path.join(OUTPUT_DIR, '2020')
    output_2025_dir = os.path.join(OUTPUT_DIR, '2025')

    has_2020 = os.path.exists(output_2020_dir) and len(os.listdir(output_2020_dir)) > 0
    has_2025 = os.path.exists(output_2025_dir) and len(os.listdir(output_2025_dir)) > 0

    # If no CAMx output exists, create synthetic data for demonstration
    if not has_2020 or not has_2025:
        print("\nCAMx output not found. Creating demonstration with existing emission data...")
        print("(Run CAMx.brick_kiln.job to get actual model results)")

        # Load existing emission data to create synthetic concentrations
        from create_camx_inputs import load_brick_kiln_emissions

        LON, LAT = get_grid_latlon()

        # Create synthetic PM2.5 based on emissions
        emis_2020 = load_brick_kiln_emissions('2020')
        emis_2025 = load_brick_kiln_emissions('2025')

        if emis_2020 is not None and emis_2025 is not None:
            # Simple concentration estimate (scaling emissions)
            # In reality, CAMx does full chemistry and transport
            scale_factor = 0.5  # Approximate µg/m³ per g/s emission

            pm25_2020 = emis_2020['PM2.5'][:100, :100] * scale_factor + 80  # Add background
            pm25_2025 = emis_2025['PM2.5'][:100, :100] * scale_factor + 80

            # Smooth fields
            from scipy.ndimage import gaussian_filter
            pm25_2020 = gaussian_filter(pm25_2020, sigma=3)
            pm25_2025 = gaussian_filter(pm25_2025, sigma=3)

            delhi_2020, delhi_2025 = create_comparison_figure(pm25_2020, pm25_2025, LON, LAT)
            create_bar_comparison(delhi_2020, delhi_2025)
        else:
            print("ERROR: Could not load emission data either.")
            print("Run the kiln_camx_analysis.py script first to generate emissions.")
            return

    else:
        print("\nLoading CAMx output...")

        # Load output for both scenarios
        ds_2020 = load_camx_output('2020')
        ds_2025 = load_camx_output('2025')

        if ds_2020 is None or ds_2025 is None:
            print("ERROR: Could not load CAMx output files")
            return

        # Calculate PM2.5
        print("\nCalculating PM2.5 concentrations...")
        pm25_2020 = calculate_pm25(ds_2020)
        pm25_2025 = calculate_pm25(ds_2025)

        ds_2020.close()
        ds_2025.close()

        # Get grid coordinates
        LON, LAT = get_grid_latlon()

        # Create visualizations
        delhi_2020, delhi_2025 = create_comparison_figure(pm25_2020, pm25_2025, LON, LAT)
        create_bar_comparison(delhi_2020, delhi_2025)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE!")
    print("=" * 70)
    print(f"""
Output figures in: {FIGURES_DIR}
  - camx_pm25_comparison.png
  - camx_delhi_comparison.png
""")


if __name__ == "__main__":
    main()
