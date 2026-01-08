#!/usr/bin/env python3
"""
Visualize CAMx Emission Files
=============================

Creates publication-quality visualizations of:
1. Gridded emission NetCDF files
2. Point source emission files
3. Comparison between 2020 and 2025 scenarios
"""

import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import netCDF4 as nc
import geopandas as gpd

# Paths
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHAPEFILE_PATH = os.path.join(PROJECT_ROOT, 'shapefiles', 'India_State_Boundary.shp')
CAMX_EMIS_DIR = os.path.join(PROJECT_ROOT, 'camx_emissions')
FIGURES_DIR = os.path.join(PROJECT_ROOT, 'figures')

# Colors
FCBK_COLOR = '#DC143C'
ZIGZAG_COLOR = '#228B22'

# Delhi location
DELHI = {'lon': 77.2090, 'lat': 28.6139}


def load_shapefiles():
    """Load UP and Delhi boundaries."""
    if not os.path.exists(SHAPEFILE_PATH):
        return None, None

    gdf = gpd.read_file(SHAPEFILE_PATH)
    if gdf.crs and gdf.crs != 'EPSG:4326':
        gdf = gdf.to_crs('EPSG:4326')

    up_gdf = gdf[gdf['State_Name'].str.contains('Uttar Pradesh', case=False, na=False)]
    delhi_gdf = gdf[gdf['State_Name'].str.contains('Delhi', case=False, na=False)]

    return up_gdf, delhi_gdf


def load_nc_emissions(filepath):
    """Load gridded emissions from NetCDF file."""
    ds = nc.Dataset(filepath)

    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    pm25 = ds.variables['PM2.5'][:]
    pm10 = ds.variables['PM10'][:]
    nox = ds.variables['NOx'][:]
    so2 = ds.variables['SO2'][:]

    ds.close()

    return {'x': x, 'y': y, 'PM2.5': pm25, 'PM10': pm10, 'NOx': nox, 'SO2': so2}


def load_point_source(filepath):
    """Load point source emissions from text file."""
    kilns = []

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('!') or line.startswith('#') or '=' in line:
                continue
            parts = line.split()
            if len(parts) >= 15:
                try:
                    kilns.append({
                        'id': int(parts[0]),
                        'lon': float(parts[1]),
                        'lat': float(parts[2]),
                        'pm25_gs': float(parts[9]),
                        'type': parts[14] if len(parts) > 14 else 'FCBK'
                    })
                except (ValueError, IndexError):
                    continue

    return kilns


def visualize_gridded_emissions():
    """Create gridded emission comparison figure."""
    print("Visualizing gridded emissions...")

    # Load data
    emis_2020 = load_nc_emissions(os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2020.nc'))
    emis_2025 = load_nc_emissions(os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2025.nc'))
    up_gdf, delhi_gdf = load_shapefiles()

    # Convert CAMx grid (km) to approximate lat/lon for visualization
    # Reference point: 83°E, 21.5°N
    ref_lon, ref_lat = 83.0, 21.5
    lon_grid = ref_lon + emis_2020['x'] / (111.0 * np.cos(np.radians(ref_lat)))
    lat_grid = ref_lat + emis_2020['y'] / 111.0
    LON, LAT = np.meshgrid(lon_grid, lat_grid)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.suptitle('Brick Kiln Gridded Emissions (CAMx Format)\nPM2.5 Emission Rate (g/s per cell)',
                fontsize=16, fontweight='bold')

    # PM2.5 2020
    ax = axes[0, 0]
    vmax = max(emis_2020['PM2.5'].max(), emis_2025['PM2.5'].max())
    im = ax.pcolormesh(LON, LAT, emis_2020['PM2.5'], cmap='YlOrRd', shading='auto', vmin=0, vmax=vmax)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax, color='blue', linewidth=2)
    ax.plot(DELHI['lon'], DELHI['lat'], 'k*', markersize=15)
    ax.set_title(f"2020 (100% FCBK)\nTotal: {emis_2020['PM2.5'].sum():.0f} g/s", fontsize=12, fontweight='bold')
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_xlim(75, 86)
    ax.set_ylim(23, 32)
    plt.colorbar(im, ax=ax, label='g/s')

    # PM2.5 2025
    ax = axes[0, 1]
    im = ax.pcolormesh(LON, LAT, emis_2025['PM2.5'], cmap='YlOrRd', shading='auto', vmin=0, vmax=vmax)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax, color='blue', linewidth=2)
    ax.plot(DELHI['lon'], DELHI['lat'], 'k*', markersize=15)
    ax.set_title(f"2025 (50% Zigzag)\nTotal: {emis_2025['PM2.5'].sum():.0f} g/s", fontsize=12, fontweight='bold')
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_xlim(75, 86)
    ax.set_ylim(23, 32)
    plt.colorbar(im, ax=ax, label='g/s')

    # Emission reduction
    ax = axes[1, 0]
    reduction = emis_2020['PM2.5'] - emis_2025['PM2.5']
    im = ax.pcolormesh(LON, LAT, reduction, cmap='Greens', shading='auto', vmin=0)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax, color='blue', linewidth=2)
    ax.plot(DELHI['lon'], DELHI['lat'], 'k*', markersize=15)
    ax.set_title(f"Emission Reduction\nTotal: {reduction.sum():.0f} g/s ({reduction.sum()/emis_2020['PM2.5'].sum()*100:.0f}%)",
                fontsize=12, fontweight='bold', color='green')
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_xlim(75, 86)
    ax.set_ylim(23, 32)
    plt.colorbar(im, ax=ax, label='g/s reduction')

    # Bar chart comparison
    ax = axes[1, 1]
    species = ['PM2.5', 'PM10', 'NOx', 'SO2']
    vals_2020 = [emis_2020[s].sum() for s in species]
    vals_2025 = [emis_2025[s].sum() for s in species]

    x = np.arange(len(species))
    width = 0.35

    bars1 = ax.bar(x - width/2, vals_2020, width, label='2020 (FCBK)', color=FCBK_COLOR, alpha=0.8)
    bars2 = ax.bar(x + width/2, vals_2025, width, label='2025 (Mixed)', color=ZIGZAG_COLOR, alpha=0.8)

    ax.set_ylabel('Total Emission Rate (g/s)')
    ax.set_xticks(x)
    ax.set_xticklabels(species)
    ax.set_title('Total Emissions by Species', fontsize=12, fontweight='bold')
    ax.legend()

    # Add reduction percentages
    for i, (v1, v2) in enumerate(zip(vals_2020, vals_2025)):
        pct = (v1 - v2) / v1 * 100
        ax.text(i, max(v1, v2) + 100, f'-{pct:.0f}%', ha='center', fontsize=10, color='green', fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_path = os.path.join(FIGURES_DIR, 'gridded_emissions.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path}")
    plt.close()


def visualize_point_sources():
    """Create point source emission comparison figure."""
    print("Visualizing point source emissions...")

    # Load data
    kilns_2020 = load_point_source(os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2020.ptsrc'))
    kilns_2025 = load_point_source(os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2025.ptsrc'))
    up_gdf, delhi_gdf = load_shapefiles()

    fig, axes = plt.subplots(1, 2, figsize=(18, 10))
    fig.suptitle('Brick Kiln Point Source Locations\n10,000 Kilns in Uttar Pradesh',
                fontsize=16, fontweight='bold')

    for ax, kilns, title, scenario in [
        (axes[0], kilns_2020, '2020: All FCBK (Red)', '2020'),
        (axes[1], kilns_2025, '2025: 50% Converted (Red=FCBK, Green=Zigzag)', '2025')
    ]:
        ax.set_facecolor('#E8F4F8')

        if up_gdf is not None:
            up_gdf.plot(ax=ax, facecolor='#E6E6FA', edgecolor='#4a4a4a', linewidth=2, alpha=0.8)
        if delhi_gdf is not None:
            delhi_gdf.plot(ax=ax, facecolor='#FFE4B5', edgecolor='#FF8C00', linewidth=3, alpha=0.9)

        # Separate by type
        fcbk = [k for k in kilns if k['type'] == 'FCBK']
        zigzag = [k for k in kilns if k['type'] == 'Zigzag']

        if fcbk:
            ax.scatter([k['lon'] for k in fcbk], [k['lat'] for k in fcbk],
                      c=FCBK_COLOR, s=5, alpha=0.5, marker='o', label=f'FCBK ({len(fcbk):,})')
        if zigzag:
            ax.scatter([k['lon'] for k in zigzag], [k['lat'] for k in zigzag],
                      c=ZIGZAG_COLOR, s=6, alpha=0.6, marker='s', label=f'Zigzag ({len(zigzag):,})')

        ax.plot(DELHI['lon'], DELHI['lat'], '*', color='gold', markersize=20,
               markeredgecolor='black', markeredgewidth=2, label='Delhi')

        ax.set_xlim(75, 86)
        ax.set_ylim(23, 32)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='lower right')

        # Add emission total
        total_pm25 = sum(k['pm25_gs'] for k in kilns)
        ax.text(0.02, 0.98, f'Total PM2.5: {total_pm25:.0f} g/s',
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output_path = os.path.join(FIGURES_DIR, 'point_source_locations.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path}")
    plt.close()


def create_summary_table():
    """Create summary statistics."""
    print("\nEmission Summary:")
    print("=" * 60)

    emis_2020 = load_nc_emissions(os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2020.nc'))
    emis_2025 = load_nc_emissions(os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2025.nc'))

    print(f"{'Species':<10} {'2020 (g/s)':<15} {'2025 (g/s)':<15} {'Reduction':<15}")
    print("-" * 60)

    for species in ['PM2.5', 'PM10', 'NOx', 'SO2']:
        v1 = emis_2020[species].sum()
        v2 = emis_2025[species].sum()
        red = v1 - v2
        pct = red / v1 * 100
        print(f"{species:<10} {v1:<15.0f} {v2:<15.0f} {red:.0f} ({pct:.0f}%)")

    print("=" * 60)

    # Convert to tonnes/year
    operating_hours = 180 * 24  # 180 days
    seconds_per_year = operating_hours * 3600

    print(f"\nAnnual Emissions (tonnes/year, {operating_hours/24:.0f} operating days):")
    print("-" * 60)

    for species in ['PM2.5', 'PM10']:
        v1 = emis_2020[species].sum() * seconds_per_year / 1e6
        v2 = emis_2025[species].sum() * seconds_per_year / 1e6
        red = v1 - v2
        print(f"{species:<10} {v1:<15,.0f} {v2:<15,.0f} {red:,.0f}")


if __name__ == "__main__":
    print("=" * 60)
    print("VISUALIZING CAMx EMISSION FILES")
    print("=" * 60)

    visualize_gridded_emissions()
    visualize_point_sources()
    create_summary_table()

    print("\n" + "=" * 60)
    print("COMPLETE!")
    print("=" * 60)
