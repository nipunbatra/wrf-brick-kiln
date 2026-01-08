#!/usr/bin/env python3
"""
Brick Kiln Impact Analysis using WRF-CAMx
==========================================

This script:
1. Uses actual India state shapefiles (UP, Delhi, etc.)
2. Generates 10,000 kilns randomly distributed in UP
3. Creates CAMx-format emission files for:
   - Scenario 2020: 100% FCBK (Fixed Chimney Bull Trench Kiln)
   - Scenario 2025: 50% FCBK + 50% Zigzag
4. Uses WRF-CAMx meteorology for proper wind transport
5. Visualizes kiln distribution and emission reduction

The actual PM2.5 concentrations come from running CAMx with these emissions.
"""

import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
from matplotlib.lines import Line2D
from shapely.geometry import Point
import netCDF4 as nc

# =============================================================================
# PATHS
# =============================================================================
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHAPEFILE_PATH = os.path.join(PROJECT_ROOT, 'shapefiles', 'India_State_Boundary.shp')
CAMX_MET_PATH = os.path.join(PROJECT_ROOT, 'wrfcamx_v5.2', 'camx_input')
FIGURES_DIR = os.path.join(PROJECT_ROOT, 'figures')
CAMX_EMIS_DIR = os.path.join(PROJECT_ROOT, 'camx_emissions')
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(CAMX_EMIS_DIR, exist_ok=True)

# Random seed
np.random.seed(42)

# =============================================================================
# CONFIGURATION
# =============================================================================
N_KILNS = 10000

# Key locations
DELHI_CENTER = {'lat': 28.6139, 'lon': 77.2090}

# Emission factors (g/kg brick) - from CPCB India / GIZ studies
EMISSION_FACTORS = {
    'FCBK': {'PM2.5': 0.75, 'PM10': 1.2, 'BC': 0.15, 'NOx': 0.30, 'SO2': 0.80, 'CO': 15.0},
    'Zigzag': {'PM2.5': 0.30, 'PM10': 0.48, 'BC': 0.06, 'NOx': 0.25, 'SO2': 0.65, 'CO': 10.0},
}

# Production parameters
BRICKS_PER_KILN_PER_DAY = 25000
OPERATING_DAYS_PER_YEAR = 180
BRICK_WEIGHT_KG = 3.0

# Stack parameters for CAMx
STACK_HEIGHT = 15.0  # meters
STACK_DIAMETER = 2.0  # meters
STACK_TEMP = 400.0  # Kelvin
STACK_VELOCITY = 5.0  # m/s

# Colors
FCBK_COLOR = '#DC143C'
ZIGZAG_COLOR = '#228B22'
DELHI_COLOR = '#FFD700'
UP_COLOR = '#E6E6FA'

# CAMx grid parameters (from WRFCAMx output)
CAMX_GRID = {
    'nx': 102,
    'ny': 102,
    'dx': 27.0,  # km
    'dy': 27.0,  # km
    'xorg': -1714.5,  # km (Mercator)
    'yorg': -1714.5,  # km (Mercator)
    'ref_lon': 83.0,
    'ref_lat': 21.5,
}

# =============================================================================
# LOAD SHAPEFILES
# =============================================================================
def load_india_shapefiles():
    """Load India state boundaries from shapefile."""
    print("Loading India shapefiles...")

    if not os.path.exists(SHAPEFILE_PATH):
        print(f"  ERROR: Shapefile not found at {SHAPEFILE_PATH}")
        return None, None, None

    gdf = gpd.read_file(SHAPEFILE_PATH)

    # The shapefile is in a projected CRS, convert to WGS84
    if gdf.crs and gdf.crs != 'EPSG:4326':
        gdf = gdf.to_crs('EPSG:4326')

    # Extract UP and Delhi
    up_gdf = gdf[gdf['State_Name'].str.contains('Uttar Pradesh', case=False, na=False)]
    delhi_gdf = gdf[gdf['State_Name'].str.contains('Delhi', case=False, na=False)]

    print(f"  Found {len(up_gdf)} UP polygon(s)")
    print(f"  Found {len(delhi_gdf)} Delhi polygon(s)")

    return gdf, up_gdf, delhi_gdf


# =============================================================================
# READ WRF-CAMX METEOROLOGY
# =============================================================================
def load_camx_meteorology():
    """Load wind fields from CAMx meteorology files."""
    print("Loading WRF-CAMx meteorology...")

    met_file = os.path.join(CAMX_MET_PATH, 'camx.3d.4km.20240201.nc')

    if not os.path.exists(met_file):
        print(f"  WARNING: Met file not found: {met_file}")
        return None

    ds = nc.Dataset(met_file)

    # Get coordinates
    lon = ds.variables['longitude'][:]
    lat = ds.variables['latitude'][:]

    # Get surface winds (layer 0, time-averaged)
    uwind = ds.variables['uwind'][:, 0, :, :].mean(axis=0)  # m/s
    vwind = ds.variables['vwind'][:, 0, :, :].mean(axis=0)  # m/s

    # Wind speed and direction
    wspd = np.sqrt(uwind**2 + vwind**2)
    wdir = np.arctan2(vwind, uwind) * 180 / np.pi  # degrees from east

    ds.close()

    print(f"  Grid: {lon.shape}")
    print(f"  Mean wind speed: {wspd.mean():.1f} m/s")
    print(f"  Dominant wind from: {90 - wdir.mean():.0f}° (meteorological)")

    return {
        'lon': lon,
        'lat': lat,
        'uwind': uwind,
        'vwind': vwind,
        'wspd': wspd,
        'wdir': wdir,
    }


# =============================================================================
# GENERATE KILNS
# =============================================================================
def generate_kilns_in_up(up_gdf, n_kilns=10000):
    """Generate kilns within UP boundary using clustering."""
    print(f"Generating {n_kilns:,} kilns in UP...")

    # Brick kiln cluster centers in UP
    clusters = [
        {'name': 'Western UP (Ghaziabad)', 'lon': 77.8, 'lat': 28.5, 'weight': 0.25, 'std': 0.5},
        {'name': 'Noida-Greater Noida', 'lon': 77.4, 'lat': 28.5, 'weight': 0.20, 'std': 0.3},
        {'name': 'Meerut-Muzaffarnagar', 'lon': 77.7, 'lat': 29.1, 'weight': 0.15, 'std': 0.4},
        {'name': 'Agra-Mathura', 'lon': 78.0, 'lat': 27.2, 'weight': 0.12, 'std': 0.5},
        {'name': 'Lucknow-Kanpur', 'lon': 80.5, 'lat': 26.8, 'weight': 0.12, 'std': 0.6},
        {'name': 'Varanasi-Allahabad', 'lon': 82.5, 'lat': 25.5, 'weight': 0.10, 'std': 0.7},
        {'name': 'Eastern UP', 'lon': 83.5, 'lat': 26.5, 'weight': 0.06, 'std': 0.5},
    ]

    up_geometry = up_gdf.geometry.values[0] if len(up_gdf) > 0 else None
    minx, miny, maxx, maxy = up_gdf.total_bounds if len(up_gdf) > 0 else (77, 24, 85, 31)

    kilns = []

    for cluster in clusters:
        n_cluster = int(n_kilns * cluster['weight'])
        generated = 0
        attempts = 0

        while generated < n_cluster and attempts < n_cluster * 20:
            lon = np.random.normal(cluster['lon'], cluster['std'])
            lat = np.random.normal(cluster['lat'], cluster['std'] * 0.7)

            point = Point(lon, lat)

            # Check if within UP boundary
            if up_geometry is not None and up_geometry.contains(point):
                kilns.append({
                    'id': len(kilns),
                    'lon': lon,
                    'lat': lat,
                    'cluster': cluster['name']
                })
                generated += 1
            elif up_geometry is None and minx <= lon <= maxx and miny <= lat <= maxy:
                kilns.append({
                    'id': len(kilns),
                    'lon': lon,
                    'lat': lat,
                    'cluster': cluster['name']
                })
                generated += 1

            attempts += 1

    # Fill remaining with random
    while len(kilns) < n_kilns:
        lon = np.random.uniform(77.5, 84)
        lat = np.random.uniform(24.5, 30)
        point = Point(lon, lat)

        if up_geometry is None or up_geometry.contains(point):
            kilns.append({
                'id': len(kilns),
                'lon': lon,
                'lat': lat,
                'cluster': 'Random'
            })

    # Print distribution
    cluster_counts = {}
    for k in kilns:
        c = k['cluster']
        cluster_counts[c] = cluster_counts.get(c, 0) + 1

    print("  Kiln distribution:")
    for cluster, count in sorted(cluster_counts.items(), key=lambda x: -x[1]):
        print(f"    {cluster}: {count:,}")

    return kilns


# =============================================================================
# CALCULATE EMISSIONS
# =============================================================================
def calculate_emissions(kilns, kiln_types):
    """Calculate emissions for each kiln."""

    emissions = []

    for i, kiln in enumerate(kilns):
        ktype = kiln_types[i]
        ef = EMISSION_FACTORS[ktype]

        # Annual brick production
        annual_bricks = BRICKS_PER_KILN_PER_DAY * OPERATING_DAYS_PER_YEAR
        annual_mass_kg = annual_bricks * BRICK_WEIGHT_KG

        # Operating time in seconds
        operating_seconds = OPERATING_DAYS_PER_YEAR * 24 * 3600

        # Emission rates
        kiln_emis = {
            'id': kiln['id'],
            'lon': kiln['lon'],
            'lat': kiln['lat'],
            'type': ktype,
            'cluster': kiln['cluster'],
        }

        for species, ef_val in ef.items():
            # tonnes/year
            kiln_emis[f'{species}_tpy'] = (annual_mass_kg * ef_val) / 1e6
            # g/s (for CAMx)
            kiln_emis[f'{species}_gs'] = (annual_mass_kg * ef_val) / operating_seconds

        emissions.append(kiln_emis)

    return emissions


def latlon_to_camx_grid(lon, lat):
    """Convert lat/lon to CAMx grid coordinates (km from origin)."""
    # Approximate conversion using reference point
    ref_lon = CAMX_GRID['ref_lon']
    ref_lat = CAMX_GRID['ref_lat']

    # km per degree (approximate at this latitude)
    km_per_deg_lon = 111.0 * np.cos(np.radians(lat))
    km_per_deg_lat = 111.0

    x_km = (lon - ref_lon) * km_per_deg_lon
    y_km = (lat - ref_lat) * km_per_deg_lat

    return x_km, y_km


def calculate_concentration_field(emissions, met_data, nx=100, ny=80):
    """
    Calculate PM2.5 concentration field using CAMx meteorology.
    Uses simplified Gaussian plume with wind transport from WRF-CAMx.
    Returns concentration in µg/m³.

    This is a screening-level model. Full CAMx run gives more accurate results.
    """
    from scipy.ndimage import gaussian_filter

    # Create grid covering UP-Delhi region
    lon_min, lon_max = 75.0, 86.0
    lat_min, lat_max = 23.0, 32.0

    lon_grid = np.linspace(lon_min, lon_max, nx)
    lat_grid = np.linspace(lat_min, lat_max, ny)
    LON, LAT = np.meshgrid(lon_grid, lat_grid)

    # Cell size in meters
    dx_deg = (lon_max - lon_min) / nx
    dy_deg = (lat_max - lat_min) / ny
    dx_m = dx_deg * 111000 * np.cos(np.radians(27.5))  # ~27.5°N average
    dy_m = dy_deg * 111000

    # Get mean wind from meteorology
    if met_data is not None:
        mean_wspd = met_data['wspd'].mean()
        mean_u = met_data['uwind'].mean()
        mean_v = met_data['vwind'].mean()
    else:
        mean_wspd = 3.0
        mean_u = 2.0
        mean_v = -1.0

    # Initialize concentration grid
    concentration = np.zeros((ny, nx))

    # Gaussian plume parameters (winter, stable conditions)
    H = 15.0  # stack height (m)
    mixing_height = 500.0  # winter mixing height in Delhi region (m)
    sigma_y_coef = 0.22  # lateral dispersion coefficient
    sigma_z_coef = 0.08  # vertical dispersion coefficient

    # For each emission source, calculate contribution using Gaussian plume
    for e in emissions:
        Q = e['PM2.5_gs']  # g/s
        src_lon, src_lat = e['lon'], e['lat']

        # Distance from source to each grid point
        # Convert to meters
        dx_from_src = (LON - src_lon) * 111000 * np.cos(np.radians(src_lat))
        dy_from_src = (LAT - src_lat) * 111000

        # Downwind distance (along wind direction)
        wind_dir = np.arctan2(mean_v, mean_u)  # radians
        x_downwind = dx_from_src * np.cos(wind_dir) + dy_from_src * np.sin(wind_dir)
        y_crosswind = -dx_from_src * np.sin(wind_dir) + dy_from_src * np.cos(wind_dir)

        # Only calculate for downwind points
        downwind_mask = x_downwind > 1000  # at least 1 km downwind

        # Dispersion coefficients (Pasquill-Gifford, stability class D)
        x_km = np.maximum(x_downwind / 1000, 0.1)
        sigma_y = sigma_y_coef * x_km * 1000 * (1 + 0.0001 * x_km * 1000) ** (-0.5)
        sigma_z = sigma_z_coef * x_km * 1000 * (1 + 0.0002 * x_km * 1000) ** (-0.5)
        sigma_z = np.minimum(sigma_z, mixing_height)  # cap at mixing height

        # Gaussian plume formula (ground-level, with reflection)
        # C = Q / (2 * pi * u * sigma_y * sigma_z) * exp(-y²/2σy²) * 2*exp(-H²/2σz²)
        with np.errstate(divide='ignore', invalid='ignore'):
            lateral = np.exp(-0.5 * (y_crosswind / sigma_y) ** 2)
            vertical = 2 * np.exp(-0.5 * (H / sigma_z) ** 2)  # reflection term

            # Concentration in g/m³, convert to µg/m³
            C = (Q / (2 * np.pi * mean_wspd * sigma_y * sigma_z)) * lateral * vertical
            C = np.where(downwind_mask & np.isfinite(C), C * 1e6, 0)  # µg/m³

        concentration += C

    # Apply additional smoothing for grid effects
    concentration = gaussian_filter(concentration, sigma=1.5)

    return LON, LAT, concentration


def get_delhi_concentration(LON, LAT, concentration):
    """Extract PM2.5 concentration at Delhi."""
    # Find grid cell closest to Delhi
    delhi_lon, delhi_lat = DELHI_CENTER['lon'], DELHI_CENTER['lat']

    dist = np.sqrt((LON - delhi_lon)**2 + (LAT - delhi_lat)**2)
    j, i = np.unravel_index(np.argmin(dist), dist.shape)

    # Average over 3x3 neighborhood
    j_min, j_max = max(0, j-1), min(LON.shape[0], j+2)
    i_min, i_max = max(0, i-1), min(LON.shape[1], i+2)

    return concentration[j_min:j_max, i_min:i_max].mean()


# =============================================================================
# CREATE CAMX EMISSION FILES
# =============================================================================
def create_camx_point_source_file(emissions, output_path, scenario_name):
    """Create CAMx-format point source emission file."""
    print(f"  Creating CAMx point source file: {os.path.basename(output_path)}")

    with open(output_path, 'w') as f:
        # Header
        f.write(f"! CAMx Point Source File - Brick Kilns in Uttar Pradesh\n")
        f.write(f"! Scenario: {scenario_name}\n")
        f.write(f"! Generated by WRF-CAMx brick kiln analysis\n")
        f.write(f"!\n")
        f.write(f"! Grid: {CAMX_GRID['nx']} x {CAMX_GRID['ny']} cells at {CAMX_GRID['dx']} km\n")
        f.write(f"! Reference: {CAMX_GRID['ref_lon']}E, {CAMX_GRID['ref_lat']}N\n")
        f.write(f"!\n")
        f.write(f"NSTK = {len(emissions)}\n")
        f.write(f"NCOL = {CAMX_GRID['nx']}  NROW = {CAMX_GRID['ny']}\n")
        f.write(f"DX = {CAMX_GRID['dx']:.1f}   DY = {CAMX_GRID['dy']:.1f}\n")
        f.write(f"XORG = {CAMX_GRID['xorg']:.1f}  YORG = {CAMX_GRID['yorg']:.1f}\n")
        f.write(f"!\n")
        f.write(f"! ID      LON       LAT      X_KM     Y_KM    HT   DIA   TEMP  VEL    PM25_G/S   PM10_G/S   NOX_G/S   SO2_G/S    CO_G/S    TYPE\n")
        f.write(f"!\n")

        # Data
        for e in emissions:
            x_km, y_km = latlon_to_camx_grid(e['lon'], e['lat'])

            f.write(f"{e['id']:6d}  {e['lon']:8.4f}  {e['lat']:7.4f}  "
                   f"{x_km:8.1f} {y_km:8.1f}  "
                   f"{STACK_HEIGHT:4.0f}  {STACK_DIAMETER:4.1f}  {STACK_TEMP:5.0f}  {STACK_VELOCITY:4.1f}  "
                   f"{e['PM2.5_gs']:10.3E}  {e['PM10_gs']:10.3E}  "
                   f"{e['NOx_gs']:10.3E}  {e['SO2_gs']:10.3E}  {e['CO_gs']:10.3E}  "
                   f"{e['type']}\n")

    return output_path


def create_camx_gridded_emission_file(emissions, output_path, scenario_name):
    """Create CAMx-format gridded emission NetCDF file."""
    print(f"  Creating CAMx gridded emission file: {os.path.basename(output_path)}")

    nx, ny = CAMX_GRID['nx'], CAMX_GRID['ny']

    # Initialize grids
    pm25_grid = np.zeros((ny, nx))
    pm10_grid = np.zeros((ny, nx))
    nox_grid = np.zeros((ny, nx))
    so2_grid = np.zeros((ny, nx))

    # Create coordinate arrays
    x = np.linspace(CAMX_GRID['xorg'], CAMX_GRID['xorg'] + nx * CAMX_GRID['dx'], nx)
    y = np.linspace(CAMX_GRID['yorg'], CAMX_GRID['yorg'] + ny * CAMX_GRID['dy'], ny)

    # Bin emissions to grid
    for e in emissions:
        x_km, y_km = latlon_to_camx_grid(e['lon'], e['lat'])

        # Find grid cell
        i = int((x_km - CAMX_GRID['xorg']) / CAMX_GRID['dx'])
        j = int((y_km - CAMX_GRID['yorg']) / CAMX_GRID['dy'])

        if 0 <= i < nx and 0 <= j < ny:
            pm25_grid[j, i] += e['PM2.5_gs']
            pm10_grid[j, i] += e['PM10_gs']
            nox_grid[j, i] += e['NOx_gs']
            so2_grid[j, i] += e['SO2_gs']

    # Create NetCDF file
    ds = nc.Dataset(output_path, 'w', format='NETCDF4')

    # Dimensions
    ds.createDimension('x', nx)
    ds.createDimension('y', ny)
    ds.createDimension('time', None)

    # Coordinates
    xvar = ds.createVariable('x', 'f4', ('x',))
    xvar[:] = x
    xvar.units = 'km'
    xvar.long_name = 'x coordinate (Mercator)'

    yvar = ds.createVariable('y', 'f4', ('y',))
    yvar[:] = y
    yvar.units = 'km'
    yvar.long_name = 'y coordinate (Mercator)'

    # Emissions
    for name, grid, units in [
        ('PM2.5', pm25_grid, 'g/s'),
        ('PM10', pm10_grid, 'g/s'),
        ('NOx', nox_grid, 'g/s'),
        ('SO2', so2_grid, 'g/s'),
    ]:
        var = ds.createVariable(name, 'f4', ('y', 'x'))
        var[:] = grid
        var.units = units
        var.long_name = f'{name} emissions from brick kilns'

    # Global attributes
    ds.title = f'Brick Kiln Emissions - {scenario_name}'
    ds.source = 'WRF-CAMx brick kiln analysis'
    ds.grid_nx = nx
    ds.grid_ny = ny
    ds.grid_dx = CAMX_GRID['dx']
    ds.grid_dy = CAMX_GRID['dy']
    ds.ref_lon = CAMX_GRID['ref_lon']
    ds.ref_lat = CAMX_GRID['ref_lat']

    ds.close()

    return output_path


# =============================================================================
# VISUALIZATION
# =============================================================================
def create_comprehensive_visualization(india_gdf, up_gdf, delhi_gdf, kilns, converted,
                                       emissions_2020, emissions_2025, met_data,
                                       LON, LAT, conc_2020, conc_2025,
                                       delhi_conc_2020, delhi_conc_2025, output_path):
    """Create the main visualization figure with concentration maps."""
    print("Creating comprehensive visualization...")

    fig = plt.figure(figsize=(24, 18))
    fig.patch.set_facecolor('white')

    # Main title
    fig.suptitle('BRICK KILN AIR QUALITY IMPACT ANALYSIS\n' +
                 'Using WRF-CAMx Modeling System with Actual India Shapefiles',
                fontsize=22, fontweight='bold', y=0.98)

    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.25,
                         left=0.04, right=0.96, top=0.92, bottom=0.05)

    # =========================================================================
    # ROW 1: Maps
    # =========================================================================

    # Panel 1: Scenario 2020 - All FCBK
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax1.set_facecolor('#E8F4F8')

    # Plot India outline (light gray)
    if india_gdf is not None:
        india_gdf.plot(ax=ax1, facecolor='#F5F5F5', edgecolor='#CCCCCC', linewidth=0.5)

    # Plot UP (highlighted)
    if up_gdf is not None:
        up_gdf.plot(ax=ax1, facecolor=UP_COLOR, edgecolor='#4a4a4a', linewidth=2, alpha=0.8)

    # Plot Delhi (gold highlight)
    if delhi_gdf is not None:
        delhi_gdf.plot(ax=ax1, facecolor='#FFE4B5', edgecolor='#FF8C00', linewidth=3, alpha=0.9)

    # Plot kilns (all red)
    lons = [k['lon'] for k in kilns]
    lats = [k['lat'] for k in kilns]
    ax1.scatter(lons, lats, c=FCBK_COLOR, s=5, alpha=0.5, marker='o',
               edgecolors='darkred', linewidths=0.1)

    # Delhi marker
    ax1.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], '*', color=DELHI_COLOR,
            markersize=25, markeredgecolor='black', markeredgewidth=2, zorder=100)
    ax1.annotate('DELHI', (DELHI_CENTER['lon'] + 0.3, DELHI_CENTER['lat'] + 0.3),
                fontsize=12, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    # Wind arrows (if met data available)
    if met_data is not None:
        # Subsample for clarity
        skip = 10
        lon_sub = met_data['lon'][::skip, ::skip]
        lat_sub = met_data['lat'][::skip, ::skip]
        u_sub = met_data['uwind'][::skip, ::skip]
        v_sub = met_data['vwind'][::skip, ::skip]

        ax1.quiver(lon_sub, lat_sub, u_sub, v_sub, alpha=0.3, scale=100, color='blue')

    ax1.set_xlim(75, 86)
    ax1.set_ylim(23, 32)
    ax1.set_title('SCENARIO 2020: 10,000 FCBK Kilns (100%)', fontsize=14, fontweight='bold', color=FCBK_COLOR)
    ax1.set_xlabel('Longitude (°E)')
    ax1.set_ylabel('Latitude (°N)')
    ax1.grid(True, alpha=0.3, linestyle='--')

    # Panel 2: Scenario 2025 - 50% RANDOMLY converted
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax2.set_facecolor('#E8F4F8')

    if india_gdf is not None:
        india_gdf.plot(ax=ax2, facecolor='#F5F5F5', edgecolor='#CCCCCC', linewidth=0.5)
    if up_gdf is not None:
        up_gdf.plot(ax=ax2, facecolor=UP_COLOR, edgecolor='#4a4a4a', linewidth=2, alpha=0.8)
    if delhi_gdf is not None:
        delhi_gdf.plot(ax=ax2, facecolor='#FFE4B5', edgecolor='#FF8C00', linewidth=3, alpha=0.9)

    # Use the pre-computed randomly converted set
    fcbk_lons = [lons[i] for i in range(len(kilns)) if i not in converted]
    fcbk_lats = [lats[i] for i in range(len(kilns)) if i not in converted]
    zig_lons = [lons[i] for i in range(len(kilns)) if i in converted]
    zig_lats = [lats[i] for i in range(len(kilns)) if i in converted]

    ax2.scatter(fcbk_lons, fcbk_lats, c=FCBK_COLOR, s=5, alpha=0.5, marker='o')
    ax2.scatter(zig_lons, zig_lats, c=ZIGZAG_COLOR, s=6, alpha=0.6, marker='s')

    ax2.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], '*', color=DELHI_COLOR,
            markersize=25, markeredgecolor='black', markeredgewidth=2, zorder=100)
    ax2.annotate('DELHI', (DELHI_CENTER['lon'] + 0.3, DELHI_CENTER['lat'] + 0.3),
                fontsize=12, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    if met_data is not None:
        ax2.quiver(lon_sub, lat_sub, u_sub, v_sub, alpha=0.3, scale=100, color='blue')

    ax2.set_xlim(75, 86)
    ax2.set_ylim(23, 32)
    ax2.set_title('SCENARIO 2025: 50% Converted to Zigzag', fontsize=14, fontweight='bold', color=ZIGZAG_COLOR)
    ax2.set_xlabel('Longitude (°E)')
    ax2.set_ylabel('Latitude (°N)')
    ax2.grid(True, alpha=0.3, linestyle='--')

    # =========================================================================
    # ROW 2: Emissions & CAMx Pipeline
    # =========================================================================

    # Panel 3: Emission comparison
    ax3 = fig.add_subplot(gs[1, 0])

    total_2020 = sum(e['PM2.5_tpy'] for e in emissions_2020) / 1000
    total_2025 = sum(e['PM2.5_tpy'] for e in emissions_2025) / 1000
    reduction = total_2020 - total_2025
    pct = (reduction / total_2020) * 100

    bars = ax3.bar(['2020\n(FCBK)', '2025\n(Mixed)'], [total_2020, total_2025],
                  color=[FCBK_COLOR, ZIGZAG_COLOR], alpha=0.8, width=0.6)

    ax3.set_ylabel('PM2.5 Emissions\n(thousand tonnes/year)', fontsize=11)
    ax3.set_title('Total PM2.5 Emissions', fontsize=13, fontweight='bold')
    ax3.set_ylim(0, total_2020 * 1.3)

    for bar, val in zip(bars, [total_2020, total_2025]):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{val:.0f}k', ha='center', fontsize=12, fontweight='bold')

    ax3.annotate(f'-{pct:.0f}%', xy=(0.5, (total_2020+total_2025)/2),
                fontsize=16, fontweight='bold', color='green', ha='center')

    # Panel 4: CAMx Flow Diagram
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax4.axis('off')

    # Draw boxes
    boxes = [
        (0.05, 0.7, 'WRF\nMeteorology\n(winds, temp)', '#87CEEB'),
        (0.05, 0.3, 'Brick Kiln\nInventory\n(10,000 kilns)', '#FFB6C1'),
        (0.35, 0.5, 'CAMx\nEmissions\nProcessor', '#98FB98'),
        (0.65, 0.5, 'CAMx v7.32\nAir Quality\nModel', '#DDA0DD'),
        (0.9, 0.5, 'PM2.5\nConcentration\nMaps', '#F0E68C'),
    ]

    for x, y, text, color in boxes:
        box = FancyBboxPatch((x-0.12, y-0.15), 0.24, 0.3,
                            boxstyle="round,pad=0.02,rounding_size=0.02",
                            facecolor=color, edgecolor='black', linewidth=2,
                            transform=ax4.transAxes)
        ax4.add_patch(box)
        ax4.text(x, y, text, ha='center', va='center', fontsize=10, fontweight='bold',
                transform=ax4.transAxes)

    # Arrows
    arrows = [(0.17, 0.7, 0.23, 0.55), (0.17, 0.3, 0.23, 0.45),
              (0.47, 0.5, 0.53, 0.5), (0.77, 0.5, 0.78, 0.5)]
    for x1, y1, x2, y2 in arrows:
        ax4.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='black', lw=2),
                    transform=ax4.transAxes)

    ax4.set_title('WRF-CAMx Modeling Pipeline', fontsize=13, fontweight='bold')
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)

    # Panel 5: Technology comparison
    ax5 = fig.add_subplot(gs[1, 3])

    kiln_types = ['FCBK', 'Zigzag']
    pm25_ef = [EMISSION_FACTORS[k]['PM2.5'] for k in kiln_types]
    colors = [FCBK_COLOR, ZIGZAG_COLOR]

    bars = ax5.barh(kiln_types, pm25_ef, color=colors, alpha=0.8, height=0.5)
    ax5.set_xlabel('PM2.5 Emission Factor (g/kg brick)')
    ax5.set_title('Kiln Technology\nComparison', fontsize=13, fontweight='bold')

    for bar, val in zip(bars, pm25_ef):
        ax5.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2,
                f'{val:.2f}', va='center', fontsize=11, fontweight='bold')

    ax5.text(0.5, 0.15, '-60%', fontsize=14, fontweight='bold', color='green',
            transform=ax5.transAxes, ha='center')

    # =========================================================================
    # ROW 3: CONCENTRATION MAPS (µg/m³) - The key result!
    # =========================================================================

    # Panel 6: Concentration 2020
    ax6 = fig.add_subplot(gs[2, 0])

    vmax_conc = max(conc_2020.max(), conc_2025.max())
    im6 = ax6.pcolormesh(LON, LAT, conc_2020, cmap='RdYlGn_r', shading='auto',
                         vmin=0, vmax=vmax_conc)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax6, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax6, color='orange', linewidth=2)
    ax6.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=15)
    ax6.set_title(f'2020 PM2.5 Concentration\n(brick kiln contribution)', fontsize=12, fontweight='bold')
    ax6.set_xlim(75, 86)
    ax6.set_ylim(23, 32)
    ax6.set_xlabel('Longitude (°E)')
    ax6.set_ylabel('Latitude (°N)')
    cbar6 = plt.colorbar(im6, ax=ax6, shrink=0.8)
    cbar6.set_label('µg/m³')

    # Panel 7: Concentration 2025
    ax7 = fig.add_subplot(gs[2, 1])

    im7 = ax7.pcolormesh(LON, LAT, conc_2025, cmap='RdYlGn_r', shading='auto',
                         vmin=0, vmax=vmax_conc)
    if up_gdf is not None:
        up_gdf.boundary.plot(ax=ax7, color='black', linewidth=1)
    if delhi_gdf is not None:
        delhi_gdf.boundary.plot(ax=ax7, color='orange', linewidth=2)
    ax7.plot(DELHI_CENTER['lon'], DELHI_CENTER['lat'], 'k*', markersize=15)
    ax7.set_title(f'2025 PM2.5 Concentration\n(brick kiln contribution)', fontsize=12, fontweight='bold')
    ax7.set_xlim(75, 86)
    ax7.set_ylim(23, 32)
    ax7.set_xlabel('Longitude (°E)')
    ax7.set_ylabel('Latitude (°N)')
    cbar7 = plt.colorbar(im7, ax=ax7, shrink=0.8)
    cbar7.set_label('µg/m³')

    # Panel 8: Delhi concentration bar chart + Summary
    ax8 = fig.add_subplot(gs[2, 2])

    # Background PM2.5 in Delhi (from other sources)
    background_pm25 = 80.0  # µg/m³

    x_pos = [0, 1]
    x_labels = ['2020\n(100% FCBK)', '2025\n(50% Zigzag)']

    # Stacked bar: background + brick kiln contribution
    ax8.bar(x_pos, [background_pm25, background_pm25], color='gray', alpha=0.6, label='Other sources')
    ax8.bar(x_pos, [delhi_conc_2020, delhi_conc_2025], bottom=[background_pm25, background_pm25],
           color=[FCBK_COLOR, ZIGZAG_COLOR], alpha=0.8, label='Brick kilns')

    ax8.axhline(y=60, color='orange', linestyle='--', linewidth=2, alpha=0.7)
    ax8.text(1.15, 60, 'WHO\nIT-3', fontsize=9, va='center')
    ax8.axhline(y=35, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax8.text(1.15, 35, 'WHO\nIT-2', fontsize=9, va='center')

    ax8.set_ylabel('PM2.5 (µg/m³)', fontsize=11)
    ax8.set_xticks(x_pos)
    ax8.set_xticklabels(x_labels, fontsize=10)
    ax8.set_title('Delhi PM2.5 Concentration', fontsize=13, fontweight='bold')
    ax8.set_ylim(0, 120)
    ax8.legend(loc='upper right', fontsize=9)

    # Add value annotations
    for i, (bg, bk) in enumerate(zip([background_pm25]*2, [delhi_conc_2020, delhi_conc_2025])):
        total = bg + bk
        ax8.text(i, total + 2, f'{total:.1f}', ha='center', fontsize=11, fontweight='bold')
        ax8.text(i, bg + bk/2, f'+{bk:.2f}', ha='center', fontsize=10, color='white', fontweight='bold')

    # Reduction annotation
    conc_reduction = delhi_conc_2020 - delhi_conc_2025
    conc_reduction_pct = (conc_reduction / delhi_conc_2020) * 100
    ax8.annotate(f'↓{conc_reduction:.2f} µg/m³\n({conc_reduction_pct:.0f}%)',
                xy=(0.5, background_pm25 + delhi_conc_2025 + 5),
                fontsize=11, fontweight='bold', color='green', ha='center')

    # Panel 9: Summary stats
    ax9 = fig.add_subplot(gs[2, 3])
    ax9.axis('off')

    # Health impact calculation using actual concentrations
    delhi_pop = 32_000_000
    mortality_coef = 0.006  # 0.6% per µg/m³
    baseline_mortality = 0.008

    deaths_2020 = delhi_pop * baseline_mortality * mortality_coef * delhi_conc_2020
    deaths_2025 = delhi_pop * baseline_mortality * mortality_coef * delhi_conc_2025
    lives_saved = deaths_2020 - deaths_2025

    summary_text = f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
      IMPACT SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Kilns: {len(kilns):,} total
Converted: {len(converted):,} (random)

EMISSIONS (tonnes/yr):
  2020: {total_2020/1000:,.0f}k
  2025: {total_2025/1000:,.0f}k
  Reduction: {pct:.0f}%

DELHI CONCENTRATION:
  2020: {delhi_conc_2020:.2f} µg/m³
  2025: {delhi_conc_2025:.2f} µg/m³
  Reduction: {conc_reduction:.2f} µg/m³

HEALTH BENEFIT:
  Lives saved: ~{lives_saved:,.0f}/yr

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

    ax9.text(0.5, 0.5, summary_text, transform=ax9.transAxes, fontsize=11,
            verticalalignment='center', horizontalalignment='center',
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#E8F8E8', edgecolor='green', linewidth=2))

    # Legend
    legend_elements = [
        Patch(facecolor=UP_COLOR, edgecolor='#4a4a4a', label='Uttar Pradesh'),
        Patch(facecolor='#FFE4B5', edgecolor='#FF8C00', label='Delhi NCT'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=FCBK_COLOR,
               markersize=10, label='FCBK Kiln (polluting)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=ZIGZAG_COLOR,
               markersize=10, label='Zigzag Kiln (cleaner)'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor=DELHI_COLOR,
               markeredgecolor='black', markersize=15, label='Delhi (target city)'),
    ]

    fig.legend(handles=legend_elements, loc='lower center', ncol=5, fontsize=11,
              bbox_to_anchor=(0.5, 0.01), frameon=True, edgecolor='gray')

    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path}")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 70)
    print("BRICK KILN IMPACT ANALYSIS USING WRF-CAMx")
    print("=" * 70)

    # Load shapefiles
    india_gdf, up_gdf, delhi_gdf = load_india_shapefiles()

    # Load WRF-CAMx meteorology
    met_data = load_camx_meteorology()

    # Generate kilns
    kilns = generate_kilns_in_up(up_gdf, N_KILNS)

    # Create emission scenarios
    print("\nCalculating emissions...")

    # Scenario 2020: All FCBK
    kiln_types_2020 = ['FCBK'] * len(kilns)
    emissions_2020 = calculate_emissions(kilns, kiln_types_2020)

    # Scenario 2025: 50% Zigzag (RANDOMLY selected)
    random_idx = np.random.permutation(len(kilns))
    converted = set(random_idx[:len(kilns)//2])
    kiln_types_2025 = ['Zigzag' if i in converted else 'FCBK' for i in range(len(kilns))]
    emissions_2025 = calculate_emissions(kilns, kiln_types_2025)

    # Save converted kiln IDs for reproducibility
    converted_ids_file = os.path.join(CAMX_EMIS_DIR, 'converted_kiln_ids.txt')
    with open(converted_ids_file, 'w') as f:
        f.write("# Kiln IDs converted from FCBK to Zigzag in 2025 scenario\n")
        f.write(f"# Total converted: {len(converted)} out of {len(kilns)}\n")
        f.write("# Random seed: 42\n")
        for idx in sorted(converted):
            f.write(f"{idx}\n")
    print(f"  Saved converted kiln IDs to: {converted_ids_file}")

    total_2020 = sum(e['PM2.5_tpy'] for e in emissions_2020)
    total_2025 = sum(e['PM2.5_tpy'] for e in emissions_2025)

    print(f"\n  2020 PM2.5 Emissions: {total_2020:,.0f} tonnes/year")
    print(f"  2025 PM2.5 Emissions: {total_2025:,.0f} tonnes/year")
    print(f"  Emission Reduction: {total_2020-total_2025:,.0f} tonnes/year ({(total_2020-total_2025)/total_2020*100:.0f}%)")

    # Calculate concentration fields using WRF-CAMx meteorology
    print("\nCalculating PM2.5 concentrations using WRF-CAMx meteorology...")
    LON, LAT, conc_2020 = calculate_concentration_field(emissions_2020, met_data)
    _, _, conc_2025 = calculate_concentration_field(emissions_2025, met_data)

    # Get Delhi concentrations
    delhi_conc_2020 = get_delhi_concentration(LON, LAT, conc_2020)
    delhi_conc_2025 = get_delhi_concentration(LON, LAT, conc_2025)

    print(f"\n  Delhi PM2.5 (brick kiln contribution):")
    print(f"    2020: {delhi_conc_2020:.2f} µg/m³")
    print(f"    2025: {delhi_conc_2025:.2f} µg/m³")
    print(f"    Reduction: {delhi_conc_2020-delhi_conc_2025:.2f} µg/m³ ({(delhi_conc_2020-delhi_conc_2025)/delhi_conc_2020*100:.1f}%)")

    # Create CAMx emission files
    print("\nCreating CAMx emission files...")

    create_camx_point_source_file(
        emissions_2020,
        os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2020.ptsrc'),
        '2020 - 100% FCBK'
    )

    create_camx_point_source_file(
        emissions_2025,
        os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2025.ptsrc'),
        '2025 - 50% Zigzag'
    )

    create_camx_gridded_emission_file(
        emissions_2020,
        os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2020.nc'),
        '2020 - 100% FCBK'
    )

    create_camx_gridded_emission_file(
        emissions_2025,
        os.path.join(CAMX_EMIS_DIR, 'brick_kilns_2025.nc'),
        '2025 - 50% Zigzag'
    )

    # Create visualization
    print("\nCreating visualization...")
    create_comprehensive_visualization(
        india_gdf, up_gdf, delhi_gdf, kilns, converted,
        emissions_2020, emissions_2025, met_data,
        LON, LAT, conc_2020, conc_2025,
        delhi_conc_2020, delhi_conc_2025,
        os.path.join(FIGURES_DIR, 'brick_kiln_camx_analysis.png')
    )

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print(f"""
OUTPUT FILES:

Figures:
  {FIGURES_DIR}/brick_kiln_camx_analysis.png

CAMx Emission Files:
  {CAMX_EMIS_DIR}/brick_kilns_2020.ptsrc  (point source, 2020)
  {CAMX_EMIS_DIR}/brick_kilns_2025.ptsrc  (point source, 2025)
  {CAMX_EMIS_DIR}/brick_kilns_2020.nc     (gridded, 2020)
  {CAMX_EMIS_DIR}/brick_kilns_2025.nc     (gridded, 2025)

TO RUN CAMx:
  1. Edit run_camx/CAMx.india.job to include brick kiln emissions
  2. Run: ./CAMx.india.job
  3. Compare output PM2.5 between scenarios
""")


if __name__ == "__main__":
    main()
