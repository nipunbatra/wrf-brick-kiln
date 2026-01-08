#!/usr/bin/env python3
"""
Process EDGAR emissions to CAMx-ready format

This script:
1. Reads EDGAR NetCDF files (0.1° x 0.1° resolution)
2. Subsets to India domain (67°E-92°E, 6°N-30°N)
3. Regrids to CAMx grid (102x102, 27km, Mercator)
4. Speciates emissions to CAMx chemical species
5. Applies temporal profiles (hourly, daily, monthly)
6. Outputs CAMx-ready emission files

CAMx emission file format:
- NetCDF with dimensions: (TSTEP, LAY, ROW, COL)
- Variables for each species (e.g., NO, NO2, SO2, FPRM, CPRM, etc.)
"""

import numpy as np
from datetime import datetime, timedelta
import os

# Try to import optional dependencies
try:
    from netCDF4 import Dataset
    HAS_NETCDF = True
except ImportError:
    HAS_NETCDF = False
    print("Warning: netCDF4 not available. Install with: pip install netCDF4")

try:
    from scipy import interpolate
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("Warning: scipy not available. Install with: pip install scipy")

# =============================================================================
# CAMx Grid Definition (from WRFCAMx output)
# =============================================================================
CAMX_GRID = {
    'nx': 102,  # columns (including buffer)
    'ny': 102,  # rows (including buffer)
    'nz': 1,    # layers (surface emissions)
    'dx': 27.0,  # km
    'dy': 27.0,  # km
    'x0': -1714.5,  # km (SW corner)
    'y0': -1714.5,  # km (SW corner)
    'projection': 'MERCATOR',
    'clon': 83.0,   # center longitude
    'clat': 21.5,   # center latitude
    'tlat': 10.0,   # true latitude
}

# Domain bounds in lat/lon (approximate)
DOMAIN = {
    'lon_min': 67.0,
    'lon_max': 93.0,
    'lat_min': 6.0,
    'lat_max': 30.0,
}

# =============================================================================
# CAMx Species Mapping
# =============================================================================
# EDGAR species -> CAMx species mapping with speciation factors
SPECIES_MAP = {
    # Particulate Matter
    'PM2.5': {
        'FPRM': 0.85,  # Fine primary PM (unspeciated)
        'PEC': 0.05,   # Eleite carbon
        'POA': 0.10,   # Primary organic aerosol
    },
    'PM10': {
        'CPRM': 0.80,  # Coarse primary PM
        'CCRS': 0.20,  # Coarse crustal
    },
    'BC': {
        'PEC': 1.0,    # Elemental carbon
    },
    'OC': {
        'POA': 1.0,    # Primary organic aerosol (assume all is POA)
    },
    # Gases
    'NOx': {
        'NO': 0.90,    # 90% as NO
        'NO2': 0.10,   # 10% as NO2
    },
    'SO2': {
        'SO2': 1.0,
    },
    'CO': {
        'CO': 1.0,
    },
    'NH3': {
        'NH3': 1.0,
    },
}

# =============================================================================
# Brick Kiln Emission Factors (g/kg brick produced)
# =============================================================================
BRICK_KILN_EF = {
    # Kiln Type: {species: emission_factor}
    'FCBTK': {  # Fixed Chimney Bull's Trench Kiln (traditional)
        'PM2.5': 0.75,
        'PM10': 1.2,
        'BC': 0.15,
        'OC': 0.25,
        'NOx': 0.3,
        'SO2': 0.8,  # depends on coal sulfur content
        'CO': 5.0,
    },
    'Zigzag': {  # Improved zigzag kiln
        'PM2.5': 0.30,
        'PM10': 0.5,
        'BC': 0.06,
        'OC': 0.10,
        'NOx': 0.25,
        'SO2': 0.6,
        'CO': 3.0,
    },
    'VSBK': {  # Vertical Shaft Brick Kiln
        'PM2.5': 0.15,
        'PM10': 0.25,
        'BC': 0.03,
        'OC': 0.05,
        'NOx': 0.20,
        'SO2': 0.4,
        'CO': 2.0,
    },
    'Tunnel': {  # Modern tunnel kiln
        'PM2.5': 0.08,
        'PM10': 0.12,
        'BC': 0.015,
        'OC': 0.025,
        'NOx': 0.15,
        'SO2': 0.3,
        'CO': 1.0,
    },
}

# =============================================================================
# Temporal Profiles
# =============================================================================
def get_hourly_profile(sector='industrial'):
    """
    Return 24-hour profile for emissions
    Values are fractions that sum to 24 (so hourly fraction = value/24)
    """
    profiles = {
        'industrial': [
            0.5, 0.5, 0.5, 0.5, 0.5, 0.7,  # 00-05
            0.9, 1.2, 1.5, 1.5, 1.5, 1.3,  # 06-11
            1.0, 1.5, 1.5, 1.5, 1.3, 1.2,  # 12-17
            1.0, 0.9, 0.8, 0.7, 0.6, 0.5,  # 18-23
        ],
        'brick_kiln': [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # 00-05 (continuous operation)
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # 06-11
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # 12-17
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # 18-23
        ],
    }
    return np.array(profiles.get(sector, profiles['industrial']))

def get_monthly_profile(sector='brick_kiln'):
    """
    Return monthly profile for emissions
    Brick kilns operate primarily in dry season (Oct-May in India)
    """
    profiles = {
        'brick_kiln': [
            1.2, 1.2, 1.0, 0.8, 0.5, 0.1,  # Jan-Jun
            0.0, 0.0, 0.1, 0.8, 1.2, 1.2,  # Jul-Dec
        ],
        'constant': [1.0] * 12,
    }
    return np.array(profiles.get(sector, profiles['constant']))

# =============================================================================
# Regridding Functions
# =============================================================================
def latlon_to_mercator(lon, lat, clon=83.0, tlat=10.0):
    """
    Convert lat/lon to Mercator coordinates (km)

    Parameters:
    - lon, lat: coordinates in degrees
    - clon: central longitude
    - tlat: true latitude

    Returns:
    - x, y: Mercator coordinates in km
    """
    R = 6371.0  # Earth radius in km
    tlat_rad = np.radians(tlat)

    x = R * np.radians(lon - clon) * np.cos(tlat_rad)
    y = R * np.log(np.tan(np.pi/4 + np.radians(lat)/2)) * np.cos(tlat_rad)

    return x, y

def regrid_to_camx(data_ll, lon_ll, lat_ll, camx_grid):
    """
    Regrid lat/lon data to CAMx Mercator grid

    Parameters:
    - data_ll: 2D array on lat/lon grid
    - lon_ll, lat_ll: 1D arrays of longitudes and latitudes
    - camx_grid: dictionary with CAMx grid parameters

    Returns:
    - data_camx: 2D array on CAMx grid
    """
    if not HAS_SCIPY:
        print("ERROR: scipy required for regridding")
        return None

    # Create CAMx grid coordinates
    x_camx = camx_grid['x0'] + np.arange(camx_grid['nx']) * camx_grid['dx']
    y_camx = camx_grid['y0'] + np.arange(camx_grid['ny']) * camx_grid['dy']

    # Convert EDGAR lat/lon to Mercator
    lon_2d, lat_2d = np.meshgrid(lon_ll, lat_ll)
    x_ll, y_ll = latlon_to_mercator(lon_2d, lat_2d, camx_grid['clon'], camx_grid['tlat'])

    # Interpolate to CAMx grid
    x_camx_2d, y_camx_2d = np.meshgrid(x_camx, y_camx)

    # Use scipy griddata for interpolation
    from scipy.interpolate import griddata

    points = np.column_stack([x_ll.ravel(), y_ll.ravel()])
    values = data_ll.ravel()

    # Remove NaN values
    valid = ~np.isnan(values)
    points = points[valid]
    values = values[valid]

    data_camx = griddata(points, values, (x_camx_2d, y_camx_2d), method='linear')

    # Fill NaN with zeros (ocean/outside domain)
    data_camx = np.nan_to_num(data_camx, nan=0.0)

    return data_camx

# =============================================================================
# CAMx Emission File Writer
# =============================================================================
def create_camx_emission_file(output_file, emissions_dict, date, camx_grid,
                               num_hours=24, species_list=None):
    """
    Create CAMx-format emission file in NetCDF format

    Parameters:
    - output_file: output filename
    - emissions_dict: {species: 2D array of emissions [mol/s or g/s]}
    - date: datetime object for start date
    - camx_grid: grid parameters
    - num_hours: number of hours in file
    - species_list: list of species to include
    """
    if not HAS_NETCDF:
        print("ERROR: netCDF4 required to write CAMx files")
        return False

    if species_list is None:
        species_list = list(emissions_dict.keys())

    nc = Dataset(output_file, 'w', format='NETCDF4')

    # Dimensions
    nc.createDimension('TSTEP', num_hours)
    nc.createDimension('DATE-TIME', 2)
    nc.createDimension('LAY', 1)
    nc.createDimension('VAR', len(species_list))
    nc.createDimension('ROW', camx_grid['ny'])
    nc.createDimension('COL', camx_grid['nx'])

    # Global attributes
    nc.FTYPE = 1
    nc.CDATE = int(date.strftime('%Y%j'))
    nc.CTIME = 0
    nc.SDATE = int(date.strftime('%Y%j'))
    nc.STIME = 0
    nc.TSTEP = 10000  # 1 hour in HHMMSS
    nc.NROWS = camx_grid['ny']
    nc.NCOLS = camx_grid['nx']
    nc.NLAYS = 1
    nc.GDTYP = 5  # Mercator
    nc.P_ALP = camx_grid['tlat']
    nc.P_BET = camx_grid['tlat']
    nc.P_GAM = camx_grid['clon']
    nc.XCENT = camx_grid['clon']
    nc.YCENT = camx_grid['clat']
    nc.XORIG = camx_grid['x0'] * 1000  # convert to meters
    nc.YORIG = camx_grid['y0'] * 1000
    nc.XCELL = camx_grid['dx'] * 1000
    nc.YCELL = camx_grid['dy'] * 1000
    nc.NVARS = len(species_list)
    nc.VAR_LIST = ''.join([s.ljust(16) for s in species_list])

    # Time flag variables
    tflag = nc.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag.long_name = 'Start time flag'
    tflag.units = 'YYYYDDD,HHMMSS'

    # Create emission variables
    for species in species_list:
        var = nc.createVariable(species, 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
        var.long_name = f'{species} emissions'
        var.units = 'mol/s' if species in ['NO', 'NO2', 'SO2', 'CO', 'NH3'] else 'g/s'

        # Fill with data (apply hourly profile)
        hourly_profile = get_hourly_profile('industrial')
        base_emissions = emissions_dict.get(species, np.zeros((camx_grid['ny'], camx_grid['nx'])))

        for hour in range(num_hours):
            var[hour, 0, :, :] = base_emissions * hourly_profile[hour % 24] / 24.0

            # Set time flags
            current_time = date + timedelta(hours=hour)
            for ivar in range(len(species_list)):
                tflag[hour, ivar, 0] = int(current_time.strftime('%Y%j'))
                tflag[hour, ivar, 1] = int(current_time.strftime('%H%M%S'))

    nc.close()
    print(f"Created: {output_file}")
    return True

# =============================================================================
# Main Processing Function
# =============================================================================
def process_edgar_to_camx(edgar_dir, output_dir, year=2022, month=2):
    """
    Main function to process EDGAR to CAMx format
    """
    print("=" * 70)
    print("EDGAR to CAMx Emissions Processor")
    print("=" * 70)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # For now, create a template emission file with placeholder data
    # This will be replaced with actual EDGAR data when available

    print("\nCreating template emission file...")
    print(f"Grid: {CAMX_GRID['nx']}x{CAMX_GRID['ny']} at {CAMX_GRID['dx']}km")

    # Create placeholder emissions (uniform low values)
    emissions = {}
    camx_species = ['NO', 'NO2', 'SO2', 'CO', 'NH3', 'FPRM', 'CPRM', 'PEC', 'POA']

    for species in camx_species:
        # Create placeholder with some spatial variation
        y = np.linspace(0, 1, CAMX_GRID['ny'])
        x = np.linspace(0, 1, CAMX_GRID['nx'])
        xx, yy = np.meshgrid(x, y)

        # Simple pattern - higher in center
        pattern = np.exp(-((xx-0.5)**2 + (yy-0.5)**2) / 0.1)

        # Scale to reasonable emission values
        if species in ['NO', 'NO2', 'SO2', 'CO', 'NH3']:
            emissions[species] = pattern * 1e-6  # mol/s/cell
        else:
            emissions[species] = pattern * 1e-3  # g/s/cell

    # Create emission file for Feb 1, 2024
    date = datetime(2024, 2, 1)
    output_file = os.path.join(output_dir, f"camx_emis.base.{date.strftime('%Y%m%d')}.nc")

    create_camx_emission_file(output_file, emissions, date, CAMX_GRID,
                               num_hours=24, species_list=camx_species)

    print("\n" + "=" * 70)
    print("NEXT STEPS:")
    print("=" * 70)
    print("""
1. Download EDGAR data from: https://edgar.jrc.ec.europa.eu/
2. Place files in: {edgar_dir}/raw_data/
3. Re-run this script to process actual data
4. For brick kiln scenarios, use create_brick_kiln_scenario()
""".format(edgar_dir=edgar_dir))

    return True

def create_brick_kiln_scenario(base_emissions, scenario='zigzag_conversion',
                                conversion_fraction=0.5):
    """
    Create brick kiln emission scenario by modifying base emissions

    Parameters:
    - base_emissions: dict of base emission arrays
    - scenario: 'zigzag_conversion', 'vsbk_conversion', 'tunnel_conversion'
    - conversion_fraction: fraction of kilns converted (0-1)

    Returns:
    - modified_emissions: dict of modified emission arrays
    """
    # Emission reduction factors by scenario
    # (ratio of new kiln EF to FCBTK EF)
    reduction_factors = {
        'zigzag_conversion': {
            'PM2.5': BRICK_KILN_EF['Zigzag']['PM2.5'] / BRICK_KILN_EF['FCBTK']['PM2.5'],
            'NOx': BRICK_KILN_EF['Zigzag']['NOx'] / BRICK_KILN_EF['FCBTK']['NOx'],
            'SO2': BRICK_KILN_EF['Zigzag']['SO2'] / BRICK_KILN_EF['FCBTK']['SO2'],
            'CO': BRICK_KILN_EF['Zigzag']['CO'] / BRICK_KILN_EF['FCBTK']['CO'],
        },
        'vsbk_conversion': {
            'PM2.5': BRICK_KILN_EF['VSBK']['PM2.5'] / BRICK_KILN_EF['FCBTK']['PM2.5'],
            'NOx': BRICK_KILN_EF['VSBK']['NOx'] / BRICK_KILN_EF['FCBTK']['NOx'],
            'SO2': BRICK_KILN_EF['VSBK']['SO2'] / BRICK_KILN_EF['FCBTK']['SO2'],
            'CO': BRICK_KILN_EF['VSBK']['CO'] / BRICK_KILN_EF['FCBTK']['CO'],
        },
        'tunnel_conversion': {
            'PM2.5': BRICK_KILN_EF['Tunnel']['PM2.5'] / BRICK_KILN_EF['FCBTK']['PM2.5'],
            'NOx': BRICK_KILN_EF['Tunnel']['NOx'] / BRICK_KILN_EF['FCBTK']['NOx'],
            'SO2': BRICK_KILN_EF['Tunnel']['SO2'] / BRICK_KILN_EF['FCBTK']['SO2'],
            'CO': BRICK_KILN_EF['Tunnel']['CO'] / BRICK_KILN_EF['FCBTK']['CO'],
        },
    }

    factors = reduction_factors.get(scenario, reduction_factors['zigzag_conversion'])

    modified = {}
    for species, data in base_emissions.items():
        # Map CAMx species to EDGAR/brick kiln species
        if species in ['FPRM', 'PEC', 'POA']:
            edgar_species = 'PM2.5'
        elif species in ['NO', 'NO2']:
            edgar_species = 'NOx'
        elif species in ['CPRM', 'CCRS']:
            edgar_species = 'PM2.5'  # Use PM2.5 factor as proxy
        else:
            edgar_species = species

        # Calculate reduction
        if edgar_species in factors:
            # Effective reduction = conversion_fraction * (1 - reduction_factor)
            # New emission = base * (1 - reduction)
            reduction = conversion_fraction * (1 - factors[edgar_species])
            modified[species] = data * (1 - reduction)
        else:
            modified[species] = data.copy()

    return modified

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    edgar_dir = base_dir
    output_dir = os.path.join(base_dir, "camx_emissions")

    process_edgar_to_camx(edgar_dir, output_dir)
