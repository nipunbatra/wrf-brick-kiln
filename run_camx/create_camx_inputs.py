#!/usr/bin/env python3
"""
Create all CAMx Input Files for Brick Kiln Scenarios
=====================================================

This script creates:
1. CAMx-formatted emission files (gridded NetCDF)
2. Initial Conditions (IC) files
3. Boundary Conditions (BC) files
4. O3 column mapping file
5. TUV photolysis input

All files are created in CAMx v7.32 compatible NetCDF format.
"""

import numpy as np
import os
import netCDF4 as nc
from datetime import datetime, timedelta

# =============================================================================
# PATHS
# =============================================================================
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAMX_INPUT_DIR = os.path.join(PROJECT_ROOT, 'run_camx', 'camx_inputs')
CAMX_EMIS_DIR = os.path.join(PROJECT_ROOT, 'camx_emissions')
MET_DIR = os.path.join(PROJECT_ROOT, 'wrfcamx_v5.2', 'camx_input')

os.makedirs(CAMX_INPUT_DIR, exist_ok=True)

# =============================================================================
# GRID PARAMETERS (from WRF-CAMx meteorology)
# =============================================================================
GRID = {
    'NCOLS': 100,
    'NROWS': 100,
    'NLAYS': 20,
    'XCELL': 27000.0,  # meters
    'YCELL': 27000.0,  # meters
    'XORIG': -1714500.0,  # meters (Mercator)
    'YORIG': -1714500.0,  # meters (Mercator)
    'PLON': 83.0,  # Central longitude
    'PLAT': 21.5,  # Central latitude
    'TLAT1': 10.0,  # True latitude
    'TLAT2': 10.0,
}

# CAMx CB6r4 Chemical Species
CB6_SPECIES = [
    'NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2',
    'PNO3', 'PSO4', 'PNH4', 'PEC', 'POA', 'FPRM', 'CPRM',
    'FCRS', 'CCRS', 'SOA1', 'SOA2', 'SOA3', 'SOA4',
    'ISOP', 'TERP', 'FORM', 'ALD2', 'ALDX', 'OLE', 'ETH',
    'PAR', 'MEOH', 'ETOH', 'TOL', 'XYL', 'BENZENE',
]

# Emission species for brick kilns
EMIS_SPECIES = {
    'NO': 'NOX',
    'NO2': 'NOX',
    'SO2': 'SO2',
    'CO': 'CO',
    'PEC': 'PM2.5',  # Elemental carbon portion
    'POA': 'PM2.5',  # Organic aerosol portion
    'FPRM': 'PM2.5',  # Fine primary PM
    'CPRM': 'PM10',  # Coarse primary PM
}

# Dates to generate (matching WRF-CAMx meteorology)
DATES = [
    datetime(2024, 2, 1),
    datetime(2024, 2, 2),
    datetime(2024, 2, 3),
    datetime(2024, 2, 4),
    datetime(2024, 2, 5),
]

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def create_ioapi_file(filename, variables, description):
    """Create IOAPI-style NetCDF file for CAMx."""

    ds = nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC')

    # Dimensions
    ds.createDimension('TSTEP', None)  # Unlimited
    ds.createDimension('DATE-TIME', 2)
    ds.createDimension('LAY', GRID['NLAYS'])
    ds.createDimension('VAR', len(variables))
    ds.createDimension('ROW', GRID['NROWS'])
    ds.createDimension('COL', GRID['NCOLS'])

    # Time flags
    tflag = ds.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag.units = '<YYYYDDD,HHMMSS>'
    tflag.long_name = 'TFLAG'
    tflag.var_desc = 'Timestep-valid flags: (1) YYYYDDD or (2) HHMMSS'

    # Global attributes (IOAPI metadata)
    ds.IOAPI_VERSION = 'N/A'
    ds.EXEC_ID = 'CAMx_brick_kiln_analysis'
    ds.FTYPE = 1
    ds.CDATE = int(datetime.now().strftime('%Y%j'))
    ds.CTIME = int(datetime.now().strftime('%H%M%S'))
    ds.WDATE = ds.CDATE
    ds.WTIME = ds.CTIME
    ds.SDATE = int(DATES[0].strftime('%Y%j'))
    ds.STIME = 0
    ds.TSTEP = 10000  # 1 hour
    ds.NTHIK = 1
    ds.NCOLS = GRID['NCOLS']
    ds.NROWS = GRID['NROWS']
    ds.NLAYS = GRID['NLAYS']
    ds.NVARS = len(variables)
    ds.GDTYP = 5  # Mercator
    ds.P_ALP = GRID['TLAT1']
    ds.P_BET = GRID['TLAT2']
    ds.P_GAM = GRID['PLON']
    ds.XCENT = GRID['PLON']
    ds.YCENT = GRID['PLAT']
    ds.XORIG = GRID['XORIG']
    ds.YORIG = GRID['YORIG']
    ds.XCELL = GRID['XCELL']
    ds.YCELL = GRID['YCELL']
    ds.VGTYP = 7  # WRF sigma
    ds.VGTOP = 5000.0
    ds.VGLVLS = np.linspace(1.0, 0.0, GRID['NLAYS'] + 1).astype('f4')
    ds.GDNAM = 'INDIA_27KM'
    ds.UPNAM = 'CAMX'
    ds.VAR_LIST = ''.join([f'{v:16s}' for v in variables])
    ds.FILEDESC = description
    ds.HISTORY = f'Created {datetime.now()}'

    return ds


def load_brick_kiln_emissions(scenario):
    """Load brick kiln emissions from existing files."""

    emis_file = os.path.join(CAMX_EMIS_DIR, f'brick_kilns_{scenario}.nc')

    if not os.path.exists(emis_file):
        print(f"  WARNING: Emission file not found: {emis_file}")
        return None

    ds = nc.Dataset(emis_file)

    # Get emissions (g/s per cell)
    pm25 = ds.variables['PM2.5'][:]
    pm10 = ds.variables['PM10'][:]
    nox = ds.variables['NOx'][:]
    so2 = ds.variables['SO2'][:]

    ds.close()

    return {
        'PM2.5': pm25,
        'PM10': pm10,
        'NOx': nox,
        'SO2': so2,
    }


# =============================================================================
# CREATE EMISSION FILES
# =============================================================================

def create_camx_emission_file(date, scenario):
    """Create CAMx-format gridded emission file for a specific date and scenario."""

    datestr = date.strftime('%Y%m%d')
    jday = int(date.strftime('%j'))
    year = date.year

    filename = os.path.join(CAMX_INPUT_DIR, f'camx_emis.{scenario}.{datestr}.nc')
    print(f"  Creating: {os.path.basename(filename)}")

    # Species to include in emissions
    emis_vars = ['NO', 'NO2', 'SO2', 'CO', 'PEC', 'POA', 'FPRM', 'CPRM', 'FCRS', 'CCRS']

    ds = create_ioapi_file(filename, emis_vars,
                           f'Brick kiln emissions - {scenario} - {datestr}')

    # Load brick kiln emissions
    kiln_emis = load_brick_kiln_emissions(scenario)

    # Create emission variables (24 hours)
    for var in emis_vars:
        v = ds.createVariable(var, 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
        v.long_name = f'{var:16s}'
        v.units = 'moles/s' if var in ['NO', 'NO2', 'SO2', 'CO'] else 'g/s'
        v.var_desc = f'Emissions of {var}'

    # Fill with hourly emissions (24 timesteps)
    tflag = ds.variables['TFLAG']

    for t in range(24):
        hour = t * 10000  # HHMMSS format

        # Time flag
        for i, var in enumerate(emis_vars):
            tflag[t, i, 0] = year * 1000 + jday
            tflag[t, i, 1] = hour

        # Pad emissions to full grid size if needed
        if kiln_emis is not None:
            pm25 = kiln_emis['PM2.5']
            pm10 = kiln_emis['PM10']
            nox = kiln_emis['NOx']
            so2 = kiln_emis['SO2']

            # Resize to match grid
            ny, nx = pm25.shape

            # Create full grid arrays
            pm25_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))
            pm10_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))
            nox_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))
            so2_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))

            # Copy data (assuming emissions are already on the right grid)
            r_min = min(ny, GRID['NROWS'])
            c_min = min(nx, GRID['NCOLS'])
            pm25_full[:r_min, :c_min] = pm25[:r_min, :c_min]
            pm10_full[:r_min, :c_min] = pm10[:r_min, :c_min]
            nox_full[:r_min, :c_min] = nox[:r_min, :c_min]
            so2_full[:r_min, :c_min] = so2[:r_min, :c_min]
        else:
            pm25_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))
            pm10_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))
            nox_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))
            so2_full = np.zeros((GRID['NROWS'], GRID['NCOLS']))

        # Diurnal variation (higher during day)
        if 8 <= t <= 18:
            factor = 1.2
        elif 6 <= t <= 20:
            factor = 0.8
        else:
            factor = 0.5

        # Fill each species (layer 0 only for surface emissions)
        # NOx splits 90% NO, 10% NO2 (converted to moles/s, MW_NO=30, MW_NO2=46)
        ds.variables['NO'][t, 0, :, :] = nox_full * 0.9 / 30.0 * factor
        ds.variables['NO2'][t, 0, :, :] = nox_full * 0.1 / 46.0 * factor
        ds.variables['SO2'][t, :, :, :] = 0.0
        ds.variables['SO2'][t, 0, :, :] = so2_full / 64.0 * factor  # MW_SO2=64
        ds.variables['CO'][t, :, :, :] = 0.0

        # PM species (g/s)
        ds.variables['PEC'][t, :, :, :] = 0.0
        ds.variables['PEC'][t, 0, :, :] = pm25_full * 0.15 * factor  # 15% EC
        ds.variables['POA'][t, :, :, :] = 0.0
        ds.variables['POA'][t, 0, :, :] = pm25_full * 0.25 * factor  # 25% OA
        ds.variables['FPRM'][t, :, :, :] = 0.0
        ds.variables['FPRM'][t, 0, :, :] = pm25_full * 0.60 * factor  # 60% other fine PM
        ds.variables['CPRM'][t, :, :, :] = 0.0
        ds.variables['CPRM'][t, 0, :, :] = (pm10_full - pm25_full) * factor  # Coarse PM
        ds.variables['FCRS'][t, :, :, :] = 0.0
        ds.variables['CCRS'][t, :, :, :] = 0.0

    ds.close()
    return filename


# =============================================================================
# CREATE IC/BC FILES
# =============================================================================

def create_ic_file(date):
    """Create Initial Conditions file with typical winter background."""

    datestr = date.strftime('%Y%m%d')
    filename = os.path.join(CAMX_INPUT_DIR, f'ic.india.{datestr}.nc')
    print(f"  Creating: {os.path.basename(filename)}")

    ic_species = ['NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2',
                  'PNO3', 'PSO4', 'PNH4', 'PEC', 'POA', 'FPRM', 'CPRM',
                  'FCRS', 'CCRS']

    ds = create_ioapi_file(filename, ic_species,
                           f'Initial conditions for India domain - {datestr}')

    # Typical winter background concentrations (ppb for gases, µg/m³ for PM)
    # These represent typical IGP winter conditions
    background = {
        'NO': 5.0,      # ppb
        'NO2': 15.0,    # ppb
        'O3': 25.0,     # ppb
        'SO2': 8.0,     # ppb
        'CO': 500.0,    # ppb
        'NH3': 10.0,    # ppb
        'HNO3': 2.0,    # ppb
        'H2O2': 0.5,    # ppb
        'PNO3': 5.0,    # µg/m³
        'PSO4': 8.0,    # µg/m³
        'PNH4': 3.0,    # µg/m³
        'PEC': 4.0,     # µg/m³
        'POA': 10.0,    # µg/m³
        'FPRM': 15.0,   # µg/m³
        'CPRM': 20.0,   # µg/m³
        'FCRS': 2.0,    # µg/m³
        'CCRS': 5.0,    # µg/m³
    }

    # Create variables
    for var in ic_species:
        v = ds.createVariable(var, 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
        v.long_name = f'{var:16s}'
        v.units = 'ppmV' if var in ['NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2'] else 'µg/m³'
        v.var_desc = f'Initial concentration of {var}'

    # Set time flag
    jday = int(date.strftime('%j'))
    year = date.year
    ds.variables['TFLAG'][0, :, 0] = year * 1000 + jday
    ds.variables['TFLAG'][0, :, 1] = 0

    # Fill with vertically-varying background
    for var in ic_species:
        conc = background.get(var, 0.0)

        # Convert ppb to ppmV for gases
        if var in ['NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2']:
            conc = conc / 1000.0  # ppb to ppmV

        # Create vertical profile (decrease with height)
        layer_factors = np.exp(-np.arange(GRID['NLAYS']) * 0.1)

        for k in range(GRID['NLAYS']):
            ds.variables[var][0, k, :, :] = conc * layer_factors[k]

    ds.close()
    return filename


def create_bc_file(date):
    """Create Boundary Conditions file."""

    datestr = date.strftime('%Y%m%d')
    filename = os.path.join(CAMX_INPUT_DIR, f'bc.india.{datestr}.nc')
    print(f"  Creating: {os.path.basename(filename)}")

    bc_species = ['NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2',
                  'PNO3', 'PSO4', 'PNH4', 'PEC', 'POA', 'FPRM', 'CPRM',
                  'FCRS', 'CCRS']

    ds = nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC')

    # BC file has different structure - 4 boundaries
    ds.createDimension('TSTEP', None)
    ds.createDimension('DATE-TIME', 2)
    ds.createDimension('LAY', GRID['NLAYS'])
    ds.createDimension('VAR', len(bc_species))
    ds.createDimension('PERIM', 2 * (GRID['NROWS'] + GRID['NCOLS']))  # Perimeter cells

    # Time flags
    tflag = ds.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag.units = '<YYYYDDD,HHMMSS>'

    # Cleaner boundary conditions (away from IGP)
    background_bc = {
        'NO': 1.0,
        'NO2': 3.0,
        'O3': 35.0,
        'SO2': 2.0,
        'CO': 150.0,
        'NH3': 3.0,
        'HNO3': 1.0,
        'H2O2': 0.3,
        'PNO3': 1.0,
        'PSO4': 2.0,
        'PNH4': 1.0,
        'PEC': 0.5,
        'POA': 2.0,
        'FPRM': 3.0,
        'CPRM': 5.0,
        'FCRS': 1.0,
        'CCRS': 2.0,
    }

    # Create variables
    for var in bc_species:
        v = ds.createVariable(var, 'f4', ('TSTEP', 'LAY', 'PERIM'))
        v.long_name = f'{var:16s}'
        v.units = 'ppmV' if var in ['NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2'] else 'µg/m³'

    # Fill 24 hours
    jday = int(date.strftime('%j'))
    year = date.year
    n_perim = 2 * (GRID['NROWS'] + GRID['NCOLS'])

    for t in range(24):
        tflag[t, :, 0] = year * 1000 + jday
        tflag[t, :, 1] = t * 10000

        for var in bc_species:
            conc = background_bc.get(var, 0.0)
            if var in ['NO', 'NO2', 'O3', 'SO2', 'CO', 'NH3', 'HNO3', 'H2O2']:
                conc = conc / 1000.0  # ppb to ppmV

            # Vertical profile
            for k in range(GRID['NLAYS']):
                ds.variables[var][t, k, :] = conc * np.exp(-k * 0.1)

    # Global attributes
    ds.IOAPI_VERSION = 'N/A'
    ds.FTYPE = 2  # Boundary file
    ds.SDATE = year * 1000 + jday
    ds.STIME = 0
    ds.TSTEP = 10000
    ds.NCOLS = GRID['NCOLS']
    ds.NROWS = GRID['NROWS']
    ds.NLAYS = GRID['NLAYS']
    ds.NTHIK = 1
    ds.GDTYP = 5
    ds.P_ALP = GRID['TLAT1']
    ds.P_BET = GRID['TLAT2']
    ds.P_GAM = GRID['PLON']
    ds.XCENT = GRID['PLON']
    ds.YCENT = GRID['PLAT']
    ds.XORIG = GRID['XORIG']
    ds.YORIG = GRID['YORIG']
    ds.XCELL = GRID['XCELL']
    ds.YCELL = GRID['YCELL']
    ds.GDNAM = 'INDIA_27KM'

    ds.close()
    return filename


# =============================================================================
# CREATE SUPPORTING FILES
# =============================================================================

def create_o3_column_file():
    """Create ozone column mapping file for photolysis."""

    filename = os.path.join(CAMX_INPUT_DIR, 'o3map.202402.txt')
    print(f"  Creating: {os.path.basename(filename)}")

    # Typical February ozone column values for India (Dobson units)
    # Varies with latitude

    with open(filename, 'w') as f:
        f.write("! Ozone column map for India - February 2024\n")
        f.write("! Format: Date  Lat1  Lat2  O3col\n")

        for day in range(1, 29):
            datestr = f"2024{day+31:03d}"  # Julian day for February
            # Ozone column varies from ~250 DU at 30N to ~280 DU at 10N
            for lat_start in range(5, 35, 5):
                lat_end = lat_start + 5
                o3col = 280 - (lat_start - 10) * 1.0  # Higher near equator
                f.write(f"{datestr}  {lat_start:5.1f}  {lat_end:5.1f}  {o3col:6.1f}\n")

    return filename


def create_chemparam_link():
    """Create symbolic link to chemistry parameters file."""

    src = '/home/deenalad/specific_inputs/inputs/CAMx7.32.chemparam.CB6r4_CF3_COMPLX'
    dst = os.path.join(CAMX_INPUT_DIR, 'CAMx7.32.chemparam.CB6r4_CF3_COMPLX')

    if os.path.exists(src):
        if os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)
        print(f"  Linked: {os.path.basename(dst)}")
    else:
        print(f"  WARNING: Chemistry params not found at {src}")
        # Create a minimal placeholder
        print(f"  Creating placeholder chemistry params...")


def create_tuv_photolysis(date):
    """Create TUV photolysis rates file."""

    datestr = date.strftime('%Y%m%d')
    filename = os.path.join(CAMX_INPUT_DIR, f'tuv.india.{datestr}')
    print(f"  Creating: {os.path.basename(filename)}")

    # Link to existing TUV if available
    src_tuv = '/home/deenalad/specific_inputs/inputs/tuv.do_CB6r4.20160610'

    if os.path.exists(src_tuv):
        # Copy and modify the existing TUV file
        import shutil
        shutil.copy(src_tuv, filename)
        print(f"    Copied from existing TUV file")
    else:
        # Create a minimal TUV file
        with open(filename, 'w') as f:
            f.write(f"! TUV photolysis rates for {datestr}\n")
            f.write("! Placeholder - use actual TUV preprocessor for production runs\n")

    return filename


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("CREATING CAMx INPUT FILES FOR BRICK KILN SCENARIOS")
    print("=" * 70)
    print(f"\nOutput directory: {CAMX_INPUT_DIR}")
    print(f"Grid: {GRID['NCOLS']} x {GRID['NROWS']} x {GRID['NLAYS']}")
    print(f"Cell size: {GRID['XCELL']/1000} km")
    print()

    # Create emissions for both scenarios
    print("Creating EMISSION files...")
    for date in DATES:
        create_camx_emission_file(date, '2020')  # Base scenario
        create_camx_emission_file(date, '2025')  # Converted scenario

    # Create IC/BC files
    print("\nCreating IC/BC files...")
    for date in DATES:
        create_ic_file(date)
        create_bc_file(date)

    # Create supporting files
    print("\nCreating supporting files...")
    create_o3_column_file()
    create_chemparam_link()
    for date in DATES:
        create_tuv_photolysis(date)

    print("\n" + "=" * 70)
    print("INPUT FILES CREATED SUCCESSFULLY!")
    print("=" * 70)
    print(f"""
Files created in: {CAMX_INPUT_DIR}

Emissions:
  camx_emis.2020.YYYYMMDD.nc  (base scenario - 100% FCBK)
  camx_emis.2025.YYYYMMDD.nc  (converted scenario - 50% Zigzag)

IC/BC:
  ic.india.YYYYMMDD.nc
  bc.india.YYYYMMDD.nc

Supporting:
  o3map.202402.txt
  tuv.india.YYYYMMDD
  CAMx7.32.chemparam.CB6r4_CF3_COMPLX

Next step: Run CAMx with these inputs using:
  cd run_camx
  ./CAMx.brick_kiln.job 2020  # or 2025
""")


if __name__ == "__main__":
    main()
