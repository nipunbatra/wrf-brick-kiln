#!/usr/bin/env python3
"""
Download EDGAR v8.1 emissions for India CAMx domain

EDGAR Database: https://edgar.jrc.ec.europa.eu/
Version: v8.1 (2024)

This script downloads gridded emissions data for:
- PM2.5, PM10, BC, OC (particulate matter)
- NOx, SO2, CO, NH3 (gaseous pollutants)
- NMVOC (for secondary organic aerosols)

Sectors relevant to brick kilns:
- NMM: Non-Metallic Minerals (includes brick production)
- IND: Industrial combustion
- RCO: Residential and other sectors
"""

import os
import sys
import urllib.request
from pathlib import Path

# EDGAR base URL for v8.1
EDGAR_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v81_FT2022_GHG"

# Species and their EDGAR codes
SPECIES = {
    'PM2.5': 'PM2.5',
    'PM10': 'PM10',
    'BC': 'BC',
    'OC': 'OC',
    'NOx': 'NOx',
    'SO2': 'SO2',
    'CO': 'CO',
    'NH3': 'NH3',
}

# Sectors
SECTORS = ['TOTAL', 'NMM', 'IND', 'RCO', 'ENE', 'TRO']

# Years to download
YEARS = [2022]  # Latest available

def get_edgar_url(species, sector, year):
    """
    Construct EDGAR download URL
    Note: Actual URLs may vary - check EDGAR website for exact paths
    """
    # EDGAR v8.1 URL pattern (example - verify on EDGAR website)
    # Format varies by product type
    url = f"https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v81_FT2022_AP/{species}_v81_FT2022_{sector}_{year}_TOTALS_emi.nc"
    return url

def download_file(url, output_path):
    """Download a file with progress indicator"""
    try:
        print(f"Downloading: {url}")
        print(f"       To: {output_path}")
        urllib.request.urlretrieve(url, output_path)
        print(f"  -> Success!")
        return True
    except Exception as e:
        print(f"  -> Failed: {e}")
        return False

def main():
    output_dir = Path(__file__).parent / "raw_data"
    output_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("EDGAR Emissions Download Script")
    print("=" * 70)
    print("""
NOTE: EDGAR data requires registration or direct download from:
      https://edgar.jrc.ec.europa.eu/

Alternative approach - download manually:
1. Go to: https://edgar.jrc.ec.europa.eu/dataset_ap81
2. Select species: PM2.5, NOx, SO2, etc.
3. Select sector: Total or specific (NMM for brick kilns)
4. Download NetCDF files
5. Place in: {output_dir}

For this project, we need emissions for:
- Domain: 67°E-92°E, 6°N-30°N (India region)
- Resolution: Will regrid to 27km CAMx grid
- Year: 2022 (or most recent available)
- Species: PM2.5, PM10, BC, OC, NOx, SO2, CO, NH3, NMVOC
""".format(output_dir=output_dir))

    # Create a list of files needed
    files_needed = []
    for species in SPECIES:
        for sector in ['TOTAL']:  # Start with totals
            for year in YEARS:
                files_needed.append({
                    'species': species,
                    'sector': sector,
                    'year': year,
                    'filename': f"EDGAR_v81_{species}_{sector}_{year}.nc"
                })

    print("\nFiles needed:")
    print("-" * 50)
    for f in files_needed:
        print(f"  - {f['filename']}")

    print(f"\nOutput directory: {output_dir}")
    print("\nTo proceed with automatic download, uncomment the download section")
    print("and ensure you have network access to EDGAR servers.")

    # Uncomment below to attempt automatic download
    # (May require authentication or different URL patterns)
    """
    for f in files_needed:
        url = get_edgar_url(f['species'], f['sector'], f['year'])
        output_file = output_dir / f['filename']
        if not output_file.exists():
            download_file(url, output_file)
        else:
            print(f"File exists: {f['filename']}")
    """

if __name__ == "__main__":
    main()
