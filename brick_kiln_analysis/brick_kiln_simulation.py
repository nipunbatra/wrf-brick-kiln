#!/usr/bin/env python3
"""
Brick Kiln Air Quality Impact Analysis
======================================

Simulates PM2.5 reduction in Delhi from converting brick kilns in UP.

Scenario:
- 2020 (Base): 10,000 brick kilns in UP, all FCBTK
- 2025 (Improved): 50% converted to Zigzag kilns
"""

import numpy as np
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# =============================================================================
# GRID DEFINITION - Matching WRFCAMx domain (India-centric)
# =============================================================================
# Domain: ~67°E to 93°E, ~6°N to 30°N
NX, NY = 102, 102
LON_MIN, LON_MAX = 67.0, 93.0
LAT_MIN, LAT_MAX = 6.0, 30.0

# Create coordinate arrays
LONS = np.linspace(LON_MIN, LON_MAX, NX)
LATS = np.linspace(LAT_MIN, LAT_MAX, NY)
LON_GRID, LAT_GRID = np.meshgrid(LONS, LATS)
DX = (LON_MAX - LON_MIN) / NX  # ~0.25 degrees
DY = (LAT_MAX - LAT_MIN) / NY

# Key locations
LOCATIONS = {
    'Delhi': {'lat': 28.6, 'lon': 77.2},
    'Lucknow': {'lat': 26.8, 'lon': 80.9},
    'Agra': {'lat': 27.2, 'lon': 78.0},
    'Kanpur': {'lat': 26.4, 'lon': 80.3},
}

# =============================================================================
# BRICK KILN PARAMETERS
# =============================================================================
BRICK_KILN = {
    'bricks_per_kiln_per_day': 25000,
    'operating_days_per_year': 180,
    'brick_weight_kg': 3.0,

    # Emission factors (g/kg brick)
    'FCBTK': {'PM2.5': 0.75, 'PM10': 1.2, 'BC': 0.15, 'NOx': 0.30, 'SO2': 0.80},
    'Zigzag': {'PM2.5': 0.30, 'PM10': 0.48, 'BC': 0.06, 'NOx': 0.25, 'SO2': 0.65},
}

# UP Brick Kiln Clusters (10,000 total)
UP_KILN_CLUSTERS = [
    {'name': 'Western UP (near Delhi)', 'lat': 28.5, 'lon': 77.8, 'kilns': 2500},
    {'name': 'Central UP', 'lat': 26.8, 'lon': 80.5, 'kilns': 2000},
    {'name': 'Agra Region', 'lat': 27.2, 'lon': 78.5, 'kilns': 1500},
    {'name': 'Kanpur Region', 'lat': 26.5, 'lon': 80.0, 'kilns': 1500},
    {'name': 'Eastern UP', 'lat': 25.5, 'lon': 82.5, 'kilns': 1500},
    {'name': 'Allahabad Region', 'lat': 25.4, 'lon': 81.8, 'kilns': 1000},
]

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def latlon_to_idx(lon, lat):
    """Convert lat/lon to grid indices"""
    i = int((lon - LON_MIN) / DX)
    j = int((lat - LAT_MIN) / DY)
    i = max(0, min(NX-1, i))
    j = max(0, min(NY-1, j))
    return i, j

def calculate_emissions(num_kilns, kiln_type='FCBTK'):
    """Calculate annual emissions (tonnes/year) and emission rate (g/s)"""
    ef = BRICK_KILN[kiln_type]
    daily_prod_kg = BRICK_KILN['bricks_per_kiln_per_day'] * BRICK_KILN['brick_weight_kg']
    annual_prod_kg = daily_prod_kg * BRICK_KILN['operating_days_per_year'] * num_kilns

    emissions = {}
    for species, ef_val in ef.items():
        annual_tonnes = (annual_prod_kg * ef_val) / 1e6
        operating_seconds = BRICK_KILN['operating_days_per_year'] * 24 * 3600
        rate_gs = (annual_tonnes * 1e6) / operating_seconds
        emissions[species] = {'tonnes_yr': annual_tonnes, 'g_s': rate_gs}
    return emissions

def create_emission_grid(clusters, kiln_type='FCBTK'):
    """Create gridded emissions from clusters"""
    species_list = list(BRICK_KILN[kiln_type].keys())
    grids = {sp: np.zeros((NY, NX)) for sp in species_list}

    for cluster in clusters:
        i, j = latlon_to_idx(cluster['lon'], cluster['lat'])
        emissions = calculate_emissions(cluster['kilns'], kiln_type)

        # Spread over 5x5 cells with Gaussian weighting
        for di in range(-2, 3):
            for dj in range(-2, 3):
                ni, nj = i + di, j + dj
                if 0 <= ni < NX and 0 <= nj < NY:
                    weight = np.exp(-(di**2 + dj**2) / 2.0) / 10.0
                    for sp in species_list:
                        grids[sp][nj, ni] += emissions[sp]['g_s'] * weight
    return grids

def simple_dispersion(emission_grid, wind_speed=3.0, mixing_height=1000):
    """Simple box model dispersion - converts g/s to µg/m³"""
    # Cell area in m²
    cell_area = (DX * 111000) * (DY * 111000)  # approx m²

    # Box volume (cell area × mixing height)
    volume = cell_area * mixing_height  # m³

    # Residence time (seconds) - assume air moves through in ~6 hours
    residence_time = 6 * 3600

    # Concentration = (emission rate × residence time) / volume
    # Plus downwind transport (simple smoothing)
    conc = (emission_grid * residence_time / volume) * 1e6  # µg/m³

    # Add atmospheric transport (smooth and spread)
    from scipy.ndimage import gaussian_filter
    conc = gaussian_filter(conc, sigma=3)

    # Add prevailing wind effect (winter: NW to SE)
    # Shift concentrations southeast
    conc_shifted = np.roll(np.roll(conc, 2, axis=0), -2, axis=1) * 0.3
    conc = conc + conc_shifted

    return conc

# =============================================================================
# MAIN ANALYSIS
# =============================================================================
def run_analysis():
    print("=" * 70)
    print("BRICK KILN AIR QUALITY IMPACT ANALYSIS")
    print("UP Brick Kilns → Delhi PM2.5")
    print("=" * 70)

    total_kilns = sum(c['kilns'] for c in UP_KILN_CLUSTERS)
    print(f"\nTotal brick kilns in UP: {total_kilns:,}")

    # ==========================================================================
    # EMISSIONS
    # ==========================================================================
    print("\n" + "-" * 50)
    print("EMISSIONS INVENTORY")
    print("-" * 50)

    # 2020: All FCBTK
    emis_2020 = calculate_emissions(total_kilns, 'FCBTK')
    print(f"\n2020 (100% FCBTK):")
    print(f"  PM2.5: {emis_2020['PM2.5']['tonnes_yr']:,.0f} tonnes/yr")
    print(f"  PM10:  {emis_2020['PM10']['tonnes_yr']:,.0f} tonnes/yr")
    print(f"  BC:    {emis_2020['BC']['tonnes_yr']:,.0f} tonnes/yr")

    # 2025: 50% Zigzag
    emis_fcbtk = calculate_emissions(total_kilns // 2, 'FCBTK')
    emis_zigzag = calculate_emissions(total_kilns // 2, 'Zigzag')
    emis_2025 = {sp: {'tonnes_yr': emis_fcbtk[sp]['tonnes_yr'] + emis_zigzag[sp]['tonnes_yr'],
                      'g_s': emis_fcbtk[sp]['g_s'] + emis_zigzag[sp]['g_s']}
                 for sp in emis_2020}

    print(f"\n2025 (50% FCBTK + 50% Zigzag):")
    print(f"  PM2.5: {emis_2025['PM2.5']['tonnes_yr']:,.0f} tonnes/yr")
    print(f"  PM10:  {emis_2025['PM10']['tonnes_yr']:,.0f} tonnes/yr")
    print(f"  BC:    {emis_2025['BC']['tonnes_yr']:,.0f} tonnes/yr")

    # Reductions
    pm25_reduction = emis_2020['PM2.5']['tonnes_yr'] - emis_2025['PM2.5']['tonnes_yr']
    pm25_pct = (pm25_reduction / emis_2020['PM2.5']['tonnes_yr']) * 100
    print(f"\nPM2.5 Emission Reduction: {pm25_reduction:,.0f} tonnes/yr ({pm25_pct:.0f}%)")

    # ==========================================================================
    # GRIDDED EMISSIONS & CONCENTRATIONS
    # ==========================================================================
    print("\n" + "-" * 50)
    print("CONCENTRATION MODELING")
    print("-" * 50)

    # Create grids
    grid_2020_fcbtk = create_emission_grid(UP_KILN_CLUSTERS, 'FCBTK')

    # 2025: split clusters
    clusters_half = [{**c, 'kilns': c['kilns'] // 2} for c in UP_KILN_CLUSTERS]
    grid_2025_fcbtk = create_emission_grid(clusters_half, 'FCBTK')
    grid_2025_zigzag = create_emission_grid(clusters_half, 'Zigzag')
    grid_2025 = {sp: grid_2025_fcbtk[sp] + grid_2025_zigzag[sp] for sp in grid_2020_fcbtk}

    # Calculate concentrations
    try:
        from scipy.ndimage import gaussian_filter
        conc_2020 = simple_dispersion(grid_2020_fcbtk['PM2.5'])
        conc_2025 = simple_dispersion(grid_2025['PM2.5'])
    except ImportError:
        # Fallback without scipy
        print("  (scipy not available - using simplified model)")
        cell_area = (DX * 111000) * (DY * 111000)
        volume = cell_area * 1000
        residence_time = 6 * 3600
        conc_2020 = (grid_2020_fcbtk['PM2.5'] * residence_time / volume) * 1e6
        conc_2025 = (grid_2025['PM2.5'] * residence_time / volume) * 1e6

    # Add background PM2.5 (Delhi winter: ~100 µg/m³ from other sources)
    background_pm25 = 80.0
    conc_2020 += background_pm25
    conc_2025 += background_pm25

    # ==========================================================================
    # DELHI ANALYSIS
    # ==========================================================================
    print("\n" + "-" * 50)
    print("DELHI PM2.5 IMPACT")
    print("-" * 50)

    delhi_i, delhi_j = latlon_to_idx(LOCATIONS['Delhi']['lon'], LOCATIONS['Delhi']['lat'])
    print(f"\nDelhi grid cell: ({delhi_i}, {delhi_j})")
    print(f"Delhi coordinates: {LOCATIONS['Delhi']['lat']}°N, {LOCATIONS['Delhi']['lon']}°E")

    # Get Delhi concentration (average over 3x3 neighborhood)
    j_min, j_max = max(0, delhi_j-1), min(NY, delhi_j+2)
    i_min, i_max = max(0, delhi_i-1), min(NX, delhi_i+2)

    delhi_pm25_2020 = conc_2020[j_min:j_max, i_min:i_max].mean()
    delhi_pm25_2025 = conc_2025[j_min:j_max, i_min:i_max].mean()

    brick_contrib_2020 = delhi_pm25_2020 - background_pm25
    brick_contrib_2025 = delhi_pm25_2025 - background_pm25

    print(f"\n  Background PM2.5: {background_pm25:.0f} µg/m³")
    print(f"\n  2020 (100% FCBTK):")
    print(f"    Brick kiln contribution: {brick_contrib_2020:.1f} µg/m³")
    print(f"    Total PM2.5: {delhi_pm25_2020:.1f} µg/m³")

    print(f"\n  2025 (50% Zigzag):")
    print(f"    Brick kiln contribution: {brick_contrib_2025:.1f} µg/m³")
    print(f"    Total PM2.5: {delhi_pm25_2025:.1f} µg/m³")

    reduction = delhi_pm25_2020 - delhi_pm25_2025
    pct_reduction = (reduction / delhi_pm25_2020) * 100 if delhi_pm25_2020 > 0 else 0
    brick_pct = ((brick_contrib_2020 - brick_contrib_2025) / brick_contrib_2020 * 100
                 if brick_contrib_2020 > 0 else 0)

    print(f"\n  ═══════════════════════════════════════")
    print(f"  PM2.5 REDUCTION IN DELHI")
    print(f"  ═══════════════════════════════════════")
    print(f"  Absolute: {reduction:.1f} µg/m³")
    print(f"  Relative (total): {pct_reduction:.1f}%")
    print(f"  Relative (brick kiln only): {brick_pct:.1f}%")

    # ==========================================================================
    # HEALTH IMPACT
    # ==========================================================================
    print("\n" + "-" * 50)
    print("HEALTH IMPACT ESTIMATE")
    print("-" * 50)

    delhi_pop = 32_000_000
    # WHO: ~6% mortality increase per 10 µg/m³ PM2.5
    mortality_coef = 0.006  # per µg/m³
    baseline_mortality = 0.008  # per year

    deaths_attr_2020 = delhi_pop * baseline_mortality * mortality_coef * brick_contrib_2020
    deaths_attr_2025 = delhi_pop * baseline_mortality * mortality_coef * brick_contrib_2025
    lives_saved = deaths_attr_2020 - deaths_attr_2025

    print(f"\n  Delhi Population: {delhi_pop:,}")
    print(f"\n  Premature deaths from brick kiln PM2.5:")
    print(f"    2020: ~{deaths_attr_2020:,.0f}/year")
    print(f"    2025: ~{deaths_attr_2025:,.0f}/year")
    print(f"\n  ★ LIVES SAVED: ~{lives_saved:,.0f}/year ★")

    # ==========================================================================
    # VISUALIZATION
    # ==========================================================================
    if HAS_MATPLOTLIB:
        print("\n" + "-" * 50)
        print("CREATING VISUALIZATIONS")
        print("-" * 50)
        create_plots(grid_2020_fcbtk, grid_2025, conc_2020, conc_2025,
                    emis_2020, emis_2025, delhi_pm25_2020, delhi_pm25_2025,
                    background_pm25, lives_saved)

    # ==========================================================================
    # FINAL SUMMARY
    # ==========================================================================
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"""
┌────────────────────────────────────────────────────────────────────┐
│           BRICK KILN CONVERSION IMPACT ON DELHI                    │
├────────────────────────────────────────────────────────────────────┤
│  Location: Uttar Pradesh, India                                    │
│  Total Kilns: {total_kilns:,}                                            │
│  Scenario: 50% conversion from FCBTK to Zigzag                     │
├────────────────────────────────────────────────────────────────────┤
│                          │   2020      │   2025      │  Change     │
├────────────────────────────────────────────────────────────────────┤
│  Kiln Mix                │ 100% FCBTK  │ 50% Zigzag  │     -       │
│  PM2.5 Emissions (t/yr)  │ {emis_2020['PM2.5']['tonnes_yr']:>9,.0f}  │ {emis_2025['PM2.5']['tonnes_yr']:>9,.0f}  │ -{pm25_reduction:,.0f}     │
│  Delhi Total PM2.5       │ {delhi_pm25_2020:>9.1f}  │ {delhi_pm25_2025:>9.1f}  │ -{reduction:.1f}       │
│  Delhi Brick Kiln PM2.5  │ {brick_contrib_2020:>9.1f}  │ {brick_contrib_2025:>9.1f}  │ -{brick_contrib_2020-brick_contrib_2025:.1f}       │
├────────────────────────────────────────────────────────────────────┤
│  ★ ESTIMATED LIVES SAVED IN DELHI: ~{lives_saved:,.0f}/year               │
└────────────────────────────────────────────────────────────────────┘
""")

    return {
        'emis_2020': emis_2020,
        'emis_2025': emis_2025,
        'delhi_pm25_2020': delhi_pm25_2020,
        'delhi_pm25_2025': delhi_pm25_2025,
        'reduction': reduction,
        'lives_saved': lives_saved,
    }

def create_plots(grid_2020, grid_2025, conc_2020, conc_2025,
                 emis_2020, emis_2025, delhi_2020, delhi_2025,
                 background, lives_saved):
    """Create visualization plots"""
    output_dir = os.path.dirname(os.path.abspath(__file__))

    # Figure 1: Emissions and Concentrations
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Brick Kiln PM2.5 Impact Analysis: UP → Delhi', fontsize=14, fontweight='bold')

    # Emission grid 2020
    ax = axes[0, 0]
    im = ax.pcolormesh(LON_GRID, LAT_GRID, grid_2020['PM2.5']*1000, cmap='YlOrRd', shading='auto')
    ax.plot(LOCATIONS['Delhi']['lon'], LOCATIONS['Delhi']['lat'], 'k*', ms=15, label='Delhi')
    ax.set_title('2020 Emissions (mg/s/cell)')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    plt.colorbar(im, ax=ax)
    ax.legend()

    # Emission grid 2025
    ax = axes[0, 1]
    im = ax.pcolormesh(LON_GRID, LAT_GRID, grid_2025['PM2.5']*1000, cmap='YlOrRd', shading='auto')
    ax.plot(LOCATIONS['Delhi']['lon'], LOCATIONS['Delhi']['lat'], 'k*', ms=15)
    ax.set_title('2025 Emissions (mg/s/cell)')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    plt.colorbar(im, ax=ax)

    # Emission reduction
    ax = axes[0, 2]
    reduction_grid = (grid_2020['PM2.5'] - grid_2025['PM2.5']) * 1000
    im = ax.pcolormesh(LON_GRID, LAT_GRID, reduction_grid, cmap='Greens', shading='auto')
    ax.plot(LOCATIONS['Delhi']['lon'], LOCATIONS['Delhi']['lat'], 'k*', ms=15)
    ax.set_title('Emission Reduction (mg/s/cell)')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    plt.colorbar(im, ax=ax)

    # Concentration 2020
    ax = axes[1, 0]
    im = ax.pcolormesh(LON_GRID, LAT_GRID, conc_2020, cmap='RdYlGn_r', shading='auto', vmin=0, vmax=200)
    ax.plot(LOCATIONS['Delhi']['lon'], LOCATIONS['Delhi']['lat'], 'k*', ms=15)
    ax.set_title(f'2020 PM2.5 Conc. (µg/m³)')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    plt.colorbar(im, ax=ax)

    # Concentration 2025
    ax = axes[1, 1]
    im = ax.pcolormesh(LON_GRID, LAT_GRID, conc_2025, cmap='RdYlGn_r', shading='auto', vmin=0, vmax=200)
    ax.plot(LOCATIONS['Delhi']['lon'], LOCATIONS['Delhi']['lat'], 'k*', ms=15)
    ax.set_title(f'2025 PM2.5 Conc. (µg/m³)')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    plt.colorbar(im, ax=ax)

    # Delhi bar chart
    ax = axes[1, 2]
    x = ['2020\n(100% FCBTK)', '2025\n(50% Zigzag)']
    brick = [delhi_2020 - background, delhi_2025 - background]
    ax.bar(x, [background, background], label='Background', color='gray', alpha=0.7)
    ax.bar(x, brick, bottom=[background, background],
           label='Brick Kilns', color=['firebrick', 'forestgreen'])
    ax.set_ylabel('PM2.5 (µg/m³)')
    ax.set_title('Delhi PM2.5 Breakdown')
    ax.legend()
    ax.set_ylim(0, 180)
    ax.axhline(60, color='orange', ls='--', alpha=0.7, label='WHO IT-3')
    for i, (b, bk) in enumerate(zip([background]*2, brick)):
        ax.text(i, b+bk+5, f'{b+bk:.0f}', ha='center', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'brick_kiln_analysis.png'), dpi=150, bbox_inches='tight')
    print(f"  Saved: brick_kiln_analysis.png")

    # Figure 2: Summary infographic
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.axis('off')

    summary = f"""
    ╔══════════════════════════════════════════════════════════════════╗
    ║          BRICK KILN CONVERSION IMPACT ON DELHI AIR QUALITY        ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║                                                                   ║
    ║   SCENARIO:                                                       ║
    ║   • Location: Uttar Pradesh, India                                ║
    ║   • 10,000 brick kilns                                           ║
    ║   • 50% converted from FCBTK → Zigzag (2020→2025)                ║
    ║                                                                   ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║                                                                   ║
    ║   EMISSION REDUCTIONS:                                            ║
    ║   • PM2.5: {emis_2020['PM2.5']['tonnes_yr']-emis_2025['PM2.5']['tonnes_yr']:,.0f} tonnes/year (30% reduction)                  ║
    ║   • BC:    {emis_2020['BC']['tonnes_yr']-emis_2025['BC']['tonnes_yr']:,.0f} tonnes/year (30% reduction)                   ║
    ║                                                                   ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║                                                                   ║
    ║   DELHI PM2.5 IMPACT:                                             ║
    ║   ┌─────────────────┬────────────┬────────────┐                  ║
    ║   │                 │    2020    │    2025    │                  ║
    ║   ├─────────────────┼────────────┼────────────┤                  ║
    ║   │ Total PM2.5     │ {delhi_2020:>6.0f} µg/m³│ {delhi_2025:>6.0f} µg/m³│                  ║
    ║   │ Brick Kiln Only │ {delhi_2020-background:>6.0f} µg/m³│ {delhi_2025-background:>6.0f} µg/m³│                  ║
    ║   └─────────────────┴────────────┴────────────┘                  ║
    ║                                                                   ║
    ║   REDUCTION: {delhi_2020-delhi_2025:.0f} µg/m³ ({(delhi_2020-delhi_2025)/delhi_2020*100:.0f}% of total)                       ║
    ║                                                                   ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║                                                                   ║
    ║   ★★★ ESTIMATED LIVES SAVED IN DELHI: ~{lives_saved:,.0f}/year ★★★         ║
    ║                                                                   ║
    ╚══════════════════════════════════════════════════════════════════╝
    """

    ax.text(0.5, 0.5, summary, transform=ax.transAxes, fontsize=11,
            verticalalignment='center', horizontalalignment='center',
            fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow'))

    plt.savefig(os.path.join(output_dir, 'summary_infographic.png'), dpi=150, bbox_inches='tight')
    print(f"  Saved: summary_infographic.png")
    plt.close('all')

if __name__ == "__main__":
    results = run_analysis()
