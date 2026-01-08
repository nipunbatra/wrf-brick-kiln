# Project Development Log

This document captures the key inputs, decisions, and thought process during the development of this brick kiln air quality impact analysis project.

---

## Project Genesis

### User Request
> "run ./test_camx.job. This project is on WRF-CAMX"

The initial goal was to run WRFCAMx (a preprocessor that converts WRF meteorological output to CAMx air quality model format) for an India domain.

---

## Technical Challenges & Solutions

### 1. WRFCAMx Configuration Issues

**Problem:** Multiple configuration mismatches between WRF output and WRFCAMx settings.

**Issues encountered and fixed:**
- **Heredoc delimiter mismatch**: Script used `ieo` but closed with `ieof`
- **Snow file format**: Expected binary Fortran unformatted, got text file with "36"
- **Date format**: Script used `YYMMDD HHMM` but code expected `YYMMDDHH`
- **Kv method**: MYJ method requires TKE variable not present in WRF output → switched to CMAQ diagnostic method
- **Grid domain mismatch**: CAMx grid was outside WRF domain → aligned to 100x100 cells at 27km resolution
- **WRF output frequency**: Script assumed 60 min but actual was 180 min (3-hourly)
- **Timezone**: Changed from 8 (Pacific) to 0 (UTC) for India domain

**Solution approach:** Iteratively ran the preprocessor, analyzed error messages, and adjusted configuration parameters until successful.

### 2. Grid Alignment

**Key insight:** For WRFCAMx to work, the CAMx grid must be completely contained within the WRF domain, and grid cell boundaries must align.

**Final configuration:**
```
Grid size: 100 x 100 cells
Resolution: 27 km x 27 km
Origin: -1714.5, -1714.5 km (Mercator projection)
Reference: 83°E, 21.5°N (center of India)
```

---

## Broader Project Goals

### User's Vision
> "scenario default: assume some emissions in a grid, scenario new: assume reduced emissions due to brick kilns converted. I have for each type of brick kilns the emission numbers."

> "I would like to use EDGAR, help! and CAMX is in src.v7.32. Get me the total pipeline."

> "assume 10,000 kilns in UP, India and in 2020 all were FCBK, then assume in 2025 50% are Zigzag. Assume some constants and show reduction in PM in Delhi."

### Interpretation
The user wanted to:
1. Use EDGAR global emissions inventory as baseline
2. Add brick kiln emissions as a point source layer
3. Compare scenarios: traditional FCBTK kilns vs cleaner Zigzag technology
4. Quantify the air quality improvement in Delhi from technology conversion

---

## Key Design Decisions

### 1. Brick Kiln Emission Factors

Based on published literature (CPCB India, GIZ studies):

| Kiln Type | PM2.5 (g/kg brick) | Relative |
|-----------|-------------------|----------|
| FCBTK (Fixed Chimney Bull Trench) | 0.75 | Baseline |
| Zigzag | 0.30 | 60% cleaner |
| VSBK (Vertical Shaft) | 0.15 | 80% cleaner |
| Tunnel Kiln | 0.08 | 89% cleaner |

**Rationale:** FCBTK is the dominant technology in India (~70% of kilns). Zigzag is the most realistic transition target as it requires minimal capital investment.

### 2. Kiln Distribution in UP

Created synthetic inventory of 10,000 kilns across 6 clusters:

| Cluster | Kilns | Distance to Delhi |
|---------|-------|-------------------|
| Western UP (Ghaziabad area) | 2,500 | ~50 km |
| Ghaziabad-Noida corridor | 2,000 | ~25 km |
| Meerut | 1,500 | ~70 km |
| Agra | 1,500 | ~200 km |
| Lucknow | 1,500 | ~500 km |
| Varanasi | 1,000 | ~800 km |

**Rationale:** Brick kilns cluster near construction markets. Western UP kilns have the greatest impact on Delhi due to proximity and prevailing wind patterns.

### 3. Dispersion Modeling Approach

Used simplified Gaussian plume model rather than full CAMx simulation for the demonstration:

```python
# Distance-based decay with wind direction weighting
concentration = emission_rate / (distance^1.5) * wind_factor * mixing_height_factor
```

**Rationale:**
- Full CAMx run requires IC/BC files, photolysis rates, and significant compute time
- Gaussian approximation captures first-order effects for policy demonstration
- Results are directionally correct even if not quantitatively precise

### 4. Health Impact Assessment

Used integrated exposure-response function:
```
Relative Risk = 1.06 per 10 µg/m³ PM2.5 increase
Delhi Population = 32 million
```

**Calculation:**
- Brick kiln contribution to Delhi PM2.5: ~0.5 µg/m³ (2020) → ~0.4 µg/m³ (2025)
- Attributable deaths: ~789/year → ~553/year
- **Lives saved: ~237/year**

---

## Results Summary

### Emission Reductions (50% Zigzag conversion)

| Pollutant | 2020 (100% FCBTK) | 2025 (50% Zigzag) | Reduction |
|-----------|-------------------|-------------------|-----------|
| PM2.5 | 101,250 t/yr | 70,875 t/yr | 30% |
| PM10 | 162,000 t/yr | 113,400 t/yr | 30% |
| Black Carbon | 20,250 t/yr | 14,175 t/yr | 30% |

### Delhi Air Quality Impact

| Metric | Value |
|--------|-------|
| Background PM2.5 | 80 µg/m³ |
| Brick kiln contribution (2020) | 0.5 µg/m³ |
| Brick kiln contribution (2025) | 0.4 µg/m³ |
| Absolute reduction | 0.2 µg/m³ |
| Relative reduction (brick kiln) | 30% |
| Estimated lives saved | ~237/year |

### Key Insight

While the absolute PM2.5 reduction in Delhi appears small (~0.2 µg/m³), this translates to significant health benefits at population scale. The 30% reduction in brick kiln emissions directly translates to 30% reduction in brick kiln-attributable mortality.

---

## Repository Structure Decisions

### What's Included

1. **CAMx v7.32 source code** (~13 MB) - Allows users to compile and run full simulations
2. **WRFCAMx v5.2** with India configuration - Ready-to-use preprocessor
3. **5 days of processed CAMx meteorology** (~420 MB) - Allows running without WRF
4. **1 day sample WRF output** (~600 MB) - For demonstration/testing
5. **EDGAR processing scripts** - Template for adding global emissions
6. **Brick kiln simulation** - Standalone analysis script

### What's NOT Included

1. **Full WRF output** (144 GB) - Too large, users must provide their own
2. **EDGAR raw data** - Must be downloaded from JRC website
3. **IC/BC files for CAMx** - Would need to be generated from global models

### Git LFS Usage

Large files tracked via Git LFS:
- `*.nc` (NetCDF meteorology files)
- `wrf_sample/*` (WRF output files)

This keeps the repository clone size manageable while preserving large data files.

---

## Commands for GitHub Setup

```bash
cd /home/nipun.batra/git/wrf-brick-kiln

# Authenticate with GitHub
gh auth login

# Create and push repository
gh repo create wrf-brick-kiln --public \
  --description "Brick kiln air quality impact analysis: WRF-CAMx modeling for Delhi PM2.5" \
  --source=. --push
```

---

## Future Improvements

1. **Real brick kiln locations** - Use satellite imagery or survey data
2. **Seasonal variation** - Kilns operate mainly Oct-Jun
3. **Full CAMx chemistry** - Include secondary PM formation
4. **Higher resolution** - Nest 9km or 3km domain over Delhi
5. **Validation** - Compare with CPCB monitoring station data
6. **Multi-city analysis** - Extend to Lucknow, Patna, Kolkata

---

## Acknowledgments

- WRF output provided from existing simulation
- CAMx v7.32 from RAMBOLL (open source)
- EDGAR emissions from JRC
- Brick kiln emission factors from CPCB India and GIZ studies

---

*Generated during project development session, January 2024*
