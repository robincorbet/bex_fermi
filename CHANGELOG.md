# Changelog

## Version 2.1 (2026-02-16)

### New Features
- **FL16Y Catalog Support**: Updated default catalog to `gll_psc_v40.fit` (16-year Fermi LAT source list)
- **Generalized Source Name Handling**: Improved regex-based pattern matching for catalog prefixes
  - Now handles all FLnY catalogs (FL8Y, FL16Y, FL20Y, future releases)
  - Handles all xFGL catalogs (1FGL through 5FGL and beyond)
  - Handles time-period catalogs (18M, 24M, etc.)
  - No need to update code for future catalog releases

### Bug Fixes
- Fixed PFILES collision when running multiple instances simultaneously
- Fixed parameter file behavior: `-file` now prompts with file values as defaults
- Fixed gtbary to use `outfile` parameter instead of modifying in-place
- Fixed catalog path validation to check MY_FERMI_DIR and FERMI_DIR
- Added comprehensive input file validation with helpful error messages
- Added cleanup of temporary light curve files (.fits and .dmp1)

### Performance Improvements
- **Major speedup**: Vectorized pweight algorithm using NumPy
  - 100-1000x faster for large datasets
  - Process 200k+ photons in seconds instead of hours

### Documentation
- Updated all examples to use gll_psc_v40.fit
- Added FL16Y catalog information
- Improved error messages for missing files

## Version 2.0 (2026-02-09)

### Initial Release
- Complete Python conversion of bex2020nd Perl script and pweight C program
- Full Fermi LAT analysis pipeline
- Interactive parameter prompting
- Automatic XML model generation with make4FGLxml
- Probability-weighted photometry
- All IRF versions supported (P6 through P8R3/v2.0)
- Batch and update modes
- Comprehensive documentation

## Fermi LAT Catalog History

- **FL16Y** (v40): 16-year catalog (2008-2024)
  - https://fermi.gsfc.nasa.gov/ssc/data/access/lat/fl16y/
- **4FGL-DR4** (v35): 14-year catalog (2008-2022)
- **4FGL-DR3** (v31): 12-year catalog
- **4FGL-DR2** (v28): 10-year catalog
- **3FGL** (v18): 4-year catalog
- **2FGL** (v8): 2-year catalog
- **1FGL** (v3): 1-year catalog

## Supported Catalog Naming Conventions

The script automatically handles source names from:
- **xFGL**: 1FGL, 2FGL, 3FGL, 4FGL, 5FGL, ...
- **FLnY**: FL8Y, FL16Y, FL20Y, ...
- **Time periods**: 18M, 24M, ...

All catalog identifiers are automatically prefixed with underscore for FTOOLS compatibility.

Example transformations:
- `4FGLJ2045.0+4202` → `_4FGLJ2045.0+4202`
- `FL16YJ2045.0+4202` → `_FL16YJ2045.0+4202`
- `18MJ2045.0+4202` → `_18MJ2045.0+4202`
