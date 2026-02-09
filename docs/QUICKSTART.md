# Quick Start Guide

## Installation

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Set up Fermi Science Tools:**
   ```bash
   # Using conda (recommended)
   conda activate fermitools
   
   # Or source the initialization script
   source $FERMI_DIR/fermi-init.sh
   ```

3. **Initialize HEASOFT:**
   ```bash
   heainit
   ```

## Quick Test

Run the test suite to verify the installation:
```bash
python3 test_standalone.py
```

All tests should pass!

## Running Modes

### Interactive Mode (DEFAULT - Just Like Original Perl Script!)

Simply run without any arguments to be prompted for all parameters:

```bash
python3 bex_fermi.py
```

You'll be prompted for:
- Bin size (seconds or days)
- Energy range (Emin, Emax in MeV)
- Whether to use probability photometry
- Catalog file (if using probability mode)
- Output file prefix (optional)
- Whether to calculate diffuse response
- IRF version (1-9)
- Whether to barycenter light curves
- Source parameter file location
- Photon file list location
- Spacecraft file location
- Aperture radius and annulus size
- Probability threshold (if using probability mode)
- Zenith, rock, and bore angle limits
- Spectral index
- Minimum solar distance

**Try the demo:**
```bash
python3 demo_prompting.py
```

This shows how the interactive prompting works without requiring real data files.

### Batch Mode (Non-Interactive)

For automated processing, use batch mode with a parameter file:

```bash
python3 bex_fermi.py -batch -file my_params.par
```

**Create a parameter file:**
```bash
# Option 1: Use the provided example
cp bex.par.example my_params.par
# Edit my_params.par as needed

# Option 2: Generate from current defaults
python3 bex_fermi.py -makepar my_params.par
```

**Parameter file format:**
```
# Comments start with #
bin_size 500
emin 100
emax 500000
use_probability true
catalog gll_psc_v35.fit
# ... etc
```

### Update Mode

Add new data points to existing light curves:

```bash
python3 bex_fermi.py -update
```

Can be combined with `-file` for batch processing:
```bash
python3 bex_fermi.py -update -batch -file my_params.par
```

## Basic Usage Examples

### 1. Apply probability weighting to existing light curve

```bash
python3 pweight_example.py 4FGL_J1826.2-1450 lc_mystar.fits photons_prob.fits
```

### 2. Use as Python module

```python
from bex_fermi import BexFermi

# Create processor
processor = BexFermi()

# Apply probability weighting
processor.pweight(
    source_name='4FGL_J1826.2-1450',
    light_curve_file='lc_mystar.fits',
    photon_file='photons_prob.fits',
    threshold=0.0  # probability threshold
)
```

### 3. Process light curve data

```python
from bex_fermi import BexFermi

processor = BexFermi()

# Convert gtbin output to formatted light curve
processor.prepare_light_curve(
    in_file='lc_fdump.txt',
    out_file='lc_formatted.txt',
    base_name='MySource'
)
```

## Key Functions

### Time Conversions
```python
from bex_fermi import BexFermi

# MJD to Fermi MET
met = BexFermi.mjd2met(58000.0)

# Fermi MET to MJD
mjd = BexFermi.met2mjd(met)
```

### File Processing

**pweight(source_name, light_curve_file, photon_file, threshold=0.0)**
- Applies probability weighting to light curve
- Replaces COUNTS with RCOUNTS (probability-weighted counts)
- Modifies light curve file in place

**prepare_light_curve(in_file, out_file, base_name)**
- Processes gtbin fdump output
- Converts times to MJD
- Calculates count-rate and exposure-based errors
- Creates formatted text output

**stitch_files(old_file, new_file, merge_file, base_name)**
- Merges old and new light curves (for update mode)
- Recalculates errors with consistent mean rate

## Configuration Options

```python
config = {
    'batch': True,              # Batch mode (no prompting)
    'debug': False,             # Enable debug logging
    'update': False,            # Update mode
    'bin_size': 500,            # Time bin size (seconds)
    'roi': 3,                   # Region of interest (degrees)
    'use_probability': True,    # Use probability weighting
    'pthreshold': 0.0,          # Probability threshold
    'emin': 100,                # Minimum energy (MeV)
    'emax': 500000,             # Maximum energy (MeV)
    'spectral_index': -2.1,     # Assumed spectral index
}

processor = BexFermi(config)
```

## Output Format

Light curve files have this format:
```
<source_name>
Time (MJD)
Rate (ph/cm^2/s)
<time_mjd> <rate> <count_rate_error> <bin_halfwidth_days> <exposure_error> <exposure>
```

Example line:
```
58000.1234567 1.23e-07 2.45e-08 0.002894 2.12e-08 1.45e+09
```

Columns:
1. Time (MJD)
2. Count rate (ph/cm²/s)
3. Count-rate based error (ph/cm²/s) - BaBar prescription
4. Bin half-width (days)
5. Exposure-based error (ph/cm²/s) - **Recommended**
6. Exposure (cm²·s)

## Common Issues

**"No module named 'astropy'"**
```bash
pip install astropy
```

**"HEASOFT environment not set up"**
```bash
heainit
```

**"Fermi Science Tools not set up"**
```bash
conda activate fermitools
```

**"Source not found in photon file"**
- Run `gtsrcprob` on photon file first
- Check source name matches exactly (case-sensitive)

## Next Steps

1. Read the full README.md for detailed documentation
2. Check test_standalone.py for algorithm examples
3. See pweight_example.py for a complete usage example
4. Review bex_fermi.py source code for implementation details

## Getting Help

Run with `-h` flag for command-line help:
```bash
python3 bex_fermi.py -h
python3 pweight_example.py -h
```

Enable debug mode for detailed logging:
```bash
python3 bex_fermi.py -debug
```

## References

- **Kerr 2011**, ApJ 732, 1, 38 - Probability weighting method
- **Corbet et al. 2022**, ApJ, doi:10.3847/1538-4357/ac6fe2 - Application to Fermi LAT
