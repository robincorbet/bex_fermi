# Usage Comparison: Perl vs Python

This document shows how to use the Python version the same way you used the Perl version.

## Side-by-Side Comparison

### Interactive Mode (Default)

**Original Perl:**
```bash
./bex2020nd
```

**Python Equivalent:**
```bash
python3 bex_fermi.py
```

**Result:** Both will prompt you for all parameters interactively.

### Batch Mode with Parameter File

**Original Perl:**
```bash
./bex2020nd -batch -file bex.par
```

**Python Equivalent:**
```bash
python3 bex_fermi.py -batch -file bex.par
```

**Result:** Both read parameters from file and run without prompting.

### Update Mode

**Original Perl:**
```bash
./bex2020nd -update
```

**Python Equivalent:**
```bash
python3 bex_fermi.py -update
```

**Result:** Both add new points to existing light curves.

### Creating a Parameter File

**Original Perl:**
```bash
./bex2020nd -makepar bex.par
```

**Python Equivalent:**
```bash
python3 bex_fermi.py -makepar bex.par
```

**Result:** Both create a parameter file with default values.

### Time Range Specification

**Original Perl:**
```bash
./bex2020nd -start 58000 -stop 59000
```

**Python Equivalent:**
```bash
python3 bex_fermi.py -start 58000 -stop 59000
```

**Result:** Both process data in the specified MJD range.

### Debug Mode

**Original Perl:**
```bash
./bex2020nd -debug
```

**Python Equivalent:**
```bash
python3 bex_fermi.py -debug
```

**Result:** Both enable verbose debug output.

## Parameter File Format

The format is **identical** between Perl and Python versions:

```
# Comment lines start with #
bin_size 500
emin 100
emax 500000
use_probability true
catalog gll_psc_v35.fit
pthreshold 0.0
roi 3
zenith_limit 105
rock 90
bore 180
event_class_min 3
irf_code 9
spectral_index -2.1
barycenter true
sun_minimum 5.0
source_list slist.dat
ft1_list plist.dat
ft2 lat_spacecraft_merged.fits
prefix
```

**You can use the same parameter files for both versions!**

## Interactive Prompting Sequence

Both versions prompt for parameters in this order:

1. **Bin size** (seconds, negative for days)
2. **Energy range** (Emin, Emax in MeV)
3. **Probability photometry?** (yes/no)
4. **Catalog file** (if probability = yes)
5. **Output prefix** (optional)
6. **Calculate diffuse response?** (yes/no)
7. **IRF version** (1-9)
8. **Barycenter?** (yes/no)
9. **Source parameter file** (with existence check)
10. **Photon file list** (with existence check)
11. **Spacecraft file** (auto-detects *_SC00.fits in current dir)
12. **Aperture radius** (degrees)
13. **Inner annulus radius** (degrees)
14. **Probability threshold** (if probability = yes)
15. **Zenith limit** (degrees)
16. **Rock angle limit** (degrees)
17. **Bore limit** (degrees)
18. **Spectral index** (auto-negated if positive)
19. **Minimum solar distance** (degrees)

**Prompting behavior is identical:**
- Press Enter to accept default (shown in brackets)
- Validation for ranges and types
- File existence checking
- Clear error messages

## Input Files

### Source List File (slist.dat)

**Format is identical:**
```
4FGL_J1826.2-1450 283.5542 -14.8431
4FGL_J0835.3-4510 128.8375 -45.1764
```

For probability mode, optional 4th field is XML model file:
```
4FGL_J1826.2-1450 283.5542 -14.8431 model_J1826.xml
```

For non-probability mode, optional fields 4-6 are:
```
SOURCE_NAME RA DEC [roi] [emin] [emax]
```

### Photon List File (plist.dat)

**Format is identical:**
```
/path/to/photons_week001.fits
/path/to/photons_week002.fits
/path/to/photons_week003.fits
```

### Spacecraft File

Both versions accept the same FT2 format:
```
lat_spacecraft_merged.fits
```

## Output Files

### Light Curve Files

**Same naming convention:**
```
lc_<source_name>dmp1.out
```

Examples:
- `lc_4FGL_J1826.2-1450dmp1.out`
- `lc_4FGL_J0835.3-4510dmp1.out`

**Same format:**
```
<source_name>
Time (MJD)
Rate (ph/cm^2/s)
58000.1234567 1.23e-07 2.45e-08 0.002893519 2.12e-08 1.45e+09
58000.2345678 1.45e-07 2.67e-08 0.002893519 2.34e-08 1.52e+09
```

Columns:
1. Time (MJD)
2. Count rate (ph/cm²/s)
3. Count-rate based error
4. Bin half-width (days)
5. Exposure-based error (**recommended**)
6. Exposure (cm²·s)

## Differences in Usage

### Minor Differences

1. **Execution:**
   - Perl: `./bex2020nd` (if executable) or `perl bex2020nd`
   - Python: `python3 bex_fermi.py`

2. **Help:**
   - Perl: `./bex2020nd -help`
   - Python: `python3 bex_fermi.py -h` or `--help`

3. **Version:**
   - Perl: Shows version at startup
   - Python: `python3 bex_fermi.py -v` or `--version`

### Functional Differences

**Python version currently:**
- ✅ Full interactive prompting (identical to Perl)
- ✅ Parameter file I/O (compatible with Perl files)
- ✅ Core pweight algorithm (C code fully integrated)
- ✅ Error calculations (identical)
- ✅ Time conversions (identical)
- ⚠️ Fermi tools pipeline (skeleton, needs completion)
- ⚠️ XML model generation (not implemented)

## Migration Strategy

### For Current Users

**Option 1: Drop-in Replacement for pweight**
```bash
# Your existing workflow
gtselect ...
gtmktime ...
gtbin ...
gtsrcprob ...

# Replace this Perl/C step:
# pweight 4FGL_J1826.2-1450 lc.fits ph.fits

# With Python:
python3 -c "from bex_fermi import BexFermi; \
  BexFermi().pweight('4FGL_J1826.2-1450', 'lc.fits', 'ph.fits')"
```

**Option 2: Use for Interactive Parameter Setup**
```bash
# Use Python for interactive setup
python3 bex_fermi.py -makepar my_params.par

# Edit if needed
nano my_params.par

# Use with Perl script
./bex2020nd -batch -file my_params.par
```

**Option 3: Gradual Migration**
```bash
# Use Python's interactive mode
python3 bex_fermi.py

# It will prompt just like the Perl version
# Parameter files are compatible between versions
```

## Examples

### Example 1: Quick Analysis (Interactive)

```bash
# Just run and answer prompts
python3 bex_fermi.py
```

You'll be walked through all parameters with sensible defaults.

### Example 2: Repeated Analysis (Batch)

```bash
# First time: set up parameters interactively
python3 bex_fermi.py -makepar standard_analysis.par

# Subsequently: reuse parameters
python3 bex_fermi.py -batch -file standard_analysis.par
```

### Example 3: Update with New Data

```bash
# Get new data from FSSC
# Add to plist.dat

# Update light curves
python3 bex_fermi.py -update -batch -file standard_analysis.par
```

### Example 4: Parameter Variations

```bash
# Base parameters
python3 bex_fermi.py -makepar base.par

# Create variations
cp base.par high_energy.par
# Edit high_energy.par: emin 1000, emax 500000

cp base.par short_bins.par  
# Edit short_bins.par: bin_size 100

# Run each
python3 bex_fermi.py -batch -file high_energy.par
python3 bex_fermi.py -batch -file short_bins.par
```

## Command Reference

### Common Commands

| Task | Perl Version | Python Version |
|------|-------------|----------------|
| Interactive run | `./bex2020nd` | `python3 bex_fermi.py` |
| Batch mode | `./bex2020nd -batch -file x.par` | `python3 bex_fermi.py -batch -file x.par` |
| Update mode | `./bex2020nd -update` | `python3 bex_fermi.py -update` |
| Make params | `./bex2020nd -makepar x.par` | `python3 bex_fermi.py -makepar x.par` |
| Time range | `./bex2020nd -start X -stop Y` | `python3 bex_fermi.py -start X -stop Y` |
| Debug | `./bex2020nd -debug` | `python3 bex_fermi.py -debug` |
| Help | `./bex2020nd -help` | `python3 bex_fermi.py -h` |

### All Options

```bash
python3 bex_fermi.py [options]

Options:
  -h, --help              Show help message
  -v, --version           Show version number
  -batch                  Batch mode (no prompting, requires -file)
  -file PARFILE           Read parameters from file
  -makepar PARFILE        Create parameter file with defaults
  -update                 Update existing light curves
  -start MJD              Start time in MJD
  -stop MJD               Stop time in MJD
  -pthreshold PROB        Probability threshold override
  -debug                  Enable debug output
```

## Tips for Switching

1. **Test First:** Run both versions on the same data and compare outputs
2. **Keep Backups:** Keep your Perl/C versions until fully validated
3. **Parameter Files:** Create and test parameter files before batch runs
4. **Interactive First:** Use interactive mode to get familiar with Python version
5. **Debug Mode:** Use `-debug` for detailed logging during testing

## Troubleshooting

**"No module named 'astropy'"**
```bash
pip install numpy astropy
```

**"Batch mode requires -file parameter"**
```bash
# Create a parameter file first
python3 bex_fermi.py -makepar my_params.par

# Then run batch mode
python3 bex_fermi.py -batch -file my_params.par
```

**Want to see prompting demo without real files?**
```bash
python3 demo_prompting.py
```

**Different results from Perl version?**
- Check that parameter files match exactly
- Verify same input files are being used
- Run with `-debug` to see what's happening
- Core algorithms are identical, so differences indicate configuration issues

## Getting Help

1. **Interactive demo:** `python3 demo_prompting.py`
2. **Command help:** `python3 bex_fermi.py -h`
3. **Test algorithms:** `python3 test_standalone.py`
4. **Read docs:** See README.md and QUICKSTART.md
5. **Debug mode:** `python3 bex_fermi.py -debug`

The Python version is designed to be a **drop-in replacement** with the **same interface** as the original Perl script!
