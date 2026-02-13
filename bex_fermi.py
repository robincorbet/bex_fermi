#!/usr/bin/env python3
"""
bex_fermi.py - Probability-weighted aperture photometry for Fermi LAT data

Python conversion of bex2020nd Perl script and pweight C code.
Based on methods described in:
- Corbet et al. 2022, ApJ, doi:10.3847/1538-4357/ac6fe2
- Kerr 2011, ApJ 732, 1, 38, doi:10.1088/0004-637X/732/1/38

Author: Converted to Python from Perl/C by Claude
Original Perl author: Robin Corbet (UMBC)
Original C author: Robin Corbet (UMBC)
Version: 2.0
Date: 2026-02-09
"""

import os
import sys
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
import logging

# Constants
VERSION = "2.0"
MJDREF = 51910.0 + 7.428703703703703e-4
FERMI_MET_ORIGIN = 51910.0  # MJD for 2001-01-01 00:00:00 TT

# Default parameters
DEFAULT_BIN_SIZE = 500
DEFAULT_ROI = 1
DEFAULT_ROI_PROBABILITY = 3
DEFAULT_ZENITH_LIMIT = 105
DEFAULT_ROCK = 90
DEFAULT_BORE = 180
DEFAULT_EVENT_CLASS_MIN = 3
DEFAULT_SOURCE_LIST = "slist.dat"
DEFAULT_FT1_LIST = "plist.dat"
DEFAULT_CATALOG = "gll_psc_v35.fit"
DEFAULT_PTHRESHOLD = 0.0
DEFAULT_FT2 = "lat_spacecraft_merged.fits"
DEFAULT_IRF_CODE = 9
DEFAULT_SUN_MINIMUM = 5.0
DEFAULT_EMIN = 100
DEFAULT_EMAX = 500000
DEFAULT_SPECTRAL_INDEX = -2.1
DEFAULT_USE_PROBABILITY = True
DEFAULT_BARYCENTER = True
DEFAULT_DODIFFUSE = False
DEFAULT_PREFIX = ""

# Model maker options
MODEL_MAKER = "make4FGLxml"
MODEL_MAKER_FLAG = ""


class BexFermi:
    """Main class for Fermi LAT aperture photometry with probability weighting."""
    
    def __init__(self, config=None):
        """Initialize with optional configuration dictionary."""
        self.config = config or {}
        self.setup_logging()
        
        # Set up parameters
        self.bin_size = self.config.get('bin_size', DEFAULT_BIN_SIZE)
        self.roi = self.config.get('roi', DEFAULT_ROI)
        self.zenith_limit = self.config.get('zenith_limit', DEFAULT_ZENITH_LIMIT)
        self.rock = self.config.get('rock', DEFAULT_ROCK)
        self.bore = self.config.get('bore', DEFAULT_BORE)
        self.event_class_min = self.config.get('event_class_min', DEFAULT_EVENT_CLASS_MIN)
        self.source_list = self.config.get('source_list', DEFAULT_SOURCE_LIST)
        self.ft1_list = self.config.get('ft1_list', DEFAULT_FT1_LIST)
        self.catalog = self.config.get('catalog', DEFAULT_CATALOG)
        self.pthreshold = self.config.get('pthreshold', DEFAULT_PTHRESHOLD)
        self.ft2 = self.config.get('ft2', DEFAULT_FT2)
        self.irf_code = self.config.get('irf_code', DEFAULT_IRF_CODE)
        self.sun_minimum = self.config.get('sun_minimum', DEFAULT_SUN_MINIMUM)
        self.emin = self.config.get('emin', DEFAULT_EMIN)
        self.emax = self.config.get('emax', DEFAULT_EMAX)
        self.spectral_index = self.config.get('spectral_index', DEFAULT_SPECTRAL_INDEX)
        self.use_probability = self.config.get('use_probability', DEFAULT_USE_PROBABILITY)
        self.barycenter = self.config.get('barycenter', DEFAULT_BARYCENTER)
        self.batch_mode = self.config.get('batch', False)
        self.update_mode = self.config.get('update', False)
        self.debug_mode = self.config.get('debug', False)
        self.start_time = self.config.get('start_time', 0)
        self.stop_time = self.config.get('stop_time', 9e99)
        self.prefix = self.config.get('prefix', DEFAULT_PREFIX)
        self.dodiffuse = self.config.get('dodiffuse', DEFAULT_DODIFFUSE)
        self.roi_inner = 0  # Will be set during prompting if needed
        
    def setup_logging(self):
        """Set up logging configuration."""
        level = logging.DEBUG if self.config.get('debug', False) else logging.INFO
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
    @staticmethod
    def mjd2met(mjd):
        """Convert MJD to Fermi MET (Mission Elapsed Time)."""
        return (mjd - FERMI_MET_ORIGIN) * 86400.0
    
    @staticmethod
    def met2mjd(met):
        """Convert Fermi MET to MJD."""
        return (met / 86400.0) + FERMI_MET_ORIGIN
    
    def pweight(self, source_name, light_curve_file, photon_file, threshold=0.0):
        """
        Python implementation of pweight.c functionality.
        
        Takes a LAT light curve file and the photon file used to create it.
        Removes COUNTS column from light curve and replaces it with a floating point
        RCOUNTS column which contains a sum of the probabilities that each photon came
        from a certain source.
        
        Args:
            source_name: Name of the source (column name in photon file)
            light_curve_file: Path to the light curve FITS file
            photon_file: Path to the photon (events) FITS file
            threshold: Probability threshold (default 0.0)
        """
        self.logger.info(f"Processing probability weights for {source_name}")
        self.logger.info(f"Light curve file: {light_curve_file}")
        self.logger.info(f"Photon file: {photon_file}")
        self.logger.info(f"Probability threshold: {threshold}")
        
        # Open FITS files
        with fits.open(photon_file) as ph_hdul, \
             fits.open(light_curve_file, mode='update') as lc_hdul:
            
            # Get photon file data
            ph_data = ph_hdul['EVENTS'].data
            nphotons = len(ph_data)
            self.logger.info(f"Number of photons: {nphotons}")
            
            # Get probability column for this source
            try:
                probability = ph_data[source_name]
            except KeyError:
                self.logger.error(f"Source {source_name} not found in photon file")
                raise
            
            # Get photon times
            ph_time = ph_data['TIME']
            
            # Get light curve data
            lc_data = lc_hdul['RATE'].data
            lc_header = lc_hdul['RATE'].header
            nlc_bins = len(lc_data)
            
            # Get light curve parameters
            tstart = lc_header['TSTART']
            tstop = lc_header['TSTOP']
            binsize = (tstop - tstart) / (nlc_bins - 1.0)
            hbs = binsize / 2.0  # half bin size
            
            # Get exposure and time from light curve
            exposure = lc_data['EXPOSURE']
            lc_time = lc_data['TIME']
            
            # Get GTIs from light curve file
            gti_data = lc_hdul['GTI'].data
            gti_start = gti_data['START']
            gti_stop = gti_data['STOP']
            ngtis = len(gti_data)
            
            # Initialize rate array
            rate = np.zeros(nlc_bins, dtype=np.float32)
            total_sum = 0.0
            
            # Cycle through photons calculating new light curve
            for i in range(nphotons):
                # Get time bin for this photon
                timebin = self._get_timebin(ph_time[i], lc_time, nlc_bins, 
                                           tstart, hbs, binsize)
                
                # Check if photon is in valid time bin and GTI
                if (0 <= timebin < nlc_bins and 
                    self._check_gti(gti_start, gti_stop, ph_time[i], ngtis)):
                    
                    if probability[i] > threshold:
                        rate[timebin] += probability[i]
                        total_sum += probability[i]
            
            self.logger.info(f"Probability sum: {total_sum}")
            
            # Delete existing COUNTS column and add RCOUNTS
            # Find the COUNTS column number
            cols = lc_hdul['RATE'].columns
            count_idx = None
            for idx, col in enumerate(cols):
                if col.name == 'COUNTS':
                    count_idx = idx
                    break
            
            if count_idx is not None:
                # Create new column definition
                new_cols = []
                for idx, col in enumerate(cols):
                    if idx == count_idx:
                        # Replace COUNTS with RCOUNTS
                        new_col = fits.Column(name='RCOUNTS', format='E', 
                                            array=rate)
                        new_cols.append(new_col)
                    else:
                        new_cols.append(col)
                
                # Create new table HDU with updated columns
                new_hdu = fits.BinTableHDU.from_columns(new_cols, 
                                                        header=lc_header,
                                                        name='RATE')
                
                # Replace the RATE extension
                lc_hdul['RATE'] = new_hdu
                
                # Write changes
                lc_hdul.flush()
                
            self.logger.info("Successfully updated light curve with probability weights")
    
    @staticmethod
    def _check_gti(gti_start, gti_stop, time, ngtis):
        """Check whether a specified time is within a good time interval."""
        for i in range(ngtis):
            if gti_start[i] <= time <= gti_stop[i]:
                return True
            # Assume GTIs are time ordered
            if time < gti_start[i]:
                return False
        return False
    
    @staticmethod
    def _get_timebin(ph_time, lc_time, nlc_bins, tstart, hbs, binsize):
        """
        Find the bin number of the light curve bin that corresponds
        to the photon time.
        """
        # Approximate estimate of bin number
        j = int(np.floor((ph_time - tstart) / binsize)) - 10
        j = max(j, 0)
        
        # Search for the correct bin
        for i in range(j, nlc_bins):
            if (lc_time[i] - hbs) <= ph_time <= (lc_time[i] + hbs):
                return i
        
        # Didn't find it
        logging.error(f"Error in get_timebin. ph_time = {ph_time:.10f}")
        return -999
    
    def prepare_light_curve(self, in_file, out_file, base_name):
        """
        Process the fdumped light curve that comes out of gtbin.
        Changes times to UT, counts to rates, and calculates two types of error bars.
        
        Args:
            in_file: Input file from gtbin fdump
            out_file: Output file with processed light curve
            base_name: Base name for the source
        """
        self.logger.info(f"Preparing light curve: {in_file} -> {out_file}")
        
        # Read the input file
        data = []
        with open(in_file, 'r') as f:
            # Skip header (4 lines)
            for _ in range(4):
                f.readline()
            
            # Read data lines
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 5:
                    data.append({
                        'i': int(fields[0]),
                        'time': float(fields[1]),
                        'counts': float(fields[2]),
                        'exposure': float(fields[3]),
                        'timedel': float(fields[4])
                    })
        
        # Calculate mean count rate
        countsum = 0.0
        expsum = 0.0
        for d in data:
            if d['exposure'] > 0.0:
                countsum += d['counts']
                expsum += d['exposure']
        
        meanrate = countsum / expsum if expsum > 0 else 0.0
        self.logger.info(f"Mean count rate: {meanrate}")
        
        # Write output file
        with open(out_file, 'w') as f:
            # Write header
            f.write(f"{base_name}\n")
            f.write("Time (MJD)\n")
            f.write("Rate (ph/cm^2/s)\n")
            
            # Process each data point
            for d in data:
                if d['exposure'] > 0.0:
                    # Calculate errors from population counts
                    popcounts = meanrate * d['exposure']
                    perr = np.sqrt(popcounts)
                    
                    # Calculate errors from sample rate (BaBar approach)
                    cerr1 = 0.5 + np.sqrt(d['counts'] + 0.25)
                    cerr2 = -0.5 + np.sqrt(d['counts'] + 0.25)
                    
                    # RMS type error
                    rmserr = np.sqrt((cerr1**2 + cerr2**2) / 2.0)
                    
                    # Calculate rates and errors
                    rate = d['counts'] / d['exposure']
                    rerr = rmserr / d['exposure']
                    prerr = perr / d['exposure']
                    
                    # Convert time to MJD
                    time_mjd = (d['time'] / 86400.0) + MJDREF
                    timedel_days = d['timedel'] / (2.0 * 86400.0)
                    
                    # Write output line
                    f.write(f"{time_mjd:.10f} {rate:.10e} {rerr:.10e} "
                           f"{timedel_days:.10f} {prerr:.10e} {d['exposure']:.10e}\n")
    
    def redo_errors(self, filename, expsum, countsum, base_name):
        """
        Recalculate exposure-based errors using mean rate from entire file.
        Used when merging files in update mode.
        """
        meanrate = countsum / expsum if expsum > 0 else 0.0
        self.logger.info(f"Recalculating errors with mean rate: {meanrate}")
        
        # Read original file
        lines = []
        with open(filename, 'r') as f:
            # Keep header (3 lines)
            for _ in range(3):
                lines.append(f.readline())
            
            # Process data lines
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 6:
                    time_mjd = float(fields[0])
                    rate = float(fields[1])
                    rerr = float(fields[2])
                    timedel = float(fields[3])
                    exposure = float(fields[5])
                    
                    # Recalculate population-based error
                    if exposure > 0:
                        popcounts = meanrate * exposure
                        perr = np.sqrt(popcounts)
                        prerr = perr / exposure
                        
                        lines.append(f"{time_mjd:.10f} {rate:.10e} {rerr:.10e} "
                                   f"{timedel:.10f} {prerr:.10e} {exposure:.10e}\n")
        
        # Write back to file
        with open(filename, 'w') as f:
            f.writelines(lines)
    
    def stitch_files(self, old_file, new_file, merge_file, base_name):
        """
        Join together old and new light curve files for update mode.
        """
        self.logger.info(f"Stitching files: {old_file} + {new_file} -> {merge_file}")
        
        # Read old file
        old_data = []
        old_header = []
        expsum = 0.0
        countsum = 0.0
        
        with open(old_file, 'r') as f:
            # Read header (3 lines)
            for _ in range(3):
                old_header.append(f.readline())
            
            # Read data
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 6:
                    old_data.append({
                        'time': float(fields[0]),
                        'rate': float(fields[1]),
                        'rerr': float(fields[2]),
                        'timedel': float(fields[3]),
                        'prerr': float(fields[4]),
                        'exposure': float(fields[5]),
                        'line': line
                    })
                    expsum += old_data[-1]['exposure']
                    countsum += old_data[-1]['rate'] * old_data[-1]['exposure']
        
        # Remove last bin (may be incomplete)
        if old_data:
            expsum -= old_data[-1]['exposure']
            countsum -= old_data[-1]['rate'] * old_data[-1]['exposure']
            last_old_time = old_data[-2]['time'] if len(old_data) > 1 else 0
            old_data = old_data[:-1]
        else:
            last_old_time = 0
        
        # Read new file
        new_data = []
        with open(new_file, 'r') as f:
            # Skip header (3 lines)
            for _ in range(3):
                f.readline()
            
            # Read data after last old time
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 6:
                    time = float(fields[0])
                    if time >= last_old_time:
                        new_data.append({
                            'time': time,
                            'rate': float(fields[1]),
                            'exposure': float(fields[5]),
                            'line': line
                        })
                        expsum += new_data[-1]['exposure']
                        countsum += new_data[-1]['rate'] * new_data[-1]['exposure']
        
        # Write merged file
        with open(merge_file, 'w') as f:
            # Write header
            f.writelines(old_header)
            
            # Write old data
            for d in old_data:
                f.write(d['line'])
            
            # Write new data
            for d in new_data:
                f.write(d['line'])
        
        # Redo errors with consistent mean rate
        self.redo_errors(merge_file, expsum, countsum, base_name)
    
    def run_fermi_tool(self, tool_name, **kwargs):
        """
        Run a Fermi Science Tool with given parameters.
        
        Args:
            tool_name: Name of the Fermi tool to run
            **kwargs: Tool parameters as keyword arguments
        """
        cmd = [tool_name]
        for key, value in kwargs.items():
            cmd.append(f"{key}={value}")
        
        self.logger.info(f"Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if result.stdout and not self.batch_mode:
                self.logger.debug(result.stdout)
            return result
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error running {tool_name}: {e}")
            if e.stderr:
                self.logger.error(e.stderr)
            raise
    
    def check_environment(self):
        """Check that required environment variables and tools are available."""
        # Check for HEADAS/HEASOFT
        if 'HEADAS' not in os.environ or 'PFILES' not in os.environ:
            self.logger.error("HEASOFT environment not set up. Run 'heainit' or equivalent.")
            sys.exit(1)
        
        # Check for Fermi Science Tools
        if 'FERMI_DIR' not in os.environ:
            self.logger.error("Fermi Science Tools not set up. Initialize fermitools.")
            sys.exit(1)
        
        self.logger.info("Environment check passed")
        
        # Create unique PFILES directory to avoid collisions
        # This is critical for running multiple instances simultaneously
        self._setup_pfiles()
    
    def _setup_pfiles(self):
        """
        Create a unique PFILES directory for this process.
        This prevents parameter file collisions when running multiple instances.
        Based on the Perl script's approach.
        """
        # Create unique pfiles directory in /tmp based on PID
        pid = os.getpid()
        self.pfiles_tmp = f"/tmp/pfiles_bex_{pid}.tmp"
        
        try:
            os.makedirs(self.pfiles_tmp, exist_ok=True)
            self.logger.info(f"Created PFILES directory: {self.pfiles_tmp}")
            
            # Get original PFILES path
            original_pfiles = os.environ.get('PFILES', '')
            
            # Copy parameter files from FERMI_DIR and HEADAS
            fermi_dir = os.environ.get('FERMI_DIR', '')
            headas = os.environ.get('HEADAS', '')
            
            if fermi_dir:
                syspfiles = os.path.join(fermi_dir, 'syspfiles')
                if os.path.exists(syspfiles):
                    subprocess.run(
                        f'cp {syspfiles}/*.par {self.pfiles_tmp}/.',
                        shell=True,
                        check=False  # Don't fail if no .par files
                    )
            
            if headas:
                syspfiles = os.path.join(headas, 'syspfiles')
                if os.path.exists(syspfiles):
                    subprocess.run(
                        f'cp {syspfiles}/f*.par {self.pfiles_tmp}/.',
                        shell=True,
                        check=False
                    )
            
            # Set PFILES to use our directory first, then system defaults
            new_pfiles = f"{self.pfiles_tmp}"
            if original_pfiles:
                # Add original PFILES as fallback (separated by :)
                new_pfiles += f":{original_pfiles}"
            
            os.environ['PFILES'] = new_pfiles
            self.logger.info(f"Set PFILES={new_pfiles}")
            
        except Exception as e:
            self.logger.warning(f"Failed to set up PFILES directory: {e}")
            # Continue anyway - it might still work
    
    def _cleanup_pfiles(self):
        """Clean up the temporary PFILES directory."""
        if hasattr(self, 'pfiles_tmp') and os.path.exists(self.pfiles_tmp):
            try:
                shutil.rmtree(self.pfiles_tmp)
                self.logger.info(f"Cleaned up PFILES directory: {self.pfiles_tmp}")
            except Exception as e:
                self.logger.warning(f"Failed to clean up PFILES directory: {e}")
    
    def get_value(self, prompt, default, minimum, maximum):
        """
        Prompt for a numerical value within specified range.
        In batch mode, returns the default value.
        
        Args:
            prompt: Prompt string
            default: Default value
            minimum: Minimum allowed value (or 'INDEF')
            maximum: Maximum allowed value (or 'INDEF')
        
        Returns:
            The validated numerical value
        """
        if self.batch_mode:
            return default
        
        while True:
            if default == 'INDEF':
                prompt_str = f"{prompt}: "
            else:
                prompt_str = f"{prompt} [default = {default}]: "
            
            response = input(prompt_str).strip()
            
            # Use default if empty
            if response == '' and default != 'INDEF':
                print(f'(Using default: "{default}")')
                return default
            
            # Try to convert to number
            try:
                value = float(response) if '.' in response or 'e' in response.lower() else int(response)
            except ValueError:
                print(f"Value ({response}) is not a number!")
                continue
            
            # Check range
            if minimum != 'INDEF' and value < minimum:
                print(f"Value ({value}) is not in the range {minimum} to {maximum}")
                continue
            if maximum != 'INDEF' and value > maximum:
                print(f"Value ({value}) is not in the range {minimum} to {maximum}")
                continue
            
            return value
    
    def get_yn(self, prompt, default):
        """
        Prompt for a yes/no response.
        In batch mode, returns the default value.
        
        Args:
            prompt: Prompt string
            default: Default value (True/False)
        
        Returns:
            Boolean value
        """
        if self.batch_mode:
            return default
        
        default_str = "y" if default else "n"
        
        while True:
            prompt_str = f"{prompt} [y/n; default = {default_str}]: "
            response = input(prompt_str).strip().upper()
            
            if response == '':
                print(f'(Using default: "{default_str}")')
                return default
            
            if response.startswith('Y'):
                return True
            elif response.startswith('N'):
                return False
            else:
                print(f'Response ({response}) is not "y" or "n"')
    
    def get_string(self, prompt, default=''):
        """
        Prompt for a string value.
        In batch mode, returns the default value.
        
        Args:
            prompt: Prompt string
            default: Default value (can be empty string for optional)
        
        Returns:
            The string value (can be empty if default is empty)
        """
        if self.batch_mode:
            return default
        
        while True:
            if default:
                prompt_str = f"{prompt} [default = {default}]: "
            else:
                # Empty default means optional - just press Enter for empty
                prompt_str = f"{prompt} [optional]: "
            
            response = input(prompt_str).strip()
            
            # If response is empty
            if response == '':
                if default or default == '':
                    # Accept empty if default exists (even if default is '')
                    if default:
                        print(f'(Using default: "{default}")')
                    return default
                else:
                    # This case shouldn't happen but handle it
                    print("A value is required")
            else:
                return response
    
    def get_file(self, prompt, default=''):
        """
        Prompt for a file name and verify it exists.
        In batch mode, returns the default value without checking.
        
        Args:
            prompt: Prompt string
            default: Default filename
        
        Returns:
            The validated filename
        """
        if self.batch_mode:
            if not Path(default).exists():
                raise FileNotFoundError(f"File not found: {default}")
            return default
        
        while True:
            if default:
                prompt_str = f"{prompt} [default = {default}]: "
            else:
                prompt_str = f"{prompt}: "
            
            response = input(prompt_str).strip()
            
            if response == '' and default:
                print(f'(Using default: "{default}")')
                filename = default
            else:
                filename = response
            
            if Path(filename).exists():
                return filename
            else:
                print(f"Can't open file: {filename}")
    
    def prompt_for_parameters(self):
        """
        Interactively prompt for all analysis parameters.
        This mimics the original Perl script's behavior.
        Uses current instance values as defaults (which may have been loaded from a parameter file).
        """
        print("\n" + "="*60)
        print("BEX FERMI - Parameter Setup")
        print("="*60 + "\n")
        
        # Bin size
        self.bin_size = self.get_value(
            "Give bin size (seconds, or -ve for days!)",
            self.bin_size, -1e10, 1e10
        )
        if self.bin_size < 0:
            self.bin_size = -self.bin_size * 86400.0
        print(f"bin size = {self.bin_size} seconds\n")
        
        # Energy range
        self.emin = self.get_value(
            "Give Emin (MeV)",
            self.emin, 0, 1e10
        )
        
        self.emax = self.get_value(
            "Give Emax (MeV)",
            self.emax, self.emin, 1e10
        )
        
        # Probability photometry
        self.use_probability = self.get_yn(
            "Use probability photometry?",
            self.use_probability
        )
        
        if self.use_probability:
            print("\n(Probability photometry will be done.")
            print("This will only work for cataloged sources.")
            print("Automatic XML file generation will be done unless overridden by source file contents")
            print("(Increasing default ROI))\n")
            self.roi = DEFAULT_ROI_PROBABILITY
            
            self.catalog = self.get_string(
                "Give source catalog file name",
                self.catalog
            )
            print(f"catalog = {self.catalog}\n")
        else:
            print("\n(Regular aperture photometry will be done)\n")
            self.roi = DEFAULT_ROI
        
        # Prefix
        self.prefix = self.get_string(
            "Give any prefix for output file name(s) [optional]",
            self.prefix
        )
        
        # Diffuse response
        self.dodiffuse = self.get_yn(
            "Calculate diffuse response?",
            self.dodiffuse
        )
        
        if self.dodiffuse:
            print("(gtdiffrsp will be called)\n")
        else:
            print("(gtdiffrsp won't be called)\n")
        
        # IRF selection
        print("\nIRF Options:")
        print("  1 = P6_V3")
        print("  2 = P6_V11")
        print("  3 = P7_V6")
        print("  4 = P7REP")
        print("  5 = P7REP/newbackground")
        print("  6 = P8")
        print("  7 = P8R3 new iso")
        print("  8 = 4FGL")
        print("  9 = v2.0")
        
        self.irf_code = self.get_value(
            "Which IRFs?",
            self.irf_code, 1, 9
        )
        
        # Barycenter
        self.barycenter = self.get_yn(
            "Barycenter light curves?",
            self.barycenter
        )
        
        if self.barycenter:
            print("(Barycenter correction. Light curves may be truncated by 60s if needed.)\n")
        else:
            print("(No barycenter correction)\n")
        
        # Source list file
        self.source_list = self.get_file(
            "Give source parameter file",
            self.source_list
        )
        
        # Photon file list
        self.ft1_list = self.get_file(
            "Give photon file name list",
            self.ft1_list
        )
        
        # Spacecraft file
        # Check for SC files in current directory
        sc_files = list(Path('.').glob('*_SC00.fits'))
        if len(sc_files) == 1:
            default_ft2 = str(sc_files[0])
            print(f"Number of SC files found = 1")
        else:
            default_ft2 = self.ft2
            if sc_files:
                print(f"Number of SC files found = {len(sc_files)}")
        
        self.ft2 = self.get_file(
            "Give spacecraft (FT2) file name",
            default_ft2
        )
        
        # ROI parameters
        self.roi = self.get_value(
            "Give aperture radius (degrees)",
            self.roi, 0, 180
        )
        
        self.roi_inner = self.get_value(
            "Give inner size of aperture annulus (degrees)",
            0, 0, self.roi
        )
        
        # Probability threshold
        if self.use_probability:
            self.pthreshold = self.get_value(
                "Give probability minimum to use",
                self.pthreshold, 0, 1
            )
        
        # Zenith limit
        self.zenith_limit = self.get_value(
            "Give Zenith limit (degrees)",
            self.zenith_limit, 0, 'INDEF'
        )
        
        # Rock angle
        self.rock = self.get_value(
            "Give rock angle limit (degrees)",
            self.rock, -90, 90
        )
        
        # Bore limit
        self.bore = self.get_value(
            "Give Bore limit (degrees)",
            self.bore, 0, 360
        )
        
        # Spectral index
        self.spectral_index = self.get_value(
            "Give spectral index",
            self.spectral_index, -100, 100
        )
        
        if self.spectral_index > 0:
            print(f"Forcing spectral index to be negative: ", end='')
            self.spectral_index = -self.spectral_index
            print(f"{self.spectral_index}")
        
        # Sun minimum
        self.sun_minimum = self.get_value(
            "Give minimum solar distance (degrees)",
            self.sun_minimum, 0, 180
        )
        
        print("\n" + "="*60)
        print("Parameter setup complete!")
        print("="*60 + "\n")
    
    def process_sources(self):
        """Main processing loop for all sources in the source list."""
        self.logger.info(f"Processing sources from: {self.source_list}")
        
        try:
            # Read source list
            sources = []
            with open(self.source_list, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        sources.append(line.split())
            
            self.logger.info(f"Found {len(sources)} sources to process")
            
            # Process each source
            for source_info in sources:
                source_name = source_info[0]
                self.logger.info(f"\n{'='*60}")
                self.logger.info(f"Processing source: {source_name}")
                self.logger.info(f"{'='*60}")
                
                try:
                    self.process_single_source(source_name, source_info)
                except Exception as e:
                    self.logger.error(f"Error processing {source_name}: {e}")
                    if self.debug_mode:
                        raise
        
        finally:
            # Clean up PFILES directory
            self._cleanup_pfiles()
    
    def process_single_source(self, source_name, source_info):
        """
        Process a single source through the full Fermi LAT pipeline.
        This follows the same workflow as the original Perl script.
        """
        print(f"\n{'='*60}")
        print(f"Processing: {source_name}")
        print(f"{'='*60}\n")
        
        # Parse source info
        ra = float(source_info[1])
        dec = float(source_info[2])
        
        # Check for model file (for probability mode)
        model_file = None
        model_given = False
        if self.use_probability and len(source_info) > 3:
            model_file = source_info[3]
            if Path(model_file).exists():
                model_given = True
                self.logger.info(f"Using model file: {model_file}")
            else:
                self.logger.warning(f"Model file {model_file} not found, will generate")
        
        # Create base name (replace + with p to avoid FTOOLS issues)
        base = source_name.replace('+', 'p')
        if self.prefix:
            base = self.prefix + base
        
        # Create temporary directory for this source
        temp_dir = tempfile.mkdtemp(prefix=f'bex_{base}_')
        self.logger.info(f"Working directory: {temp_dir}")
        
        try:
            # Define file names
            ft1_file_0 = os.path.join(temp_dir, f"{base}.eventfile0.0.fits")
            ft1_file = os.path.join(temp_dir, f"{base}.eventfile0.fits")
            temp2 = os.path.join(temp_dir, f"{base}.temp2.fits")
            temp3 = os.path.join(temp_dir, f"{base}.temp3.fits")
            eventfile2 = os.path.join(temp_dir, f"{base}.eventfile2.fits")
            lc_file = f"lc_{base}.fits"
            lc_dump = f"lc_{base}.dmp1"
            final_lc = f"lc_{base}dmp1.out"
            
            # Determine IRF settings based on irf_code
            irfs, event_class, event_type = self._get_irf_settings()
            
            # Step 1: gtselect - initial selection on full ROI, all times
            self.logger.info("Step 1: Running gtselect (initial selection)")
            self.run_fermi_tool(
                'gtselect',
                chatter=2 if not self.batch_mode else 0,
                infile=f'@{self.ft1_list}',
                outfile=ft1_file_0,
                ra=ra,
                dec=dec,
                rad=self.roi,
                tmin=0,
                tmax=0,
                emin=self.emin,
                emax=self.emax,
                zmax=180,
                evclass=event_class,
                evtype=event_type
            )
            
            # Apply annulus selection if needed
            if self.roi_inner > 0:
                self.logger.info(f"Applying annulus selection (inner radius: {self.roi_inner})")
                expr = f"circle({ra},{dec},{self.roi},RA,DEC) && !circle({ra},{dec},{self.roi_inner},RA,DEC)"
                self.run_fermi_tool('fselect', infile=ft1_file_0, outfile=ft1_file,
                                  expression=expr)
            else:
                ft1_file = ft1_file_0
            
            # Get TSTART and TSTOP from photon file
            ph_tstart = self._get_fits_keyword(ft1_file, 'TSTART')
            ph_tstop = self._get_fits_keyword(ft1_file, 'TSTOP')
            self.logger.info(f"Photon file time range: {ph_tstart} - {ph_tstop}")
            
            # Get TSTART and TSTOP from spacecraft file
            sc_tstart = self._get_fits_keyword(self.ft2, 'TSTART')
            sc_tstop = self._get_fits_keyword(self.ft2, 'TSTOP')
            self.logger.info(f"Spacecraft file time range: {sc_tstart} - {sc_tstop}")
            
            # Adjust times for barycenter correction if needed
            tstart = ph_tstart
            tstop = ph_tstop
            
            if self.barycenter:
                if tstart <= sc_tstart:
                    old_tstart = tstart
                    tstart = sc_tstart + 60.0
                    self.logger.info(f"Adjusted tstart for barycenter: {old_tstart} -> {tstart}")
                if tstop >= sc_tstop:
                    old_tstop = tstop
                    tstop = sc_tstop - 60.0
                    self.logger.info(f"Adjusted tstop for barycenter: {old_tstop} -> {tstop}")
            
            # Apply command-line time constraints
            tstart = max(self.start_time, tstart)
            tstop = min(self.stop_time, tstop)
            
            # Check for update mode
            if self.update_mode and Path(final_lc).exists():
                self.logger.info("Update mode: getting start time from existing file")
                tstart = self._get_update_start_time(final_lc)
                old_lc_backup = os.path.join(temp_dir, f"{base}.old_lc")
                shutil.copy(final_lc, old_lc_backup)
            
            # Step 2: gtselect - apply time and zenith cuts
            self.logger.info("Step 2: Running gtselect (time and zenith cuts)")
            self.run_fermi_tool(
                'gtselect',
                chatter=2 if not self.batch_mode else 0,
                infile=ft1_file,
                outfile=temp2,
                ra=ra,
                dec=dec,
                rad=self.roi,
                tmin=tstart,
                tmax=tstop,
                emin=self.emin,
                emax=self.emax,
                zmax=self.zenith_limit,
                evclass=event_class,
                evtype=event_type
            )
            
            # Step 3: gtmktime - apply GTI filters
            self.logger.info("Step 3: Running gtmktime")
            filter_expr = (
                f"(DATA_QUAL>0) && ABS(ROCK_ANGLE)<{self.rock} && (LAT_CONFIG==1) && "
                f"(angsep(RA_ZENITH,DEC_ZENITH,{ra},{dec})+{self.roi}<{self.zenith_limit}) && "
                f"(angsep({ra},{dec},RA_SUN,DEC_SUN)>{self.sun_minimum}+{self.roi}) && "
                f"(angsep({ra},{dec},RA_SCZ,DEC_SCZ)<{self.bore})"
            )
            
            self.run_fermi_tool(
                'gtmktime',
                chatter=2 if not self.batch_mode else 0,
                scfile=self.ft2,
                evfile=temp2,
                outfile=temp3,
                roicut='n',
                filter=filter_expr
            )
            
            # Steps for probability photometry
            if self.use_probability:
                # Generate model file if not provided
                if not model_given:
                    model_file = os.path.join(temp_dir, f"{base}.bex_LATxmlmodel.xml")
                    self.logger.info("Generating XML model file")
                    self._generate_model_file(temp3, model_file, ra, dec)
                
                # Step 4: gtdiffrsp (optional)
                if self.dodiffuse:
                    self.logger.info("Step 4: Running gtdiffrsp")
                    self.run_fermi_tool(
                        'gtdiffrsp',
                        chatter=2 if not self.batch_mode else 0,
                        evfile=temp3,
                        scfile=self.ft2,
                        srcmdl=model_file,
                        irfs=irfs
                    )
                
                # Step 5: gtsrcprob
                self.logger.info("Step 5: Running gtsrcprob")
                self.run_fermi_tool(
                    'gtsrcprob',
                    chatter=2 if not self.batch_mode else 0,
                    evfile=temp3,
                    outfile=eventfile2,
                    scfile=self.ft2,
                    srcmdl=model_file,
                    irfs=irfs
                )
            
            # Step 6: gtbin - create light curve
            self.logger.info("Step 6: Running gtbin")
            self.run_fermi_tool(
                'gtbin',
                chatter=2 if not self.batch_mode else 0,
                algorithm='LC',
                evfile=temp3,
                outfile=lc_file,
                scfile=self.ft2,
                tbinalg='LIN',
                tstart=tstart,
                tstop=tstop,
                dtime=self.bin_size
            )
            
            # Step 7: gtexposure
            self.logger.info("Step 7: Running gtexposure")
            
            # Handle source name prefix for catalog names starting with numbers
            src2 = source_name
            for prefix in ['1FGL', '2FGL', '3FGL', '4FGL', '18M', '24M', 'FL8Y']:
                src2 = src2.replace(prefix, f'_{prefix}')
            
            if self.use_probability:
                self.run_fermi_tool(
                    'gtexposure',
                    chatter=2 if not self.batch_mode else 0,
                    infile=lc_file,
                    scfile=self.ft2,
                    irfs=irfs,
                    srcmdl=model_file,
                    target=src2
                )
            else:
                self.run_fermi_tool(
                    'gtexposure',
                    chatter=2 if not self.batch_mode else 0,
                    infile=lc_file,
                    scfile=self.ft2,
                    irfs=irfs,
                    srcmdl='none',
                    specin=self.spectral_index
                )
            
            # Step 8: pweight - apply probability weighting
            if self.use_probability:
                self.logger.info("Step 8: Applying probability weighting")
                self.pweight(src2, lc_file, eventfile2, self.pthreshold)
            
            # Step 9: Barycenter correction
            if self.barycenter:
                self.logger.info("Step 9: Applying barycenter correction")
                self._apply_barycenter(lc_file, self.ft2, ra, dec)
            
            # Step 10: fdump - extract to text
            self.logger.info("Step 10: Running fdump")
            counts_col = 'RCOUNTS' if self.use_probability else 'COUNTS'
            self.run_fermi_tool(
                'fdump',
                prhead='no',
                infile=f'{lc_file}[1]',
                outfile=lc_dump,
                columns=f'TIME {counts_col} EXPOSURE TIMEDEL',
                pagewidth=256,
                rows='-'
            )
            
            # Step 11: prepare_light_curve - format output
            self.logger.info("Step 11: Formatting light curve")
            self.prepare_light_curve(lc_dump, final_lc, base)
            
            # Step 12: Handle update mode
            if self.update_mode and Path(old_lc_backup).exists():
                self.logger.info("Step 12: Merging with old light curve")
                merge_file = os.path.join(temp_dir, f"{base}.merged_lc")
                self.stitch_files(old_lc_backup, final_lc, merge_file, base)
                shutil.move(merge_file, final_lc)
            
            print(f"\n{'='*60}")
            print(f"SUCCESS! Light curve created: {final_lc}")
            print(f"{'='*60}\n")
            
        finally:
            # Clean up temporary directory
            if not self.debug_mode:
                shutil.rmtree(temp_dir)
                self.logger.info(f"Cleaned up temporary directory")
            else:
                self.logger.info(f"Debug mode: keeping temporary directory {temp_dir}")
    
    def _get_irf_settings(self):
        """Get IRF settings based on irf_code."""
        if self.irf_code >= 6:  # P8 and later
            return 'CALDB', 128, 3
        elif self.irf_code >= 3:  # P7
            return 'CALDB', 2, 0
        else:  # P6
            return f'P6_V{self.irf_code}_DIFFUSE', 3, 0
    
    def _get_diffuse_files(self):
        """Get diffuse model file paths and names based on irf_code."""
        fermi_dir = os.environ.get('FERMI_DIR', '')
        
        # Configuration for different IRF versions
        configs = {
            1: {  # P6_V3
                'gll_file': 'galdiffuse/gll_iem_v02.fit',
                'isotropic_file': 'galdiffuse/isotropic_iem_v02.txt',
                'galactic_name': 'gal_v02',
                'extra_galactic_name': 'eg_v02',
                'model_maker_flag': ''
            },
            2: {  # P6_V11
                'gll_file': 'galdiffuse/gll_iem_v02_P6_V11_DIFFUSE.fit',
                'isotropic_file': 'galdiffuse/isotropic_iem_v02_P6_V11_DIFFUSE.txt',
                'galactic_name': 'gal_v02',
                'extra_galactic_name': 'eg_v02',
                'model_maker_flag': ''
            },
            3: {  # P7_V6
                'gll_file': 'galdiffuse/gal_2yearp7v6_v0.fits',
                'isotropic_file': 'galdiffuse/iso_p7v6source.txt',
                'galactic_name': 'gal_2yearp7v6_v0',
                'extra_galactic_name': 'iso_p7v6source',
                'model_maker_flag': ''
            },
            4: {  # P7REP
                'gll_file': 'galdiffuse/gll_iem_v05.fits',
                'isotropic_file': 'galdiffuse/iso_source_v05.txt',
                'galactic_name': 'gll_iem_v05',
                'extra_galactic_name': 'iso_source_v05',
                'model_maker_flag': ''
            },
            5: {  # P7REP/newbackground
                'gll_file': 'galdiffuse/gll_iem_v05_rev1.fit',
                'isotropic_file': 'galdiffuse/iso_source_v05_rev1.txt',
                'galactic_name': 'gll_iem_v05_rev1',
                'extra_galactic_name': 'iso_source_v05',
                'model_maker_flag': ''
            },
            6: {  # P8
                'gll_file': 'galdiffuse/gll_iem_v06.fits',
                'isotropic_file': 'galdiffuse/iso_P8R2_SOURCE_V6_v06.txt',
                'galactic_name': 'gll_iem_v06',
                'extra_galactic_name': 'iso_source_v06',
                'model_maker_flag': 'E2CAT=True, oldNames=True'
            },
            7: {  # P8R3 new iso
                'gll_file': 'galdiffuse/gll_iem_v06.fits',
                'isotropic_file': 'galdiffuse/iso_P8R3_SOURCE_V2.txt',
                'galactic_name': 'gll_iem_v06',
                'extra_galactic_name': 'iso_P8R3_SOURCE_V2',
                'model_maker_flag': 'E2CAT=True, oldNames=True'
            },
            8: {  # 4FGL
                'gll_file': 'galdiffuse/gll_iem_v07.fits',
                'isotropic_file': 'galdiffuse/iso_P8R3_SOURCE_V2_v1.txt',
                'galactic_name': 'gll_iem_v07',
                'extra_galactic_name': 'iso_P8R3_SOURCE_V2_v1',
                'model_maker_flag': 'E2CAT=True, oldNames=True'
            },
            9: {  # v2.0 (default)
                'gll_file': 'galdiffuse/gll_iem_v07.fits',
                'isotropic_file': 'galdiffuse/iso_P8R3_SOURCE_V3_v1.txt',
                'galactic_name': 'gll_iem_v07',
                'extra_galactic_name': 'iso_P8R3_SOURCE_V3_v1',
                'model_maker_flag': 'E2CAT=True, oldNames=True, makeRegion=False'
            }
        }
        
        config = configs.get(self.irf_code, configs[9])
        
        # Prepend FERMI_DIR/refdata/fermi/ to file paths
        config['gll_file'] = os.path.join(fermi_dir, 'refdata', 'fermi', config['gll_file'])
        config['isotropic_file'] = os.path.join(fermi_dir, 'refdata', 'fermi', config['isotropic_file'])
        
        return config
    
    def _get_fits_keyword(self, filename, keyword):
        """Get a keyword value from a FITS file."""
        result = subprocess.run(
            ['fkeypar', f'{filename}[1]', keyword],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"Failed to read {keyword} from {filename}")
        
        result = subprocess.run(
            ['pget', 'fkeypar', 'value'],
            capture_output=True, text=True
        )
        return float(result.stdout.strip())
    
    def _get_update_start_time(self, lc_file):
        """Get the start time for update mode from existing light curve."""
        with open(lc_file, 'r') as f:
            # Skip header (3 lines)
            for _ in range(3):
                f.readline()
            
            # Read last two data lines
            previous_line = None
            for line in f:
                previous_line = line
            
            if previous_line:
                fields = previous_line.strip().split()
                # time - half_bin_width
                return self.mjd2met(float(fields[0]) - float(fields[3]))
        
        return self.start_time
    
    def _generate_model_file(self, evfile, model_file, ra, dec):
        """
        Generate XML model file using make4FGLxml.
        Creates a temporary Python script that calls the model maker.
        """
        self.logger.info("Generating XML model file using make4FGLxml")
        
        # Get MY_FERMI_DIR or use FERMI_DIR
        my_fermi_dir = os.environ.get('MY_FERMI_DIR', os.environ.get('FERMI_DIR', ''))
        fermi_dir = os.environ.get('FERMI_DIR', '')
        
        if not my_fermi_dir:
            self.logger.warning("MY_FERMI_DIR not set, using FERMI_DIR")
            my_fermi_dir = fermi_dir
        
        if not fermi_dir:
            raise RuntimeError("FERMI_DIR environment variable not set")
        
        # Get diffuse file configuration
        diffuse_config = self._get_diffuse_files()
        
        # Determine catalog path
        # If catalog is relative, look in MY_FERMI_DIR
        if not os.path.isabs(self.catalog):
            catalog_path = os.path.join(my_fermi_dir, self.catalog)
            # Also try current directory
            if not os.path.exists(catalog_path) and os.path.exists(self.catalog):
                catalog_path = self.catalog
        else:
            catalog_path = self.catalog
        
        if not os.path.exists(catalog_path):
            raise FileNotFoundError(
                f"Catalog file not found: {catalog_path}\n"
                f"Tried: {self.catalog} and {catalog_path}\n"
                f"Please ensure catalog file is in current directory or MY_FERMI_DIR"
            )
        
        self.logger.info(f"Using catalog: {catalog_path}")
        
        # Create temporary Python script
        py_script = tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False)
        
        try:
            # Write the Python script
            py_script.write("#!/usr/bin/env python\n")
            py_script.write("import sys\n")
            py_script.write(f"sys.path.append('{my_fermi_dir}')\n")
            py_script.write("\n")
            py_script.write("# Import the model maker\n")
            py_script.write("try:\n")
            py_script.write(f"    from {MODEL_MAKER} import *\n")
            py_script.write("except ImportError as e:\n")
            py_script.write(f"    print(f'Error: Could not import {MODEL_MAKER}')\n")
            py_script.write("    print(f'Make sure {MODEL_MAKER}.py is in MY_FERMI_DIR: {my_fermi_dir}')\n")
            py_script.write("    print(f'Import error: {{e}}')\n")
            py_script.write("    sys.exit(1)\n")
            py_script.write("\n")
            py_script.write(f"s1 = srcList('{catalog_path}', '{evfile}', '{model_file}')\n")
            
            # Build makeModel call
            gll_file = diffuse_config['gll_file']
            gal_name = diffuse_config['galactic_name']
            iso_file = diffuse_config['isotropic_file']
            eg_name = diffuse_config['extra_galactic_name']
            model_flags = diffuse_config['model_maker_flag']
            
            makemodel_call = (
                f"s1.makeModel('{gll_file}', '{gal_name}', "
                f"'{iso_file}', '{eg_name}', psForce=True"
            )
            
            if model_flags:
                makemodel_call += f", {model_flags}"
            
            makemodel_call += ")\n"
            
            py_script.write(makemodel_call)
            py_script.close()
            
            self.logger.info(f"Created model generation script: {py_script.name}")
            if self.debug_mode:
                self.logger.debug(f"Script contents:")
                with open(py_script.name, 'r') as f:
                    for line in f:
                        self.logger.debug(f"  {line.rstrip()}")
            
            # Run the script
            self.logger.info("Running make4FGLxml...")
            result = subprocess.run(
                ['python', py_script.name],
                capture_output=True,
                text=True
            )
            
            if result.returncode != 0:
                self.logger.error(f"make4FGLxml failed with return code {result.returncode}")
                if result.stderr:
                    self.logger.error(f"Error output:\n{result.stderr}")
                if result.stdout:
                    self.logger.error(f"Standard output:\n{result.stdout}")
                
                raise RuntimeError(
                    f"Failed to generate model file.\n"
                    f"Make sure {MODEL_MAKER}.py is available in MY_FERMI_DIR: {my_fermi_dir}\n"
                    f"Error: {result.stderr}"
                )
            
            if result.stdout:
                # Always show output in non-batch mode
                if not self.batch_mode:
                    print(result.stdout)
                # Log it regardless
                self.logger.debug(result.stdout)
            
            # Verify model file was created
            if not os.path.exists(model_file):
                raise RuntimeError(f"Model file was not created: {model_file}")
            
            self.logger.info(f"Successfully created model file: {model_file}")
            
        finally:
            # Clean up the temporary Python script unless in debug mode
            if os.path.exists(py_script.name):
                if self.debug_mode:
                    self.logger.debug(f"Debug mode: keeping script {py_script.name}")
                else:
                    os.unlink(py_script.name)
    
    def _apply_barycenter(self, lc_file, scfile, ra, dec):
        """Apply barycenter correction to light curve."""
        self.run_fermi_tool(
            'gtbary',
            evfile=lc_file,
            scfile=scfile,
            ra=ra,
            dec=dec,
            tcorrect='BARY'
        )
    
    def read_parameter_file(self, filename):
        """
        Read parameters from a parameter file.
        Format: parameter_name value (one per line)
        Lines starting with # are comments
        """
        self.logger.info(f"Reading parameters from: {filename}")
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and blank lines
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split(None, 1)  # Split on whitespace, max 2 parts
                
                if len(parts) == 2:
                    param, value = parts
                    param = param.lower()
                    
                    # Map parameters to attributes
                    if param == 'bin_size':
                        self.bin_size = float(value)
                    elif param == 'roi':
                        self.roi = float(value)
                    elif param == 'roi_probability':
                        # This is the default ROI for probability mode
                        pass
                    elif param == 'zenith_limit':
                        self.zenith_limit = float(value)
                    elif param == 'rock':
                        self.rock = float(value)
                    elif param == 'bore':
                        self.bore = float(value)
                    elif param == 'event_class_min':
                        self.event_class_min = int(value)
                    elif param == 'source_list':
                        self.source_list = value
                    elif param == 'ft1_list':
                        self.ft1_list = value
                    elif param == 'ft2':
                        self.ft2 = value
                    elif param == 'catalog':
                        self.catalog = value
                    elif param == 'pthreshold':
                        self.pthreshold = float(value)
                    elif param == 'irf_code':
                        self.irf_code = int(value)
                    elif param == 'sun_minimum':
                        self.sun_minimum = float(value)
                    elif param == 'emin':
                        self.emin = float(value)
                    elif param == 'emax':
                        self.emax = float(value)
                    elif param == 'spectral_index':
                        self.spectral_index = float(value)
                    elif param == 'use_probability':
                        self.use_probability = value.lower() in ('true', '1', 'yes', 'y')
                    elif param == 'barycenter':
                        self.barycenter = value.lower() in ('true', '1', 'yes', 'y')
                    elif param == 'prefix':
                        self.prefix = value if len(parts) == 2 else ''
                    else:
                        self.logger.warning(f"Unknown parameter: {param}")
        
        self.logger.info("Parameters loaded from file")
    
    def write_parameter_file(self, filename):
        """
        Write current parameters to a parameter file.
        Creates a file that can be read back with read_parameter_file().
        """
        self.logger.info(f"Writing parameters to: {filename}")
        
        with open(filename, 'w') as f:
            f.write("# BEX_FERMI parameter file\n")
            f.write(f"# Created: {__import__('datetime').datetime.now()}\n")
            f.write("#\n")
            f.write("# Format: parameter_name value\n")
            f.write("# Lines starting with # are comments\n")
            f.write("#\n\n")
            
            f.write("# Time binning\n")
            f.write(f"bin_size {self.bin_size}\n\n")
            
            f.write("# Energy range (MeV)\n")
            f.write(f"emin {self.emin}\n")
            f.write(f"emax {self.emax}\n\n")
            
            f.write("# Probability photometry\n")
            f.write(f"use_probability {self.use_probability}\n")
            if self.use_probability:
                f.write(f"catalog {self.catalog}\n")
            f.write(f"pthreshold {self.pthreshold}\n\n")
            
            f.write("# Region of interest\n")
            f.write(f"roi {self.roi}\n")
            f.write(f"roi_probability {DEFAULT_ROI_PROBABILITY}\n\n")
            
            f.write("# Selection criteria\n")
            f.write(f"zenith_limit {self.zenith_limit}\n")
            f.write(f"rock {self.rock}\n")
            f.write(f"bore {self.bore}\n")
            f.write(f"event_class_min {self.event_class_min}\n")
            f.write(f"sun_minimum {self.sun_minimum}\n\n")
            
            f.write("# IRF and spectral model\n")
            f.write(f"irf_code {self.irf_code}\n")
            f.write(f"spectral_index {self.spectral_index}\n\n")
            
            f.write("# Barycenter correction\n")
            f.write(f"barycenter {self.barycenter}\n\n")
            
            f.write("# Input files\n")
            f.write(f"source_list {self.source_list}\n")
            f.write(f"ft1_list {self.ft1_list}\n")
            f.write(f"ft2 {self.ft2}\n\n")
            
            f.write("# Output prefix (optional)\n")
            f.write(f"prefix {getattr(self, 'prefix', DEFAULT_PREFIX)}\n")
        
        print(f"Parameter file created: {filename}")
        print(f"You can now run: bex_fermi.py -batch -file {filename}")


def main():
    """Main entry point for command-line execution."""
    parser = argparse.ArgumentParser(
        description='Probability-weighted aperture photometry for Fermi LAT data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                          # Interactive mode (prompts for parameters)
  %(prog)s -batch -file bex.par     # Batch mode with parameter file
  %(prog)s -update                  # Update existing light curves
  %(prog)s -start 58000 -stop 59000 # Specify time range

Output format:
  Time (MJD), rate (ph/cm^2/s), count-rate error, bin half-width (days),
  exposure-based error, exposure (cm^2 s)
        """
    )
    
    parser.add_argument('-v', '--version', action='version', 
                       version=f'%(prog)s {VERSION}')
    parser.add_argument('-batch', action='store_true',
                       help='Run in batch mode (no prompting, requires -file)')
    parser.add_argument('-debug', action='store_true',
                       help='Enable debug output')
    parser.add_argument('-update', action='store_true',
                       help='Add points to existing light curve')
    parser.add_argument('-file', metavar='PARFILE',
                       help='Read parameters from file')
    parser.add_argument('-makepar', metavar='PARFILE',
                       help='Create parameter file with defaults')
    parser.add_argument('-start', type=float, metavar='MJD',
                       help='Start time in MJD')
    parser.add_argument('-stop', type=float, metavar='MJD',
                       help='Stop time in MJD')
    parser.add_argument('-pthreshold', type=float, metavar='PROB',
                       help=f'Probability threshold (default: {DEFAULT_PTHRESHOLD})')
    
    args = parser.parse_args()
    
    # Build configuration
    config = {
        'batch': args.batch,
        'debug': args.debug,
        'update': args.update,
    }
    
    if args.pthreshold is not None:
        config['pthreshold'] = args.pthreshold
    
    if args.start:
        config['start_time'] = BexFermi.mjd2met(args.start)
    if args.stop:
        config['stop_time'] = BexFermi.mjd2met(args.stop)
    
    # Create processor
    processor = BexFermi(config)
    
    print(f"bex_fermi.py: V{VERSION} (R. Corbet / Python conversion)")
    
    # Handle makepar mode
    if args.makepar:
        processor.write_parameter_file(args.makepar)
        return
    
    # Check environment
    processor.check_environment()
    
    # Read parameters from file if specified
    # This sets the defaults, but doesn't prevent prompting unless -batch is used
    if args.file:
        processor.read_parameter_file(args.file)
        processor.logger.info(f"Read defaults from parameter file: {args.file}")
    
    # In interactive mode (not batch), prompt for parameters
    # The prompts will use values from parameter file as defaults if -file was used
    if not args.batch:
        processor.prompt_for_parameters()
    elif not args.file:
        # Batch mode requires a parameter file
        print("\nError: Batch mode requires -file parameter")
        print("Use -file <parfile> to specify parameter file")
        print("Or run without -batch for interactive mode")
        sys.exit(1)
    
    # Process sources
    processor.process_sources()


if __name__ == '__main__':
    main()
