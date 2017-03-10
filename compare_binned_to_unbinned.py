#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09/03/17 at 4:19 PM

@author: neil

Program description here

Version 0.0.0
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
from clean_periodogram import bin_data, clean_periodogram

# =============================================================================
# Define variables
# =============================================================================
WORKSPACE = '/Astro/Projects/RayPaul_Work/SuperWASP/'

# Test data 1
TESTPATH = WORKSPACE + './CLEAN_periodogram_IDL/test.fits'

# test data 2
TESTPATH = WORKSPACE + '/Data/Elodie/ARG_54_lightcurve.fits'

# Save path
SAVEPATH = WORKSPACE + '/Data/Amps_for_binned_0.01_ARG_54.fits'

# whether to bin data
BINDATA = True
BINSIZE = 0.01


# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # Load data
    print('\n Loading data...')
    fitsrec = fits.getdata(TESTPATH, ext=1)
    time_arr = np.array(fitsrec['time'])
    data_arr = np.array(fitsrec['flux'])
    if 'eflux' in fitsrec.columns.names:
        edata_arr = np.array(fitsrec['eflux'])
    else:
        edata_arr = None
    # ----------------------------------------------------------------------
    # bin data
    if BINDATA:
        if edata_arr is None:
            time_arr, data_arr = bin_data(time_arr, data_arr, binsize=BINSIZE,
                                          log=True)
        else:
            time_arr, data_arr = bin_data(time_arr, data_arr, edata_arr,
                                          binsize=BINSIZE, log=True)
    # ----------------------------------------------------------------------
    # Run clean
    results = clean_periodogram(time_arr, data_arr, log=True, full=True)
    freqs, wfn, dft, cdft = results
    # ----------------------------------------------------------------------
    # Save to file
    table = Table()
    table['times'] = 1.0/freqs[0: len(cdft)]
    table['freqs'] = freqs[0: len(cdft)]
    table['amps'] = 2.0*abs(cdft)
    table.write(SAVEPATH, overwrite=True)

# =============================================================================
# End of code
# =============================================================================
