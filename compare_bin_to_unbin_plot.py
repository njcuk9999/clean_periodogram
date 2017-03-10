#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09/03/17 at 7:29 PM

@author: neil

Program description here

Version 0.0.0
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import matplotlib.pyplot as plt

# =============================================================================
# Define variables
# =============================================================================
WORKSPACE = '/Astro/Projects/RayPaul_Work/SuperWASP/'

# Save path
SAVEPATH = WORKSPACE + '/Data/'

FILES = ['Amps_for_unbinned_ARG_54.fits',
         'Amps_for_binned_0.1_ARG_54.fits',
         'Amps_for_binned_0.01_ARG_54.fits']

COLOURS = ['r', 'b', 'g']
LABELS = ['Unbinned', 'binsize=0.1','binsize=0.01']

# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # set up plot
    plt.close()
    fig, frames= plt.subplots(ncols=1, nrows=len(FILES))
    # loop around files
    for f in range(len(FILES)):

        frame = frames[f]
        # load data
        data = fits.getdata(SAVEPATH + FILES[f])
        # get columns
        time = np.array(data['times'])
        amps = np.array(data['amps'])
        # plot line
        frame.plot(time, amps/np.nanmax(amps), color=COLOURS[f], label=LABELS[f])
        frame.set_xscale('log')
        frame.legend(loc=1)
        frame.set_xlim(1e-3, 1e3)
        frame.grid(b=True, which='major', color='0.5', linestyle='-', alpha=0.5)
        frame.grid(b=True, which='minor', color='0.25', linestyle='--',
                   alpha=0.25)
        frame.set_ylabel('Amplitude')
        if f not in [0, len(FILES)-1]:
            frame.set_xticklabels([])

    frames[0].xaxis.tick_top()
    frames[0].set_xlabel('Time / days')
    frames[0].xaxis.set_label_position('top')

    frames[-1].set_xlabel('Time / days')

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
    plt.close()



# =============================================================================
# End of code
# =============================================================================
