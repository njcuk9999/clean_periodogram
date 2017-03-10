#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 010/03/17 at 10:19 PM

@author: neil

Program description here

Version 0.0.0
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from clean_periodogram import bin_data, clean_periodogram
import MySQLdb
import pandas
from tqdm import tqdm
from astropy.stats import LombScargle


# =============================================================================
# Define variables
# =============================================================================
WORKSPACE = '/Astro/Projects/RayPaul_Work/SuperWASP/'
# test data 2
DATAPATH = WORKSPACE + '/Data/messina_match_from_paul.fits'
# Save path
SAVEPATH = WORKSPACE + '/Data/Messina_matches/'
PLOTPATH = WORKSPACE + '/Plots/Messina_matches/'
# -----------------------------------------------------------------------------
# set database settings
HOSTNAME = 'localhost'
USERNAME = 'root'
PASSWORD = '1234'
DATABASE = 'swasp'
TABLE = 'swasp_sep16_tab'
# -----------------------------------------------------------------------------
# To use defaults (from clean_periodogram) set this to None
# Note: by design we don't use the smallest 50% of times
#       or the largest 50% of frequencies in the DFT
#       so MUST define frequencies accordingly
FREQ1 = None

# log spacing
TIME_GRID = 10**np.linspace(2*-3, 3, 100000)
FREQ2 = 1.0/TIME_GRID

# linear spacing
FREQ3 = np.linspace(1.0e-3, 2.0e3, 100000)
# -----------------------------------------------------------------------------
# define whether to bin and the bin size, and colours of each (for plot) etc
BINDATA_ARR = [False, True, True, False, True, True, False, True, True]
BINSIZE_ARR = [None, 0.1, 0.01, None, 0.1, 0.01, None, 0.1, 0.01]
FREQ_ARR = [FREQ1, FREQ1, FREQ1, FREQ2, FREQ2, FREQ2, FREQ3, FREQ3, FREQ3]

COLOURS = ['r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b']
LABELS = ['Unbinned computed freq',
          'binsize=0.1  computed freq',
          'binsize=0.01  computed freq',
          'Unbinned log spacing',
          'binsize=0.1  log spacing',
          'binsize=0.01  log spacing',
          'Unbinned linear spacing',
          'binsize=0.1  linear spacing',
          'binsize=0.01  linear spacing']
# define any set we don't want to run (using position in above arrays)
SKIP = [0]


# =============================================================================
# Define functions
# =============================================================================
def load_db(db_name, host, uname, pword):
    conn1 = MySQLdb.connect(host=host, user=uname, db=db_name,
                            connect_timeout=100000, passwd=pword)
    c1 = conn1.cursor()
    return c1, conn1


def get_list_of_objects_from_db(conn):
    # ----------------------------------------------------------------------
    # find all systemids
    print("\nGetting list of objects...")
    query1 = "SELECT CONCAT(c.systemid,'_',c.comp)"
    query1 += " AS sid FROM {0} AS c".format(TABLE)
    query1 += " where c.systemid is not null and c.systemid <> ''"
    rawdata = pandas.read_sql_query(query1, conn)
    rawsystemid = np.array(rawdata['sid'])
    # get list of unique ids (for selecting each as a seperate curve)
    sids = np.array(np.unique(rawsystemid), dtype=str)
    # return list of objects
    return sids


def get_targets_from_file():
    targetdata = fits.getdata(DATAPATH, ext=1)
    sysids = np.array(targetdata['systemid'])
    comps = np.array(targetdata['comp'])

    targets = []
    for sit in range(len(sysids)):
        targets.append('{0}_{1}'.format(sysids[sit], comps[sit]))

    return np.array(targets)


def get_data_from_db(sid, conn):
    # get data using SQL query on database
    query2 = "SELECT * FROM swasp_sep16_tab AS c"
    query2 += " WHERE CONCAT(c.systemid,'_',c.comp) = '{0}'"
    data = pandas.read_sql_query(query2.format(sid), conn)
    x, y = np.array(data['HJD']), np.array(data['FLUX2'])
    ey = np.array(data['FLUX2_ERR'])
    # sort into order (by x)
    sortx = np.argsort(x)
    x, y, ey = x[sortx], y[sortx], ey[sortx]
    # remove infinities
    m = np.isfinite(y) * np.isfinite(x)
    xm, ym, eym = x[m], y[m], ey[m]
    return xm, ym, eym


def run_clean_periodogram(time, data, edata, sid, bdata, bsize, filelist,
                          label, freq=None):
    # -----------------------------------------------------------------
    # bin data
    if bdata:
        if edata is None:
            time, data = bin_data(time, data, binsize=bsize, log=True)
        else:
            time, data = bin_data(time, data, edata, binsize=bsize, log=True)
    # -----------------------------------------------------------------
    # Run clean
    freqs, wfn, dft, cdft = clean_periodogram(time, data, freqs=freq,
                                              log=True, full=True)
    # -----------------------------------------------------------------
    # construct save name
    savename = '{0}_{1}.fits'.format(sid, label)

    # Save to file
    table = Table()
    table['times'] = 1.0 / freqs[0: len(cdft)]
    table['freqs'] = freqs[0: len(cdft)]
    table['amps'] = 2.0 * abs(cdft)
    table.write(SAVEPATH + savename, overwrite=True)
    filelist.append(SAVEPATH + savename)
    return filelist


def do_lombscargle(thedata, time):
    # commute lombscargle
    tvec, dvec, edvec = thedata

    mask = np.isfinite(tvec) & np.isfinite(dvec)
    tvec, dvec = tvec[mask], dvec[mask]
    tvec = tvec - tvec.min()
    power = LombScargle(tvec, dvec, dy=edvec).power(1.0/time)

    return 1.0/time, abs(power)


def plot_frame(framearr, skip, f, filelist, rdata):
    # skip if needed
    frame1, frame2 = framearr[f]
    if f in skip:
        frame1.text(0.5, 0.5, '{0} skipped'.format(LABELS[frow]),
                    transform=frame1.transAxes,
                    horizontalalignment='center', verticalalignment='center')
        frame2.text(0.5, 0.5, '{0} skipped'.format(LABELS[frow]),
                    transform=frame2.transAxes,
                    horizontalalignment='center', verticalalignment='center')
    else:
        # load data
        data = fits.getdata(filelist[f])
        # get columns
        time1 = np.array(data['times'])
        amps1 = np.array(data['amps'])

        # plot clean periodogram
        frame2.plot(time1, amps1, color=COLOURS[f], label=LABELS[f], zorder=2)
        frame2.set_ylabel('Amplitude')
        frame2.legend(loc=1, title='CLEAN periodogram')

        # plot lombscagle in background
        freqs2, power2 = do_lombscargle(rdata[f], time1)
        frame1.plot(1.0/freqs2, power2, color='k', zorder=1, label=LABELS[f])
        frame1.legend(loc=1, title='Lombscargle')
        frame1.set_ylabel('Power')

    # frame settings
    for frame in framearr[f]:
        frame.set_xscale('log')
        frame.set_xlim(1.0e-3, 1.0e3)
        frame.grid(b=True, which='major', color='0.5',
                   linestyle='-', alpha=0.5)
        frame.grid(b=True, which='minor', color='0.25',
                   linestyle='--', alpha=0.25)
        frame1.set_xlabel('Time / days')

    for frame in framearr[f]:
        if f == 0:
            frame.xaxis.tick_top()
            frame.set_xlabel('Time / days')
            frame.xaxis.set_label_position('top')
        elif f == len(framearr) - 1:
            frame.set_xlabel('Time / days')
        else:
            frame.set_xticklabels([])


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # load database
    # Must have database running
    # mysql -u root -p
    print("\nConnecting to database...")
    c, conn1 = load_db(DATABASE, HOSTNAME, USERNAME, PASSWORD)
    # -------------------------------------------------------------------------
    # get list of objects from database
    sid_list = get_list_of_objects_from_db(conn1)
    # -------------------------------------------------------------------------
    # load targets from DATAPATH
    targets = get_targets_from_file()

    for target in targets:
        # finding target in database list of objects
        s_id = str(target).replace(' ', '')
        if s_id not in sid_list:
            print('\n Error: {0} not in database, skipping...'.format(target))
        # Run loop for each
        print('\n\n\n Running clean periodogram loop....\n')
        files, results = [], []
        for frow in range(len(BINDATA_ARR)):
            # skip parameters
            if frow in SKIP:
                files.append(''), results.append([])
                continue
            # get this iterations parameters
            bindata, binsize = BINDATA_ARR[frow], BINSIZE_ARR[frow]
            # print details of this run
            pargs = ['\n'+ '='*50 + '\n', s_id, LABELS[frow]]
            print('{0}\n\t sid={1} {2} {0}'.format(*pargs))
            # -----------------------------------------------------------------
            # Load data from data base
            result = get_data_from_db(s_id, conn1)
            # run the clean periodogram code and save to file
            rargs = list(result) + [s_id, bindata, binsize, files, LABELS[frow]]
            files = run_clean_periodogram(*rargs, freq=FREQ_ARR[frow])
            # store the input data
            results.append(result)
        # set up plot
        print('\n\n\n Plotting graph...\n')
        plt.close()
        fig, frames = plt.subplots(ncols=2, nrows=len(files))
        fig.set_size_inches(16, 4*len(LABELS))
        # loop around files
        for frow in tqdm(range(len(files))):
            plot_frame(frames, SKIP, frow, files, results)
        # adjust the plot and save
        plt.subplots_adjust(wspace=0.2, hspace=0)
        sname = '{0}_bin_comparison_plot'.format(s_id)
        plt.savefig(PLOTPATH + sname + '.png', bbox_inches='tight')
        plt.savefig(PLOTPATH + sname + '.pdf', bbox_inches='tight')
        plt.close()

        # input('Press Enter to continue. Ctrl+c to Exit')

# =============================================================================
# End of code
# =============================================================================
