import xarray as xr
import numpy as np
import sys


def cor_ci_bootyears(ts1, ts2, seed=None, nboots=1000, conf=95):
    """ """
    ptilemin = (100.-conf)/2.
    ptilemax = conf + (100-conf)/2.

    if (ts1.size != ts2.size):
        print("The two arrays must have the same size")
        sys.exit()

    if (seed):
        np.random.seed(seed)

    samplesize = ts1.size
    ranu = np.random.uniform(0, samplesize, nboots*samplesize)
    ranu = np.floor(ranu).astype(int)

    bootdat1 = np.array(ts1[ranu])
    bootdat2 = np.array(ts2[ranu])
    bootdat1 = bootdat1.reshape([samplesize, nboots])
    bootdat2 = bootdat2.reshape([samplesize, nboots])
   
 
    bootcor = xr.corr(xr.DataArray(bootdat1), xr.DataArray(bootdat2), dim='dim_0')
    minci = np.percentile(bootcor,ptilemin)
    maxci = np.percentile(bootcor,ptilemax)

    return minci, maxci 

