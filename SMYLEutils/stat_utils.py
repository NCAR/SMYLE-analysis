import xarray as xr
import numpy as np
import sys
import cftime


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


def remove_drift(da, da_time, y1, y2):
    """
    Function to convert raw DP DataArray into anomaly DP DataArray with leadtime-dependent climatology removed.
    --Inputs--
        da:  Raw DP DataArray with dimensions (Y,L,M,...)
        da_time:  Verification time of DP DataArray (Y,L)
        y1:  Start year of climatology
        y2:  End year of climatology
        
    --Outputs--
        da_anom:  De-drifted DP DataArray
        da_climo:  Leadtime-dependent climatology
    
    Author: E. Maroon (modified by S. Yeager)
    """
    d1 = cftime.DatetimeNoLeap(y1,1,1,0,0,0)
    d2 = cftime.DatetimeNoLeap(y2,12,31,23,59,59)
    masked_period = da.where((da_time>d1) & (da_time<d2))
    da_climo = masked_period.mean('M').mean('Y')
    da_anom = da - da_climo

    return da_anom, da_climo
