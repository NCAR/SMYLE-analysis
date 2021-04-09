# routines for dealing with the calendar e.g., calculating timeseries of seasonal means
import xarray as xr
import numpy as np
import pandas as pd
from math import nan

def season_ts(ds, season, var=None):
    """ calculate timeseries of seasonal averages
    Args: ds (xarray.Dataset): dataset
          var (str): variable to calculate 
          season (str): 'DJF', 'MAM', 'JJA', 'SON'
    """
    ## set months outside of season to nan
    ds_season = ds.where(ds['time.season'] == season)

    # calculate 3month rolling mean (only middle months of season will have non-nan values)
    if (var):
        ds_season = ds_season[var].rolling(min_periods=3, center=True, time=3).mean().dropna("time", how='all')
    else:
        ds_season = ds_season.rolling(min_periods=3, center=True, time=3).mean().dropna("time", how="all")
    return ds_season

