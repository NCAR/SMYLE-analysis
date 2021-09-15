# routines for dealing with the calendar e.g., calculating timeseries of seasonal means
import xarray as xr
import numpy as np
import pandas as pd
from math import nan
import cftime
import cf_units


def season_ts(ds, season, var=None):
    """ calculate timeseries of seasonal averages
    Args: ds (xarray.Dataset): dataset
          var (str): variable to calculate 
          season (str): 'DJF', 'MAM', 'JJA', 'SON' or 'DJFM'
          
          Author: I. Simpson
    """

    if (season == 'DJFM'):
        ds_season = ds.where(
        (ds['time.month'] == 12) | (ds['time.month'] == 1) | (ds['time.month'] == 2) | (ds['time.month'] == 3))
        if (var):
            ds_season = ds_season[var].rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")
        else:
            ds_season = ds_season.rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")

    else:

        ## set months outside of season to nan
        ds_season = ds.where(ds['time.season'] == season)
    
        # calculate 3month rolling mean (only middle months of season will have non-nan values)
        if (var):
            ds_season = ds_season[var].rolling(min_periods=3, center=True, time=3).mean().dropna("time", how='all')
        else:
            ds_season = ds_season.rolling(min_periods=3, center=True, time=3).mean().dropna("time", how="all")

    return ds_season


def time_set_mid(ds, time_name, deep=False):
    """
    Return copy of ds with values of ds[time_name] replaced with midpoints of
    ds[time_name].attrs['bounds'], if bounds attribute exists.
    Except for time_name, the returned Dataset is a copy of ds2.
    The copy is deep or not depending on the argument deep.
    
    Author: K. Lindsay
    """

    ds_out = ds.copy(deep)

    if "bounds" not in ds[time_name].attrs:
        return ds_out

    tb_name = ds[time_name].attrs["bounds"]
    tb = ds[tb_name]
    bounds_dim = next(dim for dim in tb.dims if dim != time_name)

    # Use da = da.copy(data=...), in order to preserve attributes and encoding.

    # If tb is an array of datetime objects then encode time before averaging.
    # Do this because computing the mean on datetime objects with xarray fails
    # if the time span is 293 or more years.
    #     https://github.com/klindsay28/CESM2_coup_carb_cycle_JAMES/issues/7
    if tb.dtype == np.dtype("O"):
        units = "days since 0001-01-01"
        calendar = "noleap"
        tb_vals = cftime.date2num(ds[tb_name].values, units=units, calendar=calendar)
        tb_mid_decode = cftime.num2date(
            tb_vals.mean(axis=1), units=units, calendar=calendar
        )
        ds_out[time_name] = ds[time_name].copy(data=tb_mid_decode)
    else:
        ds_out[time_name] = ds[time_name].copy(data=tb.mean(bounds_dim))

    return ds_out


def time_set_midmonth(ds, time_name, deep=False):
    """
    Return copy of ds with values of ds[time_name] replaced with mid-month
    values (day=15) rather than end-month values.
    
    Author: S. Yeager
    """

    ds_out = ds.copy(deep)
    year = ds_out[time_name].dt.year
    month = ds_out[time_name].dt.month
    year = xr.where(month==1,year-1,year)
    month = xr.where(month==1,12,month-1)
    nmonths = len(month)
    newtime = [cftime.DatetimeNoLeap(year[i], month[i], 15) for i in range(nmonths)]
    ds_out[time_name] = newtime

    return ds_out


def mon_to_seas(da):
    """ Converts an Xarray DataArray containing monthly data to one containing 
    seasonal-average data, appropriately weighted with days_in_month. Time coordinate
    of output reflects (approximate) centered time value for DJF, MAM, JJA, SON
    averages.
    
    Author: S. Yeager
    """
    month_length = da.time.dt.days_in_month
    result = ((da * month_length).resample(time='QS-DEC',loffset='45D').sum(skipna=True,min_count=3) /
          month_length.resample(time='QS-DEC',loffset='45D').sum())
    return result

