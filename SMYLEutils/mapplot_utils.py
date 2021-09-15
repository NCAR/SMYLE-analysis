import matplotlib.pyplot as plt
import numpy as np
from SMYLEutils import colormap_utils as mycolors

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker

def contourmap_bothoceans_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, centrallon=0):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=centrallon))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="max", transform=ccrs.PlateCarree())

    ax.set_global()


    return ax

def map_contourf_global_subplot(fig, dat, lon, lat, ci, cmin, cmax, titlestr, leftstr, rightstr,
 nrow,ncol,subplot, proj, labels=True, showland=True, extend='neither', grid="latlon", cmap="blue2red", fontsize=15):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval:
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              nrow = number of rows in multipanel plot
              ncol = number of columns in multipanel plot
              subplot = subplot number
              proj = cartopy map projection to use
              labels = True/False (ticks and  labels are plotted if true) 
              showland = True/False (if False, fill over land)
              grid = ('latlon','camfv','pop')
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)
        mymap.set_over('pink')
        mymap.set_under('cyan')
    elif (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)
    else:
        mymap = cmap
        
    ax = fig.add_subplot(nrow,ncol,subplot, projection=proj)
    ax.set_aspect('auto')
    ax.set_title(titlestr, fontsize=fontsize)
    ax.set_title(leftstr, fontsize=fontsize,loc='left')
    ax.set_title(rightstr, fontsize=fontsize,loc='right')
    if showland:
        ax.add_feature(cfeature.COASTLINE)
    else:
        ax.add_feature(cfeature.LAND, edgecolor='black',linewidth=0.1, facecolor='grey', zorder=1)
  
    if grid=="latlon" or grid=="camfv":
        dat, lon = add_cyclic_point(dat, coord=lon)
        cntr = ax.contourf(lon, lat, dat, 
                           levels=clevs,          
                           cmap = mymap, extend=extend, transform=ccrs.PlateCarree())
        
    elif grid=="pop":
        lon, lat, dat = adjust_pop_grid(lon, lat, dat)
        cntr = ax.contourf(lon, lat, dat, 
                           levels=clevs,           
                           cmap = mymap, extend=extend, transform=ccrs.PlateCarree())
        
    else:
        raise ValueError('ERROR: unknown grid')
        
    ax.set_global()
    return ax, cntr


