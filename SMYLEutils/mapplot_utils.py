import matplotlib.pyplot as plt
import numpy as np
from SMYLEutils import colormap_utils as mycolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm

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
        #mymap.set_over('pink')
        #mymap.set_under('cyan')
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

def adjust_pop_grid(tlon,tlat,field):
    nj = tlon.shape[0]
    ni = tlon.shape[1]
    xL = int(ni/2 - 1)
    xR = int(xL + ni)

    tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360.,tlon)
    lon  = np.concatenate((tlon,tlon+360.),1)
    lon = lon[:,xL:xR]

    if ni == 320:
        lon[367:-3,0] = lon[367:-3,0]+360.
    lon = lon - 360.
    lon = np.hstack((lon,lon[:,0:1]+360.))
    if ni == 320:
        lon[367:,-1] = lon[367:,-1] - 360.

    #-- trick cartopy into doing the right thing:
    #   it gets confused when the cyclic coords are identical
    lon[:,0] = lon[:,0]-1e-8
    
    #-- periodicity
    lat  = np.concatenate((tlat,tlat),1)
    lat = lat[:,xL:xR]
    lat = np.hstack((lat,lat[:,0:1]))

    field = np.ma.concatenate((field,field),1)
    field = field[:,xL:xR]
    field = np.ma.hstack((field,field[:,0:1]))
    return lon,lat,field


def map_pcolor_global_subplot(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
                              nrow,ncol,subplot, proj, labels=True, showland=True, 
                              grid="latlon", cmap="blue2red", facecolor="white",
                              fontsize=15,centrallon=0,tricontour=False,cutoff=0.5):
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
              grid = ('latlon','camfv','camse','pop')
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        cmap = mycolors.blue2red_cmap(nlevs)
        #mymap.set_over('pink')
        #mymap.set_under('cyan')
    elif (cmap == "precip"):
        cmap = mycolors.precip_cmap(nlevs)
    elif (cmap == "blue2red_acc"):
        cmap = mycolors.blue2red_acc_cmap(clevs,cutoff)
    else:
        cmap = mpl.cm.get_cmap(cmap)
    
    norm = BoundaryNorm(clevs, ncolors=cmap.N, clip=True)
    ax = fig.add_subplot(nrow,ncol,subplot, projection=proj)
    ax.set_aspect('auto')
    ax.set_title(titlestr, fontsize=fontsize)
    ax.set_facecolor(facecolor)
    if showland:
        ax.add_feature(cfeature.COASTLINE)
    else:
        ax.add_feature(cfeature.LAND, edgecolor='black',linewidth=0.1, facecolor='grey', zorder=1)
  
    if grid=="latlon" or grid=="camfv":
        dat, lon = add_cyclic_point(dat, coord=lon)
        cntr = ax.pcolormesh(lon, lat, dat, shading='nearest',vmin=clevs.min(),vmax=clevs.max(),  cmap = cmap, norm=norm, rasterized=True, transform=ccrs.PlateCarree())
        
    elif grid=="camse":
        tri, z = get_refined_triang(lon,lat, dat)
        cntr = ax.tripcolor(tri, z, shading='flat',vmin=clevs.min(),vmax=clevs.max(), cmap = cmap, norm=norm,  rasterized=True, transform=ccrs.PlateCarree())
        
    elif grid=="pop":
        lon, lat, dat = adjust_pop_grid(lon, lat, dat)
        cntr = ax.pcolormesh(lon, lat, dat, shading='nearest',vmin=clevs.min(),vmax=clevs.max(), cmap = cmap, norm=norm, rasterized=True, transform=ccrs.PlateCarree())
        
    else:
        raise ValueError('ERROR: unknown grid')
        
    ax.set_global()
    return ax,cntr

def map_pvalsig_global_subplot(axis, pvals, lon, lat, siglvl,
                              facecolor='none', edgecolor='k',s=10,marker="."):
    """ scatterplot of 2D pvals with coordinates lon and lat on top of axis.
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              siglvl = plot dots anywhere below this significance level

    """
    lon2d,lat2d = np.meshgrid(lon, lat)
    tmplon = np.where(pvals>siglvl,lon2d,np.nan)
    tmplat = np.where(pvals>siglvl,lat2d,np.nan)
    axis.scatter(tmplon,tmplat,facecolor=facecolor, edgecolor=edgecolor,s=s,marker=marker)

    return 
