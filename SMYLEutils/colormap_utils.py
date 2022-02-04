import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap ## used to create custom colormaps
import matplotlib.colors as mcolors
import numpy as np

def blue2red_cmap(n):
    """ combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    """

    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.YlOrRd(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))

    colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap

def blue2red_acc_cmap(levs,cutoff):
    """ combine three color maps to create a diverging color map with white in the middle,
    grey scale for positive below cutoff, and red scale for positive above cutoff.
    levs = ACC levels
    cutoff = positive value to highlight (e.g., 0.5)
    """
    n = levs.size
    nneg = (np.where(levs<0,True,False)).sum()
    npos = (np.where(levs>0,True,False)).sum()
    nlow2 = (np.where(levs<-cutoff,True,False)).sum()
    nhigh2 = (np.where(levs>cutoff,True,False)).sum()
    nlow1 = nneg-nlow2
    nhigh1 = npos-nhigh2
    nwhite = n-(nneg+npos)+1

    colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg-2)))
    colors3 = plt.cm.binary(np.linspace(0,0.5, int(nhigh1)))
    colors4 = plt.cm.YlOrRd(np.linspace(0,1, int(nhigh2)))
    colorsw = np.ones((nwhite,4))
    colors = np.vstack((colors1, colorsw, colors3, colors4))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap

def precip_cmap(n):
    """ combine two existing color maps to create a diverging color map with white in the middle.
    browns for negative, blues for positive
    n = the number of contour intervals
    """
    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    colors1 = plt.cm.YlOrBr_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.GnBu(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))

    colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap

