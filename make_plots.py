'''
Run plots for drifters from these simulations.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob
import cm_pong
from matplotlib.mlab import find

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 26})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, llcrnrlat=27.01, 
        urcrnrlat=30.5, llcrnrlon=-97.8, urcrnrlon=-87.7)
# actually using psi grid here despite the name
xr = np.asanyarray(grid['xpsi'].T, order='C')
yr = np.asanyarray(grid['ypsi'].T, order='C')

## Model output ##
m = netCDF.Dataset(loc)

# Model time period to use
units = m.variables['ocean_time'].units
year = 2008
starttime = netCDF.date2num(datetime(year, 5, 1, 12, 0, 0), units)
endtime = netCDF.date2num(datetime(year, 9, 1, 12, 0, 0), units)
dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
dates = netCDF.num2date(m.variables['ocean_time'][:], units)

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity
cmap = cm_pong.salinity('YlGnBu_r', levels)
ilevels = [0,1,2,3,4,5,8] # which levels to label
ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
##

    # # Change axis and label color
    # ax.spines['bottom'].set_color('0.2')
    # ax.spines['top'].set_color('0.2')
    # ax.spines['left'].set_color('0.2')
    # ax.spines['right'].set_color('0.2')
    # ax.xaxis.label.set_color('0.2')
    # ax.yaxis.label.set_color('0.2')
    # ax.tick_params(axis='x', colors='0.2')
    # ax.tick_params(axis='y', colors='0.2')


# Loop through times that simulations were started
for t in ts:

    # Set up before plotting
    itmodel = find(t==m.variables['ocean_time'][:])[0] # index for model output at this time

    figname = 'figures/' + str(year) + '/' + lastsimdate.isoformat()[0:13] + '.png'

    # Don't redo plot
    if os.path.exists(figname):
        continue

    # Set up plot
    fig = plt.figure(figsize=(17,9))
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax, outline=False)

    # Date
    date = dates[itmodel].strftime('%Y %b %02d %H:%M')
    ax.text(0.735, 0.04, date, fontsize=24, color='0.2', transform=ax.transAxes, 
                bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

    # Plot surface salinity
    # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
    salt = np.squeeze(m.variables['salt'][itmodel,-1,1:-1,1:-1])
    mappable = ax.pcolormesh(xr, yr, salt, cmap=cmap, vmin=0, vmax=36)
    # Plot Sabine too, which gets covered by the basemap
    sabmask = ~salt[172:189,332:341].mask.astype(bool)
    sabmask[3,2] = False
    sabmask[3,3] = False
    sabmask[4,1] = False
    sabmask[4,2] = False
    sabmask[5,0] = False
    sabmask[5,1] = False
    sabmask[6,0] = False
    sabmask[4,7] = False
    sabmask[8:14,4] = False
    sabmask[15,7] = False
    sabmask[16,7] = False
    sabmask[3:5,5:7] = False
    salt[172:189,332:341] = np.ma.masked_where(~sabmask,salt[172:189,332:341])
    ax.pcolormesh(xr[172:189,332:341], yr[172:189,332:341], salt[172:189,332:341], cmap=cmap, vmin=0, vmax=36, zorder=2)

    # Colorbar in upper left corner
    cax = fig.add_axes([0.15, 0.75, 0.3, 0.03]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Surface salinity [psu]', fontsize=20)
    cb.ax.tick_params(labelsize=18) 
    cb.set_ticks(ticks)

    plt.savefig(figname, bbox_inches='tight', dpi=100)
    plt.close(fig)
