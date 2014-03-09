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
from datetime import datetime, timedelta
import glob
from cmPong import cmPong
from matplotlib.mlab import find

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 16})
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
        urcrnrlat=30.6, llcrnrlon=-97.8, urcrnrlon=-87.7)
# actually using psi grid here despite the name
xr = np.asanyarray(grid['xpsi'].T, order='C')
yr = np.asanyarray(grid['ypsi'].T, order='C')

## Model output ##
m = netCDF.Dataset(loc)

# Model time period to use
units = m.variables['ocean_time'].units
year = 2008
starttime = netCDF.date2num(datetime(year, 5, 1, 12, 0, 0), units)
endtime = netCDF.date2num(datetime(year, 10, 1, 12, 0, 0), units)
dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
dates = netCDF.num2date(m.variables['ocean_time'][:], units)

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity
cmap = cmPong.salinity('YlGnBu_r', levels)
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

    figname = 'figures/' + str(year) + '/' + dates[itmodel].isoformat()[0:13] + '.png'

    # Don't redo plot
    if os.path.exists(figname):
        continue

    # Set up plot
    fig = plt.figure(figsize=(10.24, 4.3), dpi=100)
    ax = fig.add_axes([0.065, 0.045, 0.925, 0.96])
    # ax = fig.add_subplot(111)
    ax.set_frame_on(False) # kind of like it without the box
    tracpy.plotting.background(grid=grid, ax=ax, outline=False, mers=np.arange(-97, -88))

    # Date
    date = dates[itmodel].strftime('%Y %b %02d %H:%M')
    ax.text(0.7, 0.04, date, fontsize=20, color='0.2', transform=ax.transAxes, 
                bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

    # # PONG
    # ax.text(0.8, 0.95, 'PONG.TAMU.EDU', fontsize=14, transform=ax.transAxes, color='0.3')

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

    # Mississippi river discharge rate


    # Wind for several days

    # Colorbar in upper left corner
    cax = fig.add_axes([0.085, 0.925, 0.35, 0.03]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label(r'Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=16)
    cb.ax.tick_params(labelsize=14, length=2) 
    cb.set_ticks(ticks)

    plt.savefig(figname)
    plt.close(fig)

# To make movie: ffmpeg -r 10 -pattern_type glob -i '2008-*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 movie.mp4