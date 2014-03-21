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
import bisect
from matplotlib import delaunay

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})
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
xpsi = np.asanyarray(grid['xpsi'].T, order='C')
ypsi = np.asanyarray(grid['ypsi'].T, order='C')
xr = np.asanyarray(grid['xr'].T, order='C')
yr = np.asanyarray(grid['yr'].T, order='C')

## Model output ##
m = netCDF.Dataset(loc)

# Model time period to use
units = m.variables['ocean_time'].units
year = 2008
starttime = netCDF.date2num(datetime(year, 5, 1, 0, 0, 0), units)
endtime = netCDF.date2num(datetime(year, 10, 1, 12, 0, 0), units)
dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
datesModel = netCDF.num2date(m.variables['ocean_time'][:], units)

plotdates = netCDF.num2date(ts, units)
monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,13)]

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity
cmap = cmPong.salinity('YlGnBu_r', levels)
ilevels = [0,1,2,3,4,5,8] # which levels to label
ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
##

## Wind forcing ##
w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
# Wind time period to use
unitsWind = w.variables['time'].units
datesWind = netCDF.num2date(w.variables['time'][:], unitsWind)
wdx = 25; wdy = 30 # in indices
##

## River forcing ##
r = netCDF.Dataset('/atch/raid1/zhangxq/Projects/txla_nesting6/TXLA_river_4dyes_2011.nc')
# River timing
unitsRiver = r.variables['river_time'].units
datesRiver = netCDF.num2date(r.variables['river_time'][:], unitsRiver)
tRiver = r.variables['river_time'][:]
# all of river input
Q = np.abs(r.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
# # start and end dates for river discharge
# tstartRiver = netCDF.date2num(datetime(year, 1, 1, 0, 0, 0), unitsRiver)
# tendRiver = netCDF.date2num(datetime(year+1, 1, 1, 0, 0, 0), unitsRiver)
# start and end indices in time for river discharge
itstartRiver = bisect.bisect_left(datesRiver, datetime(year, 1, 1, 0, 0, 0))
itendRiver = bisect.bisect_left(datesRiver, datetime(year+1, 1, 1, 0, 0, 0))
# ticks for months on river discharge
mticks = [bisect.bisect_left(datesRiver, monthdate) for monthdate in np.asarray(monthdates)]
mticknames = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
##

# Loop through times that simulations were started
for plotdate in plotdates:

    # Set up before plotting
    itmodel = bisect.bisect_left(datesModel, plotdate) # index for model output at this time
    itwind = bisect.bisect_left(datesWind, plotdate) # index for wind at this time
    itriver = bisect.bisect_left(datesRiver, plotdate) # index for river at this time

    figname = 'figures/' + str(year) + '/' + datesModel[itmodel].isoformat()[0:13] + '.png'

    # # Don't redo plot
    # if os.path.exists(figname):
    #     continue

    # Set up plot
    fig = plt.figure(figsize=(10.1, 4.2), dpi=100)
    ax = fig.add_axes([0.04, 0.04, 0.97, 0.88])
    ax.set_frame_on(False) # kind of like it without the box
    tracpy.plotting.background(grid=grid, ax=ax, outline=False, mers=np.arange(-97, -88), merslabels=[0, 0, 1, 0])

    # Date
    date = datesModel[itmodel].strftime('%Y %b %02d %H:%M')
    # greyfont = plt.matplotlib.font_manager.FontProperties() # grab default font properties
    # greyfont.set_color('')
    ax.text(0.77, 0.185, date, fontsize=14, color='0.2', transform=ax.transAxes, 
                bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

    # # PONG
    # ax.text(0.8, 0.95, 'pong.tamu.edu', fontsize=14, transform=ax.transAxes, color='0.3')

    # Plot surface salinity
    # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
    salt = np.squeeze(m.variables['salt'][itmodel,-1,1:-1,1:-1])
    mappable = ax.pcolormesh(xpsi, ypsi, salt, cmap=cmap, vmin=0, vmax=36)
    # # Plot Sabine too, which gets covered by the basemap
    # sabmask = ~salt[172:189,332:341].mask.astype(bool)
    # sabmask[3,2] = False
    # sabmask[3,3] = False
    # sabmask[4,1] = False
    # sabmask[4,2] = False
    # sabmask[5,0] = False
    # sabmask[5,1] = False
    # sabmask[6,0] = False
    # sabmask[4,7] = False
    # sabmask[8:14,4] = False
    # sabmask[15,7] = False
    # sabmask[16,7] = False
    # sabmask[3:5,5:7] = False
    # salt[172:189,332:341] = np.ma.masked_where(~sabmask,salt[172:189,332:341])
    # ax.pcolormesh(xpsi[172:189,332:341], ypsi[172:189,332:341], salt[172:189,332:341], cmap=cmap, vmin=0, vmax=36, zorder=2)

    # Mississippi river discharge rate
    axr = fig.add_axes([0.5, 0.05, 0.48, .13])
    axr.set_frame_on(False) # kind of like it without the box
    axr.fill_between(tRiver[itstartRiver:itriver+1], Q[itstartRiver:itriver+1], alpha=0.5, facecolor='0.4', edgecolor='0.4')
    axr.plot(tRiver[itstartRiver:itendRiver+1], Q[itstartRiver:itendRiver+1], '-', color='0.4')
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [5, 5], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [10000, 10000], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [20000, 20000], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [30000, 30000], '-', color='0.6', lw=0.5, alpha=0.5)
    # labels
    axr.text(tRiver[mticks[-3]]+16.5, 5, '0', fontsize=8, color='0.4')
    axr.text(tRiver[mticks[-3]]+16.5, 10000, '10000', fontsize=8, color='0.4')
    axr.text(tRiver[mticks[-3]]+16.5, 20000, '20000', fontsize=8, color='0.4')
    axr.text(tRiver[mticks[-3]]+15, 30000, r'30000 m$^3$s$^{-1}$', fontsize=8, color='0.4')
    axr.text(tRiver[mticks[-7]]+15, 30000, 'Mississippi discharge', fontsize=8, color='0.4')
    # ticks
    # axr.get_xaxis().set_ticklabels([])
    # axr.xaxis.set_ticks_position('bottom')
    # add ticks for each month
    # axr.set_xticks(tRiver[mticks])
    # axr.tick_params('x', width=1.5, color='0.4') # make ticks bigger
    axr.get_yaxis().set_visible(False)
    axr.get_xaxis().set_visible(False)
    # label month ticks
    for i in xrange(len(mticks)):
        axr.text(tRiver[mticks[i]], 2500, mticknames[i], fontsize=9, color='0.2')

    # Wind over the domain
    Uwind = w.variables['Uwind'][itwind,:,:]
    Vwind = w.variables['Vwind'][itwind,:,:]
    Q = ax.quiver(xr[::wdy,::wdx], yr[::wdy,::wdx], Uwind[::wdy,::wdx], Vwind[::wdy,::wdx], color='0.4', alpha=0.5)
    qk = ax.quiverkey(Q, 0.15, 0.65, 5, r'5 m$\cdot$s$^{-1}$', labelcolor='0.2')

    # Colorbar in upper left corner
    cax = fig.add_axes([0.09, 0.85, 0.35, 0.03]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label(r'Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=14, color='0.2')
    cb.ax.tick_params(labelsize=14, length=2, color='0.2', labelcolor='0.2') 
    cb.set_ticks(ticks)
    # change colorbar tick color http://stackoverflow.com/questions/9662995/matplotlib-change-title-and-colorbar-text-and-tick-colors
    cbtick = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cbtick, color='0.2')
    pdb.set_trace()

    plt.savefig(figname)
    plt.close(fig)

# To make movie: ffmpeg -r 10 -pattern_type glob -i '2008-*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 movie.mp4
