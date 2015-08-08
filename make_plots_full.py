'''
Make salinity plots for movies of the full domain.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
from datetime import datetime, timedelta
import glob
from cmPong import cmPong
from matplotlib.mlab import find
import bisect
from matplotlib import delaunay
import op
import cmocean
import matplotlib.patches as patches


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


year = 2014



def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid = tracpy.inout.readgrid(grid_filename, usebasemap=True, llcrnrlat=22.85, llcrnrlon=-97.9, urcrnrlat=30.5)
# actually using psi grid here despite the name
xpsi = np.asanyarray(grid['xpsi'].T, order='C')
ypsi = np.asanyarray(grid['ypsi'].T, order='C')
xr = np.asanyarray(grid['xr'].T, order='C')
yr = np.asanyarray(grid['yr'].T, order='C')

# to rotate wind vectors
m = netCDF.Dataset(loc)
anglev = m.variables['angle'][:]
# w = m  # this is to use the surface wind stress for wind arrows instead of wind velocities from the forcing files

## Model output ##
# if year <= 2013:
#     currents_filenames = np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_????.nc'))
# elif year == 2014:
#     currents_filenames = np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_??.nc'))
# m = netCDF.MFDataset(currents_filenames)

# Model time period to use
units = m.variables['ocean_time'].units
starttime = netCDF.date2num(datetime(year, 1, 1, 4, 0, 0), units)
if year==2014:
    endtime = netCDF.date2num(datetime(year, 9, 30, 20, 0, 0), units)
else:
    endtime = netCDF.date2num(datetime(year+1, 1, 1, 4, 0, 0), units)
dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
datesModel = netCDF.num2date(m.variables['ocean_time'][:], units)

plotdates = netCDF.num2date(ts, units)
if year == 2014:
    monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,10)]
else:
    monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,13)]

# if not os.path.exists('figures/' + str(year)):
#     os.makedirs('figures/' + str(year))

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity
cmap = cmPong.salinity(cmocean.cm.salt, levels)
# cmap = cmPong.salinity('YlGnBu_r', levels)
ilevels = [0,1,2,3,4,5,8] # which levels to label
ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
##

## Wind forcing ##

# There are multiple file locations
if year <= 2012:
    w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
elif year == 2013:
    w = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
elif year == 2014:
    w = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2014.nc')

# w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
# Wind time period to use
unitsWind = (w.variables['time'].units).replace('/','-')
datesWind = netCDF.num2date(w.variables['time'][:], unitsWind)
# datesWind = datesModel
wdx = 22; wdy = 30 # in indices
##

## River forcing ##
r1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012.nc') # use for through 2011
r2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012_2014.nc') # use for 2012-2014
# River timing
tr1 = r1.variables['river_time']
tunitsr1 = tr1.units
# interpolate times for this data file since at the 12 hours mark instead of beginning of the day
tr1 = op.resize(tr1, 0)
datesr1 = netCDF.num2date(tr1[:], tunitsr1)
tr2 = r2.variables['river_time']
datesr2 = netCDF.num2date(tr2[:], tr2.units)
# all of river input
Q1 = np.abs(r1.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
# interpolate this like for time
Q1 = op.resize(Q1, 0)
Q2 = np.abs(r2.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
# Combine river info into one dataset
iend1 = find(datesr1<datetime(2012,1,1,0,0,0))[-1] # ending index for file 1
tRiver = np.concatenate((tr1[:iend1], tr2[:]), axis=0)
datesRiver = np.concatenate((datesr1[:iend1], datesr2))
R = np.concatenate((Q1[:iend1], Q2))
r1.close(); r2.close()
# start and end indices in time for river discharge
itstartRiver = bisect.bisect_left(datesRiver, datetime(year, 1, 1, 0, 0, 0))
if year == 2014:
    itendRiver = bisect.bisect_left(datesRiver, datetime(year, 9, 30, 20, 0, 0))
else:
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

    figname = 'figures/full/' + datesModel[itmodel].isoformat()[0:13] + '.png'

    # Don't redo plot
    if os.path.exists(figname):
        continue

    # Set up plot
    fig = plt.figure(figsize=(10.1, 8.4), dpi=100)
    ax = fig.add_axes([0.06, 0.00, 0.93, 0.97])
    ax.set_frame_on(False) # kind of like it without the box
    tracpy.plotting.background(grid=grid, ax=ax, outline=False, mers=np.arange(-97, -87), merslabels=[0, 0, 1, 0], pars=np.arange(23, 32))

    # Date
    date = datesModel[itmodel].strftime('%Y %b %02d %H:%M')
    # greyfont = plt.matplotlib.font_manager.FontProperties() # grab default font properties
    # greyfont.set_color('')
    ax.text(0.77, 0.35, date, fontsize=15, color='0.2', transform=ax.transAxes, 
                bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

    # PONG
    ax.text(0.575, 0.955, 'pong.tamu.edu', fontsize=12, transform=ax.transAxes, color='0.3')

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
    axr = fig.add_axes([0.35, 0.05, 0.6, .2])
    axr.set_frame_on(False) # kind of like it without the box
    # make background rectangle so lines don't overlap
    axr.fill_between(tRiver[itstartRiver:itriver+1], R[itstartRiver:itriver+1], alpha=0.5, facecolor='0.4', edgecolor='0.4')
    axr.plot(tRiver[itstartRiver:itriver], R[itstartRiver:itriver], '-', color='0.4')
    axr.plot(tRiver[itriver:itendRiver+1], R[itriver:itendRiver+1], '-', color='0.4', alpha=0.3)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [5, 5], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [10000, 10000], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [20000, 20000], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [30000, 30000], '-', color='0.6', lw=0.5, alpha=0.5)
    # this makes sure alignment stays consistent in different years
    axr.autoscale(axis='x', tight=True) 
    axr.set_ylim(-1000,45000) 
    # labels
    axr.text(tRiver[mticks[-3]]+16.5, 5, '0', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+16.5, 10000, '10', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+16.5, 20000, '20', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+15, 30000, r'$30\times10^3$ m$^3$s$^{-1}$', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-7]]+15+50, 30000, 'Mississippi discharge', fontsize=10, color='0.2')
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
    axr.add_patch( patches.Rectangle( (0.3, 0.01), 0.7, 0.2, transform=ax.transAxes, color='white', zorder=1))    

    # Wind over the domain
    Uwind = w.variables['Uwind'][itwind,:,:]
    Vwind = w.variables['Vwind'][itwind,:,:]
    Q = ax.quiver(xr[wdy::wdy,wdx::wdx], yr[wdy::wdy,wdx::wdx], Uwind[wdy::wdy,wdx::wdx], Vwind[wdy::wdy,wdx::wdx], 
            color='k', alpha=0.1, scale=400, pivot='middle', headlength=3, headaxislength=2.8)
    qk = ax.quiverkey(Q, 0.18, 0.8, 10, r'10 m$\cdot$s$^{-1}$ wind', labelcolor='0.2', fontproperties={'size': '10'})

    # sustr = w.variables['sustr'][itwind,:,:]
    # svstr = w.variables['svstr'][itwind,:,:]
    # sustr, svstr = rot2d(op.resize(sustr,1)[1:-1,:], op.resize(svstr,0)[:,1:-1], anglev[1:-1, 1:-1])
    # Q = ax.quiver(xr[wdy+1::wdy,wdx+1::wdx], yr[wdy+1::wdy,wdx+1::wdx], sustr[wdy::wdy,wdx::wdx], svstr[wdy::wdy,wdx::wdx], 
    #         color='k', alpha=0.1, scale=4, pivot='middle', headlength=3, headaxislength=2.8)
    # qk = ax.quiverkey(Q, 0.18, 0.65, 0.1, label=r'0.1 N m$^{2}$', labelcolor='0.2', fontproperties={'size': '10'})

    # Colorbar in upper left corner
    cax = fig.add_axes([0.09, 0.915, 0.35, 0.025]) #colorbar axes
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
