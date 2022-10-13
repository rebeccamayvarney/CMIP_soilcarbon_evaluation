#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Friday 14 May 2021
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)

"""

#%%

# Analysis imports
import numpy as np
import numpy.ma as ma
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_model, open_netCDF, regrid_model
from rmv_cmip_analysis import select_time, time_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%%
fig = plt.figure(1, figsize=(60,38))
gs = gspec.GridSpec(4, 4, figure=fig, width_ratios=[1, 1, 1, 0.05], hspace=0.3, wspace=0.3)
n = 12
column_1 = 0
row_1 = 0
n_columns_1 = 4
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True 
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':58,
    'xtick.labelsize':58,
    'ytick.labelsize':58,
    'font.size':58,
}
plt.rcParams.update(params)

# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
#regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points

#%%
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

# for loop for each CMIP5 model
for model_i, a, in zip(range(n_models), range(n)):
    model = cmip6_models[model_i] # seleting the models
    print(model)

    # model land fraction
    landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)


    # model Rh
    rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/rh_Lmon_'+model+'_historical*', model)
    rh_cube = open_netCDF(rh_cube)        
    rh_cube = select_time(rh_cube, 1995, 2005)
    rh_cube = time_average(rh_cube)
    rh_cube = regrid_model(rh_cube, regrid_cube)
    rh_data = rh_cube.data*86400.*360.
    rh_data = ma.masked_where(np.logical_or(rh_data < 1e-10, rh_data > 1e8), rh_data)

    # model npp
    npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_historical*', model)
    npp_cube = open_netCDF(npp_cube)        
    npp_cube = select_time(npp_cube, 1995, 2005)
    npp_cube = time_average(npp_cube)
    npp_cube = regrid_model(npp_cube, regrid_cube)
    npp_data = npp_cube.data*86400.*360.
    npp_data = ma.masked_where(np.logical_or(npp_data < 1e-08, npp_data > 1e8), npp_data)
    

    # plotting
    print(row_1, column_1)
    ax = fig.add_subplot(gs[row_1, column_1], projection=ccrs.PlateCarree())
    #  add lat/lon grid lines to the figure
    gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
    gl.yformatter=LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
    gl.ylabels_left=True
    gl.xformatter=LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator([-90, 0, 90])
    gl.xlabels_bottom=True
    ax.coastlines()
    # set up the x and y coordination
    lat = n_lat
    lon = n_lon
    x, y = np.meshgrid(lon, lat)

    print(np.min(rh_data), np.max(rh_data))
    line = np.arange(0, 2.2, 0.01)
    diff = plt.contourf(x, y, npp_data, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
    ax.set_title(model)
    ax.set_ylim(-70,90)

    #%%
    # increase row and column 
    row_1 += 1 
    if row_1==4: 
        column_1 += 1
        row_1 = 0


#%%
# plot colourbar
ax=fig.add_subplot(gs[:,3])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='vertical').set_label(r'NPP (kg C m$^{-2}$ yr$^{-1}$)')


#%%
# Save figure
fig.savefig('figA08', bbox_inches='tight')
plt.close()