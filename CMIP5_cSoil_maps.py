#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tuesday 20th July 2021
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)

"""

#%%

# Analysis imports
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip5, open_netCDF, regrid_model
from rmv_cmip_analysis import select_time, time_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%%
# Set up subplot figure
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
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points

#%%
cmip5_models = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models = len(cmip5_models)

# for loop for each CMIP5 model
for model_i in range(n_models):
    model = cmip5_models[model_i] # seleting the models
    print(model)


    cSoil_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cSoil_Lmon_'+model+'_historical*', 'soil_carbon_content', model)
    cSoil_cube = open_netCDF(cSoil_cube)        
    cSoil_cube = select_time(cSoil_cube, 1990, 2000)
    cSoil_cube = time_average(cSoil_cube)
    cSoil_cube = regrid_model(cSoil_cube, regrid_cube)
    cSoil_data = cSoil_cube.data
    cSoil_data = ma.masked_where(np.logical_or(cSoil_data < 0, cSoil_data > 1e20), cSoil_data)

    if model=='BNU-ESM' or model=='CCSM4' or model=='CESM1-CAM5' or model=='CanESM2' or model=='IPSL-CM5A-LR' or model=='MIROC-ESM' or model=='MPI-ESM-LR' or model=='NorESM1-M':
        # model cLitter
        cLitter_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cLitter_Lmon_'+model+'_historical*', 'litter_carbon_content', model)
        cLitter_cube = open_netCDF(cLitter_cube)        
        cLitter_cube = select_time(cLitter_cube, 1990, 2000)
        cLitter_cube = time_average(cLitter_cube)
        cLitter_cube = regrid_model(cLitter_cube, regrid_cube)
        cLitter_data = cLitter_cube.data
        cLitter_data = ma.masked_where(np.logical_or(cLitter_data < 0, cLitter_data > 1e20), cLitter_data)
        soil_carbon_data = cSoil_data + cLitter_data
    else:
        soil_carbon_data = cSoil_data.copy()


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

    print(np.min(soil_carbon_data), np.max(soil_carbon_data))
    line = np.arange(1, 60, 1)
    diff = plt.contourf(x, y, soil_carbon_data, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
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
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$C_{s}$ (kg C m$^{-2}$)')


#%%
# Save figure
fig.savefig('individual_model_figures/figA07', bbox_inches='tight')
plt.close()