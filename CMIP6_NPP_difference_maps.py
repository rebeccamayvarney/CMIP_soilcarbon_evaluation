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
from netCDF4 import Dataset
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_model, open_netCDF, regrid_model, numpy_to_cube
from rmv_cmip_analysis import select_time, time_average, area_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
#regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points

#%%
npp_file_obs = Dataset('/home/links/rmv203/DATA/obs_datasets/MOD17A3_Science_NPP_mean_00_14_regridhalfdegree.nc')
npp_data_obs = npp_file_obs.variables['npp'][:]*1e-3
npp_data_obs = np.ma.masked_where(npp_data_obs<=0, npp_data_obs)

cubes = iris.load('/home/links/rmv203/DATA/obs_datasets/CARDAMOM_2001_2010_FL_RHE.nc')
for cube in cubes:
    if cube.var_name == 'longitude':
        lon = cube
    if cube.var_name == 'latitude':
        lat2 = cube
    if cube.var_name == 'Mean':
        mean_cube = cube
# Takes the latitude and longitude ‘cubes’ and makes them in to coordinates
lat_aux = iris.coords.AuxCoord(lat2.data, standard_name=lat2.name(), units=lat2.units)
lon_aux = iris.coords.AuxCoord(lon.data, standard_name=lon.name(), units=lon.units)
# Add latitude and longitude as coordinates
mean_cube.add_aux_coord(lat_aux, data_dims=(0))
mean_cube.add_aux_coord(lon_aux, data_dims=(1))
iris.util.promote_aux_coord_to_dim_coord(mean_cube, 'latitude')
iris.util.promote_aux_coord_to_dim_coord(mean_cube, 'longitude')
# regridding
rh_cube_obs = regrid_model(mean_cube, regrid_cube)
rh_data_regridded = rh_cube_obs.data
rh_data_regridded = rh_data_regridded*1e-3*365
rh_data_regridded = ma.masked_invalid(rh_data_regridded)
rh_data_regridded = np.ma.masked_where(rh_data_regridded<=0, rh_data_regridded)


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


#%%
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

rmse_array = ma.zeros([len(cmip6_models)])

# for loop for each CMIP5 model
for model_i, a, in zip(range(n_models), range(n)):
    model = cmip6_models[model_i] # seleting the models
    print(model)

    # model land fraction
    landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)

    # model npp
    npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/rh_Lmon_'+model+'_historical*', model)
    npp_cube = open_netCDF(npp_cube)        
    npp_cube = select_time(npp_cube, 1995, 2005)
    npp_cube = time_average(npp_cube)
    npp_cube = regrid_model(npp_cube, regrid_cube)
    npp_data = npp_cube.data*86400.*360.
    npp_data = ma.masked_where(np.logical_or(npp_data < 1e-8, npp_data > 1e8), npp_data)

    # DIFFERENCE
    npp_diff = npp_data - npp_data_obs

    #%% saving root mean square for that model
    rmse_cSoil = np.square(np.subtract(npp_data, rh_data_regridded))
    rmse_cSoil_cube = numpy_to_cube(rmse_cSoil, regrid_cube, 2)
    region = [0, 360, -90,  90]
    rmse_cSoil_cube = area_average(rmse_cSoil_cube, region)
    rmse_cSoil = rmse_cSoil_cube.data
    rmse_cSoil = np.sqrt(rmse_cSoil)
    print(model, rmse_cSoil)
    rmse_array[model_i] = rmse_cSoil

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

    print(np.min(npp_diff), np.max(npp_diff))
    line = np.arange(-1, 1, 0.01)
    diff = plt.contourf(x, y, npp_diff, line, cmap='bwr', extend='both', transform=ccrs.PlateCarree(central_longitude=0))
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
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$\Delta$ NPP (kg C m$^{-2}$ yr$^{-1}$)')


#%%
# Save figure
fig.savefig('fig06', bbox_inches='tight')
plt.close()


#save
np.save("rmse_npp_cmip6_array.npy", rmse_array.data)