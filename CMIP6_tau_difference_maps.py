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
from rmv_cmip_analysis import combine_netCDF_cmip6, open_netCDF, regrid_model, numpy_to_cube
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
WISE30sec_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/soilC/WISE30sec_Gustaf_orgC.nc')
WISE30sec_data = WISE30sec_cube.data
WISE30sec_data = np.ma.masked_where(WISE30sec_data>1e35, WISE30sec_data)
WISE30sec_data = np.sum(WISE30sec_data, axis=0)
WISE30sec_data = np.flip(WISE30sec_data, axis=0)
#%%
# Observational soil carbon
# ncscd
ncscd_file = Dataset('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
ncscd_data = ncscd_file.variables['soilc'][:]
n_lat_ncscd = ncscd_file.variables['lat'][:]
n_lon_ncscd = ncscd_file.variables['lon'][:]
# hwsd
hwsd_file = Dataset('/home/links/rmv203/DATA/obs_datasets/HWSD_soilc_0.5x0.5.nc')
hwsd_data = hwsd_file.variables['soilc'][:]
# merging the soil carbon observational datasets
merged_hwsd_ncscd = np.copy(hwsd_data)
merged_hwsd_ncscd[ncscd_data[:] > 0.] = ncscd_data[ncscd_data[:] > 0.]
merged_hwsd_ncscd_masked = ma.masked_where(np.logical_or(merged_hwsd_ncscd < 0, merged_hwsd_ncscd > 998), merged_hwsd_ncscd)
observational_cSoil = merged_hwsd_ncscd_masked.copy()

#%%
# CARDAMOM Heterotrophic Respiration (Rh)
cubes = iris.load('/home/links/rmv203/DATA/obs_datasets/CARDAMOM_2001_2010_FL_RHE.nc')
for cube in cubes:
    if cube.var_name == 'longitude':
        lon = cube
    if cube.var_name == 'latitude':
        lat = cube
    if cube.var_name == 'Mean':
        mean_cube = cube
# Takes the latitude and longitude ‘cubes’ and makes them in to coordinates
lat_aux = iris.coords.AuxCoord(lat.data, standard_name=lat.name(), units=lat.units)
lon_aux = iris.coords.AuxCoord(lon.data, standard_name=lon.name(), units=lon.units)
# Add latitude and longitude as coordinates
mean_cube.add_aux_coord(lat_aux, data_dims=(0))
mean_cube.add_aux_coord(lon_aux, data_dims=(1))
iris.util.promote_aux_coord_to_dim_coord(mean_cube, 'latitude')
iris.util.promote_aux_coord_to_dim_coord(mean_cube, 'longitude')
# regridding
rh_cube = regrid_model(mean_cube, regrid_cube)
rh_data_regridded = rh_cube.data
rh_data_regridded = rh_data_regridded*1e-3*365
rh_data_regridded = ma.masked_invalid(rh_data_regridded)
rh_data_regridded = np.ma.masked_where(rh_data_regridded<=0, rh_data_regridded)

# tau calculation
tau_s = observational_cSoil / rh_data_regridded
tau_s_masked = ma.masked_where(np.logical_or(tau_s < 1, tau_s > 1e4), tau_s)


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


#%%
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

rmse_array = ma.zeros([len(cmip6_models)])

# for loop for each CMIP5 model
for model_i, a, in zip(range(n_models), range(n)):
    model = cmip6_models[model_i] # seleting the models
    print(model)

    if model=='CESM2' or model=='NorESM2-LM':
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoilAbove1m_Emon_'+model+'_historical*', model)
    else:
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_historical*', model)
    cSoil_cube = open_netCDF(cSoil_cube)        
    cSoil_cube = select_time(cSoil_cube, 1950, 2000)
    cSoil_cube = time_average(cSoil_cube)
    cSoil_cube = regrid_model(cSoil_cube, regrid_cube)
    cSoil_data = cSoil_cube.data
    cSoil_data = ma.masked_where(np.logical_or(cSoil_data < 1e-8, cSoil_data > 1e8), cSoil_data)

    if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
        # model cLitter
        cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cLitter_Lmon_'+model+'_historical*', model)
        cLitter_cube = open_netCDF(cLitter_cube)        
        cLitter_cube = select_time(cLitter_cube, 1950, 2000)
        cLitter_cube = time_average(cLitter_cube)
        cLitter_cube = regrid_model(cLitter_cube, regrid_cube)
        cLitter_data = cLitter_cube.data
        cLitter_data = ma.masked_where(np.logical_or(cLitter_data < 1e-8, cLitter_data > 1e8), cLitter_data)
        soil_carbon_data = cSoil_data + cLitter_data
    else:
        soil_carbon_data = cSoil_data.copy()


    # model Rh
    rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/rh_Lmon_'+model+'_historical*', model)
    rh_cube = open_netCDF(rh_cube)        
    rh_cube = select_time(rh_cube, 1995, 2005)
    rh_cube = time_average(rh_cube)
    rh_cube = regrid_model(rh_cube, regrid_cube)
    rh_data = rh_cube.data*86400.*360.
    rh_data = ma.masked_where(np.logical_or(rh_data < 1e-8, rh_data > 1e8), rh_data)

    # Calculating Soil Turnover Time (tau_s)
    tau_s_data_historical = soil_carbon_data / rh_data
    tau_s_masked_data_historical = ma.masked_where(np.logical_or(tau_s_data_historical < 1, tau_s_data_historical > 1e4), tau_s_data_historical)

    # DIFFERENCE
    tau_diff = tau_s_masked_data_historical - tau_s_masked

    #%% saving root mean square for that model
    rmse_cSoil = np.square(np.subtract(tau_s_masked_data_historical, tau_s_masked))
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

    print(np.min(tau_diff), np.max(tau_diff))
    line = np.arange(-100, 100, 10)
    diff = plt.contourf(x, y, tau_diff, line, cmap='bwr', extend='both', transform=ccrs.PlateCarree(central_longitude=0))
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
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$\Delta \tau_\mathrm{s}$ (yr)')


#%%
# Save figure
fig.savefig('fig07', bbox_inches='tight')
plt.close()


#save
np.save("rmse_tau_cmip6_array.npy", rmse_array.data)