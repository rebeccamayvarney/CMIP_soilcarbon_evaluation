#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thursday 15th July 2021
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)

"""

#%%

# Analysis imports
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import iris
import iris.coord_categorisation

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_cmip5, open_netCDF, regrid_model
from rmv_cmip_analysis import select_time, time_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%%
# OBSERVATIONAL DATA

#%%
# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
#regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points

#%%
# Observational soil carbon
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
# hwsd
hwsd_file = Dataset('/home/links/rmv203/DATA/obs_datasets/HWSD_soilc_0.5x0.5.nc')
hwsd_data = hwsd_file.variables['soilc'][:]
# merging the soil carbon observational datasets
merged_hwsd_ncscd = np.copy(hwsd_data)
merged_hwsd_ncscd[ncscd_data[:] > 0.] = ncscd_data[ncscd_data[:] > 0.]
merged_hwsd_ncscd_masked = ma.masked_where(np.logical_or(merged_hwsd_ncscd < 0, merged_hwsd_ncscd > 998), merged_hwsd_ncscd)
#
observational_cSoil = merged_hwsd_ncscd_masked.copy()

#
cSoil_observations = observational_cSoil.copy()



#%%
# CMIP6
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'UKESM1-0-LL', 'CESM2', 'NorESM2-LM']
n_models = len(cmip6_models)

# empty array
spatial_ensemble_cSoil = ma.zeros((n_models, len(n_lat), len(n_lon)))
spatial_ensemble_cSoil_nodeep = ma.zeros((9, len(n_lat), len(n_lon)))
spatial_ensemble_npp = ma.zeros((n_models, len(n_lat), len(n_lon)))
spatial_ensemble_npp_nodeep = ma.zeros((9, len(n_lat), len(n_lon)))
spatial_ensemble_tau = ma.zeros((n_models, len(n_lat), len(n_lon)))
spatial_ensemble_tau_nodeep = ma.zeros((9, len(n_lat), len(n_lon)))

# for loop for each CMIP5 model
for model_i in range(n_models):
    model = cmip6_models[model_i] # seleting the models
    print(model)

    # model cSoil
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

    if model_i > 8:
        spatial_ensemble_cSoil[model_i, :, :] = soil_carbon_data
    else:
        spatial_ensemble_cSoil[model_i, :, :] = soil_carbon_data
        spatial_ensemble_cSoil_nodeep[model_i, :, :] = soil_carbon_data


# MAP ENSEMBLE MEAN
cSoil_ensemble_mean_cmip6 = np.nanmean(spatial_ensemble_cSoil, axis=0)
cSoil_ensemble_mean_cmip6_nodeep = np.nanmean(spatial_ensemble_cSoil_nodeep, axis=0)

# Standard deviations
cSoil_ensemble_std_cmip6 = ma.std(spatial_ensemble_cSoil, axis=0)

#
## covariance
co_var_cSoil_cmip6 = cSoil_ensemble_std_cmip6/cSoil_ensemble_mean_cmip6

# hatching the data
masked_cvar_cSoil_cmip6_1 = np.ma.masked_less(co_var_cSoil_cmip6, 1)
masked_cvar_cSoil_cmip6_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip6 < 5, co_var_cSoil_cmip6 < 0.75), co_var_cSoil_cmip6)


#%% PLOTTING MAPS

# Set up subplot figure
fig = plt.figure(1, figsize=(33,30))
gs = gspec.GridSpec(2, 2, figure=fig, width_ratios=[1, 0.1], hspace=0.1, wspace=0.1)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True 
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.color'] = 'k'
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':50,
    'xtick.labelsize':50,
    'ytick.labelsize':50,
    'font.size':50,
}
plt.rcParams.update(params)


ax = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
gl.xlabel_style = {'size':32, 'color': 'k'}
gl.ylabel_style = {'size':32, 'color': 'k'}
ax.coastlines()
# set up the x and y coordination
lat = n_lat
lon = n_lon
x, y = np.meshgrid(lon, lat)
print(np.min(cSoil_ensemble_mean_cmip6), np.max(cSoil_ensemble_mean_cmip6))
line = np.arange(0, 50, 1)
diff = plt.contourf(x, y, cSoil_ensemble_mean_cmip6, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_ylim(-70,90)
ax.set_title('CMIP6', fontweight='bold')

ax = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
gl.xlabel_style = {'size':32, 'color': 'k'}
gl.ylabel_style = {'size':32, 'color': 'k'}
ax.coastlines()
print(np.min(cSoil_observations), np.max(cSoil_observations))
line = np.arange(0, 50, 1)
diff = plt.contourf(x, y, cSoil_observations, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_ylim(-70,90)
ax.set_title('Observations', fontweight='bold')

ax=fig.add_subplot(gs[:,1])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='vertical').set_label(r'Soil Carbon (kg C m$^{-2}$)', labelpad=24)



#%%
# Save figure
fig.savefig('figures/varneyetal2022_keyfigure_v3', bbox_inches='tight')#,dpi = 300)
plt.close()