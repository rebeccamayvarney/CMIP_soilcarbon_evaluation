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
# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
#regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points

#%%
# OBSERVATIONAL DATA

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

# heterotrophic respiration (CARDAMOM)
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
# regrid the cube
rh_cube = regrid_model(mean_cube, regrid_cube)
rh_data_regridded = rh_cube.data
rh_data_regridded = rh_data_regridded*1e-3*365
rh_data_regridded = ma.masked_invalid(rh_data_regridded)
# tau calculation
tau_s = observational_cSoil / rh_data_regridded
tau_s_masked = ma.masked_where(np.logical_or(tau_s < 1, tau_s > 1e4), tau_s)


# MODIS Net Primary Production (NPP)
npp_file = Dataset('/home/links/rmv203/DATA/obs_datasets/MOD17A3_Science_NPP_mean_00_14_regridhalfdegree.nc')
npp_data = npp_file.variables['npp'][:]*1e-3
npp_data = np.ma.masked_where(npp_data<=0, npp_data)

#
cSoil_observations = observational_cSoil.copy()
npp_observations = npp_data.copy()
tau_observations = tau_s_masked.copy()



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


    # model npp
    npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_historical*', model)
    npp_cube = open_netCDF(npp_cube)        
    npp_cube = select_time(npp_cube, 1995, 2005)
    npp_cube = time_average(npp_cube)
    npp_cube = regrid_model(npp_cube, regrid_cube)
    npp_data = npp_cube.data*86400.*360.
    npp_data = ma.masked_where(np.logical_or(npp_data < 1e-8, npp_data > 1e8), npp_data)
    print(np.min(npp_data))
    if model_i > 8:
        spatial_ensemble_npp[model_i, :, :] = npp_data
    else:
        spatial_ensemble_npp[model_i, :, :] = npp_data
        spatial_ensemble_npp_nodeep[model_i, :, :] = npp_data


    # model tau_s
    rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/rh_Lmon_'+model+'_historical*', model)
    rh_cube = open_netCDF(rh_cube)        
    rh_cube = select_time(rh_cube, 1995, 2005)
    rh_cube = time_average(rh_cube)
    rh_cube = regrid_model(rh_cube, regrid_cube)
    rh_data = rh_cube.data*86400.*360.
    rh_data = ma.masked_where(np.logical_or(rh_data < 0, rh_data > 1e8), rh_data)
    tau_s = soil_carbon_data/rh_data
    tau_s_masked_data_historical = ma.masked_where(np.logical_or(tau_s < 1, tau_s > 1e4), tau_s)
    if model_i > 8:
        spatial_ensemble_tau[model_i, :, :] = tau_s_masked_data_historical
    else:
        spatial_ensemble_tau[model_i, :, :] = tau_s_masked_data_historical
        spatial_ensemble_tau_nodeep[model_i, :, :] = tau_s_masked_data_historical


# MAP ENSEMBLE MEAN
cSoil_ensemble_mean_cmip6 = np.nanmean(spatial_ensemble_cSoil, axis=0)
cSoil_ensemble_mean_cmip6_nodeep = np.nanmean(spatial_ensemble_cSoil_nodeep, axis=0)
npp_ensemble_mean_cmip6 = np.nanmean(spatial_ensemble_npp, axis=0)
npp_ensemble_mean_cmip6_nodeep = np.nanmean(spatial_ensemble_npp_nodeep, axis=0)
tau_ensemble_mean_cmip6 = np.nanmean(spatial_ensemble_tau, axis=0)
tau_ensemble_mean_cmip6_nodeep = np.nanmean(spatial_ensemble_tau_nodeep, axis=0)

# Standard deviations
cSoil_ensemble_std_cmip6 = ma.std(spatial_ensemble_cSoil, axis=0)
npp_ensemble_std_cmip6 = ma.std(spatial_ensemble_npp, axis=0)
tau_ensemble_std_cmip6 = ma.std(spatial_ensemble_tau, axis=0)


#%%
# CMIP5
cmip5_models = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models_cmip5 = len(cmip5_models)

# empty array
spatial_ensemble_cSoil_cmip5 = ma.zeros((n_models_cmip5, len(n_lat), len(n_lon)))
spatial_ensemble_npp_cmip5 = ma.zeros((n_models_cmip5, len(n_lat), len(n_lon)))
spatial_ensemble_tau_cmip5 = ma.zeros((n_models_cmip5, len(n_lat), len(n_lon)))

# for loop for each CMIP5 model
for model_j in range(n_models_cmip5):
    model = cmip5_models[model_j] # seleting the models
    print(model)

    # model cSoil
    cSoil_cube_cmip5 = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cSoil_Lmon_'+model+'_historical*', 'soil_carbon_content', model)
    cSoil_cube_cmip5 = open_netCDF(cSoil_cube_cmip5)        
    cSoil_cube_cmip5 = select_time(cSoil_cube_cmip5, 1950, 2000)
    cSoil_cube_cmip5 = time_average(cSoil_cube_cmip5)
    cSoil_cube_cmip5 = regrid_model(cSoil_cube_cmip5, regrid_cube)
    cSoil_data_cmip5 = cSoil_cube_cmip5.data
    cSoil_data_cmip5 = ma.masked_where(np.logical_or(cSoil_data_cmip5 < 1e-8, cSoil_data_cmip5 > 1e8), cSoil_data_cmip5)

    if model=='BNU-ESM' or model=='CCSM4' or model=='CESM1-CAM5' or model=='CanESM2' or model=='IPSL-CM5A-LR' or model=='MIROC-ESM' or model=='MPI-ESM-LR' or model=='NorESM1-M':
        # model cLitter
        cLitter_cube_cmip5 = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cLitter_Lmon_'+model+'_historical*', 'litter_carbon_content', model)
        cLitter_cube_cmip5 = open_netCDF(cLitter_cube_cmip5)        
        cLitter_cube_cmip5 = select_time(cLitter_cube_cmip5, 1950, 2000)
        cLitter_cube_cmip5 = time_average(cLitter_cube_cmip5)
        cLitter_cube_cmip5 = regrid_model(cLitter_cube_cmip5, regrid_cube)
        cLitter_data_cmip5 = cLitter_cube_cmip5.data
        cLitter_data_cmip5 = ma.masked_where(np.logical_or(cLitter_data_cmip5 < 1e-8, cLitter_data_cmip5 > 1e8), cLitter_data_cmip5)
        soil_carbon_data_cmip5 = cSoil_data_cmip5 + cLitter_data_cmip5
    else:
        soil_carbon_data_cmip5 = cSoil_data_cmip5.copy()

    spatial_ensemble_cSoil_cmip5[model_j, :, :] = soil_carbon_data_cmip5


    # model npp
    npp_cube_cmip5 = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/npp_Lmon_'+model+'_historical*', 'net_primary_productivity_of_carbon', model)
    npp_cube_cmip5 = open_netCDF(npp_cube_cmip5)        
    npp_cube_cmip5 = select_time(npp_cube_cmip5, 1995, 2005)
    npp_cube_cmip5 = time_average(npp_cube_cmip5)
    npp_cube_cmip5 = regrid_model(npp_cube_cmip5, regrid_cube)
    npp_data_cmip5 = npp_cube_cmip5.data*86400.*360.
    npp_data_cmip5 = ma.masked_where(np.logical_or(npp_data_cmip5 <= 0, npp_data_cmip5 > 1e8), npp_data_cmip5)
    spatial_ensemble_npp_cmip5[model_j, :, :] = npp_data_cmip5


    # model tau_s
    rh_cube_cmip5 = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/rh_Lmon_'+model+'_historical*', 'heterotrophic_respiration_carbon_flux', model)
    rh_cube_cmip5 = open_netCDF(rh_cube_cmip5)        
    rh_cube_cmip5 = select_time(rh_cube_cmip5, 1995, 2005)
    rh_cube_cmip5 = time_average(rh_cube_cmip5)
    rh_cube_cmip5 = regrid_model(rh_cube_cmip5, regrid_cube)
    rh_data_cmip5 = rh_cube_cmip5.data*86400.*360.
    tau_s_cmip5 = soil_carbon_data_cmip5/rh_data_cmip5
    tau_s_masked_data_historical_cmip5 = ma.masked_where(np.logical_or(tau_s_cmip5 < 1, tau_s_cmip5 > 1e4), tau_s_cmip5)
    spatial_ensemble_tau_cmip5[model_j, :, :] = tau_s_masked_data_historical_cmip5


# MAP ENSEMBLE MEAN
cSoil_ensemble_mean_cmip5 = np.nanmean(spatial_ensemble_cSoil_cmip5, axis=0)
npp_ensemble_mean_cmip5 = np.nanmean(spatial_ensemble_npp_cmip5, axis=0)
tau_ensemble_mean_cmip5 = np.nanmean(spatial_ensemble_tau_cmip5, axis=0)

# Standard deviations
cSoil_ensemble_std_cmip5 = ma.std(spatial_ensemble_cSoil_cmip5, axis=0)
npp_ensemble_std_cmip5 = ma.std(spatial_ensemble_npp_cmip5, axis=0)
tau_ensemble_std_cmip5 = ma.std(spatial_ensemble_tau_cmip5, axis=0)

# covariance
co_var_cSoil_cmip6 = cSoil_ensemble_std_cmip6/cSoil_ensemble_mean_cmip6
co_var_npp_cmip6 = npp_ensemble_std_cmip6/npp_ensemble_mean_cmip6
co_var_tau_cmip6 = tau_ensemble_std_cmip6/tau_ensemble_mean_cmip6

co_var_cSoil_cmip5 = cSoil_ensemble_std_cmip5/cSoil_ensemble_mean_cmip5
co_var_npp_cmip5 = npp_ensemble_std_cmip5/npp_ensemble_mean_cmip5
co_var_tau_cmip5 = tau_ensemble_std_cmip5/tau_ensemble_mean_cmip5

# hatching the data
#masked_cvar_cSoil_cmip6_1 = np.ma.masked_less(co_var_cSoil_cmip6, 1)
#masked_cvar_npp_cmip6_1 = np.ma.masked_less(co_var_npp_cmip6, 1) 
#masked_cvar_tau_cmip6_1 = np.ma.masked_less(co_var_tau_cmip6, 1) 
#
#masked_cvar_cSoil_cmip5_1 = np.ma.masked_less(co_var_cSoil_cmip5, 1)
#masked_cvar_npp_cmip5_1 = np.ma.masked_less(co_var_npp_cmip5, 1) 
#masked_cvar_tau_cmip5_1 = np.ma.masked_less(co_var_tau_cmip5, 1)

#%%
masked_cvar_cSoil_cmip6_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip6 < 5, co_var_cSoil_cmip6 < 0.75), co_var_cSoil_cmip6)
masked_cvar_npp_cmip6_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip6 < 5, co_var_npp_cmip6 < 0.75), co_var_npp_cmip6)
masked_cvar_tau_cmip6_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip6 < 5, co_var_tau_cmip6 < 0.75), co_var_tau_cmip6)

masked_cvar_cSoil_cmip5_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip5 < 5, co_var_cSoil_cmip5 < 0.75), co_var_cSoil_cmip5)
masked_cvar_npp_cmip5_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip5 < 5, co_var_npp_cmip5 < 0.75), co_var_npp_cmip5)
masked_cvar_tau_cmip5_1 = ma.masked_where(np.logical_or(cSoil_ensemble_mean_cmip5 < 5, co_var_tau_cmip5 < 0.75), co_var_tau_cmip5)


#%% PLOTTING MAPS


# Set up subplot figure
fig = plt.figure(1, figsize=(72,36))
gs = gspec.GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 0.1], hspace=0.35, wspace=0.2)
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
    'axes.labelsize':65,
    'xtick.labelsize':65,
    'ytick.labelsize':65,
    'font.size':65,
}
plt.rcParams.update(params)

##
ax = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
ax.text(0.5, 1.15, '(a) Soil carbon',transform=ax.transAxes,ha = 'center', fontweight = 'bold', fontsize=65)
ax.text(-0.45, 0.5, 'CMIP5',transform=ax.transAxes,ha = 'center', fontweight = 'bold', fontsize=65)
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
# set up the x and y coordination
lat = n_lat
lon = n_lon
x, y = np.meshgrid(lon, lat)
print(np.min(cSoil_ensemble_mean_cmip5), np.max(cSoil_ensemble_mean_cmip5))
line = np.arange(0, 50, 1)
diff = plt.contourf(x, y, cSoil_ensemble_mean_cmip5, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
plt.pcolor(x,y,masked_cvar_cSoil_cmip5_1, hatch='////////', alpha=0)
ax.set_ylim(-70,90)


ax = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
ax.text(-0.45, 0.5, 'CMIP6',transform=ax.transAxes,ha = 'center', fontweight = 'bold', fontsize=65)
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
print(np.min(cSoil_ensemble_mean_cmip6), np.max(cSoil_ensemble_mean_cmip6))
line = np.arange(0, 50, 1)
diff = plt.contourf(x, y, cSoil_ensemble_mean_cmip6, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
plt.pcolor(x,y,masked_cvar_cSoil_cmip6_1, hatch='////////', alpha=0)
ax.set_ylim(-70,90)


ax = fig.add_subplot(gs[2, 0], projection=ccrs.PlateCarree())
ax.text(-0.45, 0.5, 'Observations',transform=ax.transAxes,ha = 'center', fontweight = 'bold', fontsize=65)
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
print(np.min(cSoil_observations), np.max(cSoil_observations))
line = np.arange(0, 50, 1)
diff = plt.contourf(x, y, cSoil_observations, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_ylim(-70,90)

ax=fig.add_subplot(gs[3,0])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='horizontal').set_label(r'$C_{s}$ (kg C m$^{-2}$)')


##
ax = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
ax.text(0.5, 1.15, '(b) NPP',transform=ax.transAxes,ha = 'center', fontweight = 'bold', fontsize=65)
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
print(np.min(npp_ensemble_mean_cmip5), np.max(npp_ensemble_mean_cmip5))
line = np.arange(0, 2, 0.1)
diff = plt.contourf(x, y, npp_ensemble_mean_cmip5, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
plt.pcolor(x,y,masked_cvar_npp_cmip5_1, hatch='////////', alpha=0)
ax.set_ylim(-70,90)

ax = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
print(np.min(npp_ensemble_mean_cmip6), np.max(npp_ensemble_mean_cmip6))
line = np.arange(0, 2, 0.1)
diff = plt.contourf(x, y, npp_ensemble_mean_cmip6, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
plt.pcolor(x,y,masked_cvar_npp_cmip6_1, hatch='////////', alpha=0)
ax.set_ylim(-70,90)


ax = fig.add_subplot(gs[2, 1], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
print(np.min(npp_observations), np.max(npp_observations))
line = np.arange(0, 2, 0.1)
diff = plt.contourf(x, y, npp_observations, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_ylim(-70,90)

ax=fig.add_subplot(gs[3,1])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='horizontal').set_label(r'NPP (kg C m$^{-2}$ yr$^{-1}$)')


##
ax = fig.add_subplot(gs[0, 2], projection=ccrs.PlateCarree())
ax.text(0.5, 1.15, '(c) Soil carbon turnover time',transform=ax.transAxes, ha='center', fontweight = 'bold', fontsize=65)
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
x, y = np.meshgrid(lon, lat)
tau_ensemble_mean_cmip5 = ma.log(tau_ensemble_mean_cmip5)
print(np.min(tau_ensemble_mean_cmip5), np.max(tau_ensemble_mean_cmip5))
line = np.arange(1, 8, 0.01)
diff = plt.contourf(x, y, tau_ensemble_mean_cmip5, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
plt.pcolor(x,y,masked_cvar_tau_cmip5_1, hatch='////////', alpha=0)
ax.set_ylim(-70,90)

ax = fig.add_subplot(gs[1, 2], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
tau_ensemble_mean_cmip6 = ma.log(tau_ensemble_mean_cmip6)
print(np.min(tau_ensemble_mean_cmip6), np.max(tau_ensemble_mean_cmip6))
line = np.arange(1, 8, 0.1)
diff = plt.contourf(x, y, tau_ensemble_mean_cmip6, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
plt.pcolor(x,y,masked_cvar_tau_cmip6_1, hatch='////////', alpha=0)
ax.set_ylim(-70,90)


ax = fig.add_subplot(gs[2, 2], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60, 80])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
gl.xlabels_bottom=True
ax.coastlines()
tau_observations_log = ma.log(tau_observations)
print(np.min(tau_observations), np.max(tau_observations))
line = np.arange(1, 8, 0.1)
diff = plt.contourf(x, y, tau_observations_log, line, cmap='YlGn', alpha=0.6, extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_ylim(-70,90)

ax=fig.add_subplot(gs[3,2])
ax=plt.gca()
cbar = fig.colorbar(diff, ax, orientation='horizontal')
cbar.set_label(r'Inferred $\tau_\mathrm{s}$ (yr)')
cbar.ax.set_xticklabels(['3', '5', '7', '12', '20', '33', '55', '90', '150', '400'])


#%%
# Save figure
fig.savefig('figures/fig03', bbox_inches='tight')
plt.close()