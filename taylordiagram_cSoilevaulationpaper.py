#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:28:08 2021
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)

Edited code using: https://gist.github.com/ycopin/3342888
Thank you!!

"""

from taylor_diagram_functionclass import TaylorDiagram
from matplotlib.lines import Line2D

import numpy as np
import matplotlib.pyplot as plt

# Analysis imports
import numpy.ma as ma
from netCDF4 import Dataset
import iris
import iris.coord_categorisation

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_cmip5, open_netCDF, regrid_model, numpy_to_cube
from rmv_cmip_analysis import select_time, time_average, area_average



#%% Observations

# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
#regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points
#
WISE30sec_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/WISE30sec_Gustaf_orgC.nc')
WISE30sec_data = WISE30sec_cube.data
WISE30sec_data = np.ma.masked_where(WISE30sec_data>1e35, WISE30sec_data)
WISE30sec_data = np.sum(WISE30sec_data, axis=0)
WISE30sec_data = np.flip(WISE30sec_data, axis=0)
#
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
flatten_cSoil_obs = observational_cSoil.flatten()

#STD
soil_carbon_cube = numpy_to_cube(observational_cSoil, regrid_cube, 2)
soil_carbon_av = area_average(soil_carbon_cube, [0, 360, -90,  90])
soil_carbon_av_mean = soil_carbon_av.data

spatial_std1 = np.square(np.subtract(observational_cSoil, soil_carbon_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, regrid_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_cSoil_obs = np.sqrt(spatial_std1_data)

print('obs mean and std', soil_carbon_av_mean, spatial_std_cSoil_obs)

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
flatten_tau_obs = tau_s_masked.flatten()

#STD
tau_cube = numpy_to_cube(tau_s_masked, regrid_cube, 2)
tau_av = area_average(tau_cube, [0, 360, -90,  90])
tau_av_mean = tau_av.data

spatial_std1 = np.square(np.subtract(tau_s_masked, tau_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, regrid_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_tau_obs = np.sqrt(spatial_std1_data)


# MODIS Net Primary Production (NPP)
npp_file = Dataset('/home/links/rmv203/DATA/obs_datasets/MOD17A3_Science_NPP_mean_00_14_regridhalfdegree.nc')
npp_data_obs = npp_file.variables['npp'][:]*1e-3
npp_data_obs = np.ma.masked_where(npp_data_obs<=0, npp_data_obs)
flatten_npp_data_obs = npp_data_obs.flatten()

#STD
npp_cube = numpy_to_cube(npp_data_obs, regrid_cube, 2)
npp_av = area_average(npp_cube, [0, 360, -90,  90])
npp_av_mean = npp_av.data

spatial_std1 = np.square(np.subtract(npp_data_obs, npp_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, regrid_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_npp_obs = np.sqrt(spatial_std1_data)


# Reference std
stdrefs = dict(cSoil=spatial_std_cSoil_obs, npp=spatial_std_npp_obs, tau=spatial_std_tau_obs)


#%% CMIP6
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

cmip6_cSoil_dict = []
cmip6_npp_dict = []
cmip6_tau_dict = []

cmip6_cSoil_ensemble = np.zeros((n_models, len(n_lat), len(n_lon)))
cmip6_npp_ensemble = np.zeros((n_models, len(n_lat), len(n_lon)))
cmip6_tau_ensemble = np.zeros((n_models, len(n_lat), len(n_lon)))

for model_i,model in enumerate(cmip6_models):
#    model = cmip6_models[model_i] # seleting the models
    print(model)
    # csoil
    if model=='CESM2' or model=='NorESM2-LM':
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoilAbove1m_Emon_'+model+'_historical*', model)
    else:
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_historical*', model)
    cSoil_cube = open_netCDF(cSoil_cube)        
    cSoil_cube = select_time(cSoil_cube, 1950, 2000)
    cSoil_cube = time_average(cSoil_cube)
    cSoil_cube = regrid_model(cSoil_cube, regrid_cube)
    cSoil_data = cSoil_cube.data
    if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
        # model cLitter
        cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cLitter_Lmon_'+model+'_historical*', model)
        cLitter_cube = open_netCDF(cLitter_cube)        
        cLitter_cube = select_time(cLitter_cube, 1950, 2000)
        cLitter_cube = time_average(cLitter_cube)
        cLitter_cube = regrid_model(cLitter_cube, regrid_cube)
        cLitter_data = cLitter_cube.data
        soil_carbon_data = cSoil_data + cLitter_data
    else:
        soil_carbon_data = cSoil_data.copy()
    
    # 
    cmip6_cSoil_ensemble[model_i, :, :] = soil_carbon_data
    #
    soil_carbon_cube = numpy_to_cube(soil_carbon_data, cSoil_cube, 2)
    soil_carbon_av = area_average(soil_carbon_cube, [0, 360, -90,  90])
    soil_carbon_av_mean = soil_carbon_av.data
    
    spatial_std1 = np.square(np.subtract(soil_carbon_data, soil_carbon_av_mean))
    spatial_std1_cube = numpy_to_cube(spatial_std1, cSoil_cube, 2)
    spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
    spatial_std1_data = spatial_std1_av.data
    spatial_std_cSoil = np.sqrt(spatial_std1_data)
    #
    print('cSoil', model, 'mean', soil_carbon_av_mean, 'std', spatial_std_cSoil)
    #
    soil_carbon_data_flatten = soil_carbon_data.flatten()
    r_coeffient = ma.corrcoef(flatten_cSoil_obs, soil_carbon_data_flatten)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    r_coeffient_cSoil = r_coeffient.copy()
    
    # SOIL CARBON DICTIONARY
    cmip6_cSoil_dict.append([spatial_std_cSoil, r_coeffient_cSoil, model])
    
    
    
    # model NPP
    npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_historical*', model)
    npp_cube = open_netCDF(npp_cube)        
    npp_cube = select_time(npp_cube, 1995, 2005)
    npp_cube = time_average(npp_cube)
    npp_cube = regrid_model(npp_cube, regrid_cube)
    npp_av_mean_cube = area_average(npp_cube, [0, 360, -90,  90])
    npp_av_mean = npp_av_mean_cube.data*86400.*360.
    npp_data = npp_cube.data*86400.*360.
    cmip6_npp_ensemble[model_i, :, :] = npp_data
    
    spatial_std1 = np.square(np.subtract(npp_data, npp_av_mean))
    spatial_std1_cube = numpy_to_cube(spatial_std1, npp_cube, 2)
    spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
    spatial_std1_data = spatial_std1_av.data
    spatial_std_npp = np.sqrt(spatial_std1_data)
    #
    #print('npp', model, 'mean', npp_av_mean, 'std', spatial_std_npp)
    #
    npp_data_flatten = npp_data.flatten()
    r_coeffient = ma.corrcoef(flatten_npp_data_obs, npp_data_flatten)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    r_coeffient_npp = r_coeffient.copy()
    
    # NPP DICTIONARY
    cmip6_npp_dict.append([spatial_std_npp, r_coeffient_npp, model])
    
    
    
    # model tau_s
    rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/rh_Lmon_'+model+'_historical*', model)
    rh_cube = open_netCDF(rh_cube)        
    rh_cube = select_time(rh_cube, 1995, 2005)
    rh_cube = time_average(rh_cube)
    rh_cube = regrid_model(rh_cube, regrid_cube)
    rh_data = rh_cube.data*86400.*360.
    tau_s = soil_carbon_data/rh_data
    tau_s_masked_data_historical = ma.masked_where(np.logical_or(tau_s < 1, tau_s > 1e4), tau_s)
    #
    cmip6_tau_ensemble[model_i, :, :] = tau_s_masked_data_historical
    #
    tau_cube = numpy_to_cube(tau_s_masked_data_historical, cSoil_cube, 2)
    tau_av = area_average(tau_cube, [0, 360, -90,  90])
    tau_av_mean = tau_av.data
    #
    spatial_std1 = np.square(np.subtract(tau_s_masked_data_historical, tau_av_mean))
    spatial_std1_cube = numpy_to_cube(spatial_std1, cSoil_cube, 2)
    spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
    spatial_std1_data = spatial_std1_av.data
    spatial_std_tau = np.sqrt(spatial_std1_data)
    #
    #print('tau', model, 'mean', tau_av_mean, 'std', spatial_std_tau)
    
    tau_data_flatten = tau_s_masked_data_historical.flatten()
    r_coeffient = ma.corrcoef(flatten_tau_obs, tau_data_flatten)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    r_coeffient_tau = r_coeffient.copy()
    
    # tau DICTIONARY
    cmip6_tau_dict.append([spatial_std_tau, r_coeffient_tau, model])
    
#%%
samples = dict(cSoil = cmip6_cSoil_dict, npp = cmip6_npp_dict, tau = cmip6_tau_dict)

#%% ENSEMBLE MEANS

cmip6_cSoil_ensemble_dict = []
cmip6_npp_ensemble_dict = []
cmip6_tau_ensemble_dict = []

# cSoil
cmip6_cSoil_ensemble_mean = np.mean(cmip6_cSoil_ensemble, axis=0)
cSoil_cube = numpy_to_cube(cmip6_cSoil_ensemble_mean, cSoil_cube, 2)
cSoil_av_cube = area_average(cSoil_cube, [0, 360, -90,  90])
cSoil_av_mean = cSoil_av_cube.data
#
spatial_std1 = np.square(np.subtract(cmip6_cSoil_ensemble_mean, cSoil_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, cSoil_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_ensemble_cSoil = np.sqrt(spatial_std1_data)

print('CMIP6 ensemble mean & std', cSoil_av_mean, spatial_std_ensemble_cSoil)

cSoil_data_flatten = cmip6_cSoil_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(flatten_cSoil_obs, cSoil_data_flatten)
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
r_coeffient_cSoil = r_coeffient.copy()

#
cmip6_cSoil_ensemble_dict.append([spatial_std_ensemble_cSoil, r_coeffient_cSoil, 'CMIP6 ensemble'])
    
    
# NPP
cmip6_npp_ensemble_mean = np.mean(cmip6_npp_ensemble, axis=0)
npp_cube = numpy_to_cube(cmip6_npp_ensemble_mean, npp_cube, 2)
npp_av_cube = area_average(npp_cube, [0, 360, -90,  90])
npp_av_mean = npp_av_cube.data
#
spatial_std1 = np.square(np.subtract(cmip6_npp_ensemble_mean, npp_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, npp_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_ensemble_npp = np.sqrt(spatial_std1_data)

npp_data_flatten = cmip6_npp_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(flatten_npp_data_obs, npp_data_flatten)
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
r_coeffient_npp = r_coeffient.copy()

#
cmip6_npp_ensemble_dict.append([spatial_std_ensemble_npp, r_coeffient_npp, 'CMIP6 ensemble'])


# tau
cmip6_tau_ensemble_mean = np.mean(cmip6_tau_ensemble, axis=0)
cmip6_tau_ensemble_mean = ma.masked_where(np.logical_or(cmip6_tau_ensemble_mean < 1, cmip6_tau_ensemble_mean > 1e4), cmip6_tau_ensemble_mean)
tau_cube = numpy_to_cube(cmip6_tau_ensemble_mean, tau_cube, 2)
tau_av_cube = area_average(tau_cube, [0, 360, -90,  90])
tau_av_mean = tau_av_cube.data
#
spatial_std1 = np.square(np.subtract(cmip6_tau_ensemble_mean, tau_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, tau_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_ensemble_tau = np.sqrt(spatial_std1_data)

tau_data_flatten = cmip6_tau_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(flatten_tau_obs, tau_data_flatten)
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
r_coeffient_tau = r_coeffient.copy()

#
cmip6_tau_ensemble_dict.append([spatial_std_ensemble_tau, r_coeffient_tau, 'CMIP6 ensemble'])

#
samples_cmip6ensemble = dict(cSoil = cmip6_cSoil_ensemble_dict, npp = cmip6_npp_ensemble_dict, tau = cmip6_tau_ensemble_dict)



#%% cmip5
cmip5_models = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models = len(cmip5_models)

cmip5_cSoil_dict = []
cmip5_npp_dict = []
cmip5_tau_dict = []

cmip5_cSoil_ensemble = np.zeros((n_models, len(n_lat), len(n_lon)))
cmip5_npp_ensemble = np.zeros((n_models, len(n_lat), len(n_lon)))
cmip5_tau_ensemble = np.zeros((n_models, len(n_lat), len(n_lon)))

for model_i,model in enumerate(cmip5_models):
#    model = cmip5_models[model_i] # seleting the models
    print(model)
    # csoil
    cSoil_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cSoil_Lmon_'+model+'_historical*', 'soil_carbon_content', model)
    cSoil_cube = open_netCDF(cSoil_cube)        
    cSoil_cube = select_time(cSoil_cube, 1950, 2000)
    cSoil_cube = time_average(cSoil_cube)
    cSoil_cube = regrid_model(cSoil_cube, regrid_cube)
    cSoil_data = cSoil_cube.data
    cSoil_data = ma.masked_where(np.logical_or(cSoil_data < 1e-8, cSoil_data > 1e8), cSoil_data)
    if model=='BNU-ESM' or model=='CCSM4' or model=='CESM1-CAM5' or model=='CanESM2' or model=='IPSL-CM5A-LR' or model=='MIROC-ESM' or model=='MPI-ESM-LR' or model=='NorESM1-M':
        # model cLitter
        cLitter_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cLitter_Lmon_'+model+'_historical*', 'litter_carbon_content', model)
        cLitter_cube = open_netCDF(cLitter_cube)        
        cLitter_cube = select_time(cLitter_cube, 1950, 2000)
        cLitter_cube = time_average(cLitter_cube)
        cLitter_cube = regrid_model(cLitter_cube, regrid_cube)
        cLitter_data = cLitter_cube.data
        cLitter_data = ma.masked_where(np.logical_or(cLitter_data < 1e-8, cLitter_data > 1e8), cLitter_data)
        soil_carbon_data = cSoil_data + cLitter_data
    else:
        soil_carbon_data = cSoil_data.copy()
    #
    cmip5_cSoil_ensemble[model_i, :, :] = soil_carbon_data
    #
    soil_carbon_cube = numpy_to_cube(soil_carbon_data, cSoil_cube, 2)
    soil_carbon_av = area_average(soil_carbon_cube, [0, 360, -90,  90])
    soil_carbon_av_mean = soil_carbon_av.data
    
    spatial_std1 = np.square(np.subtract(soil_carbon_data, soil_carbon_av_mean))
    spatial_std1_cube = numpy_to_cube(spatial_std1, cSoil_cube, 2)
    spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
    spatial_std1_data = spatial_std1_av.data
    spatial_std_cSoil = np.sqrt(spatial_std1_data)
    #
    print('cSoil', model, 'mean', soil_carbon_av_mean, 'std', spatial_std_cSoil)
    #
    soil_carbon_data_flatten = soil_carbon_data.flatten()
    r_coeffient = ma.corrcoef(flatten_cSoil_obs, soil_carbon_data_flatten)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    r_coeffient_cSoil = r_coeffient.copy()
    
    # SOIL CARBON DICTIONARY
    cmip5_cSoil_dict.append([spatial_std_cSoil, r_coeffient_cSoil, model])
    
    
    
    # model NPP
    npp_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/npp_Lmon_'+model+'_historical*', 'net_primary_productivity_of_carbon', model)
    npp_cube = open_netCDF(npp_cube)        
    npp_cube = select_time(npp_cube, 1995, 2005)
    npp_cube = time_average(npp_cube)
    npp_cube = regrid_model(npp_cube, regrid_cube)
    npp_av_mean_cube = area_average(npp_cube, [0, 360, -90,  90])
    npp_av_mean = npp_av_mean_cube.data*86400.*360.
    npp_data = npp_cube.data*86400.*360.
    npp_data = ma.masked_where(np.logical_or(npp_data < 0, npp_data > 1e8), npp_data)
    #
    cmip5_npp_ensemble[model_i, :, :] = npp_data
    
    spatial_std1 = np.square(np.subtract(npp_data, npp_av_mean))
    spatial_std1_cube = numpy_to_cube(spatial_std1, npp_cube, 2)
    spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
    spatial_std1_data = spatial_std1_av.data
    spatial_std_npp = np.sqrt(spatial_std1_data)
    #
    #print('npp', model, 'mean', npp_av_mean, 'std', spatial_std_npp)
    #
    npp_data_flatten = npp_data.flatten()
    r_coeffient = ma.corrcoef(flatten_npp_data_obs, npp_data_flatten)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    r_coeffient_npp = r_coeffient.copy()
    
    # NPP DICTIONARY
    cmip5_npp_dict.append([spatial_std_npp, r_coeffient_npp, model])
    
    
    
    # model tau_s
    rh_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/rh_Lmon_'+model+'_historical*', 'heterotrophic_respiration_carbon_flux', model)
    rh_cube = open_netCDF(rh_cube)        
    rh_cube = select_time(rh_cube, 1995, 2005)
    rh_cube = time_average(rh_cube)
    rh_cube = regrid_model(rh_cube, regrid_cube)
    rh_data = rh_cube.data*86400.*360.
    rh_data = ma.masked_where(np.logical_or(rh_data < 0, rh_data > 1e8), rh_data)
    tau_s = soil_carbon_data/rh_data
    tau_s_masked_data_historical = ma.masked_where(np.logical_or(tau_s < 1, tau_s > 1e4), tau_s)
    #
    cmip5_tau_ensemble[model_i, :, :] = tau_s_masked_data_historical
    #
    tau_cube = numpy_to_cube(tau_s_masked_data_historical, cSoil_cube, 2)
    tau_av = area_average(tau_cube, [0, 360, -90,  90])
    tau_av_mean = tau_av.data
    #
    spatial_std1 = np.square(np.subtract(tau_s_masked_data_historical, tau_av_mean))
    spatial_std1_cube = numpy_to_cube(spatial_std1, cSoil_cube, 2)
    spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
    spatial_std1_data = spatial_std1_av.data
    spatial_std_tau = np.sqrt(spatial_std1_data)
    #
    #print('tau', model, 'mean', tau_av_mean, 'std', spatial_std_tau)
    
    tau_data_flatten = tau_s_masked_data_historical.flatten()
    r_coeffient = ma.corrcoef(flatten_tau_obs, tau_data_flatten)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    r_coeffient_tau = r_coeffient.copy()
    
    # tau DICTIONARY
    cmip5_tau_dict.append([spatial_std_tau, r_coeffient_tau, model])
    
#%%
samples_cmip5 = dict(cSoil = cmip5_cSoil_dict, npp = cmip5_npp_dict, tau = cmip5_tau_dict)


#%% ENSEMBLE MEANS

cmip5_cSoil_ensemble_dict = []
cmip5_npp_ensemble_dict = []
cmip5_tau_ensemble_dict = []

# cSoil
cmip5_cSoil_ensemble_mean = np.mean(cmip5_cSoil_ensemble, axis=0)
cmip5_cSoil_ensemble_mean = ma.masked_where(np.logical_or(cmip5_cSoil_ensemble_mean < 1e-8, cmip5_cSoil_ensemble_mean > 1e8), cmip5_cSoil_ensemble_mean)
cSoil_cube = numpy_to_cube(cmip5_cSoil_ensemble_mean, cSoil_cube, 2)
cSoil_av_cube = area_average(cSoil_cube, [0, 360, -90,  90])
cSoil_av_mean = cSoil_av_cube.data
#
spatial_std1 = np.square(np.subtract(cmip5_cSoil_ensemble_mean, cSoil_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, cSoil_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_ensemble_cSoil = np.sqrt(spatial_std1_data)

print('cmip5 ensemble mean & std', cSoil_av_mean, spatial_std_ensemble_cSoil)

cSoil_data_flatten = cmip5_cSoil_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(flatten_cSoil_obs, cSoil_data_flatten)
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
r_coeffient_cSoil = r_coeffient.copy()

#
cmip5_cSoil_ensemble_dict.append([spatial_std_ensemble_cSoil, r_coeffient_cSoil, 'CMIP5 Ensemble'])
    
    
# NPP
cmip5_npp_ensemble_mean = np.mean(cmip5_npp_ensemble, axis=0)
npp_cube = numpy_to_cube(cmip5_npp_ensemble_mean, npp_cube, 2)
npp_av_cube = area_average(npp_cube, [0, 360, -90,  90])
npp_av_mean = npp_av_cube.data
#
spatial_std1 = np.square(np.subtract(cmip5_npp_ensemble_mean, npp_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, npp_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_ensemble_npp = np.sqrt(spatial_std1_data)

npp_data_flatten = cmip5_npp_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(flatten_npp_data_obs, npp_data_flatten)
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
r_coeffient_npp = r_coeffient.copy()

#
cmip5_npp_ensemble_dict.append([spatial_std_ensemble_npp, r_coeffient_npp, 'CMIP5 Ensemble'])


# tau
cmip5_tau_ensemble_mean = np.mean(cmip5_tau_ensemble, axis=0)
cmip5_tau_ensemble_mean = ma.masked_where(np.logical_or(cmip5_tau_ensemble_mean < 1, cmip5_tau_ensemble_mean > 1e4), cmip5_tau_ensemble_mean)
tau_cube = numpy_to_cube(cmip5_tau_ensemble_mean, tau_cube, 2)
tau_av_cube = area_average(tau_cube, [0, 360, -90,  90])
tau_av_mean = tau_av_cube.data
#
spatial_std1 = np.square(np.subtract(cmip5_tau_ensemble_mean, tau_av_mean))
spatial_std1_cube = numpy_to_cube(spatial_std1, tau_cube, 2)
spatial_std1_av = area_average(spatial_std1_cube, [0, 360, -90,  90])
spatial_std1_data = spatial_std1_av.data
spatial_std_ensemble_tau = np.sqrt(spatial_std1_data)

tau_data_flatten = cmip5_tau_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(flatten_tau_obs, tau_data_flatten)
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
r_coeffient_tau = r_coeffient.copy()

#
cmip5_tau_ensemble_dict.append([spatial_std_ensemble_tau, r_coeffient_tau, 'CMIP5 Ensemble'])

#
samples_cmip5ensemble = dict(cSoil = cmip5_cSoil_ensemble_dict, npp = cmip5_npp_ensemble_dict, tau = cmip5_tau_ensemble_dict)


#%% PLOTTING

colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'red', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
colors_cmip5 = ['darkblue', 'dodgerblue', '#80b1d3', 'darkcyan', '#8dd3c7', 'darkseagreen', 'darkgreen', 'olive', 'gold', 'orange']


rects = dict(cSoil=131, npp=132, tau=133

fig = plt.figure(figsize=(52,30))
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':42,
    'xtick.labelsize':42,
    'ytick.labelsize':42,
    'font.size':42,
}
plt.rcParams.update(params)

for variable in ['cSoil','npp','tau']:
    
    if variable == 'cSoil':
        srange_input = (0.0, 2.6)
    if variable == 'tau':
        srange_input = (0.0, 4.2)
    if variable == 'npp':
        srange_input = (0.0, 1.8)

    dia = TaylorDiagram(stdrefs[variable], fig=fig, rect=rects[variable], label='Observations', srange=srange_input, extend = False)


    # Add samples to Taylor diagram
    for i,(stddev,corrcoef,name) in enumerate(samples[variable]):
        dia.add_sample(stddev, corrcoef, marker='o', ms=40, ls='', mfc=colors[i], mec=colors[i], label=name)
        
    for i,(stddev,corrcoef,name) in enumerate(samples_cmip5[variable]):
        dia.add_sample(stddev, corrcoef, marker='X', ms=40, ls='', mfc=colors_cmip5[i], mec=colors_cmip5[i], label=name)
        
    for i,(stddev,corrcoef,name) in enumerate(samples_cmip6ensemble[variable]):
        dia.add_sample(stddev, corrcoef, marker='o', ms=50, ls='', mfc='grey', mec='grey', label=name)
        
    for i,(stddev,corrcoef,name) in enumerate(samples_cmip5ensemble[variable]):
        dia.add_sample(stddev, corrcoef, marker='X', ms=50, ls='', mfc='grey', mec='grey', label=name)

    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5') # 5 levels
    dia.ax.clabel(contours, inline=1, fontsize=30, fmt='%.1f')

    if variable == 'cSoil':
        dia.ax.text(0.5,1.1, '(a) Soil carbon', ha='center', transform=dia.ax.transAxes, fontweight = 'bold')
        
        handels = []
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkblue', label='BNU-ESM')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='dodgerblue', label='CCSM4')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='#80b1d3', label='CanESM2')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkcyan', label='GFDL-ESM2G')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='#8dd3c7', label='GISS-E2-R')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkseagreen', label='HadGEM2-ES')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkgreen', label='IPSL-CM5A-LR')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='olive', label='MIROC-ESM')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='gold', label='MPI-ESM-LR')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='orange', label='NorESM1-LM')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='grey', label='CMIP5 ensemble mean')])
        label = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M', 'CMIP5 ensemble mean']
        leg = plt.legend(handels, label, loc='lower center', ncol=2, bbox_to_anchor=(0.6, -0.8), labelspacing=0.75, title='CMIP5')
        plt.gca().add_artist(leg)
        
        
    if variable == 'tau':
        dia.ax.text(0.5, 1.1, '(c) Soil carbon turnover time', ha='center', transform=dia.ax.transAxes, fontweight = 'bold')
        
        handels3 = []
        handels3.extend([Line2D([0,0],[0,0], linestyle='None', marker='*', markersize=50, color='k', label='Benchmark datasets')])
        label3 = ['Benchmark datasets']
        leg = plt.legend(handels3, label3, title='Observations', loc='lower center', bbox_to_anchor=(0.35, -0.392))
        plt.gca().add_artist(leg)
        
        
    if variable == 'npp':
        dia.ax.text(0.5, 1.1, '(b) NPP', ha='center', transform=dia.ax.transAxes, fontweight = 'bold')
        
        
        handels = []
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='peachpuff', label='ACCESS-ESM1-5')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='#fb8072', label='BCC-CSM2-MR')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='#80b1d3', label='CanESM5')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='dodgerblue', label='CESM2')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='red', label='CNRM-ESM2-1')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='darkcyan', label='GFDL-ESM4')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='darkgreen', label='IPSL-CM6A-LR')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='olive', label='MIROC-ES2L')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='gold', label='MPI-ESM1-2-LR')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='orange', label='NorESM2-LM')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='darkseagreen', label='UKESM1-0-LL')])
        handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='grey', label='CMIP6 ensemble mean')])
        label = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL', 'CMIP6 ensemble mean')
        leg = plt.legend(handels, label, loc='lower center', ncol=2, bbox_to_anchor=(0.6, -0.8), labelspacing=0.75, title='CMIP6')
        plt.gca().add_artist(leg)


fig.tight_layout()
plt.savefig('figures/taylor_diagram_Above1m_v4.png')
plt.close()
