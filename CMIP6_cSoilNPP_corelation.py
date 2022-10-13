#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thursday 27 May 2021
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
from rmv_cmip_analysis import combine_netCDF_cmip6, open_netCDF, regrid_model
from rmv_cmip_analysis import select_time, time_average


#%%
# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
#regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
regrid_modelcube = regrid_cube.copy()


#%%
# # loading GSDE obs cSoil
# gsde_cube = iris.load_cube('/home/links/DATA/obs_datasets/rmv203/soilC/GSDE_orgC.nc')
# n_lat = gsde_cube.coord('latitude').points
# n_lon = gsde_cube.coord('longitude').points
# gsde_data = gsde_cube.data
# gsde_data = np.ma.masked_where(gsde_data>1e35, gsde_data)
# gsde_data = np.ma.masked_where(gsde_data<0, gsde_data)
# gsde_data = np.sum(gsde_data, axis=0)
# loading WISE30sec obs cSoil
WISE30sec_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/WISE30sec_Gustaf_orgC.nc')
n_lat = WISE30sec_cube.coord('latitude').points
n_lon = WISE30sec_cube.coord('longitude').points
WISE30sec_data = WISE30sec_cube.data
WISE30sec_data = np.ma.masked_where(WISE30sec_data>1e35, WISE30sec_data)
WISE30sec_data = np.sum(WISE30sec_data, axis=0)
WISE30sec_data = np.flip(WISE30sec_data, axis=0)
print(WISE30sec_data.shape, np.min(WISE30sec_data), np.max(WISE30sec_data))
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
flat_obscSoil_data = merged_hwsd_ncscd_masked.flatten()


# Observational NPP
npp_file_obs = Dataset('/home/links/rmv203/DATA/obs_datasets/MOD17A3_Science_NPP_mean_00_14_regridhalfdegree.nc')
npp_data_obs = npp_file_obs.variables['npp'][:]*1e-3
flat_obsNPP_data = npp_data_obs.flatten()

# correlation
r_coeffient = ma.corrcoef(flat_obscSoil_data, flat_obsNPP_data)
print('r-coefficent cSoil & NPP: OBSERVATIONs', r_coeffient)


#%%
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

correlation_array = ma.zeros([len(cmip6_models)+1])
r_coeffient = r_coeffient[0]
r_coeffient = r_coeffient[1]
correlation_array[0] = r_coeffient


# for loop for each CMIP5 model
for model_i in range(n_models):
    model = cmip6_models[model_i] # seleting the models
    print(model)

    model_j = model_i+1

    # model cSoil
    if model=='CESM2' or model=='NorESM2-LM':
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoilAbove1m_Emon_'+model+'_historical*', model)
    else:
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_historical*', model)
    cSoil_cube = open_netCDF(cSoil_cube)
    n_lat = cSoil_cube.coord('latitude').points
    n_lon = cSoil_cube.coord('longitude').points
    print(len(n_lat), len(n_lon))
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

    # model npp
    npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_historical*', model)
    npp_cube = open_netCDF(npp_cube)        
    npp_cube = select_time(npp_cube, 1995, 2005)
    npp_cube = time_average(npp_cube)
    npp_cube = regrid_model(npp_cube, regrid_cube)
    npp_data = npp_cube.data*86400.*360.


    #%% Spatial Correlation
    model_cSoil_flatten = soil_carbon_data.flatten()
    model_NPP_flatten = npp_data.flatten()
    r_coeffient = ma.corrcoef(model_cSoil_flatten, model_NPP_flatten)
    print('r-coefficent cSoil & NPP:', model, r_coeffient)
    r_coeffient = r_coeffient[0]
    r_coeffient = r_coeffient[1]
    correlation_array[model_j] = r_coeffient


np.save("saved_variables/cSoilNPP_correlations_cmip6.npy", correlation_array.data)