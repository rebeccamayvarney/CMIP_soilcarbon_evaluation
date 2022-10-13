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
import iris
import iris.coord_categorisation

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, open_netCDF, regrid_model
from rmv_cmip_analysis import select_time, time_average



#%%
# Loading regrid cube
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points


#%%
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
n_models = len(cmip6_models)

# empty array
spatial_ensemble_cSoil = ma.zeros((n_models, len(n_lat), len(n_lon)))
spatial_ensemble_mrso = ma.zeros((n_models, len(n_lat), len(n_lon)))

# for loop for each CMIP5 model
for model_i in range(n_models):
    model = cmip6_models[model_i] # seleting the models
    print(model)


    # model cSoil
    cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_historical*', model)
    cSoil_cube = open_netCDF(cSoil_cube)        
    cSoil_cube = select_time(cSoil_cube, 1950, 2000)
    cSoil_cube = time_average(cSoil_cube)
    cSoil_cube = regrid_model(cSoil_cube, regrid_cube)
    cSoil_data = cSoil_cube.data
    cSoil_data = ma.masked_where(np.logical_or(cSoil_data < 1e-8, cSoil_data > 1e8), cSoil_data)

    if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM':
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


    spatial_ensemble_cSoil[model_i, :, :] = soil_carbon_data


    # model mrsos
    mrsos_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/mrsos_Lmon_'+model+'_historical*', model)
    mrsos_cube = open_netCDF(mrsos_cube)        
    mrsos_cube = select_time(mrsos_cube, 2001, 2010)
    mrsos_cube = time_average(mrsos_cube)
    mrsos_cube = regrid_model(mrsos_cube, regrid_cube)
    mrsos_data = mrsos_cube.data/100

    spatial_ensemble_mrso[model_i, :, :] = mrsos_data


# MAP ENSEMBLE MEAN
cSoil_ensemble_mean = ma.mean(spatial_ensemble_cSoil, axis=0)
mrso_ensemble_mean = ma.mean(spatial_ensemble_mrso, axis=0)

#%% Spatial Correlation
model_flatten_cSoil = cSoil_ensemble_mean.flatten()
model_flatten_mrso = mrso_ensemble_mean.flatten()
r_coeffient = ma.corrcoef(model_flatten_cSoil, model_flatten_mrso)
print('r-coefficent: CMIP6 ensemble cSoil & mrso', r_coeffient)