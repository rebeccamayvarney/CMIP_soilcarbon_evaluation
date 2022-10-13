#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 7th Jan 2022
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
from rmv_cmip_analysis import numpy_to_cube
from rmv_cmip_analysis import combine_netCDF_variable, global_total


#%% REGRID CUBE
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()

landfraction_obs = combine_netCDF_variable('/home/links/rmv203/DATA/obs_datasets/luc4c_landmask.nc', 'mask')
landfraction_obs = landfraction_obs.extract(iris.Constraint(latitude=lambda p: 60<=p<=90))

#%% WISE30sec obs cSoil
WISE30sec_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/WISE30sec_Gustaf_orgC.nc')
n_lat = WISE30sec_cube.coord('latitude').points
n_lat_WISE = np.flip(n_lat, axis=0)
n_lon_WISE = WISE30sec_cube.coord('longitude').points
WISE30sec_data = WISE30sec_cube.data
print(WISE30sec_data.shape)
#WISE30sec_data = np.ma.masked_where(WISE30sec_data>1e35, WISE30sec_data)
WISE30sec_data = np.sum(WISE30sec_data[:,:,:], axis=0)
WISE30sec_data = np.flip(WISE30sec_data, axis=0)
print(WISE30sec_data.shape, np.min(WISE30sec_data), np.max(WISE30sec_data))

WISE30sec_cube =  numpy_to_cube(WISE30sec_data, regrid_cube, 2)
WISE30sec_cube = WISE30sec_cube.extract(iris.Constraint(latitude=lambda p: 60<=p<=90))
WISE30sec_gt = global_total(WISE30sec_cube, landfrac=landfraction_obs, latlon_cons=None)
WISE30sec_gt = WISE30sec_gt.data
print('WISE30sec Global totoal', WISE30sec_gt)


#%% NCSD + HWSD
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

merged_hwsd_ncscd_cube =  numpy_to_cube(merged_hwsd_ncscd_masked, regrid_cube, 2)
merged_hwsd_ncscd_gt = global_total(merged_hwsd_ncscd_cube, landfrac=landfraction_obs, latlon_cons=None)
merged_hwsd_ncscd_gt = merged_hwsd_ncscd_gt.data
print('merged_hwsd_ncscd Global totoal', merged_hwsd_ncscd_gt)

#%%

NCSCD_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/NCSCD_orgC.nc')
n_lat = NCSCD_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = NCSCD_orgC_cube.coord('longitude').points
NCSCD_orgC_data = NCSCD_orgC_cube.data
print(NCSCD_orgC_data.shape)
NCSCD_orgC_data = np.sum(NCSCD_orgC_data[0:2,:,:], axis=0)
NCSCD_orgC_data = np.flip(NCSCD_orgC_data, axis=0)
print(NCSCD_orgC_data.shape)
# hwsd
hwsd_file = Dataset('/home/links/rmv203/DATA/obs_datasets/HWSD_soilc_0.5x0.5.nc')
hwsd_data = hwsd_file.variables['soilc'][:]
# merging the soil carbon observational datasets
merged_hwsd_ncscd_new = np.copy(hwsd_data)
merged_hwsd_ncscd_new[NCSCD_orgC_data[:] > 0.] = NCSCD_orgC_data[NCSCD_orgC_data[:] > 0.]
#merged_hwsd_ncscd_masked = ma.masked_where(np.logical_or(merged_hwsd_ncscd < 0, merged_hwsd_ncscd > 1e10), merged_hwsd_ncscd)

merged_hwsd_ncscd_new_cube =  numpy_to_cube(merged_hwsd_ncscd_new, regrid_cube, 2)
merged_hwsd_ncscd_new_gt = global_total(merged_hwsd_ncscd_new_cube, landfrac=landfraction_obs, latlon_cons=None)
merged_hwsd_ncscd_new_gt = merged_hwsd_ncscd_new_gt.data
print('merged_hwsd_ncscd Global totoal', merged_hwsd_ncscd_new_gt)

#%% S2017

S2017_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/S2017_orgC.nc')
n_lat = S2017_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = S2017_orgC_cube.coord('longitude').points
S2017_orgC_data = S2017_orgC_cube.data
#S2017_orgC_data = np.ma.masked_where(S2017_orgC_data>1e35, S2017_orgC_data)
S2017_orgC_data = np.sum(S2017_orgC_data[0:2,:,:], axis=0)
S2017_orgC_data = np.flip(S2017_orgC_data, axis=0)
print(S2017_orgC_data.shape, np.min(S2017_orgC_data), np.max(S2017_orgC_data))

S2017_orgC_cube =  numpy_to_cube(S2017_orgC_data, regrid_cube, 2)
S2017_orgC_cube = S2017_orgC_cube.extract(iris.Constraint(latitude=lambda p: 60<=p<=90))
S2017_orgC_gt = global_total(S2017_orgC_cube, landfrac=landfraction_obs, latlon_cons=None)
S2017_orgC_gt = S2017_orgC_gt.data
print('S2017_orgC Global totoal', S2017_orgC_gt)


#%% GSDE_orgC

GSDE_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/GSDE_orgC.nc')
n_lat = GSDE_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = GSDE_orgC_cube.coord('longitude').points
GSDE_orgC_data = GSDE_orgC_cube.data
GSDE_orgC_data = np.ma.masked_where(GSDE_orgC_data<0, GSDE_orgC_data)
GSDE_orgC_data = np.sum(GSDE_orgC_data[0:7,:,:], axis=0)
GSDE_orgC_data = np.flip(GSDE_orgC_data, axis=0)
#GSDE_orgC_data = ma.masked_where(np.logical_or(GSDE_orgC_data < 0, GSDE_orgC_data > 998), GSDE_orgC_data)
print(GSDE_orgC_data.shape, np.min(GSDE_orgC_data), np.max(GSDE_orgC_data))

GSDE_orgC_cube =  numpy_to_cube(GSDE_orgC_data, regrid_cube, 2)
#GSDE_orgC_cube = GSDE_orgC_cube.extract(iris.Constraint(latitude=lambda p: 60<=p<=90))
GSDE_orgC_gt = global_total(GSDE_orgC_cube, landfrac=landfraction_obs, latlon_cons=None)
GSDE_orgC_gt = GSDE_orgC_gt.data
print('GSDE_orgC Global totoal', GSDE_orgC_gt)


#%% IGBP_DIS_orgC

IGBP_DIS_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/soilC/IGBP_DIS_orgC.nc')
n_lat = IGBP_DIS_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = IGBP_DIS_orgC_cube.coord('longitude').points
IGBP_DIS_orgC_data = IGBP_DIS_orgC_cube.data
#IGBP_DIS_orgC_data = np.ma.masked_where(IGBP_DIS_orgC_data>1e35, IGBP_DIS_orgC_data)
IGBP_DIS_orgC_data = np.sum(IGBP_DIS_orgC_data[:,:,:], axis=0)
IGBP_DIS_orgC_data = np.flip(IGBP_DIS_orgC_data, axis=0)
print(IGBP_DIS_orgC_data.shape, np.min(IGBP_DIS_orgC_data), np.max(IGBP_DIS_orgC_data))

IGBP_DIS_orgC_cube =  numpy_to_cube(IGBP_DIS_orgC_data, regrid_cube, 2)
IGBP_DIS_orgC_cube = IGBP_DIS_orgC_cube.extract(iris.Constraint(latitude=lambda p: 60<=p<=90))
IGBP_DIS_orgC_gt = global_total(IGBP_DIS_orgC_cube, landfrac=landfraction_obs, latlon_cons=None)
IGBP_DIS_orgC_gt = IGBP_DIS_orgC_gt.data
print('IGBP_DIS_orgC Global totoal', IGBP_DIS_orgC_gt)