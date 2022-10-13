#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 4th Jan 2022
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)
"""

#%%

# Analysis imports
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import iris
import iris.coord_categorisation

from rmv_cmip_analysis import combine_netCDF_variable, numpy_to_cube, global_total

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%% REGRID CUBE
regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()

landfraction_obs = combine_netCDF_variable('/home/links/rmv203/DATA/obs_datasets/luc4c_landmask.nc', 'mask')

#%% WISE30sec obs cSoil
WISE30sec_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/soilC/WISE30sec_Gustaf_orgC.nc')
n_lat = WISE30sec_cube.coord('latitude').points
n_lat_WISE = np.flip(n_lat, axis=0)
n_lon_WISE = WISE30sec_cube.coord('longitude').points
WISE30sec_data = WISE30sec_cube.data
print(WISE30sec_data.shape)
WISE30sec_data = np.sum(WISE30sec_data[0:5,:,:], axis=0)
WISE30sec_data = np.flip(WISE30sec_data, axis=0)
print(WISE30sec_data.shape, np.min(WISE30sec_data), np.max(WISE30sec_data))
WISE30sec_data_flatten = WISE30sec_data.flatten()

WISE30sec_cube =  numpy_to_cube(WISE30sec_data, regrid_cube, 2)
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

merged_hwsd_ncscd_masked_flatten = merged_hwsd_ncscd_masked.flatten()
r_coeffient = ma.corrcoef(merged_hwsd_ncscd_masked_flatten, WISE30sec_data_flatten)
print('r-coefficent (WISE30sec_data_flatten):', r_coeffient)

#%%

NCSCD_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/soilC/NCSCD_orgC.nc')
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

S2017_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/soilC/S2017_orgC.nc')
n_lat = S2017_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = S2017_orgC_cube.coord('longitude').points
S2017_orgC_data = S2017_orgC_cube.data
#S2017_orgC_data = np.ma.masked_where(S2017_orgC_data>1e35, S2017_orgC_data)
S2017_orgC_data = np.sum(S2017_orgC_data[0:3,:,:], axis=0)
S2017_orgC_data = np.flip(S2017_orgC_data, axis=0)
print(S2017_orgC_data.shape, np.min(S2017_orgC_data), np.max(S2017_orgC_data))

S2017_orgC_cube =  numpy_to_cube(S2017_orgC_data, regrid_cube, 2)
S2017_orgC_gt = global_total(S2017_orgC_cube, landfrac=landfraction_obs, latlon_cons=None)
S2017_orgC_gt = S2017_orgC_gt.data
print('S2017_orgC Global totoal', S2017_orgC_gt)

S2017_orgC_data_flatten = S2017_orgC_data.flatten()
r_coeffient = ma.corrcoef(merged_hwsd_ncscd_masked_flatten, S2017_orgC_data_flatten)
print('r-coefficent (S2017_orgC_data_flatten):', r_coeffient)


#%% GSDE_orgC

GSDE_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/soilC/GSDE_orgC.nc')
n_lat = GSDE_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = GSDE_orgC_cube.coord('longitude').points
GSDE_orgC_data = GSDE_orgC_cube.data
GSDE_orgC_data = np.ma.masked_where(GSDE_orgC_data<0, GSDE_orgC_data)
GSDE_orgC_data = np.sum(GSDE_orgC_data[0:6,:,:], axis=0)
GSDE_orgC_data = np.flip(GSDE_orgC_data, axis=0)
#GSDE_orgC_data = ma.masked_where(np.logical_or(GSDE_orgC_data < 0, GSDE_orgC_data > 998), GSDE_orgC_data)
print(GSDE_orgC_data.shape, np.min(GSDE_orgC_data), np.max(GSDE_orgC_data))

GSDE_orgC_cube =  numpy_to_cube(GSDE_orgC_data, regrid_cube, 2)
GSDE_orgC_gt = global_total(GSDE_orgC_cube, landfrac=landfraction_obs, latlon_cons=None)
GSDE_orgC_gt = GSDE_orgC_gt.data
print('GSDE_orgC Global totoal', GSDE_orgC_gt)

GSDE_orgC_data_flatten = GSDE_orgC_data.flatten()
r_coeffient = ma.corrcoef(merged_hwsd_ncscd_masked_flatten, GSDE_orgC_data_flatten)
print('r-coefficent (GSDE_orgC_data_flatten):', r_coeffient)


#%% IGBP_DIS_orgC

IGBP_DIS_orgC_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/soilC/IGBP_DIS_orgC.nc')
n_lat = IGBP_DIS_orgC_cube.coord('latitude').points
n_lat = np.flip(n_lat, axis=0)
n_lon = IGBP_DIS_orgC_cube.coord('longitude').points
IGBP_DIS_orgC_data = IGBP_DIS_orgC_cube.data
#IGBP_DIS_orgC_data = np.ma.masked_where(IGBP_DIS_orgC_data>1e35, IGBP_DIS_orgC_data)
IGBP_DIS_orgC_data = np.sum(IGBP_DIS_orgC_data[:,:,:], axis=0)
IGBP_DIS_orgC_data = np.flip(IGBP_DIS_orgC_data, axis=0)
print(IGBP_DIS_orgC_data.shape, np.min(IGBP_DIS_orgC_data), np.max(IGBP_DIS_orgC_data))

IGBP_DIS_orgC_cube =  numpy_to_cube(IGBP_DIS_orgC_data, regrid_cube, 2)
IGBP_DIS_orgC_gt = global_total(IGBP_DIS_orgC_cube, landfrac=landfraction_obs, latlon_cons=None)
IGBP_DIS_orgC_gt = IGBP_DIS_orgC_gt.data
print('IGBP_DIS_orgC Global totoal', IGBP_DIS_orgC_gt)

IGBP_DIS_orgC_data_flatten = IGBP_DIS_orgC_data.flatten()
r_coeffient = ma.corrcoef(merged_hwsd_ncscd_masked_flatten, IGBP_DIS_orgC_data_flatten)
print('r-coefficent (IGBP_DIS_orgC_data_flatten):', r_coeffient)



#%% MAP PLOTS

# Set up subplot figure
fig = plt.figure(1, figsize=(24,18))
gs = gspec.GridSpec(3, 3, figure=fig, width_ratios=[1, 1, 0.05], hspace=0.3, wspace=0.3)
column = 0
row = 0
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':36,
    'xtick.labelsize':36,
    'ytick.labelsize':36,
    'font.size':36,
}
plt.rcParams.update(params)

#%%
row = 1
ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())
#  add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()
# set up the x and y coordination
lat = n_lat
lon = n_lon
x, y = np.meshgrid(lon, lat)
print(np.min(WISE30sec_data), np.max(WISE30sec_data))
line = np.arange(0, 60, 1)
diff = plt.contourf(x, y, WISE30sec_data, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_title('WISEsec30')
ax.set_ylim(-70,90)

#%%
# Observations
row = 0

ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())
# add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()
# set up the x and y coordination
lat = n_lat_ncscd
lon = n_lon_ncscd
x, y = np.meshgrid(lon, lat)
print(np.min(merged_hwsd_ncscd_masked), np.max(merged_hwsd_ncscd_masked))
line = np.arange(0, 60, 1)
diff = plt.contourf(x, y, merged_hwsd_ncscd_masked, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_title('HWSD + NCSCD')
ax.set_ylim(-70,90)

#%%
# Observations
row = 2

ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())
# add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()
# set up the x and y coordination
lat = n_lat_ncscd
lon = n_lon_ncscd
x, y = np.meshgrid(lon, lat)
print(np.min(S2017_orgC_data), np.max(S2017_orgC_data))
line = np.arange(0, 60, 1)
diff = plt.contourf(x, y, S2017_orgC_data, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_title('S2017')
ax.set_ylim(-70,90)

#%%
# Observations
row = 0
column = 1

ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())
# add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()
# set up the x and y coordination
lat = n_lat_ncscd
lon = n_lon_ncscd
x, y = np.meshgrid(lon, lat)
print(np.min(GSDE_orgC_data), np.max(GSDE_orgC_data))
line = np.arange(0, 60, 1)
diff = plt.contourf(x, y, GSDE_orgC_data, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_title('GSDE')
ax.set_ylim(-70,90)

#%%
# Observations
row = 1
column = 1

ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())
# add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()
# set up the x and y coordination
lat = n_lat_ncscd
lon = n_lon_ncscd
x, y = np.meshgrid(lon, lat)
print(np.min(IGBP_DIS_orgC_data), np.max(IGBP_DIS_orgC_data))
line = np.arange(0, 60, 1)
diff = plt.contourf(x, y, IGBP_DIS_orgC_data, line, cmap='YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
ax.set_title('IGBP DIS')
ax.set_ylim(-70,90)

#%%
# plot colourbar
ax=fig.add_subplot(gs[:,2])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$C_{s}$ (kg C m$^{-2}$)')


#%%
# Save figure
fig.savefig('figA01', bbox_inches='tight')
plt.close()

