#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Monday 24 May 2021
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
from rmv_cmip_analysis import regrid_model, select_time, time_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%%
# Set up the subplot figure
fig = plt.figure(1, figsize=(14,13))
gs = gspec.GridSpec(2, 2, figure=fig, width_ratios=[1, 0.05], hspace=0.3, wspace=0.3)
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


#%% regrid cube

regrid_cube = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/Tair_WFDEI_ann.nc')
regrid_cube = time_average(regrid_cube)
regrid_cube.coord('latitude').guess_bounds()
regrid_cube.coord('longitude').guess_bounds()
regrid_modelcube = regrid_cube
# correct lat and lon dimensions
n_lat = regrid_cube.coord('latitude').points
n_lon = regrid_cube.coord('longitude').points


#%% Observational datasets

# Soil carbon
# ncscd
ncscd_file = Dataset('/home/links/rmv203/DATA/obs_datasets/NCSCDV22_soilc_0.5x0.5.nc')
ncscd_data = ncscd_file.variables['soilc'][:]
# hwsd
hwsd_file = Dataset('/home/links/rmv203/DATA/obs_datasets/HWSD_soilc_0.5x0.5.nc')
hwsd_data = hwsd_file.variables['soilc'][:]
n_lat_ncscd = hwsd_file.variables['lat'][:]
n_lon_ncscd = hwsd_file.variables['lon'][:]
# merging the soil carbon observational datasets
merged_hwsd_ncscd = np.copy(hwsd_data)
merged_hwsd_ncscd[ncscd_data[:] > 0.] = ncscd_data[ncscd_data[:] > 0.]
merged_hwsd_ncscd_masked = ma.masked_where(np.logical_or(merged_hwsd_ncscd < 0, merged_hwsd_ncscd > 998), merged_hwsd_ncscd)
obs_Cs = merged_hwsd_ncscd_masked.copy()


# temperature
obs_temp_file = iris.load_cube('/home/links/rmv203/DATA/obs_datasets/Tair_WFDEI_ann.nc')
obs_temp_file = select_time(obs_temp_file, 2001, 2010)
obs_temp_file = time_average(obs_temp_file)
obs_temp = obs_temp_file.data
obs_temp = np.ma.masked_where(obs_temp>40, obs_temp)
current_temp = obs_temp.copy()


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
rh_cube = regrid_model(mean_cube, regrid_modelcube)
rh_data_regridded = rh_cube.data
rh_data_regridded = rh_data_regridded*1e-3*365
rh_data_regridded = ma.masked_invalid(rh_data_regridded)
rh_data_regridded = np.ma.masked_where(rh_data_regridded<=0, rh_data_regridded)
historical_observational_rh = rh_data_regridded.copy()


# MODIS Net Primary Production (NPP)
npp_file = Dataset('/home/links/rmv203/DATA/obs_datasets/MOD17A3_Science_NPP_mean_00_14_regridhalfdegree.nc')
npp_data = npp_file.variables['npp'][:]*1e-3
npp_data = np.ma.masked_where(npp_data<=0, npp_data)


#%%
# FIGURE 2a
# ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())

# # Add lat/lon grid lines to the figure
# gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
# gl.yformatter=LATITUDE_FORMATTER
# gl.ylocator = mticker.FixedLocator([-40, -20, 0, 20, 40, 60, 80])
# gl.ylabels_left=True
# gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
# gl.xformatter=LONGITUDE_FORMATTER
# gl.xlabels_bottom=True
# ax.coastlines()

# # Set up the x and y coordination
# lat = n_lat_ncscd
# lon = n_lon_ncscd
# x, y = np.meshgrid(lon, lat)

# print(np.min(obs_Cs), np.max(obs_Cs))
# line = np.arange(0, 60, 1)
# diff = plt.contourf(x, y, obs_Cs, line, cmap = 'YlGnBu', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
# #ax.text(-0.12, 1.07, 'a',transform=ax.transAxes,va = 'top',fontweight = 'bold',fontsize=34)
# ax=fig.add_subplot(gs[0,1])
# ax=plt.gca()
# fig.colorbar(diff, ax, orientation='vertical').set_label(r'$C_\mathrm{s}$ (kg C m$^{-2}$)')


# #%%
# # FIGURE 2b
# row += 1
ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())

# Add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()

# Set up the x and y coordination
lat = n_lat_ncscd
lon = n_lon_ncscd
x, y = np.meshgrid(lon, lat)

print(np.min(npp_data), np.max(npp_data))
line = np.arange(0, 2, 0.01)
diff = plt.contourf(x, y, npp_data, line, cmap = 'YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
#ax.text(-0.12, 1.07, 'b',transform=ax.transAxes,va = 'top',fontweight = 'bold', fontsize=34)
ax.set_title(r'NPP')
ax.set_ylim(-70,90)
# ax=fig.add_subplot(gs[0,1])
# ax=plt.gca()
# fig.colorbar(diff, ax, orientation='vertical').set_label(r'NPP (kg C m$^{-2}$ yr$^{-1}$)')


#%%
# FIGURE 2c
row += 1
ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())

# Add lat/lon grid lines to the figure
gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
gl.yformatter=LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
gl.ylabels_left=True
gl.xformatter=LONGITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([-90, 0, 90])
gl.xlabels_bottom=True
ax.coastlines()

# Set up the x and y coordination
lat = n_lat
lon = n_lon
x, y = np.meshgrid(lon, lat)

print(np.min(historical_observational_rh), np.max(historical_observational_rh))
line = np.arange(0, 2, 0.01)
diff = plt.contourf(x, y, historical_observational_rh, line, cmap = 'YlGn', extend='max', transform=ccrs.PlateCarree(central_longitude=0))
#ax.text(-0.12, 1.07, 'b',transform=ax.transAxes,va = 'top',fontweight = 'bold', fontsize=34)
ax.set_title(r'$R_\mathrm{h}$')
ax.set_ylim(-70,90)
# ax=fig.add_subplot(gs[1,1])
# ax=plt.gca()
# fig.colorbar(diff, ax, orientation='vertical').set_label(r'$R_\mathrm{h}$ (kg C m$^{-2}$ yr$^{-1}$)')

#%%
# plot colourbar
ax=fig.add_subplot(gs[:,1])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$C$ (kg C m$^{-2}$ yr$^{-1}$)')


#%%
fig.savefig('figA02', bbox_inches='tight')
plt.close()