#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)

"""

#%%

import numpy as np
import iris
import iris.coord_categorisation
import iris.analysis.cartography
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units


#%%
def combine_netCDF_variable(directory, variable):
    """ 
    Function combines netCDF files to create one file for entire time period considered
    (variable input - OBSERVATIONAL DATASETS)

    directory - where data sits
    variable - being considered (must be the name used in the netCDF file)
    model - name of model
    """

    # Make a list of the files in the above folder to loop through
    list_files = glob.glob(directory)
    list_files = np.array(list_files)
    newlist = np.sort(list_files)

    # Make a cubelist to add each file (cube) to
    Cubelist = iris.cube.CubeList([])

    # loop for each file in newlist
    for i in range(0, len(newlist)):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FutureWarning)
            warnings.simplefilter('ignore', UserWarning)

            # Load each file named variable as a cube
            cube = iris.load_cube(newlist[i], variable)

            # Append this cube to the cubelist
            Cubelist.append(cube)

    # matching attributes
    unify_time_units(Cubelist)
    equalise_attributes(Cubelist)

    # Concatenate each cube in cubelist together to make one data file (cube)
    new_cube = Cubelist.concatenate_cube()

    return new_cube


#%%
def combine_netCDF_model(directory, model):
    """ 
    Function combines netCDF files to create one file for entire time period considered
    (model input - STANDARD MODEL)

    directory - where data sits
    variable - being considered (must be the name used in the netCDF file)
    model - name of model
    """

    # Make a list of the files in the above folder to loop through
    list_files = glob.glob(directory)
    list_files = np.array(list_files)
    newlist = np.sort(list_files)

    # Make a cubelist to add each file (cube) to
    Cubelist = iris.cube.CubeList([])

    # loop for each file in newlist
    for i in range(0, len(newlist)):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FutureWarning)
            warnings.simplefilter('ignore', UserWarning)

            # Load cube
            cube = iris.load_cube(newlist[i])

            # Append this cube to the cubelist
            Cubelist.append(cube)

    # matching attributes
    unify_time_units(Cubelist)
    equalise_attributes(Cubelist)

    # Concatenate each cube in cubelist together to make one data file (cube)
    new_cube = Cubelist.concatenate_cube()

    return new_cube


#%%
def combine_netCDF_cmip5(directory, variable, model):
    """ 
    Function combines netCDF files to create one file for entire time period considered (CMIP5 models)

    directory - where data sits
    variable - being considered (must be the name used in the netCDF file)
    model - name of model
    """

    # Make a list of the files in the above folder to loop through
    list_files = glob.glob(directory)
    list_files = np.array(list_files)
    newlist = np.sort(list_files)

    # Make a cubelist to add each file (cube) to
    Cubelist = iris.cube.CubeList([])

    # loop for each file in newlist
    for i in range(0, len(newlist)):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FutureWarning)
            warnings.simplefilter('ignore', UserWarning)

            # Load each file named variable as a cube
            cube = iris.load_cube(newlist[i], variable)

            # CRUDE WORKAROUND TO REMOVE OVERLAPPING COORDINATE
            if model == 'HadGEM2-ES':
                if i == 9: # change to 9 if timeseries plots
                    cube = cube[0:-1]

            # Append this cube to the cubelist
            Cubelist.append(cube)

    # matching attributes
    iris.util.unify_time_units(Cubelist)
    equalise_attributes(Cubelist)
    # Concatenate each cube in cubelist together to make one data file (cube)
    new_cube = Cubelist.concatenate_cube()

    return new_cube


#%%
def combine_netCDF_cmip6(directory, model):
    """ 
    Function combines netCDF files to create one file for entire time period considered
    (CMIP6 models, not rh & cSoil)

    directory - where data sits
    variable - being considered (must be the name used in the netCDF file)
    model - name of model
    """

    # Make a list of the files in the above folder to loop through
    list_files = glob.glob(directory)
    list_files = np.array(list_files)
    newlist = np.sort(list_files)

    # Make a cubelist to add each file (cube) to
    Cubelist = iris.cube.CubeList([])

    # loop for each file in newlist
    for i in range(0, len(newlist)):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FutureWarning)
            warnings.simplefilter('ignore', UserWarning)
            
            # Load cube
            cube = iris.load_cube(newlist[i])

            # remove latitude & longitude attributes
            cube.coord('latitude').attributes = {}
            cube.coord('longitude').attributes = {}  

            # matching cube metadata
            if i == 0:
                metadata1 = cube.metadata
            else:
                cube.metadata = metadata1
            
            # creating latitude and longitude bounds
            if model=='IPSL-CM6A-LR' or model=='CNRM-ESM2-1':
                 cube.coord('latitude').guess_bounds()
                 cube.coord('longitude').guess_bounds()
         
            # CESM2 bound issue fix
            if (model=='CESM2') & (i==0):
                lat_data = cube.coord('latitude').points
                lon_data = cube.coord('longitude').points
                lat_bounds = cube.coord('latitude').bounds
                lon_bounds = cube.coord('longitude').bounds
            elif (model=='CESM2') & (i>0):
                cube.coord('latitude').points = lat_data
                cube.coord('longitude').points = lon_data
                cube.coord('latitude').bounds = lat_bounds
                cube.coord('longitude').bounds = lon_bounds
    
            # removing time attributes
            if model=='IPSL-CM6A-LR':
                cube.coord('time').attributes.pop('time_origin')
            
            # Append this cube to the cubelist
            Cubelist.append(cube)

    # matching attributes
    unify_time_units(Cubelist)
    equalise_attributes(Cubelist)

    for cube in Cubelist:
        lon_bounds = Cubelist[0].coord('longitude').bounds
        cube.coord('longitude').bounds = lon_bounds

    for i, cube in enumerate(Cubelist):
        if cube.coord('time').units == Cubelist[0].coord('time').units:
            pass
        else:
            print(i)
            
    # Concatenate each cube in cubelist together to make one data file (cube)
    new_cube = Cubelist.concatenate_cube()

    return new_cube


#%%
def open_netCDF(new_cube):
    """ Function adds time coordinates to cube """

    iris.coord_categorisation.add_year(new_cube, 'time', name='year') # add year
    iris.coord_categorisation.add_month(new_cube, 'time', name ='month') # add month
    iris.coord_categorisation.add_month(new_cube, 'time', name ='decade') # add month

    return new_cube


#%%
def define_attributes(new_cube):
    """ Function returns the time, latitude and longitude coordinates of the cube """

    time_cmip5 = new_cube.coord('time').points # Defining the time variable
    lats_cmip5 = new_cube.coord('latitude').points # Defining the lat variable
    lons_cmip5 = new_cube.coord('longitude').points # Defining the lon variable

    return time_cmip5, lats_cmip5, lons_cmip5


#%%
def select_time(new_cube, lower, upper):
    """ Function that selects the time period in years """

    sliced_cube = new_cube.extract(iris.Constraint(year=lambda y: lower<=y<=upper))

    return sliced_cube


#%%
def time_average(new_cube):
    """ Function calculates time average for cube in current time unit """

    time_average_cube = new_cube.collapsed('time', iris.analysis.MEAN)

    return time_average_cube


#%%
def annual_average(new_cube):
    """ Function calculates annual average for cube """

    annual_average_cube = new_cube.aggregated_by('year', iris.analysis.MEAN)

    return annual_average_cube


#%%
def decadal_average(new_cube):
    """ Function calculates annual average for cube """

    decadal_average_cube = new_cube.aggregated_by('decade', iris.analysis.MEAN)

    return decadal_average_cube


#%%
def numpy_to_cube(np_array, similar_cube, dimensions):
    """
    Function that converts a 1, 2, or 3 dimensional numpy array to a cube.
    (Inverse is array = cube.data)
    """

    new_cube = iris.cube.Cube.copy(similar_cube) # copy similar cube

    # time, lat, lon
    if dimensions == 3:
        new_cube.data[:,:,:] = np.nan # convert new cube entries to nan
        new_cube.data[:,:,:] = np_array # fill with numpy array data

    # lat, lon
    elif dimensions == 2:
        new_cube.data[:,:] = np.nan # convert new cube entries to nan
        new_cube.data[:,:] = np_array # fill with numpy array data

    # either time, lat or lon only
    elif dimensions == 1:
        new_cube.data[:] = np.nan # convert new cube entries to nan
        new_cube.data[:] = np_array # fill with numpy array data

    # return the numpy array, failed to convert to a cube
    else:
        print('failed to convert')
        new_cube = np_array

    return new_cube


#%%
def regrid_model(cube, regridcube):
    """ Function regrids cube to the dimensions of regridcube """

    regridcube.coord('latitude').standard_name = 'latitude'
    regridcube.coord('longitude').standard_name = 'longitude'

    model_units = cube.coord('latitude').units
    regridcube.coord('latitude').units = model_units
    regridcube.coord('longitude').units = model_units

    new_model_cube = cube.regrid(regridcube, iris.analysis.Linear())

    return new_model_cube


#%%
def area_average(cube, region):
    """
    Function to create weighted area average, by collapse a cube to a weighted area average over a specified region,
    global: region = [0, 360, -90,  90]
    """
    
    # Specify the latitudes and longitudes starting from the smallest number to largest or in latitude and longitude from south to north and east to west
    lon1, lon2, lat1, lat2 = region[0], region[1], region[2], region[3] 
    # Then intersect the data at these points
    cube = cube.intersection(longitude=(lon1, lon2),latitude=(lat1, lat2))

    #cube.coord('latitude').guess_bounds()
    #cube.coord('longitude').guess_bounds()

    #  area weighting
    weights = iris.analysis.cartography.area_weights(cube)
    # Average that area by latitude and longitudes by the weighted mean
    cube = cube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=weights)

    return cube


#%%
def area_average_obs(cube, region, model_units):
    """
    Observational dataset version - 
    Function to create weighted area average, by collapse a cube to a weighted area average over a specified region,
    global: region = [0, 360, -90,  90]
    """
    
    # Specify the latitudes and longitudes starting from the smallest number to largest or in latitude and longitude from south to north and east to west
    lon1, lon2, lat1, lat2 = region[0], region[1], region[2], region[3]

    print(cube.coord('latitude').var_name)
    print(cube.coord('latitude').units.modulus)
    cube.coord('latitude').units = model_units
    cube.coord('longitude').units = model_units
    print(cube.coord('latitude').units.modulus)

    # Then intersect the data at these points
    cube = cube.intersection(longitude=(lon1, lon2),latitude=(lat1, lat2))

    # cube.coord('latitude').guess_bounds()
    # cube.coord('longitude').guess_bounds()

    # area weighting
    weights = iris.analysis.cartography.area_weights(cube)
    # Average that area by latitude and longitudes by the weighted mean
    cube = cube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=weights)

    return cube


#%%
def area_average_landfrac(cubein, landfrac=None, latlon_cons=None):
    """
    Function to create weighted area average, using landfrac for accuracy of coastal points
    """

    # landfrac
    cube = cubein.copy()
    if landfrac is not None:
        try:
            cube.data = cube.data * (landfrac.data/100)
        except:
            landfrac = landfrac.extract(latlon_cons)
            cube.data = cube.data * (landfrac.data/100)

    # if no bounds
    if cube.coord('latitude').bounds is not None:
        pass
    else:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

    #  area weighting
    weights = iris.analysis.cartography.area_weights(cube)
    # Average that area by latitude and longitudes by the weighted mean
    cube = cube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=weights)

    return cube

#%%
def weighted_area_average_landfrac(cubein, additional_weight, landfrac=None, latlon_cons=None):
    """
    Function to create weighted area average, using landfrac for accuracy of coastal points
    """

    # landfrac
    cube = cubein.copy()
    if landfrac is not None:
        try:
            cube.data = cube.data * (landfrac.data/100)
        except:
            landfrac = landfrac.extract(latlon_cons)
            cube.data = cube.data * (landfrac.data/100)

    # if no bounds
    if cube.coord('latitude').bounds is not None:
        pass
    else:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

    #  area weighting
    weights = iris.analysis.cartography.area_weights(cube)*additional_weight
    # Average that area by latitude and longitudes by the weighted mean
    cube = cube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=weights)

    return cube


#%%
def global_total(cubein, landfrac=None, latlon_cons=None):
    '''
    Function to calculate area weighted sum of non missing data: written by Eleanor Burke, Met Office
 
    Input
    ----
    * cube: :class:`iris.cube.Cube`
        variable for which total amount needs to be calculated

    * landfrac: :class:`iris.cube.Cube`
        landfrac mask so only the land fraction of the coastal points are included in the global totals

    * latlon_cons: :class:'iris.Constraint'
          used to extract the land frac sub-region for analysis if the varaible cube and land frac are on different grids      

    Returns
    -------
    * cube_areaweight: :class:`iris.cube.Cube`
        this cube is the global total

    '''

    cube = cubein.copy()
    if landfrac is not None:
        try:
            cube.data = cube.data * landfrac.data
        except:
            landfrac = landfrac.extract(latlon_cons)
            cube.data = cube.data * landfrac.data

    if cube.coord('latitude').bounds is not None:
        pass
    else:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

    weights = iris.analysis.cartography.area_weights(cube)

    cube_areaweight = cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=weights)/1e12

    return cube_areaweight


#%%
def global_total_percentage(cubein, landfrac=None, latlon_cons=None):
    '''
    Function to calculate area weighted sum of non missing data: written by Eleanor Burke,
    edited by Rebecca Varney
 
    Input
    ----
    * cube: :class:`iris.cube.Cube`
        variable for which total amount needs to be calculated

    * landfrac: :class:`iris.cube.Cube`
        landfrac mask so only the land fraction of the coastal points are included in the global totals

    * latlon_cons: :class:'iris.Constraint'
          used to extract the land frac sub-region for analysis if the varaible cube and land frac are on different grids      

    Returns
    -------
    * cube_areaweight: :class:`iris.cube.Cube`
        this cube is the global total

    '''

    cube = cubein.copy()
    if landfrac is not None:
        try:
            cube.data = cube.data * (landfrac.data/100)
        except:
            landfrac = landfrac.extract(latlon_cons)
            cube.data = cube.data * (landfrac.data/100)

    if cube.coord('latitude').bounds is not None:
        pass
    else:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

    weights = iris.analysis.cartography.area_weights(cube)

    cube_areaweight = cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=weights)/1e12

    return cube_areaweight

