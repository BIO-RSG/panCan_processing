#!/usr/bin/env python

"""
Created on Mon Apr  9 15:18:49 2018

@author: detraceyb
"""

import numpy as np
import netCDF4 as nc4
import datetime
from pathlib import Path
import re
import sys


#==================================================================
# VARIABLES THAT CAN BE CHANGED BY THE USER

# define variable that is being extracted (CHL, PAR, SST, or RRS, case-sensitive)
# example: varlist = ["CHL", "PAR", "RRS", "SST"]
# years = range(2021, 2022) # remember the end value is excluded
# doys = range(97, 123)
# sensors = ["MODIS", "VIIRS-SNPP"]
# varlist = ["CHL", "PAR", "RRS", "SST"]


years = range(int(sys.argv[1]), int(sys.argv[1])+1)
doys = range(int(sys.argv[2]), int(sys.argv[3])+1)
sensors = [sys.argv[4]]
varlist = [sys.argv[5]]


# path to script
base_script_path = '/home/claysa/panCan_processing/'

# path to input and output
base_input_path = '/mnt/data2/claysa/'
base_output_path = '/mnt/data3/claysa/'



#==================================================================

with nc4.Dataset(base_script_path + 'data/CAN_4km.nc') as nc_grid:
    bins_out = nc_grid['bin_num'][:]

master_full_f_out = np.ones(bins_out.size, dtype=np.float32 ) * np.float32(np.nan)
master_full_i_out = np.zeros(bins_out.size, dtype=np.int16 )


def extract_bins(filein, fileout, data_datetime, bins_out, var, sensor):
    
    with nc4.Dataset(filein, mode='r') as rtgrp_in, \
         nc4.Dataset(fileout, mode='w', format='NETCDF4_CLASSIC') as rtgrp_out :
             
        bin_list = rtgrp_in['/level-3_binned_data'].variables['BinList']
        bins_in = bin_list[:]['bin_num']
        
#        define bins in data
        aind = np.clip(np.searchsorted(bins_out, bins_in), None, bins_out.size - 1)
        aind_bins_out = aind[np.flatnonzero(bins_out[aind] == bins_in)]
        aind_bins_in = np.flatnonzero(bins_out[aind] == bins_in)

#       create output variables    
        rtgrp_out.createDimension('time', None)
        rtgrp_out.createDimension('binDataDim', bins_out.size)

#        define time for output
        time_out = rtgrp_out.createVariable('time', 'i4', ('time',))
        time_out.units = 'seconds since 1970-01-01 00:00:00.000 +0000'
        time_out.calendar = 'standard'
        time_out[:] = nc4.date2num(data_datetime, units=time_out.units, calendar=time_out.calendar)

#       get data from netcdf file
        if var == 'CHL' or var == 'PAR' or var == 'SST' or var == 'CHL_POLY4' or var == 'CHL_GSM_GS':
               if var == 'CHL':
                       var = 'chlor_a'
               elif var == 'PAR':
                       var = 'par'
               elif var == 'SST':
                       var = 'sst'
               x_in_compound = rtgrp_in['/level-3_binned_data'].variables[var]
               x_in = x_in_compound[:]['sum']
               weights_in = bin_list[:]['weights']
               x_in = x_in / weights_in
#               define variable names for output data
               x_out = rtgrp_out.createVariable(var, 'f4', ('time', 'binDataDim'), fill_value=np.float32(np.nan), zlib=True)
               full_x_out = np.copy(master_full_f_out)
               full_x_out[aind_bins_out] = x_in[aind_bins_in]
               x_out[0,:] = full_x_out
               
        elif var == 'RRS':
               variables = list(rtgrp_in['/level-3_binned_data'].variables.keys())
               r = re.compile('Rrs_[0-9]{3}')
               varnames = list(filter(r.match, variables))
               for variable in varnames:
                       x_in_compound = rtgrp_in['/level-3_binned_data'].variables[variable]
                       x_in = x_in_compound[:]['sum']
                       weights_in = bin_list[:]['weights']
                       x_in = x_in / weights_in
#                      define variable names for output data
                       x_out = rtgrp_out.createVariable(variable, 'f4', ('time', 'binDataDim'), fill_value=np.float32(np.nan), zlib=True)
                       full_x_out = np.copy(master_full_f_out)
                       full_x_out[aind_bins_out] = x_in[aind_bins_in]
                       x_out[0,:] = full_x_out

#       for modis and viirs, get nobs and nscenes too
        if (sensor != "SeaWiFS"):
            nobs_in = bin_list[:]['nobs']
            nscenes_in = bin_list[:]['nscenes']
#            number of observations
            x_out = rtgrp_out.createVariable('nobs', 'i2', ('time', 'binDataDim'))
            full_x_out = np.copy(master_full_i_out)
            full_x_out[aind_bins_out] = nobs_in[aind_bins_in]
            x_out[0,:] = full_x_out
#            number of scenes        
            x_out = rtgrp_out.createVariable('nscenes', 'i2', ('time', 'binDataDim'))
            full_x_out = np.copy(master_full_i_out)
            full_x_out[aind_bins_out] = nscenes_in[aind_bins_in]
            x_out[0,:] = full_x_out


# PROCESSING LOOP
#==================================================================

for sensor in sensors:
    
    print(sensor)
    
    if sensor=="MODIS":
        sletter = "A"
    elif sensor=="VIIRS-SNPP":
        sletter = "V"
    elif sensor=="SeaWiFS":
        sletter = "S"
    
    for var in varlist:
        
        print(var)
        
        
        for year in years:
            
            print(year)
            
            
            for doy in doys:
                
                # set up file names
                data_datetime = datetime.datetime(year, 1, 1, 0, 0) + datetime.timedelta(doy - 1)
                
                pathin = base_input_path + sensor + '/' + var + '/' + '{:4d}'.format(year)
                
                if var == 'SST':
                    d = datetime.date(year,1,1) + datetime.timedelta(doy - 1)
                    month = d.strftime("%m")
                    day = d.strftime("%d")
                    if sensor == 'VIIRS-SNPP':
                        filein = 'SNPP_VIIRS.' + '{:4d}'.format(year) + month + day + '.L3b.DAY.SST.nc'
                    elif sensor == 'MODIS':
                        filein = 'AQUA_MODIS.' + '{:4d}'.format(year) + month + day + '.L3b.DAY.SST.nc'
                    else:
                        filein = sletter + '{:4d}'.format(year) + month + day + '.L3b.DAY.SST.nc'
                else:
                    if sensor == 'VIIRS-SNPP':
                        filein = sletter + '{:4d}'.format(year) + '{:03d}'.format(doy) + '.L3b_DAY_SNPP_' + var + '.nc'
                    else:
                        filein = sletter + '{:4d}'.format(year) + '{:03d}'.format(doy) + '.L3b_DAY_' + var + '.nc'
                
                if var == 'CHL':
                    pathout = base_output_path + sensor + '/CHL_OCX/PANCAN/' + '{:4d}'.format(year)
                else:
                    pathout = base_output_path + sensor + '/' + var + '/PANCAN/' + '{:4d}'.format(year)
                
                fileout = sletter + '{:4d}'.format(year) + '{:03d}'.format(doy) + '.L3b_DAY_' + var + '_panCan.nc'
                
                # check to see if pathout has been created
                checkpathout = Path(pathout)
                if(not(checkpathout.exists())):
                        checkpathout.mkdir(parents = True) 
                
                # check to see if filein exists
                checkinputfile = Path(pathin + '/' + filein)
                if(not(checkinputfile.exists())):
                        print('MISSING INPUT FILE :',checkinputfile)
                        continue
                
                # check to see if filein has already been processed
                check = Path(pathout + '/' +  fileout)
                
                if(not(check.exists())):
                    print(filein, fileout, data_datetime)   
                    pfin = pathin + '/' + filein
                    pfout = pathout + '/' + fileout
                    extract_bins(pfin, pfout, data_datetime, bins_out, var, sensor)
        
