#!/usr/bin/env python2

'''
Description:
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:
Last Modified:
License:        Apache 2.0
'''

from netCDF4 import Dataset as ncdf
import netcdftime
import datetime
from numpy import array as nparray

if __name__=="__main__":
    filename = '/home/ronald/NLeSC/SiTC/Postcodekaart_validatiedata/IDRENTHE44.nc'
    time_window = 5 # extract using time_extract +/- time_window (minutes)
    
    ncfile = ncdf(filename, 'r', format='NETCDF4')
    time_axis = ncfile.variables['time']
    tempC = ncfile.variables['TemperatureC'][:]
    time_cal = netcdftime.num2date(time_axis[:],units=time_axis.units,
                                   calendar=time_axis.calendar)
    time_steps = 24*(time_cal[-1] - time_cal[0]).days + (
        time_cal[-1] - time_cal[0]).seconds/3600
    start_time = datetime.datetime(time_cal[0].year, time_cal[0].month,
                                     time_cal[0].day, time_cal[0].hour, 00)
    # initialize current_time as start_time
    current_time = start_time
    timedelta = datetime.timedelta(minutes=60)
    timewindow = datetime.timedelta(minutes=6)
    index_array = nparray(range(0,len(time_cal)))
    time_out = []
    temp_out = []
    
    min_index = 0
    while current_time < time_cal[-1]:
        # create smaller search window
        max_index = min_index + (60*24/5)
        if not max_index < len(time_cal):
            max_index = len(time_cal)            
        time_window = time_cal[min_index:max_index]
        time_index = index_array[min_index:max_index]
        # TODO: if current_time is in time_cal use that, otherwise do...
        try:
            index_up = time_index[(time_window < (
                current_time + timedelta)) & (time_window > current_time)][0]
        except IndexError:
            index_up = []
        try:
            index_down = time_index[(time_window > (
                current_time - timedelta)) & (time_window < current_time)][-1]            
        except IndexError:
            index_down = []
        if not index_up and not index_down:
            value = None
        elif not index_up and index_down:
            value = tempC[index_down]
        elif index_up and not index_down:
            value = tempC[index_up]
        elif index_up and index_down:            
            total_length = float((time_cal[index_up] -
                                  time_cal[index_down]).seconds)
            lower_length = float((current_time-time_cal[index_down]).seconds)
            value = tempC[index_down] + (
                tempC[index_up]-tempC[index_down]) * (
                    lower_length/total_length)
        time_out.append(current_time)
        temp_out.append(value)
        # increment time
        current_time += timedelta
        print current_time, value
        # update min_index
        if index_down:
            min_index = index_down
        
    import pdb; pdb.set_trace()