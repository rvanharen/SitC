#!/usr/bin/env python2

'''
Description:    Method to time filter the temperature data in the netCDF file
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:        -
Last Modified:  -
License:        Apache 2.0
Notes:          -
'''

from netCDF4 import Dataset as ncdf
import netcdftime
import datetime
from numpy import array as nparray
from numpy import where as npwhere
from numpy import concatenate as npconcatenate
from scipy.stats import nanmean
import argparse
import os
import cPickle as pickle


class time_filter_ncfile:
    def __init__(self, filename, timedelta, timewindow, method, months, hours):
        self.filename = filename
        self.timedelta = timedelta
        self.timewindow = timewindow
        #self.season = season
        #self.timeofday = timeofday
        # TODO: don't hardcode self.months and self.hours
        # self.months = [6, 7, 8]
        # self.hours = [22, 23, 0, 1, 2, 3, 4, 5]
        #self.months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        self.months = months
        self.hours = hours
        #self.months = range(1,13)
        #self.hours = range(0,24)
        #self.hours = [22, 23, 0, 1, 2, 3, 4, 5]
        self.var = ['TemperatureC', 'DewpointC', 'PressurehPa', 'Humidity',
                    'WindSpeedKMH', 'dailyrainMM', 'HourlyPrecipMM']
        self.check_file()  # check if file exists and has nc extension
        self.read_ncfile()  # open netCDF file and read variables
        if method == 'interpolate':
            self.time_filter_ncfile()  # time filter variables
        else:
            self.time_average_ncfile()
        self.close_ncfile()  # close netCDF file

    def check_file(self):
        '''
        Function to check if file exists and has a nc extension
        '''
        if not os.path.isfile(self.filename):
            raise IOError('File ' + self.filename + ' is not a file')
        elif (not os.path.splitext(self.filename)[1][1:].lower() in ['nc']):
            raise IOError('File ' + self.filename + ' has no required ' +
                          'nc extension')

    def read_ncfile(self):
        '''
        Function to open netCDF file and read the required variables, and
        convert timeaxis to datetime object
        '''
        # open netCDF file
        self.ncfile = ncdf(self.filename, 'r', formate='NETCDF4')
        # load required variables in netCDF file
        time_axis = self.ncfile.variables['time']
        # convert time_axis to datetime object
        self.time_cal = netcdftime.num2date(time_axis[:],
                                            units=time_axis.units,
                                            calendar=time_axis.calendar)
        # extract season and time of day
        idx_extract = ([y for y,c in enumerate(self.time_cal) if
                          (c.month in self.months and c.hour in self.hours)])
        self.time_cal = self.time_cal[idx_extract]
        
        self.variables = {} # create empty dictionary
        for variable in self.var:
            # fill dictionary
            try:
                self.variables[variable] = self.ncfile.variables[
                    variable][idx_extract]
            except KeyError:
                if variable == 'TemperatureC':
                    self.variables['TemperatureC'] = (self.ncfile.variables[
                        'TemperatureF'][idx_extract] - 32)/1.8
                elif variable == 'DewpointC':
                    self.variables['DewpointC'] = (self.ncfile.variables[
                        'DewpointF'][idx_extract] - 32)/1.8
                elif variable == 'PressurehPa':
                    self.variables['PressurehPa'] = 33.8639 * self.ncfile.variables[
                        'PressureIn'][idx_extract]
                elif variable == 'dailyrainMM':
                    self.variables['dailyrainMM'] = (1./0.039370) * self.ncfile.variables[
                        'dailyrainin'][idx_extract]
                elif variable == 'HourlyPrecipMM':
                    self.variables['HourlyPrecipMM'] = (1./0.039370) * self.ncfile.variables[
                        'HourlyPrecipIn'][idx_extract]                    
                elif variable == 'WindSpeedKMH':
                    self.variables['WindSpeedKMH'] = 1.609344 * self.ncfile.variables[
                        'WindSpeedMPH'][idx_extract]                    
                else:
                    raise KeyError
                
    def time_filter_ncfile(self):
        '''
        Function to time filter the measurements
        '''
        # Define the time where the filtering needs to start
        start_time = datetime.datetime(self.time_cal[0].year,
                                       self.time_cal[0].month,
                                       self.time_cal[0].day,
                                       self.time_cal[0].hour, 00)
        # initialize current_time as start_time
        current_time = start_time
        timedelta = datetime.timedelta(minutes=int(self.timedelta))
        timewindow = datetime.timedelta(minutes=int(self.timewindow))
        index_array = nparray(range(0, len(self.time_cal)))
        self.filtered = {}
        min_index = 0
        # loop until we are at the end of the time in the netCDF file
        while current_time < self.time_cal[-1]:
            # create smaller search window to speed up the process
            max_index = min_index + (60*24/5)
            if not max_index < len(self.time_cal):
                max_index = len(self.time_cal)
            time_window = self.time_cal[min_index:max_index]
            time_index = index_array[min_index:max_index]
            values = {}
            for variable in self.var:
                # check if there is a measurement coinciding with current_time
                if current_time in time_window:
                    # get index in the temperature array
                    index = time_index[npwhere(time_window == current_time)[0][0]]
                    # extract the value in the temperature array
                    values[variable] = self.variables[variable][index]
                # if there is no measurement coinciding the current_time,
                # calculate temperature from nearby measurements within the
                # defined timewindow
                else:
                    try:
                        # index of first measurent within timewindow after
                        # current_time
                        index_up = time_index[(time_window < (
                            current_time + timedelta)) & (
                                time_window > current_time)][0]
                    except IndexError:
                        # no measurements within timewindow after current_time
                        index_up = []
                    try:
                        # index of first measurent within timewindow before
                        # current_time
                        index_down = time_index[(time_window > (
                            current_time - timedelta)) & (
                                time_window < current_time)][-1]
                    except IndexError:
                        # no measurements within timewindow before current_time
                        index_down = []
                    if not index_up and not index_down:
                        # no value is found within the timewindow
                        values[variable] = None
                    elif not index_up and index_down:
                        # use first value before current time if no value after
                        # current time is found in timewindow
                        values[variable] = self.variables[variable][index_down]
                    elif index_up and not index_down:
                        # use first value after current time if no value before
                        # current time is found in timewindow
                        values[variable] = self.variables[variable][index_up]
                    elif index_up and index_down:
                        # linear interpolation if a value before and after the
                        # current time is found within the timewindow
                        total_length = float((self.time_cal[index_up] -
                                            self.time_cal[index_down]).seconds)
                        lower_length = float((current_time -
                                            self.time_cal[index_down]).seconds)
                        values[variable] = self.variables[variable][
                            index_down] + (self.variables[variable][
                                index_up] - self.variables[variable][
                                    index_down]) * (lower_length/total_length)
                # append to output
                try:
                    self.filtered[variable].append(values[variable])
                except KeyError:
                    self.filtered[variable] = [values[variable]]
            try:
                self.filtered['datetime'].append(current_time)
            except KeyError:
                self.filtered['datetime'] = [current_time]                
            # increment time
            current_time += timedelta
            # update min_index
            if index_down:
                min_index = index_down

    def time_average_ncfile(self):
        '''
        Function to time filter the measurements
        '''
        # Define the time where the filtering needs to start
        start_time = datetime.datetime(self.time_cal[0].year,
                                       self.time_cal[0].month,
                                       self.time_cal[0].day,
                                       self.time_cal[0].hour, 00)
        # initialize current_time as start_time
        current_time = start_time
        timedelta = datetime.timedelta(minutes=int(self.timedelta))
        timewindow = datetime.timedelta(minutes=int(self.timewindow))
        index_array = nparray(range(0, len(self.time_cal)))
        self.filtered = {}  # empty dictionary
        min_index = 0
        # loop until we are at the end of the time in the netCDF file
        while current_time < self.time_cal[-1]:
            # create smaller search window to speed up the process
            max_index = min_index + (60*24/5)
            if not max_index < len(self.time_cal):
                max_index = len(self.time_cal)
            time_window = self.time_cal[min_index:max_index]
            time_index = index_array[min_index:max_index]
            values = {}  # empty dictionary
            for variable in self.var:  # loop over all variables
                # check if there is a measurement coinciding with current_time
                if current_time in time_window:
                    # get index in the temperature array
                    index = time_index[npwhere(
                        time_window == current_time)[0][0]]
                    # extract the value in the temperature array
                    values[variable] = self.variables[variable][index]
                # if there is no measurement coinciding the current_time,
                # calculate temperature from nearby measurements within the
                # defined timewindow
                else:
                    try:
                        # index of first measurent within timewindow after
                        # current_time
                        index_up = time_index[(time_window < (
                            current_time + timedelta)) & (
                                time_window > current_time)]
                    except IndexError:
                        # no measurements within timewindow after current_time
                        index_up = []
                    try:
                        # index of first measurent within timewindow before
                        # current_time
                        index_down = time_index[(time_window > (
                            current_time - timedelta)) & (
                                time_window < current_time)]
                    except IndexError:
                        # no measurements within timewindow before current_time
                        index_down = []
                    if (len(index_up) == 0 and len(index_down) == 0):
                        # no value is found within the timewindow
                        values[variable] = None
                    else:
                        index_both = npconcatenate((index_up, index_down))
                        values[variable] = nanmean(
                            self.variables[variable][index_both])
                # append to output
                try:
                    self.filtered[variable].append(values[variable])
                except KeyError:
                    self.filtered[variable] = [values[variable]]
            try:
                self.filtered['datetime'].append(current_time)
            except KeyError:
                self.filtered['datetime'] = [current_time]                
            # increment time
            current_time += timedelta
            # update min_index
            try:
                if len(index_down) > 0:
                    min_index = index_down[-1]
            except UnboundLocalError:
                pass

    def close_ncfile(self):
        '''
        Function to close the netCDF file
        '''
        self.ncfile.close()

if __name__ == "__main__":
    # define argument menu
    description = 'Time filter Wunderground netCDF data'
    parser = argparse.ArgumentParser(description=description)
    # fill argument groups
    parser.add_argument('-i', '--inputfile', help='Input netCDF file ',
                        required=True)
    parser.add_argument('-o', '--outputdir', help='Data output directory',
                        default=os.getcwd(), required=False)
    parser.add_argument('--timedelta', help='length of time step in minutes',
                        default=60, required=False)
    parser.add_argument('--timewindow', help='lenght of search window in ' +
                        'minutes (+-timewindow)', default=6, required=False)
    parser.add_argument('--method', default='interpolate',
                        help='use time averaged ' +
                        'or interpolated values',
                        choices=['interpolate', 'average'], required=False)

    # extract user entered arguments
    opts = parser.parse_args()

    opts.months = range(1,13)
    opts.hours = range(0,24)
    # time filter data
    filtered = time_filter_ncfile(opts.inputfile, opts.timedelta,
                                  opts.timewindow, opts.method, opts.months,
                                  opts.hours)
    
    # save filtered as a pickled object
    if not os.path.isdir('pickled'):
        os.mkdir('pickled')  # create pickled dir if it does not exist yet
    filename = 'pickled/' + os.path.splitext(
        os.path.basename(opts.inputfile))[0] + '.pckl'
    f = open(filename, 'w')
    pickle.dump(filtered, f)
    f.close()
