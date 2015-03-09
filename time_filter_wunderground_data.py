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
from numpy import mean as npmean
import argparse
import os


class time_filter_ncfile:
    def __init__(self, opts):
        self.filename = opts.inputfile
        self.timedelta = opts.timedelta
        self.timewindow = opts.timewindow
        self.method = opts.method
        self.check_file()  # check if file exists and has nc extension
        self.read_ncfile()  # open netCDF file and read variables
        if opts.method == 'interpolate':
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
        self.tempC = self.ncfile.variables['TemperatureC'][:]
        # convert time_axis to datetime object
        self.time_cal = netcdftime.num2date(time_axis[:],
                                            units=time_axis.units,
                                            calendar=time_axis.calendar)

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
        self.time_out = []
        self.temp_out = []
        min_index = 0
        # loop until we are at the end of the time in the netCDF file
        while current_time < self.time_cal[-1]:
            # create smaller search window to speed up the process
            max_index = min_index + (60*24/5)
            if not max_index < len(self.time_cal):
                max_index = len(self.time_cal)
            time_window = self.time_cal[min_index:max_index]
            time_index = index_array[min_index:max_index]
            # check if there is a measurement coinciding with current_time
            if current_time in time_window:
                # get index in the temperature array
                index = time_index[npwhere(time_window == current_time)[0][0]]
                # extract the value in the temperature array
                value = self.tempC[index]
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
                    value = None
                elif not index_up and index_down:
                    # use first value before current time if no value after
                    # current time is found in timewindow
                    value = self.tempC[index_down]
                elif index_up and not index_down:
                    # use first value after current time if no value before
                    # current time is found in timewindow
                    value = self.tempC[index_up]
                elif index_up and index_down:
                    # linear interpolation if a value before and after the
                    # current time is found within the timewindow
                    total_length = float((self.time_cal[index_up] -
                                          self.time_cal[index_down]).seconds)
                    lower_length = float((current_time -
                                          self.time_cal[index_down]).seconds)
                    value = self.tempC[index_down] + (
                        self.tempC[index_up] - self.tempC[index_down]) * (
                            lower_length/total_length)
            # append to output
            self.time_out.append(current_time)
            self.temp_out.append(value)
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
        self.time_out = []
        self.temp_out = []
        min_index = 0
        # loop until we are at the end of the time in the netCDF file
        while current_time < self.time_cal[-1]:
            # create smaller search window to speed up the process
            max_index = min_index + (60*24/5)
            if not max_index < len(self.time_cal):
                max_index = len(self.time_cal)
            time_window = self.time_cal[min_index:max_index]
            time_index = index_array[min_index:max_index]
            # check if there is a measurement coinciding with current_time
            if current_time in time_window:
                # get index in the temperature array
                index = time_index[npwhere(time_window == current_time)[0][0]]
                # extract the value in the temperature array
                value = self.tempC[index]
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
                    value = None
                else:
                    index_both = npconcatenate((index_up, index_down))
                    value = npmean(self.tempC[index_both])
            # append to output
            self.time_out.append(current_time)
            self.temp_out.append(value)
            # increment time
            current_time += timedelta
            # update min_index
            if len(index_down) > 0:
                min_index = index_down[-1]

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

    # time filter data
    filtered = time_filter_ncfile(opts)
