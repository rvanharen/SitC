#!/usr/bin/env python2

'''
Description:    Class to combine data of Wunderground csv data into a single
                output file
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:
Last Modified:
License:        Apache 2.0
Notes:          * User gives input directory containing csv txt files of
                  Wunderground data as input
                * User optionally specifies output directory
                * Combined netCDF (or TODO: csv) is created
'''

import csv
import glob
import os
from netCDF4 import Dataset as ncdf
from datetime import datetime
import netcdftime
from numpy import zeros
from numpy import argsort
from numpy import array as nparray
import time
from dateutil import tz
import argparse
import itertools


class process_raw_data:
    ''''
    Class to read the raw input data and combine them into a single output
    file
    '''
    def __init__(self, opts):
        # set class variables
        self.inputdir = opts.inputdir
        self.outputdir = opts.outputdir
        # define filename as basename inputdir with .nc extension
        filename = os.path.basename(self.inputdir) + '.nc'
        self.outputfile = os.path.join(self.outputdir, filename)
        # do we want to hardcode this?
        self.get_field_names()
        self.dateUTCstring = [s for s in self.field_names if "DateUTC" in s][0]
        self.field_names.append('<br>')
        # call functions
        self.combine_raw_data()
        self.write_combined_data_netcdf()

    def combine_raw_data(self):
        '''
        combine them
        into single output variable
        '''
        # get a list of all txt files in inputdir, sorted by filename
        filelist = sorted(glob.glob(os.path.join(self.inputdir, '*.txt')))
        if len(filelist) == 0:
            raise IOError('No files found in ' + self.inputdir)
        for inputfile in filelist:
            with open(inputfile, 'r') as csvin:
                reader = csv.DictReader(csvin, fieldnames=self.field_names,
                                        delimiter=',')
                try:
                    self.data
                except AttributeError:
                    reader.next()
                    try:
                        self.data = {k.strip(): [fitem(v)] for k, v in
                                     reader.next().items()}
                    except StopIteration:
                        pass
                current_row = 0
                for line in reader:
                    current_row += 1
                    if current_row == 1:  # header
                        # skip the header
                        continue
                    elif line['Time'] == '<br>':
                        # not a valid csv line, so skip
                        continue
                    else:
                        try:
                            datetime.strptime(line[self.dateUTCstring],
                                              ('%Y-%m-%d %H:%M:%S'))
                        except ValueError:
                            # Not a valid csv line, so skip
                            continue
                    for k, v in line.items():
                        if k is not None:  # skip over empty fields
                            k = k.strip()
                            self.data[k].append(fitem(v))
        # verify that everything is sorted with time
        if not self.verify_sorting():
            # sort data if needed according to time
            self.sort_data()

    def write_combined_data_netcdf(self):
        ncfile = ncdf(self.outputfile, 'w', format='NETCDF4')
        # description of the file
        ncfile.description = 'Hobby meteorologists data ' + self.inputdir
        ncfile.history = 'Created ' + time.ctime(time.time())
        # create time dimension
        timevar = ncfile.createDimension('time', None)
        # inititalize time axis
        timeaxis = zeros(len(self.data[self.dateUTCstring][1:]))
        # create time variable local time Europe/Amsterdam
        timeaxisLocal = zeros(len(self.data[self.dateUTCstring][1:]))
        # define UTC and local time-zone (hardcoded)
        from_zone = tz.gettz('UTC')
        to_zone = tz.gettz('Europe/Amsterdam')
        # convert time string to datetime object
        for idx in range(1, len(self.data[self.dateUTCstring])):
            # define time object from string
            timeObject = datetime.strptime(self.data[self.dateUTCstring][idx],
                                           '%Y-%m-%d %H:%M:%S')
            # tell timeObject that it is in UTC
            timeObject = timeObject.replace(tzinfo=from_zone)
            # time axis UTC
            timeaxis[idx-1] = netcdftime.date2num(
                timeObject, units='minutes since 2010-01-01 00:00:00',
                calendar='gregorian')
            # create timeaxis and convert to local time-zone
            timeaxisLocal[idx-1] = netcdftime.date2num(
                timeObject.astimezone(to_zone),
                units='minutes since 2010-01-01 00:00:00',
                calendar='gregorian')
        # netcdf time variable UTC
        timevar = ncfile.createVariable('time', 'i4', ('time',),
                                        zlib=True)
        timevar[:] = timeaxis
        timevar.units = 'minutes since 2010-01-01 00:00:00'
        timevar.calendar = 'gregorian'
        timevar.standard_name = 'time'
        timevar.long_name = 'time in UTC'
        # netcdf time variable local time Europe/Amseterdam
        timevar2 = ncfile.createVariable('local_time', 'i4', ('time',),
                                         zlib=True)
        timevar2[:] = timeaxisLocal
        timevar2.units = 'minutes since 2010-01-01 00:00:00'
        timevar2.calendar = 'gregorian'
        timevar2.standard_name = 'local_time'
        timevar2.long_name = 'time in local time Europe/Amsterdam'

        # create other variables in netcdf file
        for self.variable in self.field_names:
            if self.variable not in [self.dateUTCstring, 'Time', '<br>']:
                # add variables in netcdf file
                # check if variable is a string
                if not isinstance(self.data[self.variable][1], str):
                    # fill variable
                    if self.variable == 'SolarRadiationWatts/m^2':
                        variableName = 'SolarRadiation'
                    else:
                        variableName = self.variable
                    self.values = ncfile.createVariable(
                        variableName, type(self.data[self.variable][1]),
                        ('time',), zlib=True, fill_value=-999)
                else:
                    # string variables cannot have fill_value
                    self.values = ncfile.createVariable(
                        self.variable, type(self.data[self.variable][1]),
                        ('time',), zlib=True)
                # TODO: convert C->K, km/h->m/s ??
                try:  # fill variable
                    self.values[:] = self.data[self.variable][1:]
                except IndexError:
                    # for strings the syntax is slightly different
                    self.values = self.data[self.variable][1:]
                self.fill_attribute_data()

    def fill_attribute_data(self):
        '''
        Function that fills the attribute data of the netcdf file
        '''
        if self.variable == 'TemperatureC':
            self.values.units = 'C'
            self.values.standard_name = 'air_temperature'
            self.values.long_name = 'air temperature'
        elif self.variable == 'DewpointC':
            self.values.units = 'C'
            self.values.standard_name = 'dew_point_temperature'
            self.values.long_name = 'dewpoint temperature'
        elif self.variable == 'PressurehPa':
            self.values.units = 'hPa'
            self.values.long_name = 'surface pressure'
            self.values.standard_name = 'surface_air_pressure'
            pass
        elif self.variable == 'WindDirection':
            pass
        elif self.variable == 'WindDirectionDegrees':
            self.values.units = 'degrees'
            pass
        elif self.variable == 'WindSpeedKMH':
            self.values.units = 'km/h'
            self.values.standard_name = 'wind_speed'
            self.values.long_name = 'wind speed'
        elif self.variable == 'WindSpeedGustKMH':
            self.values.units = 'km/h'
            self.values.standard_name = 'wind_speed_of_gust'
            self.values.long_name = 'gust wind speed'
        elif self.variable == 'Humidity':
            self.values = ''
        elif self.variable == 'HourlyPrecipMM':
            self.values.units = 'mm/h'
            self.values.long_name = 'hourly precipitation'
        elif self.variable == 'Conditions':
            pass
        elif self.variable == 'Clouds':
            pass
        elif self.variable == 'dailyrainMM':
            self.values.units = 'mm/day'
            self.values.long_name = 'daily precipitation'
        elif self.variable == 'SoftwareType':
            pass
        elif self.variable == 'SolarRadiationWatts/m^2':
            self.values.units = 'Watts/m2'
            self.values.long_name = 'solar radiation'
        else:
            raise Exception('Unkown field name ' + self.variable)

    def write_combined_data_csv(self):
        '''
        Function to write the output to a csv file
        '''
        pass

    def verify_sorting(self):
        '''
        Function to verify that the data is sorted according to the time axis
        defined by self.data['DateUtC'][1:]
        '''
        if not all(earlier <= later for earlier, later in
                   itertools.izip(self.data[self.dateUTCstring][1:],
                                  self.data[self.dateUTCstring][2:])):
            return False
        else:
            return True

    def sort_data(self):
        '''
        Function to sort the data according to the time axis defined by
        self.data['DateUTC'][1:]
        '''
        idx_sort = argsort(self.data[self.dateUTCstring][1:])
        for field_name in self.field_names:
            if field_name is not self.dateUTCstring:
                self.data[field_name][1:] = nparray(
                    self.data[field_name][1:])[idx_sort].tolist()

    def get_field_names(self):
        # get a list of all txt files in inputdir, sorted by filename
        filelist = sorted(glob.glob(os.path.join(self.inputdir, '*.txt')))
        for inputfile in filelist:
            with open(inputfile, 'r') as csvin:
                reader = csv.DictReader(csvin, delimiter=',')
                # loader header information
                self.field_names = {k: [v] for k, v in
                                    reader.next().items()}
                self.field_names = self.field_names[None][0][:]
                try:
                    reader.next()
                    # first txt file with data in it found
                    # use field_names from this file
                    break
                except StopIteration:
                    pass


def fitem(item):
    item = item.strip()
    try:
        item = float(item)
    except ValueError:
        pass
    return item


if __name__ == "__main__":
    # define argument menu
    description = 'Combine csv files weather underground in one output file'
    parser = argparse.ArgumentParser(description=description)
    # fill argument groups
    parser.add_argument('-i', '--inputdir', help='Data input directory ' +
                        'containing txt files', required=True)
    parser.add_argument('-o', '--outputdir', help='Data output directory',
                        default=os.getcwd(), required=False)
    # extract user entered arguments
    opts = parser.parse_args()
    # process data
    process_raw_data(opts)
