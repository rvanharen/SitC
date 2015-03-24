#!/usr/bin/env python2

'''
Description:
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:        -
Last Modified:  -
License:        Apache 2.0
Notes:          -
'''

import zipfile
import csv
import StringIO
from combine_wunderground_data import fitem
from numpy import vstack
from numpy import zeros
from numpy import mean as nmean
from numpy import unique
from numpy import delete as npdelete
from numpy import array as nparray
from numpy import bincount
import os
from datetime import datetime
from time_filter_wunderground_data import time_filter_ncfile
import utils

class load_reference_data:
    def __init__(self, filename):
        self.filename = filename
        self.load_file()
        self.process_reference_data()
        
    def load_file(self):
        '''
        function description
        '''
        # load the zip file
        zipf = zipfile.ZipFile(self.filename)
        # name of csv name in zip file
        txtname = os.path.splitext(os.path.basename(self.filename))[0] + '.txt'
        # read the data in the txt file
        data = StringIO.StringIO(zipf.read(txtname))
        reader = csv.reader(data)
        # file content is ignored while start_data==False
        start_data = False
        # loop through all rows of the txt file
        for row in reader:
            if not start_data:
                if '# STN' in row:  # look for header definition in csv file
                    header = [item.strip() for item in row]
                    # found the header
                    # set start_data = True ==>> use content from here
                    start_data = True
            else:
                # we found the header definition, data starts from here
                if len(row) > 0:
                    # strip data in row and convert to int
                    # empty fiels are filled with -999
                    data_row = [int(item.strip()) if item.strip() else -999 for
                                item in row if item]
                    # create array with output data
                    try:
                        csvdata = vstack((csvdata, data_row))
                    except UnboundLocalError:
                        csvdata = data_row
        # create a dictionary from the header and output data
        self.csvdata = dict(zip(header, csvdata.T))

    def process_reference_data(self):
        '''
        process the reference csv data
        '''
        # Convert time to datetime.datetime object
        dtime = [datetime.strptime(str(item), ('%Y%m%d')) for item in
                 self.csvdata['YYYYMMDD']]
        # rain is (-1 for <0.05 mm), set to 0
        self.csvdata['RH'] = [0 if item == -1 else
                                  item for item in self.csvdata['RH']]
        # process all variables that need to be divided by 10
        for variable in ['T10', 'T', 'RH', 'FF']:
            # T10: temperature at 10 cm height, divide by 10 to convert to degC
            # T: temperature at 1.50 m height, divide by 10 to convert to degC
            # RH: rain
            # FF: wind speed
            self.csvdata[variable] = [round(0.1 * item, 1) if item != -999 else
                                item for item in self.csvdata[variable]]
        # SWD
        self.csvdata['Q'] = [round(float(10000 * item)/3600, 5) if
                             item != -999 else item for item in
                             self.csvdata['Q']]
        # api
        c = 0.95  # constant factor
        # add empty array for api to self.csvdata dictionary
        self.csvdata['api'] = zeros(len(self.csvdata['RH']))
        self.csvdata['api'][0] = c * self.csvdata['RH'][0]  # set first element
        for index, item in enumerate(self.csvdata['RH']):
            self.csvdata['api'][index] = c * item + self.csvdata['api'][index-1]

class UHI:
    '''
    class description
    '''
    def __init__(self, reference_data, wund_data):
        self.reference_data = reference_data
        self.wund_data = wund_data
        # calculate wind components reference data
        U, V = utils.wind_components(reference_data['FF'],
                                     reference_data['DD'])
        # find all Nones in self.wund_data['TemperatureC'] list
        idx_nones = [i for i,j in enumerate(self.wund_data['TemperatureC']) if
                     j is None]
        # remove idx_nones for all keys in dictionary
        for key in self.wund_data.keys():
            self.wund_data[key] = npdelete(nparray(self.wund_data[key]),
                                           idx_nones)
        # get time of day (hour) for each datetime object
        hoursx = [item.hour for item in self.wund_data['datetime']]
        # find all unique hours and its indices
        uniqueID, uniqueInd, uniqueCounts = unique(hoursx, return_inverse=True,
                                                   return_counts=True)
        temp = [float(item) for item in self.wund_data['TemperatureC']]  # ugly hack
        mtemp = bincount(uniqueInd, weights = temp) / uniqueCounts
        # TODO: check if mtemp starts from 00 or from another time
        import pdb; pdb.set_trace()

if __name__ == "__main__":
    filename = 'schiphol/uurgeg_240_2011-2020.zip'
    reference_data = load_reference_data(filename)
    filtered_wund_data = time_filter_ncfile('INOORDHO63.nc', 60, 30, 'average')
    UHI(reference_data.csvdata, filtered_wund_data.filtered)
    import pdb; pdb.set_trace()
