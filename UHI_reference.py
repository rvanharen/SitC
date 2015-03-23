#!/usr/bin/env python2

import zipfile
import csv
import StringIO
from combine_wunderground_data import fitem
from numpy import vstack
import os

class load_reference_data:
    def __init__(self, filename):
        self.filename = filename
        self.load_file()
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
        import pdb; pdb.set_trace()


if __name__=="__main__":
    filename = 'schiphol/uurgeg_240_2011-2020.zip'
    load_reference_data(filename)