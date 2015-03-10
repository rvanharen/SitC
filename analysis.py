#!/usr/bin/env python2


'''
Description:    
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:        -
Last Modified:  -
License:        Apache 2.0
Notes:          -
'''

from time_filter_wunderground_data import time_filter_ncfile
import datetime
from numpy import array as nparray

class analysis:
    def __init__(self, time_object, variable, time_of_day_object):
        '''
        description
        '''
        self.time_object = time_object
        self.variable = variable
        self.time_of_day_object = time_of_day_object
        self.extract_time()
        
    def extract_time(self):
        '''
        description
        '''
        if not isinstance(self.time_of_day_object, datetime.time):
            raise TypeError('time_of_day_object is not of type datetime.time')
        time_of_day = [x.time() for x in self.time_object]
        # index of items at time time_of_ay
        indices = [i for i, x in enumerate(time_of_day) if x == self.time_of_day_object]
        # extract variables
        self.variable = nparray(self.variable)[indices].tolist()
        self.time_object = nparray(self.time_object)[indices].tolist()
        
    def get_results(self):
        '''
        description
        '''
        return self.time_object, self.variable

if __name__ == "__main__":
    # define some test variables
    filename = '/home/ronald/NLeSC/SitC/github/SitC/IDRENTHE44.nc'
    timedelta = 60
    timewindow = 6
    method = 'average'
    time_out = datetime.time(13,0) # 13:00
    # filter ncfile 
    filtered = time_filter_ncfile(filename, timedelta, timewindow, method)
    # get an instance of the class
    analysis = analysis(filtered.time_out, filtered.temp_out, time_out)
    # extract time
    analysis.extract_time()
    # get results
    time_object, variable = analysis.get_results()
    
    # plot 
    from pylab import *
    scatter(time_object, variable)
    show()
    import pdb; pdb.set_trace()