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
from numpy import concatenate as npconcatenate
from collections import namedtuple


class analysis:
    def __init__(self, time_object, variable, time_of_day_object):
        '''
        class the analyze data with the following methods:
        * extract_time(self)
        * extract_seasn(self, start_season, end_season=None)
        * get_results(self)
        '''
        self.time_object = time_object
        self.variable = variable
        self.time_of_day_object = time_of_day_object
        self.extract_time()

    def extract_time(self):
        '''
        extract data on a specific time of the day
        '''
        if not isinstance(self.time_of_day_object, datetime.time):
            raise TypeError('time_of_day_object is not of type datetime.time')
        time_of_day = [x.time() for x in self.time_object]
        # index of items at time time_of_ay
        indices = [i for i, x in enumerate(time_of_day) if x == self.time_of_day_object]
        # extract variables
        self.variable = nparray(self.variable)[indices].tolist()
        self.time_object = nparray(self.time_object)[indices].tolist()

    def extract_season(self, start_season, end_season=None):
        '''
        extract measurements within a season between start_season and
        end_season. If end_season is before start_season, it is assumed that 
        end_season is in the following year.
        '''
        if end_season:
            # both start_season and end_season are defined
            # extract values between start_season and end_season
            if not (isinstance(start_season, date_season) and
                    isinstance(end_season, date_season)):
                raise TypeError('start_season and end_season should be of ' +
                    'type date_season')
            else:
                # initialize indices as an empty list
                indices = []
                # find all yrs in the dataset
                yrs = set([x.year for x in self.time_object])
                # check if the season is within one year or ends in the
                # following year
                if datetime.datetime(
                    2000, end_season.month, end_season.day
                    ) >= datetime.datetime(
                        2000, start_season.month, start_season.day):
                    # end_season is after begin_season in the same year
                    for yr in yrs:
                        indices = npconcatenate((indices, [
                            i for i, x in enumerate(self.time_object) if
                            (x>=datetime.datetime(
                                yr, start_season.month, start_season.day)) and
                            (x<=datetime.datetime(yr, end_season.month,
                                                end_season.day))]))
                else:
                    # end_season is in the beginning of the following year
                    for yr in yrs:
                        indices = npconcatenate((indices, [
                            i for i, x in enumerate(self.time_object) if
                            (x>=datetime.datetime(
                                yr, start_season.month, start_season.day)) and
                            (x<=datetime.datetime(yr+1, end_season.month,
                                                end_season.day))]))
        else:
            # only start_season is defined
            # extract values on this date only
            if not (isinstance(start_season, date_season)):
                raise TypeError('start_season should be of ' +
                    'type date_season')
            else:
                indices = [i for i, x in enumerate(self.time_object) if
                           (x.month==start_season.month and
                            x.day==start_season.day)]
        # extract variables
        self.variable = nparray(self.variable)[indices.tolist()].tolist()
        self.time_object = nparray(self.time_object)[indices.tolist()].tolist()


    def get_results(self):
        '''
        return the results of an instance of the analysis class
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
    
    date_season = namedtuple('date_season', ['month', 'day'])
    dkey_start = date_season(month=10, day=1)
    dkey_end = date_season(month=3, day=1)
    analysis.extract_season(dkey_start,dkey_end)
    time_object2, variable2 = analysis.get_results()
    
    # plot 
    from pylab import *
    scatter(time_object, variable)
    show()
    import pdb; pdb.set_trace()