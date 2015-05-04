#!/usr/bin/env python2

import glob
import os
from netCDF4 import Dataset as ncdf
import netcdftime
import datetime

class extend:
    def __init__(self, ncdir):
        self.ncdir = ncdir
        self.getfilenames()
        extend = [self.startdate(self.ncfiles[c]) for c,j in enumerate(self.ncfiles)]
        files_to_ext = [self.ncfiles[c] for c,j in enumerate(extend) if j==True]
        self.stationids = [os.path.splitext(os.path.basename(c))[0] for c in
                           files_to_ext]
        for station in self.stationids:       
            outdir = '/home/ronald/NLeSC/SitC/github/SitC/Download'
            command = './wunderground_getdata.py -o ' + outdir + ' -b 2008 -k -s ' + station
            os.system(command)
            self.remove_existing_ncfile(station)
            self.combine_files_to_ncdf(station)
            
    def getfilenames(self):
        '''
        get a list of all netcdf files already downloaded
        '''
        self.ncfiles = glob.glob(self.ncdir + '/*.nc')


    def startdate(self, filename):
        '''
        return true if first date in netcdf file is before 2010,
        else return false
        '''
        # open netCDF file
        self.ncfile = ncdf(filename, 'r', format='NETCDF4')
        # load required variables in netCDF file
        time_axis = self.ncfile.variables['time']
        # convert time_axis to datetime object
        self.time_cal = netcdftime.num2date(time_axis[0],
                                            units=time_axis.units,
                                            calendar=time_axis.calendar)
        self.ncfile.close()
        if (datetime.datetime(2009,12,1) <= self.time_cal <= datetime.datetime(2010,1,3)):
            return True
        else:
            return False

    def remove_existing_ncfile(self, station):
        '''
        check for existing ncfile, remove if found
        '''
        if os.path.exists(os.path.join(self.ncdir,  station + '.nc')):
            os.remove(os.path.join(self.ncdir, station + '.nc'))
    
    def combine_files_to_ncdf(self, station):
        '''
        combine all files to a single netcdf file -> calls other script
        '''
        inputdir = './Download/' + station
        command = './combine_wunderground_data.py -i ' + inputdir + ' -o ' + self.ncdir
        os.system(command)

    
if __name__=="__main__":
    extend('ncfiles')