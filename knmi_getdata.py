#!/usr/bin/env python2

'''
Description:
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:        -
Last Modified:  -
License:        Apache 2.0
Notes:          -
'''

from lxml.html import parse
import csv
import urllib2
from lxml import html
import numbers
import json
import os
import utils
from numpy import vstack
import argparse


class get_knmi_reference_data:
    '''
    description
    '''
    def __init__(self, opts):
        #self.outputdir = opts.outputdir
        self.csvfile = opts.csvfile
        self.outputdir = opts.outputdir
        self.keep = opts.keep
        self.check_output_dir()
        if len(opts.stationid)==0:
            self.get_station_ids()
        else:
            self.stationdids = [opts.stationid]
        self.download_station_data()
        self.get_station_locations()

    def get_station_ids(self):
        '''
        get all stationids from the KNMI website
        '''
        self.url = 'http://www.knmi.nl/klimatologie/uurgegevens/'
        page = parse(self.url)
        # get list of ids
        rows = page.xpath(".//tbody/@id")
        #self.stationids = [int(stationid[3:]) for stationid in rows]
        self.stationids = [str(stationid) for stationid in rows]
        

    def download_station_data(self):
        page = parse(self.url)
        for stationid in self.stationids:
            print stationid
            relpaths = page.xpath(".//tbody[@id='" + stationid + "']/tr/td/span/a/@href")
            for path in relpaths:
                fullpath = os.path.join(self.url, path)
                request = urllib2.urlopen(fullpath)
                filename = os.path.basename(path)
                outputfile = os.path.join(self.outputdir, filename)
                if self.keep:
                    if os.path.exists(outputfile):
                        # check if filesize is not null
                        if os.path.getsize(outputfile) > 0:
                            # file exists and is not null, continue next iteration
                            continue
                        else:
                            # file exists but is null, so remove and redownload
                            os.remove(outputfile)
                elif os.path.exists(outputfile):
                    os.remove(outputfile)
                #save
                output = open(outputfile, "w")
                output.write(request.read())
                output.close()

    def get_station_locations(self):
        # get station names for stationids
        url = 'http://www.knmi.nl/klimatologie/metadata/stationslijst.html'
        page = parse(url)
        url_metadata = page.xpath(".//table/tr/td/a/@href")
        station_name_id = [c.text for c in page.xpath(".//table/tr/td/a")]
        station_id = [s.split()[0] for s in station_name_id]
        station_names = [" ".join(s.split()[1:]) for s in station_name_id]
        for idx, stationid in enumerate(station_id):
            station_url = os.path.join(os.path.split(url)[0],
                                       url_metadata[idx])
            page = parse(station_url)
            rows = [c.text for c in page.xpath(".//table/tr/td")]
            idx_position = rows.index('Positie:') + 1
            idx_startdate = rows.index('Startdatum:') + 1
            lat, lon = rows[idx_position].encode('UTF-8').replace(
                '\xc2\xb0','').replace(' N.B. ', ',').replace(
                    'O.L.','').strip().split(',')
            lat,lon = self.latlon_conversion(lat,lon)
            try:
                dataout = vstack((dataout,
                                 [station_id[idx], station_names[idx],
                                  lat, lon, station_url]))
            except NameError:
                dataout = [station_id[idx], station_names[idx],
                           lat, lon, station_url]
        header = ['station_id', 'station_name','latitude', 'longitude', 'url']
        dataout = vstack((header, dataout))
        # write to csv file
        utils.write_csvfile(self.csvfile, dataout)
        
        # get station locations
        pass

    def latlon_conversion(self, lat, lon):
        '''
        conversion of GPS position to lat/lon decimals
            example string for lat and lon input: "52 11'"
        '''
        # latitude conversion
        latd = lat.replace("'","").split()
        lat = float(latd[0]) + float(latd[1])/60
        # longitude conversion
        lond = lon.replace("'","").split()
        lon = float(lond[0]) + float(lond[1])/60
        return lat,lon

    def check_output_dir(self):
        '''
        check if outputdir exists and create if not
        '''
        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)

if __name__ == "__main__":
    # define argument menu
    description = 'Get data KNMI reference stations'
    parser = argparse.ArgumentParser(description=description)
    # fill argument groups
    parser.add_argument('-o', '--outputdir', help='Data output directory',
                        default=os.path.join(os.getcwd(),'KNMI'),
                        required=False)
    parser.add_argument('-s', '--stationid', help='Station id',
                        default='', required=False, action='store')
    parser.add_argument('-c', '--csvfile', help='CSV data file',
                        required=True, action='store')
    parser.add_argument('-k', '--keep', help='Keep downloaded files',
                        required=False, action='store_true')
    parser.add_argument('-l', '--log', help='Log level',
                        choices=utils.LOG_LEVELS_LIST,
                        default=utils.DEFAULT_LOG_LEVEL)
    # extract user entered arguments
    opts = parser.parse_args()
    # define logger
    logname = os.path.basename(__file__) + '.log'
    logger = utils.start_logging(filename=logname, level=opts.log)
    # process data
    get_knmi_reference_data(opts)
