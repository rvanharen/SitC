#!/usr/bin/env python2

'''
Description:
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:        -
Last Modified:  -
License:        Apache 2.0
Notes:          -
'''

from datetime import date
from datetime import timedelta
import urllib2
import htmllib
import cStringIO
import formatter
import argparse
import os
import sys
from lxml import html
import numbers
import json
import csv
from combine_wunderground_data import fitem
import utils

class get_wundergrond_data:
    def __init__(self, opts):
        self.outputdir = opts.outputdir
        self.stationid = opts.stationid
        self.csvfile = opts.csvfile
        self.keep = opts.keep
        self.startdate = date(opts.startyear, 1, 1)
        self.enddate = date(opts.endyear, 12, 31)
        if not any([opts.stationid, self.csvfile]):
            raise IOError('stationid or csv file with stationids should ' +
                          'be specified')
        if self.csvfile:
            self.load_csvfile()
            # TODO: add check if required dict keys exist
            if not opts.stationid:
                stationids = self.csvdata['Station ID']
            else:
                self.stationids = [opts.stationid]
            for self.stationid in stationids:
                print self.stationid
                try:
                    station_index = self.csvdata['Station ID'].index(self.stationid)
                    self.zipcode = self.csvdata['zipcode'][station_index]
                    lon = self.csvdata['lon'][station_index]
                    lat = self.csvdata['lat'][station_index]
                    height = self.csvdata['height'][station_index]
                    # create a dictionary for the location
                    location_items = ['lat', 'lon', 'height']
                    location = dict(zip(location_items, [lon, lat, height]))
                except ValueError:
                    print "stationid not found in csvfile"
                    # get location and zipcode if not in csvfile
                    self.location = self.get_station_location(self.stationid)
                    self.zipcode = self.get_station_zipcode(self.location)
                self.outputdir = os.path.join(opts.outputdir, self.stationid)
                if not os.path.exists(self.outputdir):
                    os.makedirs(self.outputdir)
                self.get_data()

    def get_data(self):
        '''
        Download data from Weather Underground website for a given stationid
            , a startyar, and an endyear. The html file is parsed and written
            as csv to a separate txt file for each day.
        '''
        for td in utils.progressbar(range(0, (self.enddate - self.startdate).days +
                                    1), "Downloading: ", 60):
            # increase the date by 1 day for the next download
            current_date = self.startdate + timedelta(days=td)
            # set download url
            url = 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID=' + \
                self.stationid + '&day=' + str(current_date.day) + '&year=' + \
                str(current_date.year) + '&month=' + \
                str(current_date.month) + '&format=1'
            # define outputfile
            outputfile = self.stationid + '_' + str(current_date.year) \
                + str(current_date.month).zfill(2) + \
                str(current_date.day).zfill(2) + '.txt'
            # check if we want to keep previous downloaded files
            if self.keep:
                if os.path.exists(os.path.join(self.outputdir, outputfile)):
                    continue  # file exists, continue with next iteration
            elif os.path.exists(os.path.join(self.outputdir, outputfile)):
                os.remove(os.path.join(self.outputdir, outputfile))
            # open outputfile
            with open(os.path.join(self.outputdir, outputfile),
                      'wb') as outfile:
                # open and read the url
                handler = urllib2.urlopen(url)
                content = handler.read()
                # convert spaces to non-breaking spaces
                content = content.replace(' ', '&nbsp;')
                # Removing all the HTML tags from the file
                outstream = cStringIO.StringIO()
                parser = htmllib.HTMLParser(
                    formatter.AbstractFormatter(
                        formatter.DumbWriter(outstream)))
                parser.feed(content)
                # convert spaces back to regular whitespace (' ')
                content = outstream.getvalue().replace('\xa0', ' ')
                # write output
                outfile.write(content)
                # close handler and outstream
                outstream.close()
                handler.close()

    def get_station_location(self, stationid):
        '''
        get the location of a Wunderground stationid
        '''
        # set url to get the location from
        url = 'http://dutch.wunderground.com/personal-weather-station/dashboard?ID=' + stationid
        # open and read url
        handler = urllib2.urlopen(url)
        content = handler.read()
        # find the correct html tag that has the location info in it
        tree = html.fromstring(content).find_class('subheading')
        # get the string of the location
        if len(tree) == 1:
            raw_location = str(tree[0].text_content())
        else:
            raise IOError('Cannot parse location from html file')
        # remove anything non-numeric from the string and create a list
        location_list = [float(s) for s in raw_location.split() if
                         utils.is_number(s)]
        location_items = ['lat', 'lon', 'height']
        # create a dictionary for the location
        location = dict(zip(location_items, location_list))
        # check if latitude and longitude are not zero
        if int(location['lat']) == 0 or int(location['lon']) == 0:
            raise ValueError('Could not extract a valid location for ' +
                             'stationid: ' + stationid)
        return location

    def get_station_zipcode(self, location):
        '''
        get zipcode for a given location
        location['lat'] gives latitude of the location
        location['lon'] gives the longitude of the location
        '''
        # google maps api url
        url = 'https://maps.googleapis.com/maps/api/geocode/json?latlng=' + \
            str(location['lat']) + ',' + str(location['lon'])
        # open url
        handler = urllib2.urlopen(url)
        # load json
        js = json.load(handler)
        # extract the address_component
        address_components = js['results'][0]['address_components']
        # extract the zipcode from the address component
        zipcode = [address_components[x]['long_name'] for x in
                   range(0, len(address_components)) if
                   address_components[x]['types'][0] == 'postal_code'][0]
        # return the zipcode
        return zipcode.encode('utf-8')

    def load_csvfile(self):
        '''
        load data csvfile
        '''
        with open(self.csvfile, 'r') as csvin:
            reader = csv.DictReader(csvin, delimiter=',')
            try:
                self.csvdata
            except AttributeError:
                reader.next()
                try:
                    self.csvdata = {k.strip(): [fitem(v)] for k, v in
                                 reader.next().items()}
                except StopIteration:
                    pass
            current_row = 0
            for line in reader:
                current_row += 1
                if current_row == 1:  # header
                    # skip the header
                    continue
                for k, v in line.items():
                    if k is not None:  # skip over empty fields
                        k = k.strip()
                        self.csvdata[k].append(fitem(v))                 


if __name__ == "__main__":
    # define argument menu
    description = 'Combine csv files weather underground in one output file'
    parser = argparse.ArgumentParser(description=description)
    # fill argument groups
    parser.add_argument('-o', '--outputdir', help='Data output directory',
                        default=os.getcwd(), required=False)
    parser.add_argument('-b', '--startyear', help='Start year',
                        default=2010, required=False)
    parser.add_argument('-e', '--endyear', help='End year',
                        default=date.today().year, required=False)
    parser.add_argument('-s', '--stationid', help='Station id',
                        default='', required=False, action='store')
    parser.add_argument('-c', '--csvfile', help='CSV data file',
                        required=False, action='store')
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
    get_wundergrond_data(opts)
