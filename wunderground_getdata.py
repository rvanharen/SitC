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

class get_wundergrond_data:
    def __init__(self,opts):
        self.outputdir = opts.outputdir
        self.stationid = opts.s__stationid        
        self.startdate = date(opts.startyear, 1, 1)
        self.enddate = date(opts.endyear, 12, 31)
        self.location = self.get_station_location(self.stationid)
        self.get_station_zipcode(self.location)
        self.get_data()

    def get_data(self):
        #for td in range(0,(self.enddate-self.startdate).days +1):
        for td in progressbar(range(0,(self.enddate-self.startdate).days +1),
                              "Downloading: ", 60):            
            current_date = self.startdate + timedelta(days=td)
            url = 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID=' + \
            self.stationid + '&day=' + str(current_date.day) +'&year=' + \
            str(current_date.year) + '&month=' + str(current_date.month) + '&format=1'
            outputfile = self.stationid + '_' + str(current_date.year) \
                + str(current_date.month).zfill(2) + str(current_date.day).zfill(2) + '.txt'
            with open(os.path.join(self.outputdir, outputfile), 'wb') as outfile:
                handler = urllib2.urlopen(url)
                content = handler.read()
                # convert spaces to non-breaking spaces
                content = content.replace(' ', '&nbsp;')

                # Removing all the HTML tags from the file
                outstream = cStringIO.StringIO()
                parser = htmllib.HTMLParser(
                    formatter.AbstractFormatter(formatter.DumbWriter(outstream)))
                parser.feed(content)
                content = outstream.getvalue().replace('\xa0', ' ')
                outstream.close()
                outfile.write(content)
                handler.close()

    def get_station_location(self, stationid):
        '''
        get the location of a Wunderground stationid
        '''
        url = 'http://dutch.wunderground.com/personal-weather-station/dashboard?ID=' + stationid
        handler = urllib2.urlopen(url)
        content = handler.read()
        tree = html.fromstring(content).find_class('subheading')
        if len(tree)==1:
            raw_location = str(tree[0].text_content())
        else:
            raise IOError('Cannot parse location from html file')
        location_list = [float(s) for s in raw_location.split() if is_number(s)]
        location_items = ['lat', 'lon', 'height']
        location = dict(zip(location_items, location_list))
        if int(location['lat'])==0 or int(location['lon'])==0:
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
                   range(0,len(address_components)) if
                   address_components[x]['types'][0]=='postal_code'][0]
        # return the zipcode
        return zipcode.encode('utf-8')

def progressbar(it, prefix = "", size = 60):
    '''
    progressbar for a loop
    '''
    count = len(it)
    def _show(_i):
        x = int(size*_i/count)
        sys.stdout.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), _i, count))
        sys.stdout.flush()
    _show(0)
    for i, item in enumerate(it):
        yield item
        _show(i+1)
    sys.stdout.write("\n")
    sys.stdout.flush()
    
def is_number(s):
    '''
    check if the value in the string is a number
    '''
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False
    
if __name__ == "__main__":
    # define argument menu
    description = 'Combine csv files weather underground in one output file'
    parser = argparse.ArgumentParser(description=description)
    # fill argument groups
    parser.add_argument('-o', '--outputdir', help='Data output directory',
                        default=os.getcwd(), required=True)
    parser.add_argument('-b', '--startyear', help='Start year',
                        default=2010, required=False)
    parser.add_argument('-e', '--endyear', help='End year',
                        default=date.today().year, required=False)
    parser.add_argument('-s' '--stationid', help='Station id',
                        default='IGELDERL5', required=False, action='store')
    # extract user entered arguments
    opts = parser.parse_args()
    # process data
    get_wundergrond_data(opts)