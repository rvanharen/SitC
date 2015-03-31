#!/usr/bin/env python2

'''
Description:    Shared functionality:
                    * start_logging(filename, level=DEFAULT_LOG_LEVEL)
                    * progressbar(it, prefix="", size=60)
                    * is_number(s)
                    * wind_components(wind_speed, wind_direction)
                    * ismember(a, b)
Author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
Created:
Last Modified:
License:        Apache 2.0
Notes:
'''

import logging
import sys
from numpy import sin as npsin
from numpy import cos as npcos
from numpy import radians as npradians
import csv
from math import radians, cos, sin, asin, sqrt
from osgeo import osr

# define global LOG variables
DEFAULT_LOG_LEVEL = 'debug'
LOG_LEVELS = {'debug': logging.DEBUG,
              'info': logging.INFO,
              'warning': logging.WARNING,
              'error': logging.ERROR,
              'critical': logging.CRITICAL}
LOG_LEVELS_LIST = LOG_LEVELS.keys()
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
DATE_FORMAT = "%Y/%m/%d/%H:%M:%S"

def start_logging(filename, level=DEFAULT_LOG_LEVEL):
    "Start logging with given filename and level."
    try:
        logger
    except UnboundLocalError:
        logger = logging.getLogger()
    else:  # wish there was a logger.close()
        for handler in logger.handlers[:]:  # make a copy of the list
            logger.removeHandler(handler)
    logger.setLevel(LOG_LEVELS[level])
    formatter = logging.Formatter(LOG_FORMAT, datefmt=DATE_FORMAT)
    fh = logging.FileHandler(filename)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def progressbar(it, prefix="", size=60):
    '''
    progressbar for a loop
    '''
    count = len(it)

    def _show(_i):
        x = int(size*_i/count)
        sys.stdout.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x),
                                               _i, count))
        sys.stdout.flush()
    _show(0)
    for i, item in enumerate(it):
        yield item
        _show(i+1)
    sys.stdout.write("\n")
    sys.stdout.flush()

def is_number(s):
    '''
    check if the value in the string is a number and return True or False
    '''
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

def wind_components(wind_speed, wind_direction):
    '''
    return U and V wind components from wind speed and 
    wind direction (in degrees)
    '''
    U = wind_speed * npsin(npradians(wind_direction)) * -1
    V = wind_speed * npcos(npradians(wind_direction)) * -1
    return U, V

def ismember(a, b):
    '''
    return items in a that are also in b
    '''
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    # None can be replaced by any other "not in b" value
    return [bind.get(itm, None) for itm in a]

def write_csvfile(csvfile, data_out):
    # TODO: check if csv file exists already
    with open(csvfile, 'w') as fp:
        a = csv.writer(fp, delimiter=',')
        a.writerows(data_out)

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

def ReprojectCoords(coords,src_srs,tgt_srs):
    ''' Reproject a list of x,y coordinates.

        @type geom:     C{tuple/list}
        @param geom:    List of [[x,y],...[x,y]] coordinates
        @type src_srs:  C{osr.SpatialReference}
        @param src_srs: OSR SpatialReference object
        @type tgt_srs:  C{osr.SpatialReference}
        @param tgt_srs: OSR SpatialReference object
        @rtype:         C{tuple/list}
        @return:        List of transformed [[x,y],...[x,y]] coordinates

        tgt_srs and src_srs can be defined in a way similar to:
            tgt_srs=osr.SpatialReference()
            tgt_srs.ImportFromEPSG(28992)
            src_srs=osr.SpatialReference()
            src_srs.ImportFromEPSG(4326)
    '''
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    x,y,z = transform.TransformPoint(coords[0],coords[1])
    return x, y, z 
