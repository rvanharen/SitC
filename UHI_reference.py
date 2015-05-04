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
from numpy import nanmax as nmax
from numpy import unique
from numpy import delete as npdelete
from numpy import array as nparray
from numpy import bincount
from numpy import percentile as nppercentile
from numpy import newaxis as npnewaxis
from numpy import meshgrid as npmeshgrid
from numpy import concatenate
from numpy import arange as nparange
from numpy import where as npwhere
from numpy import isnan
import os
from datetime import datetime
from datetime import date
from datetime import timedelta
from time_filter_wunderground_data import time_filter_ncfile
import utils
import operator
from netCDF4 import Dataset as ncdf
from osgeo import osr, gdal
from math import sqrt
import argparse
import cPickle as pickle
from numpy import corrcoef as npcor
from time import time
import glob
import random
from pylab import scatter, show, plot, savefig, clf, title, xlim,  ylim
import matplotlib.pyplot as plt
import matplotlib
from numpy import *

#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}

#matplotlib.rc('font', **font)
matplotlib.rcParams.update({'font.size': 22})

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
                    try:
                        data_row = [int(item.strip()) if item.strip()
                                    else -999 for item in row if item]
                    except ValueError:
                        import pdb; pdb.set_trace()
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
        ## Convert time to datetime.datetime object
        # hours should be 0-23 instead of 1-24
        self.csvdata['HH'] = [item if item!=24 else 0 for item in
                              self.csvdata['HH']]
        # combine date and hour into one string
        dstring = [str(self.csvdata['YYYYMMDD'][idx]) +
                   str(self.csvdata['HH'][idx]).zfill(2) for idx in
                   range(0,len(self.csvdata['YYYYMMDD']))]
        # create datetime object
        self.csvdata['datetime'] = [datetime.strptime(
            str(item), ('%Y%m%d%H')) for item in dstring]
        # Correct conversion of the datestring
        # the date of the night -> HH=24 is HH=0 on the next day!
        self.csvdata['datetime'] = [c+ timedelta(days=1) if c.hour==0 else
                                    c for c in self.csvdata['datetime']]
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
            self.csvdata['api'][index] = c * item + self.csvdata['api'][
                index-1]


class calculate_UHI:
    '''
    calculate UHI and UHI 95th percentile for a Wunderground station
    '''
    def __init__(self, reference_data, wund_data, stationid, opts,
                 weights=None):
        self.reference_data = reference_data
        self.wund_data = wund_data
        self.stationid = stationid
        self.months = opts.months
        self.starttime = opts.starttime
        self.timemax = opts.timemax
        self.lengthwindow = opts.lengthwindow
        self.drydays = opts.drydays
        self.weights = weights
        self.remove_idx_nones()
        if not self.weights:
            self.results =[self.find_shared_datetime_temperature(
                reference_data)]
            tru = self.results[0]['tru']
        else:
            self.results = [self.find_shared_datetime_temperature(
                reference_data[idx]) for idx in range(0,len(reference_data))]
            tru = [(sum([(weights[idx]*self.results[idx]['tru'][
                dtime]) for idx in range(0, len(
                    reference_data))])/sum(weights)) for dtime in range(
                        0, len(self.results[0]['dtime']))]
        if len(self.results[0]['turn'])==0:
            # no times shared or no dry days (or rain data=-999)
            self.corrcoef = 0
            return
        self.UHI = nparray(self.results[0]['turn']) - nparray(tru)
        self.calculate_daily_UHI()
        if hasattr(self, 'UHI_av'):
            self.UHI95 = nppercentile(self.UHI_av, 95)
            self.UHI50 = nppercentile(self.UHI_av, 50)
        self.corrcoef = npcor(self.results[0]['turn'], tru)[0][1]
        print "correlation coefficient: " + str(self.corrcoef)
        # calculate UHI by subtracting reference temperature from station
        # measurement for each timestep
        if not os.path.isfile('figs/' + stationid + '.png'):
            plot(range(0,len(self.results[0]['turn'])),
                 self.results[0]['turn'], 'b')
            plot(range(0,len(tru)), tru, 'r')
            savefig('figs/' + stationid + '.png', bbox_inches='tight')
            clf()
        #show()

    def remove_idx_nones(self):
        '''
        find all Nones in self.wund_data['TemperatureC'] list
        '''
        idx_nones = [i for i,j in enumerate(self.wund_data['TemperatureC']) if
                     j is None or j<0 or j>45]
        # remove idx_nones for all keys in dictionary
        for key in self.wund_data.keys():
            self.wund_data[key] = npdelete(nparray(self.wund_data[key]),
                                           idx_nones)

    def find_shared_datetime_temperature(self, ref_data):
        '''
        find station temperature and reference temperature on shared datetimes
        '''
        # load datetime objects for reference data and Wunderground data
        ref_datetime = ref_data['datetime']
        # TODO: move this to creation of pckl file
        ref_datetime = [c + timedelta(days=1) if c.hour==0 else c for c in ref_datetime]
        wund_datetime = self.wund_data['datetime']
        # find shared datetimes between two reference/Wunderground data
        shared_datetime_1 = utils.ismember2(ref_datetime, wund_datetime)
        shared_datetime_2 = utils.ismember2(wund_datetime, ref_datetime)
        # extract temperature in both sets for shared datetimes
        tru = [ref_data['T'][i] for i in shared_datetime_2
               if i is not None]
        # datetime object of shared datetimes
        dtime = [ref_datetime[i] for i in shared_datetime_2] # if i is not None]
        turn = [self.wund_data['TemperatureC'][i] for i in shared_datetime_1
                if i is not None]
        # remove measurements for datetimes where temperature is not defined
        # in Wunderground data
        valid = [idx for idx,c in enumerate(turn) if not isnan(c)]
        turn, tru = nparray(turn)[valid], nparray(tru)[valid]
        dtime = nparray(dtime)[valid]
        #         
        turn = [turn[j] for j,c in enumerate(dtime) if (c.month in self.months and c.year>=2007)]
        tru = [tru[j] for j,c in enumerate(dtime) if (c.month in self.months and c.year>=2007)]
        dtime = [c for c in dtime if (c.month in self.months and c.year>=2007)]
        # use dry days only
        if self.drydays:
            # Extract hourly rain at each extracted time
            shared_datetime1 = utils.ismember2(dtime, ref_datetime)
            refrain = [ref_data['RH'][i] for i in shared_datetime1 if i is
                       not None]
            shared_datetime2 = utils.ismember2(dtime, wund_datetime)
            #ccc = [ref_datetime[i] for i in shared_datetime1]
            #ddd = [wund_datetime[i] for i in shared_datetime2]
    
            wundrain = [self.wund_data['HourlyPrecipMM'][i] for i in
                        shared_datetime2 if i is not None]
            # require no rain both at reference station (refrain) and 
            # weather underground station (wundrain)
            refrain = [int(round(i*10)) if not isnan(i) else NaN for i in refrain]
            wundrain = [int(round(i*10)) if not isnan(i) else NaN for i in wundrain]
            idx_dry = [c for c,j in enumerate
                       ([all(c) for c in
                         zip([True if i==0 else False for i in refrain],
                             [True if i==0 else False for i in wundrain])])
                         if j==True]
            turn = nparray(turn)[idx_dry]
            tru = nparray(tru)[idx_dry]
            dtime = nparray(dtime)[idx_dry]
        # create a dictionary containing the results and return the results
        results = dict(zip(['turn', 'tru', 'dtime'], [turn, tru, dtime]))
        return results

    def calculate_daily_UHI(self):
        '''
        calculate daily UHI
        '''
        # set startdate as the first available date
        startdate = datetime(self.results[0]['dtime'][0].year,
                             self.results[0]['dtime'][0].month,
                             self.results[0]['dtime'][0].day, 0, 0)
        # set enddate as the last available date
        enddate = datetime(self.results[0]['dtime'][-1].year,
                           self.results[0]['dtime'][-1].month,
                           self.results[0]['dtime'][-1].day, 0, 0)
        # loop through all days in between
        for td in range(0, (enddate - startdate).days):
            # increase the date by 1 day for the next download
            start_window = startdate + timedelta(days=td, hours=self.starttime)
            if start_window.month in self.months:
                end_window = start_window + timedelta(hours=self.lengthwindow)
                # average UHI over time interval
                within = [idx for idx, date in enumerate(
                    self.results[0]['dtime']) if (start_window
                                                  <= date <= end_window)]
                # require at least self.lengthwindow UHI times per time interval
                if (len(within) >= self.lengthwindow) or (self.drydays and
                                                          len(within)>=1):
                    if self.timemax:
                        UHI_av = [nmax(self.UHI[within])]
                    else:
                        UHI_av = [nmean(self.UHI[within])]
                    if not isnan(UHI_av[0]):
                        try:
                            self.UHI_av = concatenate((self.UHI_av, UHI_av))
                        except AttributeError:
                            self.UHI_av = UHI_av

def load_csv_data(csvfile):
    '''
    load data csvfile
    '''
    logger.info('Load stationdata from csv file')
    with open(csvfile, 'r') as csvin:
        reader = csv.DictReader(csvin, delimiter=',')
        try:
            csvdata
        except UnboundLocalError:
            reader.next()
            try:
                csvdata = {k.strip(): [fitem(v)] for k, v in
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
                    csvdata[k].append(fitem(v))
    return csvdata

def find_zipcode_map(lon_in, lat_in):
    '''
    return closest [ufrac, greenfrac, inwfrac] for input lon/lat
    '''
    ## convert lon_in, lat_in to EPSG:28992 spatial reference system
    # target spatial reference system
    tgt_srs=osr.SpatialReference()
    tgt_srs.ImportFromEPSG(28992)
    # source spatial reference system
    src_srs=osr.SpatialReference()
    src_srs.ImportFromEPSG(4326)
    # transformed lon_in and lat_in coordinates
    lon_in_t, lat_in_t = utils.ReprojectCoords(
        [lon_in, lat_in], src_srs, tgt_srs)
    # open geotif
    dataset = gdal.Open("zipcode/buurt_bev_dichth.tif", gdal.GA_ReadOnly)
    # create lon/lat arrays of geotif
    gt = dataset.GetGeoTransform()
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    lon = gt[0] + gt[1]*nparange(0,cols)
    lat = gt[3] + gt[5]*nparange(0,rows)
    # extract window surrounding point
    lon_window = lon[(lon >= lon_in_t - 50) & (lon <= lon_in_t + 50)]
    lat_window = lat[(lat >= lat_in_t - 50) & (lat <= lat_in_t + 50)]
    # create meshgrid
    lonx, latx = npmeshgrid(lon_window,lat_window)
    # reshape to one dimensional arrays
    lonx = lonx.reshape(-1)
    latx = latx.reshape(-1)
    # calculate distance to each point in the surrounding window
    distance = [sqrt((lon_in_t-lonx[idx])**2 + (lat_in_t-latx[idx])**2) for idx
                in range(0,len(lonx))]
    # find index of closest reference station to wunderground station
    min_index, min_value = min(enumerate(distance), key=operator.itemgetter(1))
    lon_sel, lat_sel = lonx[min_index], latx[min_index]
    # indices of gridpoint
    lat_idx = lat[:].tolist().index(lat_sel)
    lon_idx = lon[:].tolist().index(lon_sel)
    # extract gridpoint
    inwfrac = dataset.GetRasterBand(1).ReadAsArray()[lat_idx,lon_idx]
    dataset = None  # close geotiff
    # green fraction
    dataset = gdal.Open("zipcode/conv_green.tif", gdal.GA_ReadOnly)
    greenfrac = dataset.GetRasterBand(1).ReadAsArray()[lat_idx,lon_idx]
    dataset = None  # close geotiff
    # urban fraction
    dataset = gdal.Open("zipcode/conv_water.tif", gdal.GA_ReadOnly)
    ufrac = dataset.GetRasterBand(1).ReadAsArray()[lat_idx,lon_idx]
    dataset = None  # close geotiff
    # return 
    #import pdb; pdb.set_trace()
    return [ufrac, greenfrac, inwfrac]

def create_list_of_stations():
    '''
    Return a randomly ordered list of stationids from the netCDF files in the
    subdirectory ncfiles. Remove stationids that have measurement issues.
    '''
    from numpy import hstack
    stationfiles = glob.glob('ncfiles/*')  # get a list of all station files
    # list of stations with data issues that cause an exception
    issuelist = ['IDAVISVA3', 'IDRENTHE13', 'IDRENTH15', 'IGELDERL145',
                 'IFRIESLA43', 'IGRONING120', 'IGRONING103', 'INOORDBR139',
                 'IDRENTHE50', 'INOORDBR138', 'INOORDBR129', 'IGELDERL153',
                 'INOORDBR140', 'IZUIDHOL181', 'ILIMBURG72', 'IFRIESLA51',
                 'IZUIDHOL107', 'IGELDERL142', 'IUTRECHT23', 'IUTRECHT125',
                 'IUTRECHT96', 'IUTRECHT24', 'ISOUTHHO16', 'INOORDHO179',
                 'INORTHHO20', 'IZEELAND55', 'IZEELAND25', 'IOVERIJS111',
                 'IZEELAND4']
    issue_measurements = ['I90578852', 'IDRENTHE28', 'IDRENTHE52',
                          'IFRIESLAN2', 'IFLEVOLA25', 'IFRIESLA22',
                          'IFRIESLA71', 'IGELDER29', 'IGELDER5',
                          'IGELDER70', 'IGELDER75', 'IGRONING41',
                          'IGRONING56', 'IHAARLEM1', 'ILIMBURG29',
                          'ILIMBURG36', 'ILIMBURG47', 'INHOOSTZ1',
                          'INOORDBR48', 'INOORDHO102', 'INOORDHO108',
                          'INOORDHO070', 'INORTHH08', 'IOVERIJS50',
                          'IUTRECHT56', 'IUTRECHT63', 'IZEELAND36',
                          'IZHLEIDE2', 'IZHVLAAR1', 'IZUIDHOL104',
                          'IZUIDHOL78',
                          'IFLEVOLA16', 'INBROOSE2', 'INORTHHO3', 'ILIMBURG40']
    # combine issue lists
    issuelist = hstack((issuelist, issue_measurements))
    # create stationids from stationfiles for stationids not in issuelist
    stationids = [os.path.splitext(os.path.basename(c))[0] for c in
                  stationfiles if os.path.splitext(os.path.basename(c))[0]
                  not in issuelist]
    random.shuffle(stationids)  # randomize the order of the list
    return stationids

def main(opts):
    # load csv file KNMI reference stations
    reference_stations = load_csv_data(opts.knmifile)
    # station_ids were converted to floats when reading csvdata, convert to int
    reference_stations['station_id'] = [int(c) for c in
                                        reference_stations['station_id']]
    # load csv file wunderground stations
    wunderground_stations = load_csv_data(opts.wundfile)
    stationids = create_list_of_stations()
    # load reference station data
    reference_data = []
    for reference_station in reference_stations['station_id']:
        print str(reference_station)
        filenameKNMI = 'pickled/KNMI_uurgeg_' + str(reference_station) + '.pckl'
        reference_data.append(calculate_load_knmi_data(filenameKNMI,
                                                       str(reference_station)))
    i=0       
    for stationid in stationids:
        # find index of wunderground station
        try:
            index_wunderground = wunderground_stations['Station ID'].index(
                stationid)
        except ValueError:
            continue
        # if stationtype argument is given, only stationids with the
        # (lowercase converted) stationtype measurement devices are considered
        try:
            if opts.stationtype:
                # require specific station type
                if not True in [c in  wunderground_stations
                                ['Station Type'][index_wunderground].lower()
                                for c in opts.stationtype]:
                    raise Exception
        except Exception:
            continue
        # calculate distance of wunderground station to all KNMI reference stations
        distance = [utils.haversine(
            wunderground_stations['lon'][index_wunderground],
            wunderground_stations['lat'][index_wunderground],
            reference_stations['longitude'][idx],
            reference_stations['latitude'][idx]) for idx in range(
                0,len(reference_stations['longitude']))]
        # calculate or load (if already calculated before) the filtered
        # Wunderground data
        filtered_wund_data = calculate_load_wund_data(stationid)            
        if opts.interpolate:
            # inverse distance weighted interpolateion
            # en.wikipedia.org/wiki/inverse_distance_weighting
            power = 2
            weights = [1/(c**power) for c in distance]
            UHI = calculate_UHI(reference_data,
                        filtered_wund_data.filtered, stationid, opts,
                        weights)
        else:
            # use value at closest reference station instead of interpolated
            # find index of closest reference station to wunderground station
            min_index, min_value = min(enumerate(distance),
                                    key=operator.itemgetter(1))
            UHI = calculate_UHI(reference_data[min_index],
                                filtered_wund_data.filtered, stationid,
                                opts)
        # Require a correlation of at least 0.7 between reference temperature
        # serie and station temperature serie. In addition, UHI.UHI95 should
        # exist
        if UHI.corrcoef < 0.85 or not hasattr(UHI, 'UHI95'):
            continue  # someting must be wrong, skip the station
        UHI_station = [wunderground_stations['lat'][index_wunderground],
            wunderground_stations['lon'][index_wunderground],
            UHI.UHI95, UHI.UHI50, len(UHI.UHI_av)]
        # create a list of station names
        try:
            stationlist = vstack((stationlist, stationid))
        except UnboundLocalError:
            stationlist = stationid
        zipcode = find_zipcode_map(
            wunderground_stations['lon'][index_wunderground],
            wunderground_stations['lat'][index_wunderground])
        # combine UHI_station and zipcode
        UHIzip_station = concatenate((UHI_station, zipcode))
        try:
            UHIzip = vstack((UHIzip, UHIzip_station))
        except NameError:
            UHIzip = UHIzip_station
        i += 1
        print i
        if i>300:
            break
    # require at least 360 days of data
    import pdb; pdb.set_trace()
    UHIdata = nparray([UHIzip[idx,:] for idx,c in
                       enumerate(UHIzip[:,4]) if c>len(opts.months)*30*4])
    stationlist = [stationlist[idx,:] for idx,c in
                       enumerate(UHIzip[:,4]) if c>len(opts.months)*30*4]
    # create spatial scatter plot of UHI values
    plot_scatter_spatial(UHIdata[:,1], UHIdata[:,0], UHIdata[:,2],
                         'spatial.png')
    # fit statistical model
    fit, recons = fit_statistical_model(UHIdata, 'reconst.png')
    import pdb; pdb.set_trace()
    plot_scatter_spatial(UHIdata[:,1], UHIdata[:,0], recons,
                         'spatial_recons.png')

def fit_statistical_model(UHIdata, filename):
    '''
    Fit a statistical model to the computed UHI values using the zipcode map.
    Statistical model has the form
        UHI = a*ufrac + b*inw + c*greenfrac + d
    '''
    import scipy.optimize as optimize
    from numpy import ones, hstack
    # data to fit the statistical model
    ydata = UHIdata[:,2]
    xdata = UHIdata[:,6:]
    # bevd should at least be >0
    ydata = ydata[xdata[:,1]>0]
    xdata = xdata[xdata[:,1]>0]
    # scale variables to correct units
    xdata[:,0] = xdata[:,0] / 100
    
    #xdata[:,1] = xdata[:,1] / 10
    xxdata = range(1,10000)
    yydata = ydata - -0.04225098 * xdata[:,1]
    aaa = [1.965*c**0.138 for c in xxdata]
    bbb = [6.16036848*c**(-0.01301373) for c in xxdata]
    plt.figure()
    plot(xxdata, aaa, 'r')
    scatter(xdata[:,1], yydata)
    plot(xxdata, bbb, 'b')
    
    # add an array of ones to optionally fit a constant
    xdata = hstack((xdata,ones((len(xdata[:,1]),1))))
    # initial guess of fit
    x0 = nparray([0.0, 0.0, 0.0])
    # fit model 1
    fit = optimize.leastsq(func, x0, args=(xdata, ydata))[0]
    # create reconstruction from model 1
    recons = fit[0] * xdata[:,0] + fit[1] * xdata[:,1] + fit[2] * xdata[:,2]
    print fit
    # fit model 2
    x0 = nparray([0.0, 0.0, 0.0])
    fit = optimize.leastsq(func2, x0, args=(xdata, ydata))[0]
    #from leastsq_bounds import leastsq_bounds
    #fit2 =leastsq_bounds(func2, x0, bounds=[[-inf,inf],[-inf,inf],[0,1]],boundsweight=100,args=(xdata, ydata))[0]

    # create reconstruction from model 2
    recons2 = fit[0]*xdata[:,0] + fit[1]*(xdata[:,1]**fit[2])
    print fit
    #from leastsq_bounds import leastsq_bounds
    #fit2 =leastsq_bounds(func2, x0, bounds=[[-inf,inf],[-inf,inf],[0,1]],boundsweight=100,args=(xdata, ydata))[0]

    # create reconstruction from model 2
    plt.figure()
    plt.plot(xdata[:,0], ydata, 'o')
    # calc the trendline
    z = polyfit(xdata[:,0], ydata, 1)
    p = poly1d(z)
    # the line equation:
    trendline = str('y=%.6fx+(%.6f)'%(z[0],z[1]))
    plt.plot(xdata[:,0], p(xdata[:,0]),'r--', label=trendline)
    plt.legend()
    plt.xlabel('green fraction [%]')
    plt.ylabel('UHI')
    plt.savefig('trend.png', bbox_inches='tight')
    plt.figure()
    #xlim(0,5)
    #ylim(0,5)
    plt.scatter(recons, ydata)
    plt.scatter(recons2, ydata, color='r')
    plt.title('Correlation: ' + str(npcor(recons, UHIdata[:,2])[0][1]))
    print npcor(recons, UHIdata[:,2])[0][1], npcor(recons2, UHIdata[:,2])[0][1]
    plt.xlabel('reconstructed')
    plt.ylabel('observed')
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    plt.savefig('reconst.png', bbox_inches='tight')
    return fit, recons

def calculate_load_wund_data(stationid):
    '''
    Calculate or load Wunderground station data:
        pickled file exists -> load
        pickled file doesn't exist -> calculate
    '''
    if not os.path.isfile('pickled/' + stationid + '.pckl'):
        filtered_wund_data = time_filter_ncfile('ncfiles/' + stationid + '.nc',
                                                60, 10, 'average', range(1,13),
                                                range(0,24))
        if not os.path.isdir('pickled'):
            # create pickled dir if it does not exist yet
            os.mkdir('pickled')
        # save filtered as a pickled object            
        filename = 'pickled/' + stationid +'.pckl'
        f = open(filename, 'w')
        pickle.dump(filtered_wund_data, f, -1)
        f.close()
    else:
        # load data from pickled object
        f = open('pickled/' + stationid + '.pckl')
        filtered_wund_data = pickle.load(f)
        f.close()
    return filtered_wund_data

def calculate_load_knmi_data(filenameKNMI, reference_station):
    '''
    Calculate or load KNMI reference data:
        pickled file exists -> load
        pickled file doesn't exist -> calculate
    '''
    # add check for filesize 0
    if not os.path.isfile(filenameKNMI):
        # save filtered as a pickled object
        if not os.path.isdir('pickled'):
            os.mkdir('pickled')  # create pickled dir if it does not exist yet
        # generate filename of closest reference station
        filename = 'KNMI/uurgeg_' + reference_station + '_2001-2010.zip'
        filename2 = 'KNMI/uurgeg_' + reference_station + '_2011-2020.zip'
        # load reference data from csv files in zipfiles
        dict1 = load_reference_data(filename).csvdata
        dict2 = load_reference_data(filename2).csvdata
        # combine dicts
        reference_data = dict((k, concatenate((dict1.get(k), dict2.get(k))))
                            for k in set(dict1.keys() + dict2.keys()))
        # save reference data as a pickled object
        f = open(filenameKNMI, 'w')
        pickle.dump(reference_data, f, -1)
        f.close()                        
    else:
        # load pickled object
        f = open(filenameKNMI)
        reference_data = pickle.load(f)
        f.close()
    return reference_data

def plot_scatter_spatial(lon, lat, var, filename):
    '''
    Scatter plot of spatial point data using the colormap RdYlBu_r for the
    domain of the Netherlands. The plot is saved to the filename that is given
    as argument to the function.
    '''
    from mpl_toolkits.basemap import Basemap
    cm = plt.cm.get_cmap('RdYlBu_r')
    plt.clf()  # close all open plots (if any)
    fig = plt.figure()  # create empty plot
    fig.set_size_inches(18.5, 10.5)
    # create discrete colorbar
    m = Basemap(projection='merc', lat_0 = 52.5, lon_0 = 5.5, resolution = 'i',
                area_thresh = 1000.0, llcrnrlon = 3, llcrnrlat = 50.8,
                urcrnrlon = 7.2, urcrnrlat = 54)
    # compute map projection coordinates for lat/lon grid.
    x, y = m(lon,lat)
    m.drawcoastlines()
    m.drawcountries()
    sc = m.scatter(x, y, c=var, s=200, marker='o', cmap=cm) #, vmin=0.5, vmax=6.5)
    plt.colorbar(sc)
    plt.savefig(filename, bbox_inches='tight')

def func(params, xdata, ydata):
    '''
    The function whose square is to be minimised.
    params ... list of parameters tuned to minimise function.
    Further arguments:
    xdata ... design matrix for a linear model.
    ydata ... observed data.
    '''
    from numpy import dot as npdot
    return (ydata - npdot(xdata, params))

def func2(params, xdata, ydata):
    return ydata - ((params[0]*xdata[:,0]) + (params[1]*(xdata[:,1]**params[2]))) #+ params[3])

def func3(params, xdata, ydata):
    return ydata - ((params[0]*(xdata[:,1]**params[1]))) #+ params[3])

if __name__ == "__main__":
    # define logger
    logname = os.path.basename(__file__) + '.log'
#    logger = utils.start_logging(filename=logname, level=opts.log)
    logger = utils.start_logging(filename=logname, level='info')
    # define argument menu
    description = 'Time filter Wunderground netCDF data'
    parser = argparse.ArgumentParser(description=description)
    # fill argument groups
    parser.add_argument('-w', '--wundfile', help='Wunderground csv file ' +
                        '[default: wunderground_stations.csv]',
                        default='wunderground_stations.csv', required=False)
    parser.add_argument('-k', '--knmifile', help='KNMI csv file [default: ' +
                        'knmi_reference_data.csv]',
                        default='knmi_reference_data.csv', required=False)
    parser.add_argument('-i', '--interpolate', help='Distance weighted ' +
                        'interpolation of KNMI reference data instead of ' +
                        'nearest reference station', required=False,
                        action='store_true')
    parser.add_argument('-t', '--timemax', help='Use maximum of time window ' +
                        'instead of average', required=False,
                        action='store_true')                        
    parser.add_argument('-s', '--stationtype', help='Require a certain ' +
                        'instrument for the Wunderground station',
                        required=False, type=str, nargs='+')
    parser.add_argument('-m', '--months', required=True, type=int, nargs='+',
                        help='month numbers (1-12) separated by space used to ' +
                        'calculate UHI')
    parser.add_argument('-b', '--starttime', required=False, type=int,
                        default=22, help='start time (hour) [default = 22]')
    parser.add_argument('-l', '--lengthwindow', required=False, type=int,
                        default=7,
                        help='Length of time window (hours) [default = 7]')
    parser.add_argument('-d', '--drydays', required=False, action='store_true',
                        help='Use only dry time intervals')    
    # extract user entered arguments
    opts = parser.parse_args()
    # main function
    main(opts)