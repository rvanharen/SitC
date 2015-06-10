#!/usr/bin/env python2

from pylab import *
from numpy import array as nparray
from numpy import nanmax as nmax
from UHI_reference import load_csv_data
from numpy import percentile as nppercentile
from numpy import nanmean as nmean
from datetime import datetime
from datetime import timedelta
from dateutil import tz
import utils
from numpy import isnan
import os
import shutil

def main(station):
    basename = '/home/ronald/Documenten/wunderground_gj/csv/'
    if station == "apeldoorn":
        filename = 'Apeldoornuuranalyse2008.csv'
    if station == "ijsselmuiden":
        filename = 'ijsselmuidenuuranalyse.csv'
    if station == "losser":
        filename = 'Losseruurgem.csv'
    if station == 'heerhugowaard':
        filename = "Heerhugowaardanalyse_fixed.csv"
    if station == 'leiden':
        filename = 'Leidenuuranalyse4.csv'
    if station == 'haarlem':
        filename = 'Haarlemuuranalyse.csv'
    if station == 'ijsselmonde':
        filename = 'Ijsselmondeuuranalyse_corrected.csv'
    if station == 'delft':
        filename = 'Delftuuranalyse_fixed.csv'
    if station == 'denhaag':
        filename = 'Den Haaganalyse.csv'
    if station == 'rotterdam':
        filename = 'Rotterdamuuranalyse2.csv'
    if station == 'damwoude':
        filename = "Damwoudeuuranalyse_fixed.csv"
    if station == 'purmerend':
        filename = "PurmerendanalyseGJ_fixed.csv"
    if station == 'wageningen':
        filename = 'Wageningenuuranalyse2.csv'
    if station == 'heemskerk':
        filename = 'Heemskerkuuranalyse3_fixed.csv'
    if station == 'doornenburg':
        filename = 'Doornenburguuranalyse.csv'
    fileloc = basename + filename
    # load data from csv file
    dtime, tw, tr = load_data(fileloc, station)
    #dtime = dtime[:44000]
    #tw = tw[:44000]
    #tr = tr[:44000]
    # remove nans
    dtime, tw, tr = remove_nans(dtime, tw, tr)
    dtime_utc_knmi_fixed, tw_fixed, tr_fixed = convert_wund_timezone_to_utc(
        dtime, tw,tr, station)
    # calculate daily max uhi for original dataset
    uniquetime, uhi_daily = calculate_uhi_distribution(dtime, tw, tr)
    # calculate daily max uhi for fixed dataset
    uniquetime_fixed, uhi_daily_fixed = calculate_uhi_distribution(
        dtime_utc_knmi_fixed, tw_fixed, tr_fixed)
    # plot daily cycle
    plot_daily_cycle(dtime, tr, tw, dtime_utc_knmi_fixed, tr_fixed, tw_fixed, station)
    # print_results
    print_results(uhi_daily, uhi_daily_fixed)
    import pdb; pdb.set_trace()
def load_data(fileloc, station):
    '''
    load data from csv file
    '''
    aa = load_csv_data(fileloc)
    tw = aa['Tw']
    tr = aa['Tr']

    # convert variables in csv file to floats
    #if station == in ['']:
    #    tw = [(0.1*float(str(c).replace(',','.').replace('#N/B','nan'))-32)/1.8 if c else nan for c in tw[:-1]]
    #else:
    tw = [float(str(c).replace(',','.').replace('#N/B','nan').replace('--','nan')) if c else nan for c in tw[:-1]]
    if station in ['leiden', 'haarlem', 'damwoude', 'wageningen']:
        tr = [float(str(c).replace(',','.').replace('#N/B','nan')) if c else nan for c in tr[:-1]]
    else:
        tr = [0.1*float(str(c).replace(',','.').replace('#N/B','nan')) if c else nan for c in tr[:-1]]
    # strings of date and hour
    try:
        datestring = [str(int(day)) for day in aa['YYYYMMDD'][:-1]]
    except KeyError:
        yy = [str(int(yr)) for yr in aa['year'][:-1]]
        mm = [str(int(mn)).zfill(2) for mn in aa['month'][:-1]]
        dd = [str(int(dy)).zfill(2) for dy in aa['day'][:-1]]
        datestring = [yy[idx]+mm[idx]+dd[idx] for idx,c in enumerate(yy)] 
    hourstring = [str(int(hour)).zfill(2) for hour in aa['HH'][:-1]]
    # set HH=24 to HH=0 -> need to add 1 day to datetime object!
    if '24' in hourstring:
        hourstring = [item if item!='24' else '0' for item in
                                    hourstring]
        add_day = True
    else:
        add_day = False
    # combine string in datestring
    dstring = [datestring[idx] + hourstring[idx] for idx in range(0,len(datestring))]
    # convert to datetime object
    dtime = [datetime.strptime(
                str(item), ('%Y%m%d%H')) for item in dstring]
    # now add 1 day to datetime object where we change HH=24 to HH=0
    if add_day:
        # the date of the night -> HH=24 is HH=0 on the next day!
        dtime = [c+ timedelta(days=1) if c.hour==0 else
                                            c for c in dtime]
    # sort time axis for plotting
    idx = argsort(dtime)
    tw = nparray(tw)[idx]
    tr = nparray(tr)[idx]
    dtime= nparray(dtime)[idx]
    return dtime, tw, tr

def remove_nans(dtime, tw, tr):
    dtime = nparray(dtime)[~isnan(tw)]
    tr = nparray(tr)[~isnan(tw)]
    tw = nparray(tw)[~isnan(tw)]
    dtime = dtime[~isnan(tr)]
    tw = tw[~isnan(tr)]    
    tr = tr[~isnan(tr)]
    return dtime, tw, tr

def convert_wund_timezone_to_utc(dtime, tw, tr, station):
    '''
    assume wunderground data is in local time -> convert to UTC
    and find the datetime objects shared with the KNMI UTC data
    '''
    # timezone conversion
    from_zone = tz.gettz('Europe/Amsterdam')
    to_zone = tz.gettz('UTC')
    # tell datetime object that it is in local time
    dtime_local = [d.replace(tzinfo=from_zone) for d in dtime]
    dtime_utc_knmi = [d.replace(tzinfo=to_zone) for d in dtime]  # KNMI data already in UTC
    dtime_utc_wund = [d.astimezone(to_zone) for d in dtime_local]
    print station
    if station in ['ijsselmonde', 'losser', 'leiden', 'delft', 'rotterdam', 'heemskerk']:
        # seems to be 1 hour of when converting to UTC
        print "add 1 hour"
        dtime_utc_wund = [c + timedelta(hours=1) for c in dtime_utc_wund]
        #pass
    # find shared datetime
    shared_datetime_1 = utils.ismember2(dtime_utc_knmi, dtime_utc_wund)
    shared_datetime_2 = utils.ismember2(dtime_utc_wund, dtime_utc_knmi)
    # fix temperature tw
    tw_fixed = [tw[i] for i in shared_datetime_1 if i is not None]
    tr_fixed = [tr[i] for i in shared_datetime_2 if i is not None]
    dtime_utc_knmi_fixed = [dtime_utc_knmi[i] for i in shared_datetime_2 if i is not None]
    return dtime_utc_knmi_fixed, tw_fixed, tr_fixed
    
def calculate_uhi_distribution(dtime, tw, tr):
    '''
    calculate daily max uhi
    '''
    # calculate uhi
    uhi = nparray(tw) - nparray(tr)
    # extract date for each timestep
    dt= [step.date() for step in dtime]
    uniquetime = sort(list(set(dt)))
    uhi_daily = [nmax([uhi[i] for i, j in enumerate(dt) if
                       j==uniquetime[k] ]) for k in range(0,len(uniquetime))]
    return uniquetime, uhi_daily
    
def plot_daily_cycle(dtime, tr, tw, dtime_utc_knmi_fixed, tr_fixed, tw_fixed, station):
    directory = os.path.join('gj', station)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    for month in range(1,13):
        try:
            fig = figure()
            fig.subplots_adjust(hspace=0.5)
            subplot(2,1,1)
            # find first day in June for original time axis
            dmonth = [d.month for d in dtime]
            idx = dmonth.index(month) # first day in June
            plot(dtime[idx:idx+440], tr[idx:idx+440],)
            plot(dtime[idx:idx+440], tw[idx:idx+440], 'r')
            title('original time-axis')
            subplot(2,1,2)
            # find first day in June for corrected time axis
            dmonth = [d.month for d in dtime_utc_knmi_fixed]
            idx = dmonth.index(month) # first day in June
            plot(dtime_utc_knmi_fixed[idx:idx+440], tr_fixed[idx:idx+440])
            plot(dtime_utc_knmi_fixed[idx:idx+440], tw_fixed[idx:idx+440], 'r')
            title('time-axis corrected')
            plt.savefig(directory + '/dailycycle_' + station + '_' + str(month) + '.png')
            plt.close()
        except Exception:
            pass
        
    dailycycle_tw = zeros(23)
    dailycycle_tr = zeros(23)
    dailycycle_tw_fixed = zeros(23)
    dailycycle_tr_fixed = zeros(23)
    dhour = [d.hour for d in dtime]
    dhour2 = [d.hour for d in dtime_utc_knmi_fixed]
    for hr in range(1,24):
        indices = [i for i, x in enumerate(dhour) if x == hr]
        indices2 = [i for i, x in enumerate(dhour2) if x == hr]
        dailycycle_tw[hr-1] = nmean(nparray(tw)[indices])
        dailycycle_tr[hr-1] = nmean(nparray(tr)[indices])
        dailycycle_tw_fixed[hr-1] = nmean(nparray(tw_fixed)[indices2])
        dailycycle_tr_fixed[hr-1] = nmean(nparray(tr_fixed)[indices2])
    fig = figure()
    fig.subplots_adjust(hspace=0.5)
    subplot(2,1,1)
    plot(range(1,24), dailycycle_tw, 'r')
    plot(range(1,24), dailycycle_tr, 'b')
    title('original time-axis')
    subplot(2,1,2)
    plot(range(1,24), dailycycle_tw_fixed, 'r')
    plot(range(1,24), dailycycle_tr_fixed, 'b')
    title('corrected time-axis')
    show()
    
def print_results(uhi_daily, uhi_daily_fixed):
    print "orig"
    print "50th percentile: " + str(nppercentile(uhi_daily, 50))
    print "95th percentile: " + str(nppercentile(uhi_daily, 95))
    print "fixed"
    print "50th percentile: " + str(nppercentile(uhi_daily_fixed, 50))
    print "95th percentile: " + str(nppercentile(uhi_daily_fixed, 95))

if __name__=="__main__":
    stations = ['ijsselmonde', 'apeldoorn', 'losser', 'leiden', 'heerhugowaard',
                'apeldoorn', 'ijsselmuiden']
    stations = ['doornenburg']
    for station in stations:
        print '\n' + station
        main(station)
    
    # correct original: leiden, heerugowaard
    # fixed: apeldoorn, ijsselmuiden, 
    # ijsselmonde: to fix local->UTC and add 1 hour ? (see may plot)
    # (losser?)