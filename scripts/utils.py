""" utils.py
    
    Random functions used by the various scripts and notebooks in this repo
    Initial code written by Alek Petty and Marco Bagnardi (02/01/2020) 

    Python dependencies:
        See below for the relevant module imports. More information on installation is given in the README file.

    Update history:
        02/01/2020: Version 1.
"""

import matplotlib, sys
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
import h5py
from glob import glob
#import seaborn as sns
from tqdm import tqdm
from astropy.time import Time
from datetime import datetime
import pyproj

def get_atl07_data_beam(ATL07, beamStr):
    """ Gran ATL07 data and put in pandas dataframer
    
    Args:
        ATL07 (hdf5 file): ATL07 file
        beamStr (str): ATLAS beam 
    
    Returns:
        pd_data: pandas dataframe

    Caveats:
        need to add data gap condtion but not sure how..

    """

    # grab location data and along-track distance
    lon_this = ATL07[beamStr + '/sea_ice_segments/longitude'][:]
    lat_this = ATL07[beamStr + '/sea_ice_segments/latitude'][:]
    dist_this = ATL07[beamStr + '/sea_ice_segments/seg_dist_x'][:]
    # Height
    height_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_height'][:]
    # Seg Type
    seg_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_type'][:]
    # SSH_Flag
    ssh_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_ssh_flag'][:]
    # segment length
    seglength_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_length_seg'][:]
    # ice_conc
    ice_conc = ATL07[beamStr + '/sea_ice_segments/stats/ice_conc'][:]
    # delta time
    delta_time_this=ATL07[beamStr + '/sea_ice_segments/delta_time'][:]

    # #Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL07['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # subtract leap seconds offset (will change to 38 soon, watch out)
    gpsseconds = atlas_epoch + delta_time_this - 37.

    ## Use astropy to convert GPS time to UTC time
    #tiso=Time(gps_seconds,format='gps').utc.datetime

     # seg stats
    seg_h_mean_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_mean_h'][:]
    seg_h_width_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_w'][:]
    seg_surf_err_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_surface_error_est'][:]

    # Assign to pandas dataframe
    pd_data = pd.DataFrame({'lons': lon_this, 'lats': lat_this, 'along_dist':dist_this, 'seg_length': seglength_this, 
            'seg_type':seg_this, 'height':height_this, 'ssh_flag':ssh_this, 'seg_h_mean':seg_h_mean_this, 
            'seg_h_width':seg_h_width_this, 'seg_surf_err':seg_surf_err_this, 'gpsseconds':gpsseconds, 'ice_conc':ice_conc})
    

    print('seg limits:', np.amin(pd_data['seg_length'].values), np.amax(pd_data['seg_length'].values))
    print('lon limits:', np.amin(pd_data['lons'].values), np.amax(pd_data['lons'].values))
    print('lat limits:', np.amin(pd_data['lats'].values), np.amax(pd_data['lats'].values))

    # Filter data
    pd_data= pd_data.dropna() 
    pd_data = pd_data[pd_data['seg_length']<200]
    pd_data = pd_data[pd_data['seg_length']>2]
    pd_data = pd_data[pd_data['height']<1e6]
    pd_data = pd_data[pd_data['seg_h_mean']<1e6]

    return pd_data

def get_atl07_data_beam_extra(ATL07, beamStr):
    """ Gran ATL07 data and put in pandas dataframe, add some extra info for diagnostics
    
    Args:
        ATL07 (hdf5 file): ATL07 file
        beamStr (str): ATLAS beam 
    
    Returns:
        pd_data: pandas dataframe

    Caveats:
        need to add data gap condtion but not sure how..

    """

    # grab location data and along-track distance
    lon_this = ATL07[beamStr + '/sea_ice_segments/longitude'][:]
    lat_this = ATL07[beamStr + '/sea_ice_segments/latitude'][:]
    dist_this = ATL07[beamStr + '/sea_ice_segments/seg_dist_x'][:]
    # Height
    height_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_height'][:]
    # Seg Type
    seg_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_type'][:]
    # SSH_Flag
    ssh_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_ssh_flag'][:]
    # segment length
    seglength_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_length_seg'][:]
    # photon rate
    photon_rate = ATL07[beamStr + '/sea_ice_segments/stats/photon_rate'][:]
    # background rate
    bground_rate = ATL07[beamStr + '/sea_ice_segments/stats/backgr_r_200'][:]

    # ice_conc
    ice_conc = ATL07[beamStr + '/sea_ice_segments/stats/ice_conc'][:]
    # delta time
    delta_time_this=ATL07[beamStr + '/sea_ice_segments/delta_time'][:]

    # #Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL07['/ancillary_data/atlas_sdp_gps_epoch'][0] 

    # subtract leap seconds offset (will change to 38 soon, watch out)
    gpsseconds = atlas_epoch + delta_time_this - 37.

    ## Use astropy to convert GPS time to UTC time
    #tiso=Time(gps_seconds,format='gps').utc.datetime

     # seg stats
    seg_h_mean_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_mean_h'][:]
    seg_h_width_this = ATL07[beamStr + '/sea_ice_segments/stats/hist_w'][:]
    seg_surf_err_this = ATL07[beamStr + '/sea_ice_segments/heights/height_segment_surface_error_est'][:]

    # Assign to pandas dataframe
    pd_data = pd.DataFrame({'lons': lon_this, 'lats': lat_this, 'along_dist':dist_this, 'seg_length': seglength_this, 
            'seg_type':seg_this, 'height':height_this, 'ssh_flag':ssh_this, 'seg_h_mean':seg_h_mean_this, 
            'seg_h_width':seg_h_width_this, 'seg_surf_err':seg_surf_err_this, 'gpsseconds':gpsseconds, 'ice_conc':ice_conc, 
            'photon_rate':photon_rate, 'bground_rate':bground_rate})
    

    print('seg limits:', np.amin(pd_data['seg_length'].values), np.amax(pd_data['seg_length'].values))
    print('lon limits:', np.amin(pd_data['lons'].values), np.amax(pd_data['lons'].values))
    print('lat limits:', np.amin(pd_data['lats'].values), np.amax(pd_data['lats'].values))

    # Filter data
    pd_data= pd_data.dropna() 
    pd_data = pd_data[pd_data['seg_length']<200]
    pd_data = pd_data[pd_data['seg_length']>2]
    pd_data = pd_data[pd_data['height']<1e6]
    pd_data = pd_data[pd_data['seg_h_mean']<1e6]

    return pd_data


def calc_leadheight_thresh(data, idxT, height_var='height', seg_h_width_var='seg_h_width', seg_h_mean_var='seg_h_mean', seg_surf_err_var='seg_surf_err',percentile_thresh=20):
    """ Calculate the local lead height
    
    Args:
        data (pandas dataframe): original ATL07 data
        idxT (list): section index
        percentile_thresh (int): snow density variable of choosing
    
    Returns:
        lead_htmax: lead height

    Caveats:
        need to add data gap condtion but not sure how..

    """

    height_thresh = np.percentile(data.iloc[idxT][height_var].values, percentile_thresh)
    #print(height_thresh)
    
    idx_emin=np.where(data.iloc[idxT][seg_h_width_var].values<0.13)
    #print(idx_emin)
    if np.size(idx_emin)>0:
        emin=np.nanmin(data.iloc[idxT][seg_h_mean_var].values[idx_emin])
    else:
        emin=0.13
    

    height_error = (2.*data.iloc[idxT][seg_surf_err_var]) + emin
    height_error=height_error[np.isfinite(height_error)]
    height_error_max=np.nanmax(height_error)
    lead_htmax = max(height_error_max, height_thresh)
    #print('lead htmax:', lead_htmax, height_error_max, height_thresh)
    return lead_htmax


def calc_chords(leads, pd_data, version_str, minsize=3, minlength=60., maxdiff=300, max_chord_length=50000, footprint_size=11):
    """ Finding chord lengths
    
     Args:
        leads (numpy array): array of along track distance with lead segments set to negative numbers
        pd_data (dataframe): pandas dataframe of ATL07 data
        minsize (int): minimum number of segments for a valid chord
        minlength (int): minimum length for a valid chord
        maxdiff (int): maximum segment gap within a chord
        max_chord_length (int): maximum chord length (set to nan if higher than this)
        footprint_size (int): footprint size, add this to the chord length estimate (half each side)

    To do:
        Since switching to pandas it makes sense to use groupby to group the floes instead of this ad-hoc attempt.
    
    """

    splitdata=np.split(leads, np.where(leads< -100)[0])

    splitdata_lon=np.split(pd_data.lons, np.where(leads< -100)[0])
    splitdata_lat=np.split(pd_data.lats, np.where(leads< -100)[0])
    splitdata_time=np.split(pd_data.gpsseconds, np.where(leads< -100)[0])


    lengths=[]
    groups=[]
    lons=[]
    lats=[]
    time=[]

    #print('first group from granule', splitdata[0])
    print('processing chords...')
    for x in np.arange(np.size(splitdata)):
        # start at index 1 as index 0 is the lead
        # If it's a chord (i.e. positive values)
        if np.mean(splitdata[x][1:])>-900:
            # If greater that minimum number of segments
            if np.size(splitdata[x][1:]) >= minsize:
                #If greater than minimum length 
                 if (splitdata[x][-1]-splitdata[x][1]) >= minlength:
                    #If no segment gap greater than max gap 
                    if (np.max(np.diff(splitdata[x][1:]))<maxdiff):
                        
                        groups.append(splitdata[x][1:])
                        lengths.append(splitdata[x][-1]-splitdata[x][1]+footprint_size)
                        lons.append(np.mean(splitdata_lon[x][1:]))
                        lats.append(np.mean(splitdata_lat[x][1:]))
                        time.append(np.mean(splitdata_time[x][1:]))
                    else:
                        print('dropped floe group as max segment separation = ', np.max(np.diff(splitdata[x][1:])))
                        
    #print('last group from granule', splitdata[x])

    pd_data_chords = pd.DataFrame({'lengths': lengths, 'lons': lons, 'lats':lats, 'time': time})
    #print(pd_data_chords)

    # nan chord lenghts > 50 km
    pd_data_chords.loc[pd_data_chords['lengths']>max_chord_length, ['lengths']] = np.nan

    return pd_data_chords, groups



def get_chords_2versions(data, km_unit=True, lon_var = 'lons', lat_var='lats', alongtrack_var='along_dist', sshvar = 'ssh_flag', segtype_var = 'seg_type', height_var='height', percentileThresh=20, percentile_section_size=10000):
    """ Calculate two versions of chord lengths from ATL07

    Use the SSH flag and also a variable height filter to seperate and calculate chord lengths
    
    Args:
        data (pandas dataframe): data
        percentileThresh (int): local perncetile height threshold
        percentile_section_size (int): size of sections for calculating local height threshold (in meters)
    
    Returns:
        version 1 and version 2 chord length data
    
    Caveats:
        should probably split into version 1 and 2
    """

    if km_unit:
        data=data.copy()
        # convert along track to meters
        data[alongtrack_var] = data[alongtrack_var].multiply(1000.)

    
    lead_indexv1=np.copy(data[alongtrack_var]) 
    lead_indexv2=np.copy(data[alongtrack_var]) 
    #print('along min/max:', np.amin(lead_indexv2), np.amax(lead_indexv2))

    lead_indexv1[np.where(data[sshvar]>0.5)]=-999
    
    # loop over 10 km sections to do the percentile height /seg type classification as in ATL10 processing
    for x in range(int(data[alongtrack_var].iloc[0]), np.int(data[alongtrack_var].iloc[-1]), percentile_section_size):
        #print(x, percentile_section_size)
        #print('x:', x)
        idxT=np.where((lead_indexv2>x)&(lead_indexv2<(x+percentile_section_size)))[0]
        #print('idx:', idxT)
        if (np.size(idxT)>10):
            # require at least 10 segments

            height_thresh = calc_leadheight_thresh(data, idxT, percentile_thresh=percentileThresh)
            #print(np.percentile(data[height_var].iloc[idxT], percentileThresh))
            idxM=np.where((data[segtype_var].iloc[idxT]>1.5 ) & (data[segtype_var].iloc[idxT]<5.5 ) & (data[height_var].iloc[idxT]<height_thresh))
            #print(idxT[0][0]+idxM)
            lead_indexv2[idxT[0]+idxM]=-999

        else:
            continue
    
    # Find raw floe lengths
    pddata_v1, groupsv1 = calc_chords(lead_indexv1, data, 'v1')
    pddata_v2, groupsv2 = calc_chords(lead_indexv2, data, 'v2')
            
    return pddata_v1, pddata_v2, lead_indexv1, lead_indexv2, groupsv1, groupsv2


def reset_matplotlib():
    """ Reset matplotlib to a common default. """
    
    # Force agg backend.
    plt.switch_backend('agg')
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    plt.rcParams['ytick.major.size'] = 2
    plt.rcParams['axes.linewidth'] = .6
    plt.rcParams['lines.linewidth'] = .6
    plt.rcParams['patch.linewidth'] = .6
    plt.rcParams['ytick.labelsize']=8
    plt.rcParams['xtick.labelsize']=8
    plt.rcParams['legend.fontsize']=9
    plt.rcParams['font.size']=9
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def get_beam_config(ATL07file):
    """ Get the beam string configuration for the spacecraft orientation of the given file"""
    
    orientation_flag=ATL07file['orbit_info']['sc_orient'][0]
    print('orientation flag:', orientation_flag)
    
    if (orientation_flag==0):
        print('Backward orientation')
        beamStrs=['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
                
    elif (orientation_flag==1):
        print('Forward orientation')
        beamStrs=['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
        
    elif (orientation_flag==2):
        print('Transitioning, do not use for science!')
    
    return beamStrs

def create_grid_nsidc_arctic(epsg_string='3411', nx=304, ny=448, leftx=-3837500, dxRes=25000, uppery=5837500, dyRes=25000):
    """ Use pyproj to generate a grid covering the given domain (defined by the projection and the corner lat/lons)"""

    crs = pyproj.CRS.from_string("epsg:"+epsg_string)
    p=pyproj.Proj(crs)
    
    print(dxRes, dyRes)

    x=leftx+dxRes*np.indices((ny,nx),np.float32)[1]
    y=uppery-dxRes*np.indices((ny,nx),np.float32)[0]
    lons, lats = p(x, y, inverse=True)
    
    return x, y, lats, lons, p


def create_grid_nsidc_antarctic(epsg_string='3412', nx=316, ny=332, leftx=-3937500, dxRes=25000, uppery=4337500, dyRes=25000):
    """ Use pyproj to generate a grid covering the given domain (defined by the projection and the corner lat/lons)"""

    crs = pyproj.CRS.from_string("epsg:"+epsg_string)
    p=pyproj.Proj(crs)
    
    print(dxRes, dyRes)

    x=leftx+dxRes*np.indices((ny,nx),np.float32)[1]
    y=uppery-dxRes*np.indices((ny,nx),np.float32)[0]
    lons, lats = p(x, y, inverse=True)
    
    return x, y, lats, lons, p


def bin_data_nsidc(xpts, ypts, var, xptsG, yptsG, dx):
    """ Bin data using numpy histogram 2d
    
    Adapted for the NSIDC grid which has orgin in the top left corner.

    """
    xptsG2=np.flipud(xptsG)
    yptsG2=np.flipud(yptsG)

    xbins=xptsG2[0]-(dx/2)
    ybins=yptsG2[:, 0]-(dx/2)
    # add on one more bin in each direction
    xbins2=np.append(xbins, xbins[-1]+dx)
    ybins2=np.append(ybins, ybins[-1]+dx)
    print('binning..')
    print(xbins2.shape)
    print(ybins2.shape)
    z, _, _ = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2), weights=var)
    counts, _, _ = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2))

    varG = z / counts
    
    # Need to flip the arrayback as we did this to make the y axis go from negative to positive, then we need to transpose because of weirdness with histogram2d
    varG=np.flipud(varG.T)
    counts=np.flipud(counts.T)

    return varG, counts

def get_region_mask(ancDataPath, mplot, xypts_return=0):
    header = 300
    datatype='uint8'
    file_mask = ancDataPath+'/region_n.msk'
    
    #1 Non-regional ocean  
    #2 Sea of Okhotsk 
    #3 Bering Sea  
    #4 Hudson Bay 
    #5 Baffin Bay/Davis Strait/Labrador Sea    
    #6 Greenland Sea   Bellingshausen 
    #7 Kara and Barents Seas

    #8 - Arctic Ocean
    #9 - Canadian Archipelago
    #10 - Gulf of St Lawrence
    #11 - Land

    fd = open(file_mask, 'rb')
    region_mask = np.fromfile(file=fd, dtype=datatype)
    region_mask = np.reshape(region_mask[header:], [448, 304])

    if (xypts_return==1):
        mask_latf = open(ancDataPath+'/psn25lats_v3.dat', 'rb')
        mask_lonf = open(ancDataPath+'/psn25lons_v3.dat', 'rb')
        lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
        lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

        xpts, ypts = mplot(lons_mask, lats_mask)

        return region_mask, xpts, ypts
    else:
        return region_mask

def get_region_mask_s(ancDataPath, mplot, xypts_return=0):
    header = 300
    datatype='uint8'
    file_mask = ancDataPath+'/region_s.msk'
    

    #2 Weddell Sea
    #3 Indian Ocean  
    #4 Pacific Ocean
    #5 Ross Sea
    #6 Bellingshausen Amundsen Sea

    #11 - Land
    #12 - Coast

    fd = open(file_mask, 'rb')
    region_mask = np.fromfile(file=fd, dtype=datatype)
    region_mask = np.reshape(region_mask[header:], [332, 316])

    if (xypts_return==1):
        mask_latf = open(ancDataPath+'/pss25lats_v3.dat', 'rb')
        mask_lonf = open(ancDataPath+'/pss25lons_v3.dat', 'rb')
        lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [332, 316])
        lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [332, 316])

        xpts, ypts = mplot(lons_mask, lats_mask)

        return region_mask, xpts, ypts
    else:
        return region_mask
