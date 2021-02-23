""" calc_chords.py
    
    Processing chord length estimates with ATL07
    Initial code written by Alek Petty (02/01/2020)
    
    Note this uses Basemap for plotting, but the other scripts have been converted to Cartopy (as I slowly make the transition).

    Input:
        ATL07 data

    Output:
        Raw chord length estimates in a netcdf file

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

import utils as ut

ut.reset_matplotlib()


def gen_raw_chord_netcdf(pd_data, out_path, extra_str):
    """ convert Pandas dataframe to xarray dataset to easily save as netcdf file"""

    chord_data = xr.Dataset.from_dataframe(pd_data)

    #chord_data['date-time'].attrs={'units':'datetime64', 'coordinates':'along-track chord','long_name':'mean time of chord acquisition'}
    chord_data['lats'].attrs={'units':'degrees North', 'coordinates':'along-track chord','long_name':'mean chord latitude'}
    chord_data['lons'].attrs={'units':'degrees East', 'coordinates':'along-track chord','long_name':'mean chord longitude'}
    chord_data['beam_number'].attrs={'units':'beam number', 'coordinates':'along-track chord','long_name':'ATLAS beam number (1 3 5 = strong beams)'}
    chord_data['lengths'].attrs={'units':'meters', 'coordinates':'along-track chord','long_name':'along-track chord length from a given ATLAS beam. Data ordered by ascending beam number and then time.'}
              

    chord_data.attrs={'Description':'Along-track raw chord length estimates from ICESat-2 ('+extra_str+' algorithm).',
        'Contact':'Alek Petty (alek.a.petty@nasa.gov)',
        'Reference': 'Petty, A. A., Bagnardi, M., Kurtz, N., Tilling, R., Fons, S., Armitage, T., et al. (2021). Assessment of ICESat‐2 sea ice surface classification with Sentinel‐2 imagery: implications for freeboard and new estimates of lead and floe geometry. Earth and Space Science, 8, e2020EA001491. https://doi.org/10.1029/2020EA001491',
        'Conventions':'CF-1.8',
        'featureType':'point',
        'creation date':"Created " + datetime.today().strftime("%d/%m/%y"),
        'projection':'Data referenced to WGS84 (EPSG:4326)'}

    chord_data.to_netcdf(out_path+'.nc')



##---- Config --------

release='rel003'
#ATL07path = '/sea_ice_pso/aivanoff/atl0710_testrun_rel003_code20200420_segtype2-5_mod_redo/'
ATL07path = '/cooler/I2-ASAS/'+release
figPath='../figures/'+release+'/'
savePath='../data/'+release+'/'
beams = [1, 3, 5] # 1 3 5 should always be strong beams!
print(str(len(beams))+'beams')

version='001'


#----different hemispheres/time periods of analysis
period=3

if (period==2):
    # Arctic start of winter 2018 rel002/003 crossover
    mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)

    hem='nh'
    start_date='20181101'
    end_date = '20190430'
    dateStr=start_date+'-'+end_date
    filesAll = glob(ATL07path+'/ATL07-01/data/ATL07-01_*.h5')
    date_range = pd.date_range(start_date,end_date, freq='d').strftime("%Y%m%d").tolist()
    
    ATL07files = [s for s in filesAll for sub in date_range if sub in s ]

elif (period==3):
    # Arctic start of winter 2018 rel002/003 crossover
    mapProj = Basemap(epsg=3412,resolution='l', llcrnrlon=225, llcrnrlat=-45, urcrnrlon=42, urcrnrlat=-44)

    hem='sh'
    start_date='20190501'
    end_date = '20190930'
    dateStr=start_date+'-'+end_date
    filesAll = glob(ATL07path+'/ATL07-02/data/ATL07-02_*.h5')
    date_range = pd.date_range(start_date,end_date, freq='d').strftime("%Y%m%d").tolist()

    ATL07files = [s for s in filesAll for sub in date_range if sub in s ]



pddata_v1_list=[]
pddata_v2_list=[]


# Loop over files

for h5file in tqdm(ATL07files, leave=False):
    
    input_filename = os.path.basename(h5file)[:-3]
    print(input_filename)
    ATL07 = h5py.File(h5file, 'r')
    
    beamStrs=ut.get_beam_config(ATL07)
    
    # Loop over beams
    for beamNum in beams:
        beamStr=beamStrs[beamNum-1]

        try:
            pd_data = ut.get_atl07_data_beam(ATL07, beamStr)

            if (len(pd_data)>1000):

                # Find chords 
                try:
                    pddata_v1_granule, pddata_v2_granule, _, _, _, _ = ut.get_chords_2versions(pd_data, km_unit=False )
                except:
                    print('error calculating chords') 

                # add beam number info to dataframe (to enable beam filtering)
                pddata_v1_granule['beam_number']=str(beamNum)
                pddata_v2_granule['beam_number']=str(beamNum)
                pddata_v1_list.append(pddata_v1_granule)
                pddata_v2_list.append(pddata_v2_granule)

            else:
                print('Not enough beam data to calculate chord lengths')
        
        except:
           print('No data in ', h5file, ' ', beamStr) 
    
    ATL07.close()



pddata_v1_all = pd.concat(pddata_v1_list,axis=0)
pddata_v2_all = pd.concat(pddata_v2_list,axis=0)

# arrange in order of beam number
pddata_v1_all=pddata_v1_all.sort_values(['beam_number', 'time'], ascending=[True, True])
pddata_v2_all=pddata_v2_all.sort_values(['beam_number', 'time'], ascending=[True, True])

# Reset row indexing
pddata_v1_all=pddata_v1_all.reset_index(drop=True)
pddata_v2_all=pddata_v2_all.reset_index(drop=True)

gen_raw_chord_netcdf(pddata_v1_all, savePath+'IS2_raw_chord_lengths_v1_'+dateStr+'_'+hem+'_'+release[-3:]+'_001', 'v1')
gen_raw_chord_netcdf(pddata_v2_all, savePath+'IS2_raw_chord_lengths_v2_'+dateStr+'_'+hem+'_'+release[-3:]+'_001', 'v2')


print('Total number of sections: ', len(pddata_v1_all.shape))
#print('Total number of sections v2: ', len(lengths07v2))

plot=True
if (plot==True):

    xptsv1, yptsv1=mapProj(pddata_v1_all.lons, pddata_v1_all.lats)
    xptsv2, yptsv2=mapProj(pddata_v2_all.lons, pddata_v2_all.lats)

    xpts=[xptsv1, xptsv2]
    ypts=[yptsv1, yptsv2]

    data=[pddata_v1_all['lengths'], pddata_v2_all['lengths']]

    means=['%0d' %np.nanmean(dataT) for dataT in data]

    minval=[0, 0]
    maxval=[10000, 10000]
    labels=['(a) v1 (ssh)', '(b) v2 (20%)']
    cbarlabels=['Chord length (m)', 'Chord length (m)']
    cmaps=[cm.viridis,  cm.viridis, cm.RdBu_r]


    if (hem=='nh'):
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7, 5))
    else:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7, 3.5))

    for i in range(len(data)):
        ax=axs.flatten()[i]
        plt.sca(ax)
        
        im1=plt.hexbin(xpts[i], ypts[i], C=data[i], vmin=minval[i], vmax=maxval[i], gridsize=(300, 300), cmap=cmaps[i], zorder=2, rasterized=True)
        
        mapProj.drawcoastlines(linewidth=0.25, zorder=5)
        mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
        mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=80, linewidth = 0.25, zorder=10)
        mapProj.fillcontinents(color='0.95',lake_color='0.9', zorder=3)

        if (hem=='nh'):
            cax1 = fig.add_axes([0.35+(i*0.48), 0.13, 0.1, 0.03])
        else:
            cax1 = fig.add_axes([0.14+(i*0.33), 0.48, 0.1, 0.03])

        cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='both')
        #cbar.set_ticks(np.linspace(minval[i], maxval[i], 6) )
        cbar.set_label(cbarlabels[i], labelpad=3)
        cbar.set_ticks([minval[i], maxval[i]])

        ax.annotate(labels[i], xy=(0.97, 0.95), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

        ax.annotate('mean: '+means[i]+' m', xy=(0.97, 0.91), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')


        if (i==0):
            ax.annotate(' '+dateStr, xy=(0.02, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

    plt.tight_layout()
    fig.savefig(figPath+'/section_chordlength_version_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'3raw.png', dpi=300)


