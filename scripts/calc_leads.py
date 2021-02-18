""" calc_leads.py
    
    Processing lead fraction estimates with ATL07
    Initial code written by Alek Petty (02/01/2020)
    
    Input:
        ATL07 data

    Output:
        Lead fraction estimates in a netcdf file

    Python dependencies:
        See below for the relevant module imports. More information on installation is given in the README file.

    Update history:
        02/01/2020: Version 1.
    
"""

import matplotlib, sys
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap

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
import matplotlib.colorbar as mcbar
import cartopy.feature as cfeature

import common_functions as cF

cF.reset_matplotlib()


def get_sections_radioleads(pd_dataT, section_size=10000):
    maxDist=np.int(pd_dataT.along_dist.iloc[-1])
    #print('seg dist of last seg', maxDist, 'max seg dist', np.nanmax(dist))
    lons_sections=[]
    lats_sections=[]
    lead_fractions=[]
    lead_fractions2=[]
    radio_fractions=[]
    radio_ratios=[]
    lead_radio_fractions=[]
    ice_conc_sections=[]
    time_sections=[]

    for x in range(int(pd_dataT.along_dist.iloc[0]), maxDist, section_size):
        idx = where((pd_dataT.along_dist>=x)& (pd_dataT.along_dist<x+section_size))
        #print(np.nanmean(dist[idx]))
        #total_segs=len(pd_dataT.ssh.iloc[idx])
        seglength_section = pd_dataT.iloc[idx].seg_length
        seglength_total = np.sum(pd_dataT.iloc[idx].seg_length)

        lons_sections.append(np.nanmean(pd_dataT.iloc[idx].lons))
        lats_sections.append(np.nanmean(pd_dataT.iloc[idx].lats))
        ice_conc_sections.append(np.nanmean(pd_dataT.iloc[idx].ice_conc))
        time_sections.append(Time(np.nanmean(pd_dataT.iloc[idx].gpsseconds),format='gps').utc.datetime64)
        #print(where(ssh[idx]==1)[0])
        #ssh_segs=len(where(pd_dataT.ssh.iloc[idx]==1)[0])
        
        #print('sum all seg length:', np.sum(seglength_total))
        
        if (len(seglength_section)>1000.):


            seglength_ssh = pd_dataT.iloc[idx].query('ssh_flag >0.5').seg_length
            seglength_radio = pd_dataT.iloc[idx].query('seg_type >1.5').seg_length
            seglength_radiodark = pd_dataT.iloc[idx].query('seg_type >5.5').seg_length
            seglength_radiospecular = pd_dataT.iloc[idx].query('seg_type >1.5 and seg_type < 5.5').seg_length

            height_thresh = ut.calc_leadheight_thresh(pd_dataT, idx)
            #print(str(height_thresh))
            seglength_ssh2 = pd_dataT.iloc[idx].query('seg_type >1.5 and seg_type < 5.5 and height < '+str(height_thresh)).seg_length

            radio_fraction=np.sum(seglength_radio)/seglength_total
            radio_ratio=np.sum(seglength_radiospecular)/np.sum(seglength_radio)
            lead_radio_fraction=np.sum(seglength_ssh)/np.sum(seglength_radiospecular)
            lead_fraction=np.sum(seglength_ssh)/seglength_total
            lead_fraction2=np.sum(seglength_ssh2)/seglength_total

            radio_fractions.append(radio_fraction)
            radio_ratios.append(radio_ratio)
            lead_radio_fractions.append(lead_radio_fraction)
            lead_fractions.append(lead_fraction)
            lead_fractions2.append(lead_fraction2)

        else:
            #print('no segs')
            
            radio_fractions.append(np.nan)
            radio_ratios.append(np.nan)
            lead_radio_fractions.append(np.nan)
            lead_fractions.append(np.nan)
            lead_fractions2.append(np.nan)
            
    return lons_sections, lats_sections, radio_fractions, radio_ratios, lead_radio_fractions, lead_fractions, lead_fractions2, ice_conc_sections, time_sections


# Map projection
#mapProj = Basemap(projection='npstere',boundinglat=60,lon_0=0, resolution='l' , round=False)
mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)

## ICESat-2 ATL10 data directory
release='rel003'
#ATL07path = '/sea_ice_pso/aivanoff/atl0710_testrun_rel003_code20200420_segtype2-5_mod_redo/'
ATL07path = '/cooler/I2-ASAS/'+release
figPath='./figures/'+release+'/'
savePath='/sea_ice_pso/aapetty/classification_data/'
beams = [1, 3, 5] # 1 3 5 should always be strong beams!
print(str(len(beams))+'beams')
#beams = ['gt1l']

version='v2'
period=14



if (period==14):
    # Arctic start of winter 2018 rel002/003 crossover
    mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)

    hem='nh'
    start_date='20191001'
    end_date = '20191102'
    dateStr=start_date+'-'+end_date
    filesAll = glob(ATL07path+'/ATL07-01/data/ATL07-01_*.h5')
    date_range = pd.date_range(start_date,end_date, freq='d').strftime("%Y%m%d").tolist()
    
    ATL07files = [s for s in filesAll for sub in date_range if sub in s ]

elif (period==15):
    # Arctic start of winter 2018 rel002/003 crossover
    mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)

    hem='nh'
    start_date='20181101'
    end_date = '20190430'
    dateStr=start_date+'-'+end_date
    filesAll = glob(ATL07path+'/ATL07-01/data/ATL07-01_*.h5')
    date_range = pd.date_range(start_date,end_date, freq='d').strftime("%Y%m%d").tolist()
    
    ATL07files = [s for s in filesAll for sub in date_range if sub in s ]

elif (period==16):
    # Arctic start of winter 2018 rel002/003 crossover
    mapProj = Basemap(epsg=3412,resolution='l', llcrnrlon=225, llcrnrlat=-45, urcrnrlon=42, urcrnrlat=-44)

    hem='sh'
    start_date='20190501'
    end_date = '20190930'
    dateStr=start_date+'-'+end_date
    filesAll = glob(ATL07path+'/ATL07-02/data/ATL07-02_*.h5')
    date_range = pd.date_range(start_date,end_date, freq='d').strftime("%Y%m%d").tolist()
    
    ATL07files = [s for s in filesAll for sub in date_range if sub in s ]


print(ATL07files)
print('Number of files:', len(ATL07files))

lons_sections = np.array([])
lats_sections = np.array([])
lead_fractions = np.array([])
lead_fractions2 = np.array([])
radio_fractions = np.array([])
radio_ratios = np.array([])
lead_specular_fractions = np.array([])
ice_conc_sections = np.array([])
time_sections = np.array([])

def gen_netcdf(pd_data, out_path, extra_str):
    """ convert Pandas dataframe to xarray dataset to easily save as netcdf file"""

    chord_data = xr.Dataset.from_dataframe(pd_data)

    chord_data['date-time'].attrs={'units':'datetime64', 'coordinates':'along-track 10 km swath section','long_name':'mean datetime of swath section acquisition'}
    chord_data['lats'].attrs={'units':'degrees North', 'coordinates':'along-track 10 km swath section','long_name':'mean chord latitude'}
    chord_data['lons'].attrs={'units':'degrees East', 'coordinates':'along-track 10 km swath section','long_name':'mean chord longitude'}
    chord_data['lead_fraction_v1'].attrs={'units':'fraction (0 to 1)', 'coordinates':'along-track 10 km swath section','long_name':'Linear lead fraction within a 10 km section from the three strong ATLAS beams (version 1 algorithm)'}
    chord_data['lead_fraction_v2'].attrs={'units':'fraction (0 to 1)', 'coordinates':'along-track 10 km swath section','long_name':'Linear lead fraction within a 10 km section from the three strong ATLAS beams (version 2 algorithm)'}
              

    chord_data.attrs={'Description':'Along-track lead fraction estimates from ICESat-2. Mean lead fraction calculated from the three strong beams in 10 km along-track swath scetions using two SSH metrics (v1/v2).',
        'Contact':'Alek Petty (alek.a.petty@nasa.gov)',
        'Reference': 'Petty, A. A., Bagnardi, M., Kurtz, N., Tilling, R., Fons, S., Armitage, T., et al. (2021). Assessment of ICESat‐2 sea ice surface classification with Sentinel‐2 imagery: implications for freeboard and new estimates of lead and floe geometry. Earth and Space Science, 8, e2020EA001491. https://doi.org/10.1029/2020EA001491',
        'Conventions':'CF-1.8',
        'featureType':'point',
        'creation date':"Created " + datetime.today().strftime("%d/%m/%y"),
        'projection':'Data referenced to WGS84 (EPSG:4326)'}

    chord_data.to_netcdf(out_path+'.nc') 
                

# Loop over ATL07 files
for h5file in tqdm(ATL07files, leave=False):
    
    input_filename = os.path.basename(h5file)[:-3]
    print(input_filename)
    ATL07 = h5py.File(h5file, 'r')
    
    beamStrs=rF.get_beam_config(ATL07)
    
    # Loop over beams
    beam_dataframes=[]
    for beamNum in beams:
        beamStr=beamStrs[beamNum-1]
        print(beamStr)
        try:
            # Coordinates and height
            beam_dataframes.append(ut.get_data_beam(ATL07, beamStr))
        except:
            print('No beam data in ', h5file, ' ', beamStr)
    
    try:
        # Combine data frames across beams
        pd_data=pd.concat(beam_dataframes)
        pd_data=pd_data.sort_values(by=['along_dist'])
    
        #print(pd_data.head(20))

        ATL07.close()
        
        if (len(pd_data)>1000):
            # Floe processing
            lons_section, lats_section, radio_fraction, radio_ratio, lead_specular_fraction, lead_fraction, lead_fraction2, ice_conc_section = get_sections_radioleads(pd_data)

            lons_sections = np.hstack([lons_sections, lons_section])
            lats_sections = np.hstack([lats_sections, lats_section])
            radio_fractions = np.hstack([radio_fractions, radio_fraction])
            radio_ratios = np.hstack([radio_ratios, radio_ratio])
            lead_specular_fractions = np.hstack([lead_specular_fractions, lead_specular_fraction])
            lead_fractions = np.hstack([lead_fractions, lead_fraction])
            lead_fractions2 = np.hstack([lead_fractions2, lead_fraction2])
            ice_conc_sections = np.hstack([ice_conc_sections, ice_conc_section])
            time_sections = np.hstack([time_sections, time_section])
        else:
            print('Dataframe too small')

    except:
        print('No data in file')



#print('Total number of ICESat-2 segments: ', len(lons))
#lons_sections.dump(savePath+'section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'lons_sections')
#lats_sections.dump(savePath+'section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'lats_sections')
#radio_fractions.dump(savePath+'section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'radio_fractions')
#radio_ratios.dump(savePath+'section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'radio_ratios')
#lead_specular_fractions.dump(savePath+'section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'lead_specular_fractions')
#lead_fractions.dump(savePath+'section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'leadfractions')
#lead_fractions2.dump(savePath+'section_leadtypefraction2_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'leadfractions2')
#ice_conc_sections.dump(savePath+'section_leadtypefraction2_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'ice_conc_sections')


pd_data = pd.DataFrame({'lead_fraction_v1': lead_fractions,'lead_fraction_v2': lead_fractions2, 'lons': lons_sections, 'lats':lats_sections, 'date-time': time_sections})


data=[radio_fractions*100., radio_ratios*100, lead_specular_fractions*100., lead_fractions*100.]
means=['%.01f' %np.nanmean(dataT) for dataT in data]
#stds=['%02f' %np.nanstd(dataT) for dataT in data]

minval=[0, 0, 0, 0]
maxval=[5, 100, 100, 5]
labels=['radio leads:all segments (%)', 'specular:radio leads (%)', 'ssh:specular leads (%)', 'ssh leads:all segments (%)']
cmaps=[cm.YlOrRd,  cm.RdYlBu_r, cm.plasma_r, cm.YlOrRd]

xpts_sections, ypts_sections=mapProj(lons_sections, lats_sections)
    

if (hem=='nh'):
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(11, 5))
    plt.subplots_adjust(bottom=0.01, left=0.02, top=0.95, right=0.98)
    for i in range(len(data)):
        ax=axs.flatten()[i]
        sca(ax)
        
        im1=hexbin(xpts_sections, ypts_sections, C=data[i], vmin=minval[i], vmax=maxval[i], gridsize=2000, cmap=cmaps[i], zorder=2, rasterized=True)
        
        mapProj.drawcoastlines(linewidth=0.25, zorder=5)
        mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
        mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=80, linewidth = 0.25, zorder=10)
        mapProj.fillcontinents(color='0.95',lake_color='0.9', zorder=3)


        #cax1 = fig.add_axes([0.13+(i*0.25), 0.09, 0.1, 0.03])
        cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.01,shrink=0.9)
        cbar=fig.colorbar(im1,cax=cax,extend='both',**kw)

        #cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='both')
        cbar.set_ticks(np.linspace(minval[i], maxval[i], 5) )
        cbar.set_label(labels[i], labelpad=3)
        #cbar.set_ticks(np.arange(vmin, vmax+0.001, 0.4))

        ax.annotate('mean: '+means[i]+' %', xy=(0.97, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

        if (i==0):
            ax.annotate(' '+dateStr, xy=(0.02, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

    #plt.tight_layout()
    fig.savefig(figPath+'/section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'4.png', dpi=300)


if (hem=='sh'):
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(11, 3.1))
    plt.subplots_adjust(bottom=0.02, left=0.02, top=0.95, right=0.98)
    for i in range(len(data)):
        ax=axs.flatten()[i]
        sca(ax)
        
        im1=hexbin(xpts_sections, ypts_sections, C=data[i], vmin=minval[i], vmax=maxval[i], gridsize=2000, cmap=cmaps[i], zorder=2, rasterized=True)
        
        mapProj.drawcoastlines(linewidth=0.25, zorder=5)
        mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
        mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=80, linewidth = 0.25, zorder=10)
        mapProj.fillcontinents(color='0.95',lake_color='0.9', zorder=3)

        #cax1 = fig.add_axes([0.13+(i*0.25), 0.09, 0.1, 0.03])
        cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.01,shrink=0.9)
        cbar=fig.colorbar(im1,cax=cax,extend='both',**kw)

        cbar.set_ticks(np.linspace(minval[i], maxval[i], 5) )
        cbar.set_label(labels[i], labelpad=3)
        #cbar.set_ticks(np.arange(vmin, vmax+0.001, 0.4))

        ax.annotate('mean: '+means[i]+' %', xy=(0.97, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

        if (i==0):
            ax.annotate(' '+dateStr, xy=(0.02, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

    fig.savefig(figPath+'/section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'4.png', dpi=300)




plot_section_map=False
if plot_section_map:
    
    xpts_sections, ypts_sections=mapProj(lons_sections, lats_sections)
    
    vmin=0
    vmax=1.2

    if (hem=='nh'):
        fig=plt.figure(figsize=(4., 6))
    else:
        fig=plt.figure(figsize=(4., 4))
    ax1=gca()
    im1=hexbin(xpts_sections, ypts_sections, C=lead_fractions*100, gridsize=2000, vmin=vmin, vmax=vmax, cmap=cm.viridis, zorder=2, rasterized=True)
    #im1=mapProj.scatter(df.xpts, df.ypts, c=df.ssh, s=100, cmap=cm.viridis, zorder=2, rasterized=True)
    mapProj.drawcoastlines(linewidth=0.25, zorder=5)
    mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
    mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=80, linewidth = 0.25, zorder=10)
    mapProj.fillcontinents(color='0.95',lake_color='0.9', zorder=3)
    #im11 = mapProj.contour(xptsc , yptsc, iceconc,levels=0.15, colors='k', linewidths=0.8, zorder=5, alpha=1)


    if ((period<3) | (period>4)):
        cax1 = fig.add_axes([0.65, 0.95, 0.3, 0.03])
    else:
        cax1 = fig.add_axes([0.35, 0.47, 0.35, 0.04])
        
    cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='max')
    cbar.set_label('Lead fraction (%)', labelpad=3)
    cbar.set_ticks(np.arange(vmin, vmax+0.001, 0.4))
    

    ax1.annotate(' '+dateStr, xy=(0.02, 0.97), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

    subplots_adjust(bottom=0.12, left=0.02, top = 0.98, right=0.98)

    plt.tight_layout()
    fig.savefig(figPath+'/section_leadtypefraction_ATL07'+release+hem+dateStr+str(len(beams))+'beams'+version+'p'+str(period)+'test.png', dpi=300)
    #fig.savefig(figPath+'/map1_ATL07'+varStr+dateStr+'ShotData.pdf', dpi=400)
