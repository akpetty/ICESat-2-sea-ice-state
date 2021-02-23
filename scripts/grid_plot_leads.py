""" grid_leads.py
    
    Gridding and plotting the lead fraction estimates
    Initial code written by Alek Petty (02/01/2020)
    
    Input:
        Swath section lead fraction estimates

    Output:
        Gridded lead fraction estimates in a netcdf file, gridded plots

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
import matplotlib.colorbar as mcbar
import pyproj
import utils as ut
import netCDF4 as nc4
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ut.reset_matplotlib()

def gen_netcdf(savePathT, xpts2d, ypts2d, lons2d, lats2d, lead_fraction1, lead_fraction2, epsg):
    f = nc4.Dataset(savePathT+'IS2_gridded_lead_fraction_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.nc','w', format='NETCDF4') #'w' stands for write

    print ('dimensions:', lons2d.shape[0], lons2d.shape[1])
    f.createDimension('x', lons2d.shape[0])
    f.createDimension('y', lons2d.shape[1])

    longitude = f.createVariable('longitude', 'f4', ('x', 'y'))
    latitude = f.createVariable('latitude', 'f4', ('x', 'y'))  

    xgrid = f.createVariable('xgrid', 'f4', ('x', 'y'))
    ygrid = f.createVariable('ygrid', 'f4', ('x', 'y')) 

    # Assume all three main variables
    lead_fraction_v1 = f.createVariable('lead_fractionv1', 'f4', ('x', 'y'))
    lead_fraction_v2 = f.createVariable('lead_fractionv2', 'f4', ('x', 'y'))

    longitude[:] = lons2d #The "[:]" at the end of the variable instance is necessary
    latitude[:] = lats2d

    xgrid[:] = xpts2d #The "[:]" at the end of the variable instance is necessary
    ygrid[:] = ypts2d

    lead_fraction_v1[:] = np.around(lead_fraction1, decimals=4)
    lead_fraction_v2[:] = np.around(lead_fraction2, decimals=4)
    
    longitude.units = 'degrees East'
    longitude.long_name='longitude'

    xgrid.units = 'meters'
    xgrid.long_name='projection grid x values'
    xgrid:description = "center values of projection grid in x direction" 

    ygrid.units = 'meters'
    ygrid.long_name='projection grid y values'
    ygrid:description = "center values of projection grid in y direction" 

    latitude.units = 'degrees North'
    latitude.long_name='latitude'

    lead_fraction_v1.description = 'Mean lead fraction (v1 lead methodology) within 25 km x 25 km grid-cells.'
    lead_fraction_v1.units = '%'
    lead_fraction_v1.long_name ='Lead fraction (v1 method)'

    lead_fraction_v2.description = 'Mean lead fraction (v2 lead methodology) within 25 km x 25 km grid-cells.'
    lead_fraction_v2.units = '%'
    lead_fraction_v2.long_name ='Lead fraction (v2 method)'

    #Add global attributes
    #Add global attributes
    f.description = "Gridded (swath) lead fraction estimates from ICESat-2. Gridded data calculated from 10 km swath section data derived from combining data from three strong beams. Two sea surface metrics are used to drive lead fraction, both variables are included in this dataset (v1 and v2). Data gridded onto the NSIDC Polar Stereographic grid using a simple binning procedure (EPSG:"+epsg+")."
    f.contact = "Alek Petty (alek.a.petty@nasa.gov)"
    f.reference = "'Petty, A. A., Bagnardi, M., Kurtz, N., Tilling, R., Fons, S., Armitage, T., et al. (2021). Assessment of ICESat‐2 sea ice surface classification with Sentinel‐2 imagery: implications for freeboard and new estimates of lead and floe geometry. Earth and Space Science, 8, e2020EA001491. https://doi.org/10.1029/2020EA001491'"
    f.history = "Created " + datetime.today().strftime("%d/%m/%y")
    #f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

    f.close()

## ICESat-2 ATL10 data directory
release='rel003'

figPath='../figures/'+release+'/maps/'
savePath='../data/'+release+'/'
anc_data_path='/sea_ice_pso/aapetty/raw_data/OTHER/'
beams = [1, 3, 5] # 1 3 5 should always be strong beams!
print(str(len(beams))+'beams')

period=2
dx=25000


if (period==2):  
    #mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)
    epsg='3411'
    mapProj = pyproj.Proj("+init=EPSG:"+epsg)
    hem='nh'
    start_date='20181101'
    end_date = '20190430'
    dateStr=start_date+'-'+end_date
    #gridbins=(175, 350)
    xptsG, yptsG, latG, lonG, proj= ut.create_grid_nsidc_arctic()
    region_mask, xptsR, yptsR=ut.get_region_mask(anc_data_path, mapProj, xypts_return=1)

elif (period==3):
    # Arctic start of winter 2018 rel002/003 crossover
    #mapProj = Basemap(epsg=3412,resolution='l', llcrnrlon=225, llcrnrlat=-45, urcrnrlon=42, urcrnrlat=-44)
    epsg='3412'
    mapProj = pyproj.Proj("+init=EPSG:"+epsg)
    hem='sh'
    start_date='20190501'
    end_date = '20190930'
    dateStr=start_date+'-'+end_date
    #gridbins=(300, 300)
    xptsG, yptsG, latG, lonG, proj= ut.create_grid_nsidc_antarctic()
    region_mask, xptsR, yptsR=ut.get_region_mask_s(anc_data_path, mapProj, xypts_return=1)


leads = xr.open_dataset(savePath+'IS2_swath_section_lead_stats_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.nc')


leads = leads.dropna(dim='index')


xv1, yv1 = mapProj(leads.lons.values, leads.lats.values)
leadsG, _=ut.bin_data_nsidc(xv1, yv1, leads.lead_fraction_v1.values, xptsG, yptsG, dx)

xv2, yv2 = mapProj(leads.lons.values, leads.lats.values)
leadsG2, _=ut.bin_data_nsidc(xv2, yv2, leads.lead_fraction_v2.values, xptsG, yptsG, dx)

leadsG[region_mask>10]=np.nan
leadsG2[region_mask>10]=np.nan

# Convert to percent
leadsG=leadsG*100.
leadsG2=leadsG2*100.

gen_netcdf(savePath, xptsG, yptsG, lonG, latG, leadsG,leadsG2, epsg)

data=[leadsG, leadsG2, (leadsG2-leadsG)]

means=['%.02f' %np.nanmean(dataT) for dataT in data]


#fig=plt.figure(figsize=(10, 5))
#ax=plt.gca()
#im1 = ax.imshow(leadsG)

#plt.tight_layout()
#fig.savefig(figPath+'/IS2_gridded_section_lead_stats_'+dateStr+'_'+hem+'_'+release[-3:]+'_001_test.png', dpi=300)

minval=[0, 0, -1]
maxval=[6, 6, 1]
if (hem=='nh'):
    labels=['(a) v1 (ssh)', '(b) v2 (20%)', '(c) difference (b - a)']
else:
    labels=['(d) v1 (ssh)', '(e) v2 (20%)', '(f) difference (e - d)']

cbarlabels=['Lead fraction (%)', 'Lead fraction (%)', 'lead fraction (%)']
cmaps=[plt.cm.YlOrRd,  plt.cm.YlOrRd, plt.cm.RdBu_r]

if (hem=='nh'):
    fig, axs = plt.subplots(1, 3, figsize=(8, 5), subplot_kw=dict(projection=ccrs.NorthPolarStereo(central_longitude=-45)))
    plt.subplots_adjust(bottom=0.01, left=0.02, top=0.95, right=0.98)
else:
   fig, axs = plt.subplots(1, 3, figsize=(8, 3.1), subplot_kw=dict(projection=ccrs.SouthPolarStereo(central_longitude=0)))
   plt.subplots_adjust(bottom=0.02, left=0.02, top=0.95, right=0.98)

for i in range(len(data)):
    ax=axs.flatten()[i]
    plt.sca(ax)
    
    im1 = ax.pcolormesh(lonG , latG, data[i], cmap=cmaps[i], vmin=minval[i], vmax=maxval[i], transform=ccrs.PlateCarree(), edgecolors='None', zorder=1, rasterized=True)

    ax.coastlines(zorder=3)
    ax.gridlines(draw_labels=False,
              linewidth=0.22, color='gray', alpha=0.5, linestyle='--',zorder=3)
    ax.add_feature(cfeature.LAND, facecolor='0.95', zorder=2)

    if (hem=='nh'):
        ax.set_extent([-179.5, 179.5, 50, 90], ccrs.PlateCarree())
        ax.set_xlim([np.percentile(xptsG, 10), np.percentile(xptsG, 90)])
    else:
        ax.set_extent([-179, 179, -90, -50], ccrs.PlateCarree())

    cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.01,shrink=0.9)
    cb=fig.colorbar(im1,cax=cax,extend='both',**kw)
    #cb.set_label(varStr+' ('+units_lab+')',size=8)
    #ax.set_title(varStr+' '+date_string+month_string+extra)
    cb.set_ticks(np.linspace(minval[i], maxval[i], 5) )
    cb.set_label(cbarlabels[i], labelpad=3)
    #im1=hexbin(xpts_sections, ypts_sections, C=data[i], vmin=minval[i], vmax=maxval[i], gridsize=gridbins, cmap=cmaps[i], zorder=2, rasterized=True)
    

    if (hem=='nh'):
        ax.annotate(labels[i], xy=(0.98, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')
        #ax.annotate('mean: '+means[i]+' m', xy=(0.97, 0.87), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')
    else:
        ax.annotate(labels[i], xy=(0.98, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

    ax.annotate('mean: '+means[i]+' %', xy=(0.97, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')


    if (i==0):
        ax.annotate(dateStr, xy=(0.02, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

#plt.tight_layout()
fig.savefig(figPath+'/IS2_gridded_section_lead_stats_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.png', dpi=300)

