""" grid_plot_chords.py
    
    Gridding and plotting the raw chord length estimates
    Initial code written by Alek Petty (02/01/2020)
    
    Input:
        Raw (along-track) chord length estimates

    Output:
        Gridded chord length estimates in a netcdf file, gridded plots

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
from cartopy.util import add_cyclic_point

ut.reset_matplotlib()

def gen_netcdf(savePathT, xpts2d, ypts2d, lons2d, lats2d, chord_lengthv1, chord_lengthv2, epsg):
    f = nc4.Dataset(savePathT+'IS2_gridded_chord_lengths_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.nc','w', format='NETCDF4') #'w' stands for write

    print ('dimensions:', lons2d.shape[0], lons2d.shape[1])
    f.createDimension('x', lons2d.shape[0])
    f.createDimension('y', lons2d.shape[1])

    longitude = f.createVariable('longitude', 'f4', ('x', 'y'))
    latitude = f.createVariable('latitude', 'f4', ('x', 'y'))  

    xgrid = f.createVariable('xgrid', 'f4', ('x', 'y'))
    ygrid = f.createVariable('ygrid', 'f4', ('x', 'y')) 

    # Assume all three main variables
    chord_length_v1 = f.createVariable('chord_length_v1', 'f4', ('x', 'y'))
    chord_length_v2 = f.createVariable('chord_length_v2', 'f4', ('x', 'y'))

    longitude[:] = lons2d #The "[:]" at the end of the variable instance is necessary
    latitude[:] = lats2d

    xgrid[:] = xpts2d #The "[:]" at the end of the variable instance is necessary
    ygrid[:] = ypts2d

    chord_length_v1[:] = np.around(chord_lengthv1, decimals=4)
    chord_length_v2[:] = np.around(chord_lengthv2, decimals=4)
    
    longitude.units = 'degrees East'
    longitude.long_name='longitude'

    latitude.units = 'degrees North'
    latitude.long_name='latitude'

    xgrid.units = 'meters'
    xgrid.long_name='projection grid x values'
    xgrid:description = "center values of projection grid in x direction" 

    ygrid.units = 'meters'
    ygrid.long_name='projection grid y values'
    ygrid:description = "center values of projection grid in y direction" 


    chord_length_v1.description = 'Mean chord length (v1 ssh algorithm) within 25 km x 25 km grid-cells.'
    chord_length_v1.units = 'm'
    chord_length_v1.long_name ='Chord length (v1)'

    chord_length_v2.description = 'Mean chord length (v2 ssh algorithm) within 25 km x 25 km grid-cells.'
    chord_length_v2.units = 'm'
    chord_length_v2.long_name ='Chord length (v2)'

    #Add global attributes
    f.description = "Gridded chord length estimates from ICESat-2. Gridded chord length calculated from individual chord length estimates derived independently from the three strong beams. Two sea surface metrics are used to drive chord length, both variables are included in this dataset (v1 and v2). Data gridded onto the NSIDC Polar Stereographic grid using a simple binning procedure (EPSG:"+epsg+")."
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

if (period==1):
    # Arctic start of winter 2018 rel002/003 crossover
    epsg='3411'
    mapProj = pyproj.Proj("+init=EPSG:"+epsg)

    hem='nh'
    start_date='20181101'
    end_date = '20181101'
    dateStr=start_date+'-'+end_date
    xptsG, yptsG, latG, lonG, proj= ut.create_grid_nsidc_arctic()
    region_mask, xptsR, yptsR=ut.get_region_mask(anc_data_path, mapProj, xypts_return=1)


if (period==2):  
    #mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)
    mplot = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=44., urcrnrlon=105, urcrnrlat=41.37)
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
    mplot = Basemap(epsg=3412,resolution='l', llcrnrlon=225, llcrnrlat=-45, urcrnrlon=42, urcrnrlat=-44)
    epsg='3412'
    mapProj = pyproj.Proj("+init=EPSG:"+epsg)
    hem='sh'
    start_date='20190501'
    end_date = '20190930'
    dateStr=start_date+'-'+end_date
    #gridbins=(300, 300)
    xptsG, yptsG, latG, lonG, proj= ut.create_grid_nsidc_antarctic()
    region_mask, xptsR, yptsR=ut.get_region_mask_s(anc_data_path, mapProj, xypts_return=1)



chords_v1 = xr.open_dataset(savePath+'IS2_raw_chord_lengths_v1_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.nc')
chords_v2 = xr.open_dataset(savePath+'IS2_raw_chord_lengths_v2_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.nc')

chords_v1 = chords_v1.dropna(dim='index')
chords_v2 = chords_v2.dropna(dim='index')

xv1, yv1 = mapProj(chords_v1.lons.values, chords_v1.lats.values)
chordsG, _=ut.bin_data_nsidc(xv1, yv1, chords_v1.lengths.values, xptsG, yptsG, dx)

xv2, yv2 = mapProj(chords_v2.lons.values, chords_v2.lats.values)
chordsG2, _=ut.bin_data_nsidc(xv2, yv2, chords_v2.lengths.values, xptsG, yptsG, dx)

chordsG[region_mask>10]=np.nan
chordsG2[region_mask>10]=np.nan

# Convert to km
chordsG=chordsG*0.001
chordsG2=chordsG2*0.001

gen_netcdf(savePath, xptsG, yptsG, lonG, latG, chordsG,chordsG2, epsg)


data=[chordsG, chordsG2, (chordsG2-chordsG)]

means=['%.02f' %np.nanmean(dataT) for dataT in data]


#fig=plt.figure(figsize=(10, 5))
#ax=plt.gca()
#im1 = ax.imshow(chordsG)
#plt.tight_layout()
#fig.savefig(figPath+'/IS2_gridded_chord_lengths_v1_'+dateStr+'_'+hem+'_'+release[-3:]+'_001_test.png', dpi=300)


# Use the basemap projection/library for plotting as having problems with cartopy
xptsGM , yptsGM = mplot(lonG, latG)

minval=[0, 0, -1]
maxval=[16, 16, 1]
if (hem=='nh'):
    labels=['(a) v1 (ssh)', '(b) v2 (20%)', '(c) difference (b - a)']
else:
    labels=['(d) v1 (ssh)', '(e) v2 (20%)', '(f) difference (e - d)']

cbarlabels=['Chord length (km)', 'Chord length (km)', 'difference (km)']
cmaps=[plt.cm.viridis,  plt.cm.viridis, plt.cm.RdBu_r]

if (hem=='nh'):
    #fig, axs = plt.subplots(1, 3, figsize=(8, 5), subplot_kw=dict(projection=ccrs.NorthPolarStereo(central_longitude=-45)))
    fig, axs = plt.subplots(1, 3, figsize=(8, 5))
    plt.subplots_adjust(bottom=0.01, left=0.02, top=0.95, right=0.98)
else:
   fig, axs = plt.subplots(1, 3, figsize=(8, 3.1), subplot_kw=dict(projection=ccrs.SouthPolarStereo(central_longitude=0)))
   plt.subplots_adjust(bottom=0.02, left=0.02, top=0.95, right=0.98)

for i in range(len(data)):
    ax=axs.flatten()[i]
    plt.sca(ax)
    
    im1 = mplot.pcolormesh(xptsGM , yptsGM, data[i], cmap=cmaps[i], vmin=minval[i], vmax=maxval[i], edgecolors='None', zorder=1, rasterized=True)

    mplot.drawcoastlines(linewidth=0.25, zorder=5)
    mplot.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
    mplot.drawmeridians(np.arange(-180.,180.,30.), latmax=80, linewidth = 0.25, zorder=10)
    mplot.fillcontinents(color='0.95',lake_color='0.9', zorder=3)

    cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.01,shrink=0.9)
    cb=fig.colorbar(im1,cax=cax,extend='both',**kw)
    cb.set_ticks(np.linspace(minval[i], maxval[i], 5) )
    cb.set_label(cbarlabels[i], labelpad=3)

    if (hem=='nh'):
        ax.annotate(labels[i], xy=(0.98, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')
    else:
        ax.annotate(labels[i], xy=(0.98, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

    ax.annotate('mean: '+means[i]+' km', xy=(0.97, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

    if (i==0):
        ax.annotate(dateStr, xy=(0.02, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

fig.savefig(figPath+'/IS2_gridded_chord_lengths_v1_'+dateStr+'_'+hem+'_'+release[-3:]+'_001.png', dpi=300)

