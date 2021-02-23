

### Import necessary Python modules

import cartopy.crs as ccrs
#import cv2
import gdal
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
from pyproj import Proj, transform
import osr
from scipy import ndimage, misc
import sys
from glob import glob
import os
from matplotlib.ticker import MaxNLocator
import numpy.ma as ma
import utils as ut
import argparse

# ## 1. Data ingestion and preparation

################### Input data file names ########################################

relStr='rel003'
data_path = '/cooler/I2-ASAS/'
s2_data_path = '/sea_ice_pso/aapetty/raw_data/'
figPath='../figures/'+relStr+'/profiles/'
savePath='../data/'+relStr+'/profiles/'
bufferRatio=60 # Bigger number smaler width (try 100 for zoom in)

# Function to parse command line arguments
def parse():
    """ Parse command line input arguments. """
    parser = argparse.ArgumentParser(description='S2 plot arguments')
    parser.add_argument('-e', action='store', default='', dest='example',
                        help='Arctic or Antarctic',
                        type=str)
    parser.add_argument('-b', action='store', default='', dest='beam_str',
                        help='beam str (e.g. gt1r)',
                        type=str)
    inps = parser.parse_args()
    return inps

# Read command line input parameters
inps = parse()
# Check for missing command line input parameters
if inps.example is None:
    print('Missing hemisphere. Use option -e')
    sys.exit(1)
if inps.beam_str is None:
    print('Missing beam. Use option -b')
    sys.exit(1)

example = inps.example
beam_str = inps.beam_str

print(example, beam_str)

if example=='Arctic':
    
    beam_name=beam_str
    date='20190526_NH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/RT_T14XMQ_20190525T230121_B04.tif'
    # ATLAS DATA
    # -1 indicates we're using the latest version of the data (i.e. the higher number)
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/data/ATL07-01_20190526*_08820301_'+relStr[-3:]+'_'+version+'.h5')[-1]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-01/data/ATL10-01_20190526*_08820301_'+relStr[-3:]+'_'+version+'.h5')[-1]
    
    # 1st 50 km (good)
    Slim = 80.15
    Nlim = 80.57

    # 2nd 50 km (2 leads)
    #Slim = 80.57
    #Nlim = 81
    
    # paper draft
    #Slim = 80.73
    #Nlim = 80.93

    # focussing on the big lead in gt21 (need to also manually change the size of the top line when doing so!)
    #Slim = 80.6
    #Nlim = 80.66

else:
    print('Antarctic example')
    beam_name=beam_str
    date='20190317_SH'
    version='*'
    # Sentinel-2 surface reflectance single-band geoTIFF file
    imageFile = s2_data_path+'/RT_T20CPE_20190317T122209_B04.tif'
    # ATLAS DATA
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-02/data/ATL07-02_20190317*_12070201_'+relStr[-3:]+'_'+version+'.h5')[0]
    ATL10_dataFile = glob(data_path+relStr+'/ATL10-02/data/ATL10-02_20190317*12070201_'+relStr[-3:]+'_'+version+'.h5')[0]
    
    # paper draft (first 50 km roughly)
    #Slim = -73
    #Nlim = -72.6

    # include a little bit more cloud
    Slim = -73
    Nlim = -72.57

    #Slim = -73
    #Nlim = -72.5
    
try:
    geoTIFF = gdal.Open(imageFile)
except IOError:
    sys.exit('geoTIFF file is not a valid file')

### Check if ICESat-2 ATLAS file is valid
#print(IS2_dataFile)
try:
    ATL07file = h5py.File(ATL07_dataFile, 'r')
except IOError:
    sys.exit('not a valid file')


### Check if ICESat-2 ATLAS file is valid
#print(IS2_dataFile)
try:
    ATL10file = h5py.File(ATL10_dataFile, 'r')
except IOError:
    sys.exit('not a valid file')


### Read geoTIFF file and print information

# Image size
geoTIFF_size_x = geoTIFF.RasterXSize # Size of x
geoTIFF_size_y = geoTIFF.RasterYSize # Size of y

# Geographic information (geocoded extent, projection, etc.)
geoTIFF_geo_trans = geoTIFF.GetGeoTransform()   # UL corner x-coordinate, W-E pixel size, rotation (0 if N up),
                                                # UL corner y-coordinate, rotation (0 if N up), N-S pixel size
    
geoTIFF_proj = geoTIFF.GetProjection() # Geographic projection information
geoTIFF_inproj = osr.SpatialReference()
geoTIFF_inproj.ImportFromWkt(geoTIFF_proj)

# Convert WKT projection information to a CartoPy projection
geoTIFF_projcp = geoTIFF_inproj.GetAuthorityCode('PROJCS')
#geoTIFF_projection = ccrs.epsg(geoTIFF_projcp)

# Extract edge coordinates using half resolution cell to get coordinates of pixel center
geoTIFF_x_min = geoTIFF_geo_trans[0] + geoTIFF_geo_trans[1]*0.5
geoTIFF_x_max = geoTIFF_geo_trans[0] + geoTIFF_size_x*geoTIFF_geo_trans[1] - geoTIFF_geo_trans[1]*0.5
geoTIFF_y_min = geoTIFF_geo_trans[3] + geoTIFF_size_y*geoTIFF_geo_trans[5] + geoTIFF_geo_trans[5]*0.5
geoTIFF_y_max = geoTIFF_geo_trans[3] - geoTIFF_geo_trans[5]*0.5

geoTIFF_extent = [geoTIFF_x_min, geoTIFF_x_max, geoTIFF_y_min, geoTIFF_y_max]

# Convert geoTIFF image to data array
surfReflectance = geoTIFF.ReadAsArray()

# Close dataset
geoTIFF = None

# Print basic information
print('Image size: ', geoTIFF_size_x, geoTIFF_size_x)
print('Image projection EPGS code: ', geoTIFF_projcp)
print('Xmin: ', geoTIFF_extent[0])
print('Xmax: ', geoTIFF_extent[1])
print('Ymin: ', geoTIFF_extent[2])
print('Ymax: ', geoTIFF_extent[3])


sc_orient = ATL07file['orbit_info/sc_orient'][0]
if (sc_orient==0):
    #backward
    if (beam_name[-1]=='l'):
        beamStrength='strong'
    else:
        beamStrength='weak'
elif (sc_orient==1):
    #forward
    if (beam_name[-1]=='l'):
        beamStrength='weak'
    else:
        beamStrength='strong'
else:
    print('no clear spaceraft orientation')


## Read ATL07 parameters of interest
sc_orient = ATL07file['orbit_info/sc_orient'][:]
seg_delta_time = ATL07file[beam_name + '/sea_ice_segments/delta_time'][:]
seg_id = ATL07file[beam_name + '/sea_ice_segments/height_segment_id'][:] 
seg_lat = ATL07file[beam_name + '/sea_ice_segments/latitude'][:]
seg_lon = ATL07file[beam_name + '/sea_ice_segments/longitude'][:]
seg_track_dist = ATL07file[beam_name + '/sea_ice_segments/seg_dist_x'][:]
seg_conf = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_confidence'][:]
seg_height = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_height'][:]
seg_length = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_length_seg'][:]
seg_quality_binary = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_quality'][:]
seg_quality_flag = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_fit_quality_flag'][:]
seg_rms = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_rms'][:]
seg_ssh_flag = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_ssh_flag'][:]
seg_surf_err = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_surface_error_est'][:]
seg_type = ATL07file[beam_name + '/sea_ice_segments/heights/height_segment_type'][:]
seg_gaussian = ATL07file[beam_name +'/sea_ice_segments/heights/height_segment_w_gaussian'][:]
seg_backgr_rate = ATL07file[beam_name + '/sea_ice_segments/stats/backgr_r_200'][:]
seg_exmax_u1 = ATL07file[beam_name + '/sea_ice_segments/stats/exmax_mean_1'][:]
seg_exmax_u2 = ATL07file[beam_name + '/sea_ice_segments/stats/exmax_mean_2'][:]
seg_exmax_mix = ATL07file[beam_name + '/sea_ice_segments/stats/exmax_mix'][:]
seg_exmax_s1 = ATL07file[beam_name + '/sea_ice_segments/stats/exmax_stdev_1'][:]
seg_exmax_s2 = ATL07file[beam_name + '/sea_ice_segments/stats/exmax_stdev_2'][:]
seg_photon_rate = ATL07file[beam_name + '/sea_ice_segments/stats/photon_rate'][:]
seg_h_coarse_mn = ATL07file[beam_name + '/sea_ice_segments/stats/height_coarse_mn'][:]
seg_h_coarse_std = ATL07file[beam_name + '/sea_ice_segments/stats/height_coarse_stdev'][:]
seg_h_mean = ATL07file[beam_name + '/sea_ice_segments/stats/hist_mean_h'][:]
seg_h_width = ATL07file[beam_name + '/sea_ice_segments/stats/hist_w'][:]
seg_solar_elev = ATL07file[beam_name + '/sea_ice_segments/geolocation/solar_elevation'][:]

# #Add this value to delta time parameters to compute full gps_seconds
atlas_epoch=ATL07file['/ancillary_data/atlas_sdp_gps_epoch'][0] 

# subtract leap seconds offset (will change to 38 soon, watch out)
gps_seconds = atlas_epoch + seg_delta_time - 37

# Close HDF5 file
ATL07file = None
# Create Pandas data frame with all variables
ATL07dF = pd.DataFrame({'DeltaTime':seg_delta_time,
                      'SegmentID':seg_id,
                      'lats':seg_lat,
                      'lons':seg_lon,
                      'along_dist':seg_track_dist,
                      'HConfid':seg_conf,
                      'height':seg_height,
                      'seg_length':seg_length,
                      'QualityBinary':seg_quality_binary,
                      'QualityFlag':seg_quality_flag,
                      'RMS':seg_rms,
                      'ssh_flag':seg_ssh_flag,
                      'seg_surf_err':seg_surf_err,
                      'seg_type':seg_type,
                      'WGaussian':seg_gaussian,
                      'BackgrRate':seg_backgr_rate,
                      'ExMean1':seg_exmax_u1,
                      'ExMean2':seg_exmax_u2,
                      'ExMix':seg_exmax_mix,
                      'ExStd1':seg_exmax_s1,
                      'ExStd2':seg_exmax_s2,
                      'PhRate':seg_photon_rate,
                      'HCoarseMn':seg_h_coarse_mn,
                      'HCoarseStd':seg_h_coarse_std,
                      'SolarElev':seg_solar_elev,
                      'seg_h_mean':seg_h_mean,
                      'seg_h_width':seg_h_width,
                      'gpsseconds':gps_seconds
                     })

    
# Transform WGS84 Latitude and Longitude to UTM X and Y
ATL_inProj = Proj(init='epsg:4326') # ATLAS data projection
out_proj = 'epsg:' + str(geoTIFF_projcp) # String with geoTIFF projection epsg code
ATL_outProj = Proj(init=out_proj) # GeotTIFF projection

seg_x_utm, seg_y_utm = transform(ATL_inProj, ATL_outProj, seg_lon, seg_lat) # Transform coordinates

# Add new coorinates to dataframe
ATL07dF['UTM_X'] = seg_x_utm
ATL07dF['UTM_Y'] = seg_y_utm

# Replace no-data value 3.4028235e+38 with Numpy NaN in data frame
ATL07dF = ATL07dF.replace(np.max(ATL07dF.Height), np.nan)

# Crop ICESat-2 data frame to imagery extent
ATL07dF_crop = ATL07dF[(ATL07dF['UTM_X'] > geoTIFF_extent[0]) & (ATL07dF['UTM_X'] < geoTIFF_extent[1]) & 
                   (ATL07dF['UTM_Y'] > geoTIFF_extent[2]) & (ATL07dF['UTM_Y'] < geoTIFF_extent[3])]

# Reset data frame index value
ATL07dF_crop = ATL07dF_crop.reset_index(drop=True)

# Print dataframe size 
print('Full dataset n. of columns: ', ATL07dF.shape[1])
print('Full dataset n. of segments: ', ATL07dF.shape[0])
print('Cropped dataset n. of columns: ', ATL07dF_crop.shape[1])
print('Cropped dataset n. of segments: ', ATL07dF_crop.shape[0])


## Read ATL10 parameters of interest
freeboard=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['beam_fb_height'][:]
ssh_flag=ATL10file[beam_name]['freeboard_beam_segment']['height_segments']['height_segment_ssh_flag'][:]
seg_type_flag=ATL10file[beam_name]['freeboard_beam_segment']['height_segments']['height_segment_type'][:]
seg_length=ATL10file[beam_name]['freeboard_beam_segment']['height_segments']['height_segment_length_seg'][:]
freeboard_confidence=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['beam_fb_confidence'][:]
#seg_height=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['height_segment_height'][:]
height_segment_id=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['height_segment_id'][:]
lons=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['lons'][:]
lats=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['lats'][:]
deltaTime=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['delta_time'][:]-ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['delta_time'][0]
#seg_length = ATL10file[beam_name + '/sea_ice_segments/heights/height_segment_length_seg'][:]
#seg_quality = ATL10file[beam_name + '/sea_ice_segments/heights/height_segment_quality'][:]
#seg_rms = ATL10file[beam_name + '/sea_ice_segments/heights/height_segment_rms'][:]

segDist=ATL10file[beam_name]['freeboard_beam_segment']['beam_freeboard']['seg_dist_x'][:] #-f1[beam]['freeboard_beam_segment']['beam_freeboard']['seg_dist_x'][0]


lead_loc=ATL10file[beam_name]['leads']['ssh_ndx'][:]
    
#lead_id = height_segment_id[lead_loc]
lead_flag = np.zeros((np.size(height_segment_id)))
# ADD SIZE
lead_flag[lead_loc]=1


# Create Pandas data frame with all variables
ATL10dF = pd.DataFrame({'freeboard':freeboard, 'lons':lons, 'lats':lats, 'delta_time':deltaTime,
                      'along_dist':segDist, 'height_segment_id': height_segment_id, 'lead_flag':lead_flag, 
                        'ssh_flag':ssh_flag, 'seg_type_flag':seg_type_flag, 'seg_length': seg_length})


    
seg_x_utm, seg_y_utm = transform(ATL_inProj, ATL_outProj, lons, lats) # Transform coordinates

# Add new coorinates to dataframe
ATL10dF['UTM_X'] = seg_x_utm
ATL10dF['UTM_Y'] = seg_y_utm

# Replace no-data value 3.4028235e+38 with Numpy NaN in data frame
#ATL10dF = ATL10dF.replace(np.max(ATL10dF.freeboard), np.nan)
# dF = dF[(dF['freeboard']<maxFreeboard)]
# Crop ICESat-2 data frame to imagery extent
ATL10dF_crop = ATL10dF[(ATL10dF['UTM_X'] > geoTIFF_extent[0]) & (ATL10dF['UTM_X'] < geoTIFF_extent[1]) & 
                   (ATL10dF['UTM_Y'] > geoTIFF_extent[2]) & (ATL10dF['UTM_Y'] < geoTIFF_extent[3])]

# Reset data frame index value
ATL10dF_crop = ATL10dF_crop.reset_index(drop=True)

# Print dataframe size 
print('Full dataset n. of columns: ', ATL10dF.shape[1])
print('Full dataset n. of segments: ', ATL10dF.shape[0])
print('Cropped dataset n. of columns: ', ATL10dF_crop.shape[1])
print('Cropped dataset n. of segments: ', ATL10dF_crop.shape[0])


### Extract Sentinel-2 surface reflectance values at ICESat-2 segment locations

# Function to extract imagery pixel values at locations of ICESat-2 data
def get_value_at_point(array_from, pos):
    """ Extract raster value at given position from coordinates """
    samp_x = int((pos[0] - geoTIFF_geo_trans[0]) / geoTIFF_geo_trans[1])
    samp_y = int((pos[1] - geoTIFF_geo_trans[3]) / geoTIFF_geo_trans[5])

    return array_from[samp_y-1, samp_x-1]

# Preallocate empty array
dF_length = ATL07dF_crop.shape[0]
band_value = np.empty([dF_length, 1])

# Extract imagery pixel values at locations of ICESat-2 data
for i in range(0, dF_length):
    band_value[i] = get_value_at_point(surfReflectance, (ATL07dF_crop['UTM_X'][i], ATL07dF_crop['UTM_Y'][i]))

# Add new column to data frame
ATL07dF_crop['BandValue'] = band_value

### Plot Sentinel-2 image with ICESat-2 beam overlaid

# Select ICESat-2 parameter to plot
# 'DeltaTime', 'SegmentID', 'lats', 'lons', 'along_dist', HConfid', 'height', 'seg_length', 
# 'Quality', 'RMS', 'ssh_flag'', 'seg_surf_err', 'seg_type', 'WGaussian', 'BackgrRate', 'ExMean1', 'ExMean2', 'ExMix',
# 'ExStd1', 'ExStd2', 'PhRate', 'HCoarseMn', 'HCoarseStd', 'SolarElev','UTM_X', 'UTM_Y'
IS2_param = 'ssh_flag'

# Make plot
#get_ipython().run_line_magic('matplotlib', 'inline')

#subplot_kw = dict(projection=geoTIFF_projection)
#fig, ax = plt.subplots(figsize=(10, 10))

# Plot image as grayscale
#im = ax.imshow(surfReflectance, extent=geoTIFF_extent, origin='upper', cmap='gray')

# Overlay ICESat-2 data
#plt.scatter(ATL07dF_crop['UTM_X'], ATL07dF_crop['UTM_Y'], s=1, c=ATL07dF_crop[IS2_param], cmap='viridis')
#plt.colorbar()


# ## 2. Area of Interest selection

### Subset datasets to limits
sub_dF07 = ATL07dF_crop[(ATL07dF_crop['lats'] < Nlim) & (ATL07dF_crop['lats'] > Slim) ]

# Reset data frame index value
sub_dF07 = sub_dF07.reset_index(drop=True)

# Print dataframe size 
print('AOI n. of columns: ', sub_dF07.shape[1], " * added surface refelctance value.")
print('AOI n. of points: ', sub_dF07.shape[0])


### Subset datasets to limits
sub_dF10 = ATL10dF_crop[(ATL10dF_crop['lats'] < Nlim) & (ATL10dF_crop['lats'] > Slim) ]

# Reset data frame index value
sub_dF10 = sub_dF10.reset_index(drop=True)

# Print dataframe size 
print('AOI n. of columns: ', sub_dF10.shape[1], " * added surface refelctance value.")
print('AOI n. of points: ', sub_dF10.shape[0])

print('start lat/lon of subset', sub_dF07['lats'].iloc[0], sub_dF07['lons'].iloc[0])
print('end lat/lon of subset', sub_dF07['lats'].iloc[-1], sub_dF07['lons'].iloc[-1])

print('start lat/lon of ATL10 subset', sub_dF10['lats'].iloc[0], sub_dF10['lons'].iloc[0])
print('end lat/lon of ATL10 subset', sub_dF10['lats'].iloc[-1], sub_dF10['lons'].iloc[-1])

print('start lat/lon of profile across entire image', ATL07dF_crop['lats'].iloc[0], ATL07dF_crop['lons'].iloc[0])
print('end lat/lon of profile across entire image', ATL07dF_crop['lats'].iloc[-1], ATL07dF_crop['lons'].iloc[-1])


print('length of subset', str(int(0.001*(sub_dF07['along_dist'].iloc[-1] - sub_dF07['along_dist'].iloc[0])))+' km')
print('length of ATL10 subset', str(int(0.001*(sub_dF10['along_dist'].iloc[-1] - sub_dF10['along_dist'].iloc[0])))+' km')
print('length of profile across entire image', str(int(0.001*(ATL07dF_crop['along_dist'].iloc[-1] - ATL07dF_crop['along_dist'].iloc[0])))+' km')


### Determine ICESat-2 ground track inclination in Sentinel-2 image

# Trigonometry calculations
a = sub_dF07.UTM_Y.iloc[0] - sub_dF07.UTM_Y.iloc[-1]
b = sub_dF07.UTM_X.iloc[0] - sub_dF07.UTM_X.iloc[-1]
c = np.sqrt((a**2) + (b**2))

cosA = (b**2 + c**2 - a**2) / (2*b*c)
inclin = np.arccos(cosA)
inclin_deg = np.degrees(inclin)

# Adjust size of buffer zone from ICESat-2 ground track depending on length of subset
buffer = c/bufferRatio

deltaX = buffer * np.sin(inclin)
deltaY = buffer * np.cos(inclin)

# Compute coordinates of corners of AOI

# Descending satellite pass
if sub_dF07.UTM_Y.iloc[0] > sub_dF07.UTM_Y.iloc[-1]:
    rev = True # set flag to determine how to apply rotation
    
    UL_x = sub_dF07.UTM_X.iloc[0] - deltaX
    UL_y = sub_dF07.UTM_Y.iloc[0] + deltaY

    UR_x = sub_dF07.UTM_X.iloc[0] + deltaX
    UR_y = sub_dF07.UTM_Y.iloc[0] - deltaY

    LL_x = sub_dF07.UTM_X.iloc[-1] - deltaX
    LL_y = sub_dF07.UTM_Y.iloc[-1] + deltaY

    LR_x = sub_dF07.UTM_X.iloc[-1] + deltaX
    LR_y = sub_dF07.UTM_Y.iloc[-1] - deltaY

# Ascending satellite pass    
if sub_dF07.UTM_Y.iloc[0] < sub_dF07.UTM_Y.iloc[-1]:
    rev = False # set flag to determine how to apply rotation
    
    UL_x = sub_dF07.UTM_X.iloc[0] - deltaX
    UL_y = sub_dF07.UTM_Y.iloc[0] - deltaY

    UR_x = sub_dF07.UTM_X.iloc[0] + deltaX
    UR_y = sub_dF07.UTM_Y.iloc[0] + deltaY

    LL_x = sub_dF07.UTM_X.iloc[-1] - deltaX
    LL_y = sub_dF07.UTM_Y.iloc[-1] - deltaY

    LR_x = sub_dF07.UTM_X.iloc[-1] + deltaX
    LR_y = sub_dF07.UTM_Y.iloc[-1] + deltaY

# Print AOI information
print('ICESat-2 inclination angle in imagery coordinate system: ', inclin_deg)
print('AOI coordinates:')
print('UL: ', UL_x, UL_y)
print('UR: ', UR_x, UR_y)
print('LR: ', LR_x, LR_y)
print('LL: ', LL_x, LL_y)



### Make new plot of Sentinel-2 image with ICESat-2 beam overlaid (AOI only)

# Make plot
#get_ipython().run_line_magic('matplotlib', 'inline')

#subplot_kw = dict(projection=geoTIFF_projection)
fig, ax = plt.subplots(figsize=(10, 10))

# Plot image as grayscale
im = ax.imshow(surfReflectance, extent=geoTIFF_extent, origin='upper', cmap='gray')

# Overlay ICESat-2 data
plt.scatter(sub_dF07['UTM_X'], sub_dF07['UTM_Y'], s=1, c=sub_dF07[IS2_param], cmap='viridis')
plt.plot([UL_x, LL_x, LR_x, UR_x, UL_x], [UL_y, LL_y, LR_y, UR_y, UL_y], c='r')
plt.colorbar()
plt.savefig(figpath+'overlay'+date+beam_name+relStr+str(Slim)+str(Nlim)+'.png', dpi=300)
### Create shapefile with extent of AOI to be used for image cropping

# Define function
def createShapefile(filename, LL_x, LL_y, UL_x, UL_y, UR_x, UR_y, LR_x, LR_y):
    import ogr

    driver = ogr.GetDriverByName('ESRI Shapefile')

    datasource = driver.CreateDataSource(filename)
    layer = datasource.CreateLayer('layerName',geom_type=ogr.wkbPolygon)

    outline = ogr.Geometry(type=ogr.wkbLinearRing)
    outline.AddPoint(LL_x, LL_y)
    outline.AddPoint(UL_x, UL_y)
    outline.AddPoint(UR_x, UR_y)
    outline.AddPoint(LR_x, LR_y)
    outline.AddPoint(LL_x, LL_y)
    polygon = ogr.Geometry(type=ogr.wkbPolygon)
    polygon.AddGeometry(outline)

    #create feature object with polygon geometry type from layer object:
    feature = ogr.Feature( layer.GetLayerDefn() )
    feature.SetGeometry(polygon)
    layer.CreateFeature(feature)
   
    return

# Call function to create AOI polygon shapefile
createShapefile(figpath+'AOIPolygon.shp', LL_x, LL_y, UL_x, UL_y, UR_x, UR_y, LR_x, LR_y)


# ### This part runs GDAL utilities from the command line. Must be adapted to local environment

### Run gdalwarp to crop Sentinel-2 image to AOI

# Provide full path to gdalwarp on local system
#path_to_gdalwarp = '/Users/aapetty/opt/anaconda3/envs/py37/bin/gdalwarp'
#path_to_gdalwarp = '/att/opt/other/centos/anaconda3/envs/earthml/bin/gdalwarp'
path_to_gdalwarp='/home/aapetty/.conda/envs/py37/bin/gdalwarp'
# Remove pre-existing files if needed
#get_ipython().system(" rm $figpath'AOI.'*")

shapefile=figpath+'AOIPolygon.shp'
outpath=figpath+'AOI.tif'
# Run gdalwarp
#os.system(' $path_to_gdalwarp -cutline $shapefile -crop_to_cutline -dstalpha $imageFile $outpath')

# trying the python way

ds = gdal.Warp(figpath+'AOI.tif',
               imageFile,
               format='GTiff',
               cutlineDSName= figpath+'AOIPolygon.shp',
               cropToCutline = True,
               dstAlpha = True)
ds = None
#gdal.Warp()


### Convert GEOTiff AOI file to JPG for subsequent use

# Options to convert GEOTiff image to JPEG
options_list = ['-ot Byte','-of JPEG','-b 1','-scale'] 
options_string = " ".join(options_list)

# Use gdal_translate in Python
gdal.Translate(figpath+'AOI.jpg',
               figpath+'AOI.tif',
               options=options_string)

# Read new JPEG image
#image = cv2.imread('AOI.tif')
import matplotlib.pyplot as plt
image = plt.imread(figpath+'AOI.jpg')


# Rotate image by inclination angle according to satellite pass
# Descending
if rev is True:
    rot_im = ndimage.rotate(image, 180-inclin_deg, reshape=True, order=1)
# Ascending
if rev is False:
    rot_im = ndimage.rotate(image, inclin_deg, reshape=True, order=1)

### Trim edges to remove image padding
padNum=12
crop = np.delete(rot_im,np.where(~rot_im.any(axis=0))[0], axis=1)
crop = np.delete(crop,np.where(~crop.any(axis=1))[0], axis=0)
crop = np.delete(crop, range(padNum), 0)
crop = np.delete(crop, range(crop.shape[0]-padNum,crop.shape[0]), 0)
crop = np.delete(crop, range(padNum), 1)
crop = np.delete(crop, range(crop.shape[1]-padNum,crop.shape[1]), 1)

#if you want to just load this and use
crop.dump(figpath+'cropped_image'+date+beam_name+relStr+str(Slim)+str(Nlim))
np.savetxt(figpath+'cropped_image'+date+beam_name+relStr+str(Slim)+str(Nlim)+'.txt', crop)
print(crop.shape)
#import numpy.ma as ma
#crop=ma.masked_where(crop<1, crop)


# ## 3. Generate final plot

### Prepare data for plotting

# Filp image if necessary
if rev is False:
    crop = np.flipud(np.fliplr(crop))

# Convert distances from meters to kilometers with respect to first segment in AOI
sub_dF07['along_dist'] = sub_dF07['along_dist'] - sub_dF07['along_dist'].iloc[0]
sub_dF07['along_dist'] = sub_dF07['along_dist'].multiply(0.001)

# Convert distances from meters to kilometers with respect to first segment in AOI
sub_dF10['along_dist'] = sub_dF10['along_dist'] - sub_dF10['along_dist'].iloc[0]
sub_dF10['along_dist'] = sub_dF10['along_dist'].multiply(0.001)

# Convert background rate to MHz
sub_dF07['BackgrRate'] = sub_dF07['BackgrRate'].multiply(0.000001)



atl07_segs=np.size(sub_dF07['ssh_flag'].values)
atl07_segs


ssh07_segs=np.size(np.where(sub_dF07['ssh_flag'] == 1))
ssh07_segs

specular_segs=np.size(np.where((sub_dF07['seg_type'] > 1.5)&(sub_dF07['seg_type'] < 5.5)))
specular_segs

darklead_segs=np.size(np.where(sub_dF07['seg_type'] > 5.5))
darklead_segs

cloud_segs=np.size(np.where(sub_dF07['seg_type'] < 0.5))
cloud_segs

atl10_segs=np.size(sub_dF10['ssh_flag'].values)
atl10_segs

ssh10_segs=np.size(np.where(sub_dF10['ssh_flag'] > 1.5))
ssh10_segs

mean_fb=np.round(np.mean(sub_dF10['freeboard']), decimals=3)
mean_fb



mafreeboards=ma.masked_where(sub_dF10['freeboard']<0.06, sub_dF10['freeboard'])


# Find the segments that are in ATL07 but not ATL10 
atlmask = np.isin(sub_dF07['SegmentID'].values, sub_dF10['height_segment_id'].values,  invert=True)
atlmaskint=atlmask.astype(int)
np.size(atlmask), np.size(np.nonzero(atlmaskint)[0])

#ssh ids in ATL10
ssh_ids=sub_dF07['SegmentID'].values[np.where(sub_dF07['ssh_flag'].values>0.5)]
ssh_inatl10 = np.isin(ssh_ids, sub_dF10['height_segment_id'].values)

#ssh ids
radio_ids=sub_dF07['SegmentID'].values[np.where(sub_dF07['seg_type'].values<1.5)]
radio_inatl10 = np.isin(radio_ids, sub_dF10['height_segment_id'].values,  invert=True)
radio_inatl10

sub_dF10_adj = sub_dF10[(sub_dF10['ssh_flag'] < 0.5) & (sub_dF10['seg_type_flag']==1 ) & (sub_dF10['freeboard']>0.04)]
#
mean_fb_adj=np.round(np.mean(sub_dF10_adj['freeboard']), decimals=3)


# In[47]:

shadingSize=30
cmap_ssh = colors.ListedColormap(['red', 'yellow'])
# Generate plot with white background
#%matplotlib qt
#%matplotlib inline
fig= plt.figure(figsize=(10,8))
#fig.patch.set_facecolor('xkcd:white')

#x0 = (sub_dF07['along_dist'] - sub_dF07['seg_length']/2000).values
#x1 = (sub_dF07['along_dist'] + sub_dF07['seg_length']/2000).values
#ssh = sub_dF07['ssh_flag''].values

# ATL07 heights ssh label
ax1 = plt.subplot2grid((5, 1), (0, 0), rowspan=1)


for index, row in sub_dF07.iterrows():
    #x0 = row['along_dist'] - row['seg_length']/2000
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if row['ssh_flag'] > 0.1:
        plt.axvspan(x0,x1, facecolor='k', alpha=0.2)
        
plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=sub_dF07['height'], cmap='hot', s=2, zorder=2)
plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
ax1.xaxis.set_ticklabels([])
plt.ylim(-0.5, 3)
ax1.annotate('Grey shading is ssh_flag',xy=(0.02, 0.87), xycoords='axes fraction')
plt.ylabel('ATL07 Height')
# ATL07 heights (in atl10)
ax2 = plt.subplot2grid((5, 1), (1, 0), rowspan=1)
#for x in range(np.size(x0)):
    #print(x)
#for m in range(np.size(atlmask)):
#    if atlmask[m]==False:
#        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if row['ssh_flag'] > 0.1:
        plt.axvspan(x0,x1, facecolor='k', alpha=0.2)
        
plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=atlmaskint, vmin=0, vmax=1, cmap='viridis_r', s=2, zorder=2)
plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
plt.ylim(-0.5, 3)
ax2.xaxis.set_ticklabels([])
plt.ylabel('ATL07 Height')
ax2.annotate('Purple markers not in ATL10',xy=(0.02, 0.87), xycoords='axes fraction')

ax3 = plt.subplot2grid((5, 1), (2, 0), rowspan=1)
#for x in range(np.size(x0)):
    #print(x)
#for m in range(np.size(atlmask)):
#    if atlmask[m]==False:
#        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if row['ssh_flag'] > 0.1:
        plt.axvspan(x0,x1, facecolor='k', alpha=0.2)


cmap = plt.cm.get_cmap('gnuplot', 7)    # 11 discrete colors

imq=plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=sub_dF07['QualityFlag'], cmap=cmap, vmin=-1.5, vmax=5.5, s=2, zorder=2)
plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
plt.ylim(-0.5, 3)
ax3.xaxis.set_ticklabels([])
plt.ylabel('ATL07 Height')
percent=np.round(np.size(np.where(sub_dF07['QualityFlag']>4))/np.size(sub_dF07['QualityFlag'])*100, decimals=2)

ax3.annotate('Quality flag ('+str(percent)+'% > 4)',xy=(0.02, 0.87), xycoords='axes fraction')

cax = fig.add_axes([0.925, 0.41, 0.022, 0.2])
cbar = plt.colorbar(imq,cax=cax, orientation='vertical', extend='both', use_gridspec=True)
#cbar.set_label('Quality flag', labelpad=28, rotation=0)
xticks1 = np.linspace(-1, 5, 7)
cbar.set_ticks(xticks1)
cbar.set_ticklabels(['inv', '0', '1', '2', '3', '4', '5'])
#plt.clim(-1.5,5.5)

ax4 = plt.subplot2grid((5, 1), (3, 0), rowspan=1)
#mask = np.isin(sub_DF['height_segment_id'].values, sub_DF10['height_segment_id'].values,  invert=True)
#for m in range(np.size(atlmask)):
#    if atlmask[m]==False:
#        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)
# Add lead_flag=1 shading
for index, row in sub_dF10.iterrows():
    # CHANGE TO LEAD LENGTH
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if row['ssh_flag'] > 0.1:
        plt.axvspan(x0,x1, facecolor='k', alpha=0.2)

ax4.annotate('Grey shading is leads',xy=(0.02, 0.87), xycoords='axes fraction')
plt.scatter(sub_dF10['along_dist'], sub_dF10['freeboard'], c=sub_dF10['freeboard'], cmap='hot', s=2, zorder=2)
#ax3.xaxis.set_ticklabels([])
plt.ylabel('Freeboard')
plt.xlabel('Along track distance (km)')
plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
plt.ylim(-0.5, 3)
ax4.annotate('Mean = '+str(mean_fb)+' m',xy=(0.8, 0.84), xycoords='axes fraction')

ax5 = plt.subplot2grid((5, 1), (4, 0), rowspan=1)
#mask = np.isin(sub_DF['height_segment_id'].values, sub_DF10['height_segment_id'].values,  invert=True)
#for m in range(np.size(atlmask)):
#    if atlmask[m]==False:
#        plt.axvspan(x0[m],x1[m], facecolor='k', alpha=0.1)
# Add lead_flag=1 shading
for index, row in sub_dF10.iterrows():
    # CHANGE TO LEAD LENGTH
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if row['ssh_flag'] > 0.1:
        plt.axvspan(x0,x1, facecolor='k', alpha=0.2)

ax5.annotate('Adjusted freeboard',xy=(0.02, 0.87), xycoords='axes fraction')
plt.scatter(sub_dF10_adj['along_dist'], sub_dF10_adj['freeboard'], c=sub_dF10_adj['freeboard'], cmap='hot', s=2, zorder=2)
#ax3.xaxis.set_ticklabels([])
plt.ylabel('Freeboard')
plt.xlabel('Along track distance (km)')
plt.xlim(np.amin(sub_dF10['along_dist']), np.amax(sub_dF10['along_dist']))
plt.ylim(-0.5, 3)
ax5.annotate('Mean = '+str(mean_fb_adj)+' m',xy=(0.8, 0.84), xycoords='axes fraction')


plt.subplots_adjust(left=0.06, right=0.92, top=0.96, bottom=0.08, wspace=0, hspace=0)
#plt.tight_layout()
### Save plot to file
plt.savefig(figpath+'FINAL_plot'+date+beam_name+relStr+str(Slim)+str(Nlim)+'diagnostic.png', dpi=300)


percent=np.round(np.size(np.where(sub_dF07['QualityFlag']==5))/np.size(sub_dF07['QualityFlag'])*100, decimals=2)
percent


#get_ipython().run_line_magic('matplotlib', 'inline')
#fig= plt.figure(figsize=(13,5))
#f, (ax1, ax2) = plt.subplots(2, 1, sharey=True)
f, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (11, 4), sharex=False, sharey=True)
plt.sca(ax1)
fig.patch.set_facecolor('xkcd:white')
plt.scatter(sub_dF07['along_dist'].values[atlmask], sub_dF07['height'].values[atlmask], color='r', s=2, alpha=0.3)
#plt.scatter(sub_DF07['along_dist'].values[~atlmask], sub_DF07['elev'].values[~atlmask], color='b', s=2, alpha=0.3)
plt.annotate('Segments not in ATL10', xy=(0.02, 0.92), xycoords='axes fraction')
plt.xlabel('Along track distance (km)')
plt.ylabel('Heights (m)')
plt.sca(ax2)
plt.boxplot([sub_dF07['height'].values[~np.isnan(sub_dF07['height'].values)], sub_dF07['height'].values[~atlmask][~np.isnan(sub_dF07['height'].values[~atlmask])], sub_dF07['height'].values[atlmask][~np.isnan(sub_dF07['height'].values[atlmask])]])
ax2.set_xticklabels(['ATL07 Heights', 'In ATL10', 'Not in ATL10'])
plt.savefig(figpath+'FINAL_plot'+date+beam_name+relStr+str(Slim)+str(Nlim)+'diagnostic_boxplot.png', dpi=300)


# Filter data
#sub_dF07= sub_dF07.dropna() 
#sub_dF07 = sub_dF07[sub_dF07['seg_length']<200]
#sub_dF07 = sub_dF07[sub_dF07['seg_length']>2]
#sub_dF07 = sub_dF07[sub_dF07['height']<1e6]
#sub_dF07 = sub_dF07[sub_dF07['seg_h_mean']<1e6]


pddata_v1, pddata_v2, lead_indexv1, lead_indexv2, floe_groups07v1, floe_groups07v2 = ut.get_chords_2versions(sub_dF07, km_unit=True)

#lonsv1T, latsv1T, lengths07v1, lonsv2T, latsv2T, lengths07v2, floe_groups07v1, floe_groups07v2, lead_indexv1, lead_indexv2

# Print output
floe_lengthv1_mean=int(np.mean(pddata_v1.lengths))
floe_lengthv1_std=int(np.std(pddata_v1.lengths))

floe_lengthv2_mean=int(np.mean(pddata_v2.lengths))
floe_lengthv2_std=int(np.std(pddata_v2.lengths))

floe_lengthv1_mean_str='%.02f' %(floe_lengthv1_mean*0.001)
floe_lengthv2_mean_str='%.02f' %(floe_lengthv2_mean*0.001)

print('Chord length stats')
print('v1:', floe_lengthv1_mean, floe_lengthv1_std)
print('v2:', floe_lengthv2_mean, floe_lengthv2_std)


# Lead fraction calculations
# Version 1 - use the SSH flag
seglength_ssh=sub_dF07['seg_length'][sub_dF07['ssh_flag']>0.5]

#print(idx_ssh)
#seglength_ssh = seglength_total[idx_ssh]
#seglength_total 
lead_fraction=np.sum(seglength_ssh)/np.sum(sub_dF07['seg_length'])
lead_fraction_str='%.02f' %(lead_fraction*100)
lead_fraction_str

# Version 2

seglength_sshv2=sub_dF07['seg_length'][lead_indexv2<0]
lead_fraction2=np.sum(seglength_sshv2)/np.sum(sub_dF07['seg_length'])
lead_fraction2_str='%.02f' %(lead_fraction2*100)
lead_fraction2_str

print('Lead fraction')
print('v1:', lead_fraction_str)
print('v2:', lead_fraction2_str)



print('plotting final profile...')
shadingSize=30

fig= plt.figure(figsize=(10,6))

# Plot 1
fig= plt.figure(figsize=(9,9.5))
fig.patch.set_facecolor('xkcd:white')

ax1 = plt.subplot2grid((14, 1), (0, 0), rowspan=2)
plt.imshow(crop, extent=[sub_dF07.along_dist.iloc[0], sub_dF07.along_dist.iloc[-1], -buffer/1000, +buffer/1000], aspect='auto', cmap='gray', vmin=0, vmax=255)
plt.scatter(sub_dF07['along_dist'][sub_dF07['ssh_flag']<0.5],np.zeros(len(sub_dF07['along_dist'][sub_dF07['ssh_flag']<0.5])), c='r', s=buffer/800)
plt.scatter(sub_dF07['along_dist'][sub_dF07['ssh_flag']>0.5],np.zeros(len(sub_dF07['along_dist'][sub_dF07['ssh_flag']>0.5])), c='y', s=buffer/400, zorder=2)

#plt.scatter(sub_dF07['along_dist'],np.zeros(len(sub_dF07['along_dist'])), c=sub_dF07['ssh_flag''], s=buffer/10000, cmap=cmap_ssh)

ax1.xaxis.set_ticklabels([])
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
plt.xlim(ax5.get_xlim())
plt.autoscale(enable=True, axis='x', tight=True)
#plt.suptitle(date, color='k', fontsize=12)
ax1.annotate(date+' '+relStr+' '+beam_name+' ('+beamStrength+')  N(ATL07) segments = '+str(atl07_segs)+' N(ATL10) = '+str(atl10_segs),xy=(0.02, 1.05), xycoords='axes fraction')

# Plot 2
ax2 = plt.subplot2grid((14, 1), (2, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if (row['ssh_flag'] > 0.1):
        plt.axvspan(x0,x1, facecolor='k', alpha=0.1)

plt.scatter(sub_dF07['along_dist'],sub_dF07['BandValue'], c=sub_dF07['BandValue'], s=1, cmap='gray', vmin = 0, vmax = 1.5, zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax2.xaxis.set_ticklabels([])
ax2.set_yticks([0, 1])
plt.ylim(-0.2,1.4)
plt.ylabel('Surface reflectance', color='k', fontsize=12)
ax2.tick_params(axis='y', colors='k', labelsize=12)
plt.xlim(ax5.get_xlim())
# Plot 3


# Plot 4
ax3 = plt.subplot2grid((14, 1), (4, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # Specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    # Dark lead
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['height'], c=sub_dF07['height'], s=1, vmin=0, vmax=3, cmap='hot', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax3.xaxis.set_ticklabels([])
ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
plt.ylabel('Height (m)', color='k', fontsize=12)
ax3.tick_params(axis='y', colors='k', labelsize=12)
plt.ylim(-0.2,4)
plt.xlim(ax5.get_xlim())

ax3.annotate('N='+str(specular_segs)+' (specular segs)',xy=(0.02, 0.87), color='b', xycoords='axes fraction')
ax3.annotate('N='+str(darklead_segs)+' (dark lead segs)',xy=(0.02, 0.75), color='r', xycoords='axes fraction')
ax3.annotate('N='+str(cloud_segs)+' (cloud segs)',xy=(0.02, 0.63), color='y', xycoords='axes fraction')
ax3.annotate('N='+str(ssh07_segs)+' (ATL07 candidate ssh segs)',xy=(0.02, 0.51), xycoords='axes fraction')
# Plot 5
ax4 = plt.subplot2grid((14, 1), (6, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    #x0 = row['along_dist'] - row['seg_length']/2000
    #x1 = row['along_dist'] + row['seg_length']/2000
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['PhRate'], c=sub_dF07['PhRate'], s=1, cmap='cool', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax4.xaxis.set_ticklabels([])

plt.ylabel('Photon Rate (ph/shot)', color='k', fontsize=12)
plt.xlim(ax5.get_xlim())

ax4.tick_params(axis='y', colors='k', labelsize=12)
ax4.yaxis.set_major_locator(MaxNLocator(5))
# Plot 6
ax5 = plt.subplot2grid((14, 1), (8, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['BackgrRate'], c=sub_dF07['BackgrRate'], s=1, cmap='summer', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax5.xaxis.set_ticklabels([])
ax5.yaxis.set_label_position("right")
ax5.yaxis.tick_right()
plt.ylabel('Background Rate (MHz)', color='k', fontsize=12)
plt.xlabel('Along track distance (km)', color='k', fontsize=16)
ax5.tick_params(axis='y', colors='k', labelsize=12)
ax5.tick_params(axis='x', colors='k', labelsize=12)
ax5.yaxis.set_major_locator(MaxNLocator(integer=True))
### Plot 7 (ICESat-2 background rate)
ax6 = plt.subplot2grid((14, 1), (10, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # Specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    # Dark lead
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.5)
        

fi=0
for f1 in floe_groups07v1:
    #idx0=np.where(sub_dF07['SegmentID']==floe_groups[f][0])
    #idx1=np.where(sub_dF07['SegmentID']==floe_groups[f][-1])
    x0=f1[0]/1000.
    x1=f1[-1]/1000.
    #print(x0, x1)
    if (fi % 2) == 0: 
        plt.hlines(y=3.1, xmin=x0, xmax=x1, linewidths=3, color='r', alpha=0.5)
    else:
        plt.hlines(y=3.2, xmin=x0, xmax=x1, linewidths=3, color='r', alpha=0.5)
    fi+=1
    
fi=0
for f2 in floe_groups07v2:
    #idx0=np.where(sub_dF10['height_segment_id']==floe_groups10[f][0])
    #idx1=np.where(sub_dF10['height_segment_id']==floe_groups10[f][-1])
    #x0=sub_dF10['along_dist'].iloc[idx0].item()
    #x1=sub_dF10['along_dist'].iloc[idx1].item()
    #
    x0=f2[0]/1000.
    x1=f2[-1]/1000.
    #print(x0, x1)
    if (fi % 2) == 0: 
        plt.hlines(y=2.8, xmin=x0, xmax=x1, linewidths=3, color='m', alpha=0.5)
    else:
        plt.hlines(y=2.9, xmin=x0, xmax=x1, linewidths=3, color='m', alpha=0.5)
    fi+=1

ax6.annotate(r'Chords',xy=(0.01, 0.88), color='r',xycoords='axes fraction')

ax6.annotate(floe_lengthv1_mean_str+r' km  (C$_l^{v1}$)',xy=(0.01, 0.56), color='r',xycoords='axes fraction')
ax6.annotate(floe_lengthv2_mean_str+r' km  (C$_l^{v2}$)',xy=(0.01, 0.41), color='m',xycoords='axes fraction')
#ax6.annotate(r'C$_l^{v3}$: '+str(floe_lengthv3_mean)+' m',xy=(0.02, 0.26), co∂or='m',xycoords='axes fraction')
#ax6.annotate(r'C$_l^{v4}$: '+str(floe_lengthv4_mean)+' m',xy=(0.02, 0.11), color='m',xycoords='axes fraction')

ax6.annotate(lead_fraction_str+r' %  (L$_f^{v1}$)',xy=(0.01, 0.26), color='k',xycoords='axes fraction')
ax6.annotate(lead_fraction2_str+r' %  (L$_f^{v2}$)',xy=(0.01, 0.11), color='k',xycoords='axes fraction')

plt.ylim(-0.2,4)  
plt.xlim(ax5.get_xlim())
ax6.yaxis.set_ticklabels([])
ax6.xaxis.set_ticklabels([])

### Plot 7 (ICESat-2 background rate)
ax7 = plt.subplot2grid((14, 1), (12, 0), rowspan=2)


for index, row in sub_dF10.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    #if ((row['ssh_flag'] > 1.1) & (row['ssh_flag'] < 2.2)):
    if (row['ssh_flag'] > 0.5):
        plt.axvspan(x0,x1, facecolor='k', alpha=0.1)
        
plt.scatter(sub_dF10['along_dist'],sub_dF10['freeboard'], c=sub_dF10['freeboard'], vmin=0, vmax=3, s=1, cmap='hot', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
plt.ylabel('freeboard (m)', color='k', fontsize=12)
ax7.yaxis.set_label_position("right")
ax7.yaxis.tick_right()
plt.xlabel('Along track distance (km)', color='k', fontsize=12)
ax7.tick_params(axis='y', colors='k', labelsize=12)
ax7.tick_params(axis='x', colors='k', labelsize=12)
ax7.annotate('N='+str(ssh10_segs)+' (ATL10 ssh segs)',xy=(0.02, 0.87), xycoords='axes fraction')
#ax7.annotate('Mean = '+str(mean_fb)+' m / mean > 0.05 = '+str(mean_fb2)+' m',xy=(0.6, 0.84), xycoords='axes fraction')
 
plt.xlim(ax5.get_xlim())
plt.ylim(-0.2,3.8)

plt.subplots_adjust(left=0.07, right=0.94, top=0.96, bottom=0.06, wspace=0, hspace=0)
#plt.tight_layout()
### Save plot to file
plt.savefig(figpath+'/examples/FINAL_plot'+date+beam_name+relStr+str(Slim)+str(Nlim)+'withFloev5.png', dpi=300)



