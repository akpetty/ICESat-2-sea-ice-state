""" plot_atl0710_stats.py
    
    Script to generate plots using Sentinel-2 and ICESat-2 sea ice data
    Initial code written by Alek Petty (02/01/2020) based on earlier scripts by Marco Bagnardi
    
    Input:
        ATL07 data

    Output:
        Profiles of ATL07 data

    Python dependencies:
        See below for the relevant module imports. More information on installation is given in the README file.

    Update history:
        02/01/2020: Version 1.
    
"""

### Import necessary Python modules
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
from pyproj import Proj, transform
from scipy import ndimage, misc
import sys
from glob import glob
import utils as ut

# Use some default matplotlib sizes
ut.reset_matplotlib()

# ## 1. Data ingestion and preparation

################### Input data file names ########################################

relStr='rel003'
data_path = '/cooler/I2-ASAS/'
s2_data_path = '/sea_ice_pso/aapetty/raw_data/'
figPath='../figures/'+relStr+'/profiles/'
savePath='../data/'+relStr+'/profiles/'

beamNum=3
example='rgt'

if example=='waves':
    print('Arctic example')
    beam_name='gt1l'
    date='20190323_NH'

    # ATLAS DATA
    # Take the last one in the case of multiple sub-versions of the granule
    rgt_num='1294'
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/data/ATL07-01_20190323035632_'+rgt_num+'0201_'+relStr[-3:]+'_*.h5')[-1]
    #ATL10_dataFile = glob('/Users/aapetty/DATA/ICESat2/'+relStr+'/ATL10/ATL10-01_20190526*_08820301_002_01.h5')[0]
    Slim = 78.8
    Nlim = 79.7
    Emin = 0
    Emax = 90

if example=='rgt':
    print('Arctic example')
    beam_name='gt1l'
    date='20190323_NH'
    #606, 626, 646, 666, 686, 706
    rgt_num = '1046'
    # ATLAS DATA
    # Take the last one in the case of multiple sub-versions of the granule
    ATL07_dataFile = glob(data_path+relStr+'/ATL07-01/data/ATL07-01_202009*'+rgt_num+'*_'+relStr[-3:]+'_*.h5')[-1]
    #ATL10_dataFile = glob('/Users/aapetty/DATA/ICESat2/'+relStr+'/ATL10/ATL10-01_20190526*_08820301_002_01.h5')[0]
    Slim = 76.8
    Nlim = 88.8
    Emin = -180
    Emax = 180

# Grab file
ATL07 = h5py.File(ATL07_dataFile, 'r')

# Determine appropriate ground track index
beamStrs=ut.get_beam_config(ATL07)
beamStr=beamStrs[beamNum-1]

# Grab data
pd_data = ut.get_atl07_data_beam_extra(ATL07, beamStr)

print('OG Min lon:', np.amin(pd_data['lons']), 'max lon:', np.amax(pd_data['lons']))
print('OG Min lat:', np.amin(pd_data['lats']), 'max lat:', np.amax(pd_data['lats']))

### Subset datasets to limits
sub_dF07 = pd_data[(pd_data['lats'] < Nlim) & (pd_data['lats'] > Slim) &(pd_data['lons'] < Emax) & (pd_data['lons'] > Emin) ]

print('Subset min lon:', np.amin(pd_data['lons']), 'max lon:', np.amax(pd_data['lons']))
print('Subset min lat:', np.amin(pd_data['lats']), 'max lat:', np.amax(pd_data['lats']))


# Reset data frame index value
sub_dF07 = sub_dF07.reset_index(drop=True)
print('Size of subset: ', sub_dF07.shape[0])

# Convert distances from meters to kilometers with respect to first segment in AOI
sub_dF07['along_dist'] = sub_dF07['along_dist'] - sub_dF07['along_dist'].iloc[0]
sub_dF07['along_dist'] = sub_dF07['along_dist'].multiply(0.001)

# Convert background rate to MHz
sub_dF07['bground_rate'] = sub_dF07['bground_rate'].multiply(0.000001)


pddata_v1, pddata_v2, lead_indexv1, lead_indexv2, floe_groups07v1, floe_groups07v2 = ut.get_chords_2versions(sub_dF07, km_unit=True)


# Get some basic stats for output
atl07_segs=np.size(sub_dF07['height'].values)
atl07_segs

ssh07_segs=np.size(np.where(sub_dF07['ssh_flag'] > 0.5))
print('Number of candidate ssh segments:', ssh07_segs)

specular_segs=np.size(np.where((sub_dF07['seg_type'] > 1.5)&(sub_dF07['seg_type'] < 5.5)))
print('Number of specular lead segments:', specular_segs)

darklead_segs=np.size(np.where(sub_dF07['seg_type'] > 5.5))
print('Number of dark lead segments:', darklead_segs)

floe_lengthv1_mean=int(np.mean(pddata_v1.lengths))
floe_lengthv1_std=int(np.std(pddata_v1.lengths))

floe_lengthv2_mean=int(np.mean(pddata_v2.lengths))
floe_lengthv2_std=int(np.std(pddata_v2.lengths))

floe_lengthv1_mean_str='%.02f' %(floe_lengthv1_mean*0.001)
floe_lengthv2_mean_str='%.02f' %(floe_lengthv2_mean*0.001)

print('Chord length stats')
print('v1:', floe_lengthv1_mean, floe_lengthv1_std)
print('v2:', floe_lengthv2_mean, floe_lengthv2_std)

seglength_sshv1=sub_dF07['seg_length'][lead_indexv1<0]
lead_fraction=np.sum(seglength_sshv1)/np.sum(sub_dF07['seg_length'])
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


###------- Generate final plot

# Create binary colormap for SSH Flag
cmap_ssh = colors.ListedColormap(['red', 'yellow'])

# Width of ssh stripes
shadingSize=50

# Plot 1
fig= plt.figure(figsize=(9,5.5))
fig.patch.set_facecolor('xkcd:white')


# Plot 1
ax1 = plt.subplot2grid((7, 1), (0, 0), rowspan=2)

ax1.axhline(y=0, alpha=0.7)
plt.scatter(sub_dF07['along_dist'], sub_dF07['height'], c=sub_dF07['height'], cmap='viridis', s=2, zorder=2)
#plt.plot(sub_dF07['along_dist'],sub_dF07['height'], '.', markersize=2, color='r', label='all segs')

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    #if ((row['SSHFlag'] > 0.1) & (row['SSHFlag'] < 1.2)):
    if (row['ssh_flag'] > 0.1):
        plt.axvspan(x0,x1, facecolor='k', alpha=0.3)

ax1.annotate('Gray v lines = candidate ssh segments, N='+str(ssh07_segs),xy=(0.02, 1.28), xycoords='axes fraction')
ax1.annotate('Blue v lines = specular segments, N='+str(specular_segs),xy=(0.02, 1.16), color='b', xycoords='axes fraction')
ax1.annotate('Red v lines = dark lead segments, N='+str(darklead_segs),xy=(0.02, 1.04), color='r',xycoords='axes fraction')
ax1.annotate(date+' RGT:'+rgt_num+' '+relStr+' '+beamStr+' (beam #'+str(beamNum)+')  ATL07 segments = '+str(atl07_segs)+' ('+str(Slim)+'N - '+str(Nlim)+'N)',xy=(0.02, 1.44), xycoords='axes fraction')


#plt.plot(sub_dF07[sub_dF07['seg_type'] > 1.5]['along_dist'],sub_dF07[sub_dF07['seg_type'] > 1.5]['height'], 'x', markersize=2, color='b', label='possile lead segs')
#plt.plot(sub_dF07[sub_dF07['ssh_flag'] > 0.5]['along_dist'],sub_dF07[sub_dF07['ssh_flag'] > 0.5]['height'], '.', markersize=2, color='r', label='ssh segs')
#plt.plot(xnew,height_ints, '-',  color='b', zorder=2, label='50 km mean of lead segs')



#ax1.legend(loc=1, ncol=3,handlelength=1.5, labelspacing=0.2 , numpoints=1, frameon=False)
#plt.scatter(sub_dF07['along_track'],sub_dF07['Height'], c=sub_dF07['Height'], s=1, cmap='hot', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax1.xaxis.set_ticklabels([])
plt.ylabel('Height (m)', color='k')
ax1.tick_params(axis='y', colors='k')
#plt.ylim(-0.2,3)


# Plot 2
ax2 = plt.subplot2grid((7, 1), (2, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    #x0 = row['along_track'] - row['Length']/2000
    #x1 = row['along_track'] + row['Length']/2000
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.1)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.1)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['photon_rate'], c=sub_dF07['photon_rate'], s=1, cmap='Wistia_r', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax2.xaxis.set_ticklabels([])
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel('Ph Rate (/shot)', color='k')
ax2.tick_params(axis='y', colors='k')


# Plot 2
ax3 = plt.subplot2grid((7, 1), (4, 0), rowspan=2)

for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.1)
    # specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.1)
            
plt.scatter(sub_dF07['along_dist'],sub_dF07['bground_rate'], c=sub_dF07['bground_rate'], s=1, cmap='summer', zorder=2)
plt.autoscale(enable=True, axis='x', tight=True)
ax3.xaxis.set_ticklabels([])
plt.ylabel('Bgr Rate (MHz)', color='k')

ax3.tick_params(axis='y', colors='k')
ax3.tick_params(axis='x', colors='k')

# Plot 3
ax4 = plt.subplot2grid((7, 1), (6, 0), rowspan=1)


for index, row in sub_dF07.iterrows():
    x0 = row['along_dist'] - shadingSize/2000
    x1 = row['along_dist'] + shadingSize/2000
    
    if (row['ssh_flag'] > 0.1):
        plt.axvspan(x0,x1, facecolor='k', alpha=0.3)

    # cloud
    if row['seg_type'] < 0.5:
        plt.axvspan(x0,x1, facecolor='y', alpha=0.5)
    # Specular lead
    elif ((row['seg_type'] > 1.5) & (row['seg_type'] < 5.5)):
        plt.axvspan(x0,x1, facecolor='b', alpha=0.1)
    # Dark lead
    elif row['seg_type'] > 5.5:
        plt.axvspan(x0,x1, facecolor='r', alpha=0.1)
        

fi=0
for f1 in floe_groups07v1:
    #idx0=np.where(sub_dF07['SegmentID']==floe_groups[f][0])
    #idx1=np.where(sub_dF07['SegmentID']==floe_groups[f][-1])
    x0=f1[0]/1000.
    x1=f1[-1]/1000.
    #print(x0, x1)
    if (fi % 2) == 0: 
        plt.hlines(y=3.1, xmin=x0, xmax=x1, linewidths=3, color='c')
    else:
        plt.hlines(y=3.3, xmin=x0, xmax=x1, linewidths=3, color='c')
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
        plt.hlines(y=2.2, xmin=x0, xmax=x1, linewidths=3, color='m')
    else:
        plt.hlines(y=2.4, xmin=x0, xmax=x1, linewidths=3, color='m')
    fi+=1

#ax4.annotate(r'Chords',xy=(0.01, 0.88), color='r',xycoords='axes fraction')

ax4.annotate(floe_lengthv1_mean_str+r' km',xy=(1.01, 0.74), color='c',xycoords='axes fraction')
ax4.annotate(lead_fraction_str+r' %',xy=(1.01, 0.56), color='c',xycoords='axes fraction')
#ax6.annotate(r'C$_l^{v3}$: '+str(floe_lengthv3_mean)+' m',xy=(0.02, 0.26), co∂or='m',xycoords='axes fraction')
#ax6.annotate(r'C$_l^{v4}$: '+str(floe_lengthv4_mean)+' m',xy=(0.02, 0.11), color='m',xycoords='axes fraction')



ax4.annotate(floe_lengthv2_mean_str+r' km',xy=(1.01, 0.32), color='m',xycoords='axes fraction')
ax4.annotate(lead_fraction2_str+r' %',xy=(1.01, 0.1), color='m',xycoords='axes fraction')

plt.ylim(1.5,4)  
plt.autoscale(enable=True, axis='x', tight=True)
plt.xlim(ax3.get_xlim())
ax4.yaxis.set_ticklabels([])
plt.xlabel('Along track distance (km)', color='k')

plt.subplots_adjust(left=0.09, right=0.92, top=0.85, bottom=0.09, wspace=0, hspace=0)
#plt.tight_layout()
### Save plot to file
plt.savefig(figPath+'FINAL_plot'+date+'_rgt'+rgt_num+'_'+beam_name+'_'+relStr+'_'+str(Slim)+'N-'+str(Nlim)+'N.png', dpi=300)

