from PyQt5 import QtWidgets, QtGui, QtCore
from astropy.io import fits
import astropy.units as u
from pyGCS import *
from GCSgui import *
from scipy import ndimage
#import sunpy
import sunpy.map
import sys
from skimage import exposure

sys.path.append('/Users/kaycd1/STEREO_Mass/IDLport') 
from secchi_prep import secchi_prep
from wispr_prep import wispr_prep
from lasco_prep import c2_prep, c3_prep
from solohi_prep import solohi_fits2grid

import matplotlib.image as mpimg

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)

# Manually list the files for now

fnameA1 = 'fits/20120712_162400_d4c2A.fts'
fnameA2 = 'fits/20120712_172400_d4c2A.fts' 

#fnameB1 = '/Users/kaycd1/wombat/fits/C2_22800178.fts'
#fnameB2 = '/Users/kaycd1/wombat/fits/C2_22800186.fts'

#fnameC1 = '/Users/kaycd1/wombat/fits/C3_32048310.fts'
#fnameC2 = '/Users/kaycd1/wombat/fits/C3_32048311.fts'

#fnameB1 = 'fits/20120712_162400_d4c2B.fts'
#fnameB2 = 'fits/20120712_172400_d4c2B.fts' 

#fnameA1 = 'fits/20120713_004901_s4h1A.fts'
#fnameA2 = 'fits/20120713_072901_s4h1A.fts'

fnameB1 = 'fits/psp_L2_wispr_20210121T055621_V1_1211.fits'
fnameB2 = 'fits/psp_L2_wispr_20210121T062821_V1_1211.fits'
    
#fnameB1 = None
#fnameB2 = None

#fnameC1 = None
#fnameC2 = None

# solohi
file1A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T053609_V02.fits'
file2A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-2ft_20220329T053809_V02.fits'
file3A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-3fg_20220329T054157_V02.fits'
file4A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T055357_V02.fits'

file1B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T052409_V02.fits'
file2B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-2ft_20220329T052609_V02.fits'
file3B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-3fg_20220329T051757_V02.fits'
file4B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T052957_V02.fits'


fnameC1 = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T052409_V02.fits'
fnameC2 = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T053609_V02.fits'


allFiles = []
if (fnameA1 != None) & (fnameA2 != None):
    allFiles.append([fnameA1, fnameA2])
if (fnameB1 != None) & (fnameB2 != None):
    allFiles.append([fnameB1, fnameB2])
#if (fnameC1 != None) & (fnameC2 != None):
#    allFiles.append([fnameC1, fnameC2])
allFiles.append([[file1B, file2B, file3B, file4B], [file1A, file2A, file3A, file4A]])

diffMaps =  []  

for aPair in allFiles:
    if ('wispr' not in aPair[0]) & ('C2' not in aPair[0]) & ('C3' not in aPair[0]) & ('solo' not in aPair[0]) & (len(aPair[0])!=4):
        ims, hdrs = secchi_prep([aPair[0], aPair[1]])
        if hdrs[1]['DETECTOR'] not in ['HI1', 'HI2']:
            diff = ims[1] - ims[0]
        else:
            # gets unhappy if give it HI in brightness units instead of counts
            # gotta take out the extreme values so full data doesn't get dumped in
            # tiny part of color scale
            diff = ims[1] - ims[0]
            diff = diff / np.median(np.abs(diff))
            dperc = 10
            diff[np.where(diff < np.percentile(diff,dperc))] = np.percentile(diff,dperc)
            diff[np.where(diff > np.percentile(diff,100-dperc))] = np.percentile(diff,100-dperc)
            diff = exposure.equalize_hist(diff)
    elif 'wispr' in aPair[0]:
        ims, hdrs = wispr_prep([aPair[0], aPair[1]], straylightOff=True)
        diff = ims[1] - ims[0]
        # treat detector as int not string so fix that
        for hdr in hdrs:
            hdr['detector'] = str(hdr['detector'])
    elif 'C2' in aPair[0]:
        ims, hdrs = c2_prep(aPair)
        diff = ims[1] - ims[0]
        for hdr in hdrs:
            hdr['obsrvtry'] = 'SOHO'
    elif 'C3' in aPair[0]:
        ims, hdrs = c3_prep(aPair)
        diff = ims[1] - ims[0]
        for hdr in hdrs:
            hdr['obsrvtry'] = 'SOHO'
    elif len(aPair[0]) == 4:
        im0, hdr0 = solohi_fits2grid(aPair[0])
        im1, hdr1 = solohi_fits2grid(aPair[1])
        diff = im1 - im0
        hdrs = [hdr0, hdr1]
    else:
        with fits.open(aPair[0]) as hdulist:
            im1  = hdulist[0].data
            hdr1 = hdulist[0].header
        with fits.open(aPair[1]) as hdulist:
            im2  = hdulist[0].data
            hdr2 = hdulist[0].header
        diff = im2 - im1
        hdrs = [hdr1, hdr2]
        
    
        #for key in hdrs[0]:
        #    print (key, hdrs[0][key])
        
    rd = sunpy.map.Map(diff, hdrs[1])
    diffMaps.append(rd)


# Completely necessary things
'''print (diffMaps[0].data.shape)
img = mpimg.imread('/Users/kaycd1/Downloads/BC1.png')
r, g, b = img[:,:,0], img[:,:,1], img[:,:,2]
gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
hdr = diffMaps[0].meta
diffMaps[0] = sunpy.map.Map(gray[:,:], hdr)

img = mpimg.imread('/Users/kaycd1/Downloads/BC2.png')
r, g, b = img[:,:,0], img[:,:,1], img[:,:,2]
gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
hdr = diffMaps[0].meta
diffMaps[1] = sunpy.map.Map(gray[:,:], hdr)'''

# Number of points in the wireframe
# [nleg, ncircle, naxis]
#ns =[3,10,31]      
ns =[5,20,50]      
#ns =[3,6,21]      


# Pass everything to the GUI -------------------------------------------|
runGCSgui(diffMaps, ns=ns)

