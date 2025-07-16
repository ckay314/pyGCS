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

fnameB1 = 'fits/20120712_162400_d4c2B.fts'
fnameB2 = 'fits/20120712_172400_d4c2B.fts' 

fnameC1 = 'fits/20120713_004901_s4h1A.fts'
fnameC2 = 'fits/20120713_072901_s4h1A.fts'


#fnameB1 = None
#fnameB2 = None

#fnameC1 = None
#fnameC2 = None



allFiles = []
if (fnameA1 != None) & (fnameA2 != None):
    allFiles.append([fnameA1, fnameA2])
if (fnameB1 != None) & (fnameB2 != None):
    allFiles.append([fnameB1, fnameB2])
if (fnameC1 != None) & (fnameC2 != None):
    allFiles.append([fnameC1, fnameC2])

diffMaps =  []  

for aPair in allFiles:
   
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

