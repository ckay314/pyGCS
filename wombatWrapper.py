#from PyQt5 import QtWidgets, QtGui, QtCore
from astropy.io import fits
import astropy.units as u
from pyGCS import *
from GCSgui import *
from scipy import ndimage
import sunpy.map
import sys

sys.path.append('/Users/kaycd1/STEREO_Mass/IDLport') 
from secchi_prep import secchi_prep
from wispr_prep import wispr_prep
from lasco_prep import c2_prep, c3_prep
from solohi_prep import solohi_fits2grid

#import matplotlib.image as mpimg

from wombatGUI import releaseTheWombat

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)


fnameA1 = 'fits/20120712_162400_d4c2A.fts'
fnameA2 = 'fits/20120712_172400_d4c2A.fts' 

fnameB1 = 'fits/20120712_162400_d4c2B.fts'
fnameB2 = 'fits/20120712_172400_d4c2B.fts' 


#fnameB1 = '/Users/kaycd1/wombat/fits/C2_22800178.fts'
#fnameB2 = '/Users/kaycd1/wombat/fits/C2_22800186.fts'

#fnameC1 = '/Users/kaycd1/wombat/fits/C3_32048310.fts'
#fnameC2 = '/Users/kaycd1/wombat/fits/C3_32048311.fts'

aPair = [fnameA1, fnameA2]
ims, hdrs = secchi_prep([aPair[0], aPair[1]])
if hdrs[1]['DETECTOR'] not in ['HI1', 'HI2']:
     diff = ims[1] - ims[0]
     diff = sunpy.map.Map(diff, hdrs[1])

bPair = [fnameB1, fnameB2]
ims2, hdrs2 = secchi_prep([bPair[0], bPair[1]])
if hdrs[1]['DETECTOR'] not in ['HI1', 'HI2']:
     diff2 = ims2[1] - ims2[0]
     diff2 = sunpy.map.Map(diff2, hdrs2[1])


releaseTheWombat([[[diff],[hdrs[1]]], [[diff2],[hdrs2[1]]]], nWFs=3)