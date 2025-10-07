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


# STEREO
fnameA1 = 'fits/20120712_162400_d4c2A.fts'
fnameA2 = 'fits/20120712_172400_d4c2A.fts' 

#fnameB1 = 'fits/20120712_162400_d4c2B.fts'
#fnameB2 = 'fits/20120712_172400_d4c2B.fts' 

fnameB1 = 'fits/20120713_004901_s4h1A.fts'
fnameB2 = 'fits/20120713_072901_s4h1A.fts'

# LASCO
#fnameB1 = '/Users/kaycd1/wombat/fits/C2_22800178.fts'
#fnameB2 = '/Users/kaycd1/wombat/fits/C2_22800186.fts'

#fnameB1 = '/Users/kaycd1/wombat/fits/C3_32048310.fts'
#fnameB2 = '/Users/kaycd1/wombat/fits/C3_32048311.fts'

# PSP     
#fnameB1 = 'fits/psp_L2_wispr_20210121T055621_V1_1211.fits'
#fnameB2 = 'fits/psp_L2_wispr_20210121T062821_V1_1211.fits'


# SoloHI
file1A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T053609_V02.fits'
file2A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-2ft_20220329T053809_V02.fits'
file3A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-3fg_20220329T054157_V02.fits'
file4A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T055357_V02.fits'

file1B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T052409_V02.fits'
file2B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-2ft_20220329T052609_V02.fits'
file3B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-3fg_20220329T051757_V02.fits'
file4B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T052957_V02.fits'

#fnameC1 = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T052957_V02.fits'
#fnameC2 = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T055357_V02.fits'


# STEREO
aPair = [fnameA1, fnameA2]
ims, hdrs = secchi_prep([aPair[0], aPair[1]])
diff = ims[1] - ims[0]
diff = sunpy.map.Map(diff, hdrs[1])

'''bPair = [fnameB1, fnameB2]
ims2, hdrs2 = secchi_prep([bPair[0], bPair[1]])
diff2 = ims2[1] - ims2[0]
diff2 = sunpy.map.Map(diff2, hdrs2[1])'''
     
# LASCO 
'''bPair = [fnameB1, fnameB2]
ims2, hdrs2 = c3_prep([bPair[0], bPair[1]])
if hdrs2[1]['DETECTOR'] not in ['HI1', 'HI2']:
     diff2 = ims2[1] - ims2[0]
     hdrs2[1]['OBSRVTRY'] = 'SOHO'
     diff2 = sunpy.map.Map(diff2, hdrs2[1])'''
     
# PSP
'''ims2, hdrs2 = wispr_prep([fnameB1, fnameB2], straylightOff=True)
diff2 = ims2[1] - ims2[0]
# detector as int not string so fix that
for hdr in hdrs2:
    hdr['detector'] = str(hdr['detector'])
diff2 = sunpy.map.Map(diff2, hdrs2[1])'''

# SolO Quad
aPair = [[file1B, file2B, file3B, file4B], [file1A, file2A, file3A, file4A]]
im0, hdr0 = solohi_fits2grid(aPair[0])
im1, hdr1 = solohi_fits2grid(aPair[1])
diff2 = im1 - im0
hdrs2 = [hdr0, hdr1] 
diff2 = sunpy.map.Map(diff2, hdrs2[1])

# SolO single (or anything fully prepped)
'''aPair = [fnameC1, fnameC2]
with fits.open(aPair[0]) as hdulist:
    im1  = hdulist[0].data
    hdr1 = hdulist[0].header
with fits.open(aPair[1]) as hdulist:
    im2  = hdulist[0].data
    hdr2 = hdulist[0].header
diff2 = im2 - im1
hdrs2 = [hdr1, hdr2]
diff2 = sunpy.map.Map(diff2, hdrs2[1])'''



releaseTheWombat([[[diff],[hdrs[1]]], [[diff2],[hdrs2[1]]]], nWFs=4)