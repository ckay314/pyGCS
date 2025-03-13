from PyQt5 import QtWidgets, QtGui, QtCore
from astropy.io import fits
import astropy.units as u
from pyGCS import *
from GCSgui import *
from scipy import ndimage
#import sunpy
import sunpy.map
import sys

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

fnameC1 = 'fits/L20120712_162400.fts'
fnameC2 = 'fits/L20120712_171200.fts'


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

def reclip(aMap, OGmap):
    myCent = [aMap.reference_pixel.x.to_value(), aMap.reference_pixel.y.to_value()]
    OGdim = OGmap.dimensions
    OGx = OGdim.x.to_value()
    OGy = OGdim.y.to_value()
    hwx = OGx / 2 
    hwy = OGy / 2
    myData = aMap.data
    ix1, ix2 = int(myCent[0] - hwx), int(myCent[0] + hwx)
    iy1, iy2 = int(myCent[1] - hwy), int(myCent[1] + hwy)
    reclipData = myData[iy1:iy2, ix1:ix2]
    # Sunpy submap routine only changes crpix/naxis so following that
    aMap.meta['crpix1'] = (aMap.meta['crpix1'] - ix1) 
    aMap.meta['crpix2'] = (aMap.meta['crpix2'] - iy1) 
    aMap.meta['naxis1']  = OGmap.meta['naxis1']
    aMap.meta['naxis2']  = OGmap.meta['naxis2']
    # IF have issues try updating map (bottom_left_coord, dimensions, reference_pixel, top_right_coord)
    # or meta (crval, crota/pc_matrix, crvalA, crpixA, xcen, ycen, pcA )
    reclipMap = sunpy.map.Map(reclipData, aMap.meta)
    return reclipMap
    
    
for aPair in allFiles:
    my_map1 = sunpy.map.Map(aPair[0])
    my_map2 = sunpy.map.Map(aPair[1])
    if my_map1.dimensions != my_map2.dimensions:
        sys.exit('Dimension mismatch between image and base for '+ aPair[0] + ' and ' + aPair[1])
    
    flData = my_map1.data.astype(np.float32)
    my_map1F = sunpy.map.Map(flData, my_map1.meta)
    if 'crota' in my_map1.meta:
        crota = my_map1.meta['crota']
    elif 'crota1' in my_map1.meta:
        crota = my_map1.meta['crota1']
    else:
        crota = 0
     
    if np.abs(crota) > 0.001: 
        # Rotation occurs about crpix (ref pix) which is neither the Sun nor image center
        # After rot sun is centered at crpix-crval/cdelt in fits units (index from 1)
        # and exact match to various sunpy methods of obtaining (accounting for index from 0)
        # If recenter then crpix is in image center (w/indexing diff)
        # Clip defaults to true but also seems to do nothing so not explicitly including
        my_map1FR = my_map1F.rotate(angle=crota*u.deg, missing=0, recenter=True) 
        
        # Not passing an angle to rotate is equiv to rot by angle =crota*u.deg

        # Check the dimensions bc apparently clip does nothing
        if (my_map1.dimensions[0] != my_map1FR.dimensions[0]) or (my_map1.dimensions[1] != my_map1FR.dimensions[1]):
            my_map1FR = reclip(my_map1FR, my_map1)
    else:
        my_map1FR = my_map1F

    flData2 = my_map2.data.astype(np.float32)
    if 'crota' in my_map2.meta:
        crota2 = my_map2.meta['crota']
    elif 'crota1' in my_map2.meta:
        crota2 = my_map2.meta['crota1']
    else:
        crota2 = 0
    my_map2F = sunpy.map.Map(flData2, my_map2.meta)
    if np.abs(crota2) > 0.001: 
        my_map2FR = my_map2F.rotate(missing=0, clip=True)
        if (my_map2.dimensions[0] != my_map2FR.dimensions[0]) or (my_map2.dimensions[1] != my_map2FR.dimensions[1]):
           my_map2FR = reclip(my_map2FR, my_map2)
    else:
        my_map2FR = my_map2F
    
    rd = sunpy.map.Map(my_map2FR - my_map1FR.quantity)
    diffMaps.append(rd)

# Number of points in the wireframe
# [nleg, ncircle, naxis]
#ns =[3,10,31]      
ns =[5,20,50]      


# Pass everything to the GUI -------------------------------------------|
runGCSgui(diffMaps, ns=ns)

