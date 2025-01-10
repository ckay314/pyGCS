from PyQt5 import QtWidgets, QtGui, QtCore
from astropy.io import fits
from pyGCS import *
from GCSgui import *
from scipy import ndimage
#import sunpy
import sunpy.map
import sys

# Manually list the files for now
#fnameB1 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_162400_14c2A.fts'
#fnameB2 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_183900_14c2A.fts'


fnameA1 = 'fits/20120712_162400_d4c2A.fts'
fnameA2 = 'fits/20120712_185400_d4c2A.fts' 

fnameB1 = 'fits/20120712_162400_d4c2B.fts'
fnameB2 = 'fits/20120712_185400_d4c2B.fts' 


#fnameB1 = None
#fnameB2 = None


fnameC1 = 'fits/C3_20120712_163005.fts'
fnameC2 = 'fits/C3_20120712_171805.fts'

allFiles = []
if (fnameA1 != None) & (fnameA2 != None):
    allFiles.append([fnameA1, fnameA2])
if (fnameB1 != None) & (fnameB2 != None):
    allFiles.append([fnameB1, fnameB2])
if (fnameC1 != None) & (fnameC2 != None):
    allFiles.append([fnameC1, fnameC2])

diffMaps =  []  

def reclip(aMap, OGdim):
    myCent = [aMap.reference_pixel.x.to_value(), aMap.reference_pixel.y.to_value()]
    OGx = OGdim.x.to_value()
    OGy = OGdim.y.to_value()
    hwx = OGx / 2 
    hwy = OGy / 2
    myData = aMap.data
    ix1, ix2 = int(myCent[0] - hwx), int(myCent[0] + hwx)
    iy1, iy2 = int(myCent[1] - hwy), int(myCent[1] + hwy)
    reclipData = myData[iy1:iy2, ix1:ix2]
    aMap.meta['crpix1'] = aMap.meta['crpix1'] - ix1
    aMap.meta['crpix2'] = aMap.meta['crpix2'] - iy1
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
        my_map1FR = my_map1F.rotate(missing=0, clip=True)
        # Check the dimensions bc apparently clip does nothing
        if (my_map1.dimensions[0] != my_map1FR.dimensions[0]) or (my_map1.dimensions[1] != my_map1FR.dimensions[1]):
            my_map1FR = reclip(my_map1FR, my_map1.dimensions)
    else:
        my_map1FR = my_map1F

    flData2 = my_map2.data.astype(np.float32)
    if 'crota' in my_map1.meta:
        crota2 = my_map1.meta['crota']
    elif 'crota1' in my_map1.meta:
        crota2 = my_map1.meta['crota1']
    else:
        crota2 = 0
    my_map2F = sunpy.map.Map(flData2, my_map2.meta)
    if np.abs(crota2) > 0.001: 
        my_map2FR = my_map2F.rotate(missing=0, clip=True)
        if (my_map2.dimensions[0] != my_map2FR.dimensions[0]) or (my_map2.dimensions[1] != my_map2FR.dimensions[1]):
           my_map2FR = reclip(my_map2FR, my_map2.dimensions)
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

