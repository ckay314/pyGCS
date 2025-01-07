from PyQt5 import QtWidgets, QtGui, QtCore
from astropy.io import fits
from pyGCS import *
from GCSgui import *
#import sunpy
import sunpy.map
#from sunpy.coordinates.ephemeris import get_horizons_coord
#import datetime

# Manually list the files for now
fnameA1 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_162400_14c2A.fts'
fnameA2 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_183900_14c2A.fts'

fnameB1 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_162400_14c2A.fts'
fnameB2 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_202400_14c2A.fts'

fnameC1 = None
fnameC2 = None

allFiles = []
if (fnameA1 != None) & (fnameA2 != None):
    allFiles.append([fnameA1, fnameA2])
if (fnameB1 != None) & (fnameB2 != None):
    allFiles.append([fnameB1, fnameB2])
if (fnameC1 != None) & (fnameC2 != None):
    allFiles.append([fnameC1, fnameC2])

diffMaps =  []  
for aPair in allFiles:
    my_map1 = sunpy.map.Map(aPair[0])
    my_map2 = sunpy.map.Map(aPair[1])
    rd = sunpy.map.Map(my_map2 - my_map1.quantity)
    diffMaps.append(rd)

# Number of points in the wireframe
# [nleg, ncircle, naxis]
#ns =[3,10,31]      
ns =[5,20,50]      


# Pass everything to the GUI -------------------------------------------|
runGCSgui(diffMaps, ns=ns)

