from PyQt5 import QtWidgets, QtGui, QtCore

import numpy as np
import pyqtgraph as pg
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from pyGCS import *
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.append('/Users/kaycd1/STEREO_Mass/IDLport') 
from wcs_funs import fitshead2wcs, wcs_get_pixel


global nSats
global CMElat, CMElon, CMEtilt, height, k, ang, satpos



# |------ Things to potentially improve ------|
# Need to test with COR1/ C1
# Need to adapt to HI
# Unpack loop calling pts2proj
# Make image scl mode for each panel so can use diff ones at same time?
# Add toggle to do height on each panel -> velo measurements


occultDict = {'STEREO_SECCHI_COR2':[3,14], 'STEREO_SECCHI_COR1':[1.1,4], 'SOHO_LASCO_C1':[1.1,3], 'SOHO_LASCO_C2':[2,6], 'SOHO_LASCO_C3':[3.7,32], 'STEREO_SECCHI_HI1':[15,80], 'STEREO_SECCHI_HI2':[80,215]} # all in Rsun


def pts2projOLD(pts_in, obs, scale, center=[0,0], occultR=None):
        # Take in a list of points and an observer location and project into pixel coordinates
        # Expect pts as [lat, lon, r] in [deg, deg, x] where x doesn't matter as long as consistent
        # In theory can be in any 3D sphere system as long as consistent across all pts (inc obs) 
        # Scale should either be a single value or [scalex, scaley] !!! yes this is oppo of lat/lon input 
        #   order but this is how CK's brain works
        # occultR is an angle, either arcsec/deg to match scale
    
        # Useful constants 
        rad2arcsec = 206265
        dtor = np.pi / 180.
    
        # check inputs and reformat as 2D array (even if single pt)
        if isinstance(pts_in, list):
            pts_in = np.array(pts_in)
        
        if len(pts_in.shape) == 1:
            pts_in = np.array([pts_in])
        
        # check if scale is single value or array of two
        if isinstance(scale, float) or isinstance(scale, int):
            # single value, set both to it
            sclx = scale
            scly = scale
        else:
            sclx = scale[0]
            scly = scale[1]
        
        # convert all the pts to radians
        pts_lats = pts_in[:,0]*dtor 
        pts_lons = pts_in[:,1]*dtor 
        pts_rs   = pts_in[:,2]

        # convert obs location
        obs_lat = obs[0]*dtor
        obs_lon = obs[1]*dtor
        obs_r   = obs[2]
    
        # define dLon var as short hand
        dLon = pts_lons - obs_lon
    
        # Convert from Stony heliographic (or something similar) to heliocentric cartesian
        # this is with x to right, y up, z toward obs
        x = pts_rs  * np.cos(pts_lats) * np.sin(dLon)
        y = pts_rs * (np.sin(pts_lats) * np.cos(obs_lat) - np.cos(pts_lats)*np.cos(dLon)*np.sin(obs_lat))
        z = pts_rs * (np.sin(pts_lats) * np.sin(obs_lat) + np.cos(pts_lats)*np.cos(dLon)*np.cos(obs_lat))
        
        # Convert to projected vals
        d = np.sqrt(x**2 +  y**2 + (obs_r-z)**2)
        dthetax = np.arctan2(x, obs_r - z) * rad2arcsec / sclx 
        thetax  = dthetax + center[0]
        dthetay = np.arcsin(y/d)* rad2arcsec / scly 
        thetay  = dthetay + center[1]
        
        # Check if we want to throw out the points that would be behind the occulter
        
        if occultR: 
            dProj = np.sqrt(dthetax**2 +  dthetay**2)
            outs = []
            for i in range(len(d)):
                if (dProj[i] > occultR/sclx) or (z[i] > 0):
                    outs.append([thetax[i], thetay[i]])
            outs = np.array(outs)
        else:
            # repackage as array of [ [pixX1, pixY1], [pixX2, pixY2], ...]
            outs = np.array([thetax, thetay]).transpose()
        return outs

def pts2proj(pts_in, obs, scale, mywcs, center=[0,0], occultR=None):
        # Take in a list of points and an observer location and project into pixel coordinates
        # Expect pts as [lat, lon, r] in [deg, deg, x] where x doesn't matter as long as consistent
        # In theory can be in any 3D sphere system as long as consistent across all pts (inc obs) 
        # Scale should either be a single value or [scalex, scaley] !!! yes this is oppo of lat/lon input 
        #   order but this is how CK's brain works
        # occultR is an angle, either arcsec/deg to match scale
    
        # Useful constants 
        rad2arcsec = 206265
        dtor = np.pi / 180.
    
        # check inputs and reformat as 2D array (even if single pt)
        if isinstance(pts_in, list):
            pts_in = np.array(pts_in)
        
        if len(pts_in.shape) == 1:
            pts_in = np.array([pts_in])
        
        # check if scale is single value or array of two
        if isinstance(scale, float) or isinstance(scale, int):
            # single value, set both to it
            sclx = scale
            scly = scale
        else:
            sclx = scale[0]
            scly = scale[1]
        
        # convert all the pts to radians
        pts_lats = pts_in[:,0]*dtor 
        pts_lons = pts_in[:,1]*dtor 
        pts_rs   = pts_in[:,2]

        # convert obs location
        obs_lat = obs[0]*dtor
        obs_lon = obs[1]*dtor
        obs_r   = obs[2]
    
        # define dLon var as short hand
        dLon = pts_lons - obs_lon
    
        # Convert from Stony heliographic (or something similar) to heliocentric cartesian
        # this is with x to right, y up, z toward obs
        x = pts_rs  * np.cos(pts_lats) * np.sin(dLon)
        y = pts_rs * (np.sin(pts_lats) * np.cos(obs_lat) - np.cos(pts_lats)*np.cos(dLon)*np.sin(obs_lat))
        z = pts_rs * (np.sin(pts_lats) * np.sin(obs_lat) + np.cos(pts_lats)*np.cos(dLon)*np.cos(obs_lat))
        
        # Convert to projected vals
        d = np.sqrt(x**2 +  y**2 + (obs_r-z)**2)
        if mywcs['cunit'][0] == 'arcsec':
            rad2unit = rad2arcsec
        elif mywcs['cunit'][0] == 'deg':
            rad2unit = 180. / np.pi
        dthetax = np.arctan2(x, obs_r - z) * rad2unit 
        dthetay = np.arcsin(y/d)* rad2unit 
        
        coord = wcs_get_pixel(mywcs, [dthetax, dthetay], doQuick=False)
        
        # Check if we want to throw out the points that would be behind the occulter
        outs = np.array([coord[0,:], coord[1,:]]).transpose()
        thetax, thetay = coord[0,:], coord[1,:]
        
        if occultR: 
            dProj = np.sqrt(dthetax**2 +  dthetay**2)
            outs = []
            for i in range(len(d)):
                if (dProj[i] > occultR) or (z[i] > 0):
                    outs.append([thetax[i], thetay[i]])
            outs = np.array(outs)
        else:
            # repackage as array of [ [pixX1, pixY1], [pixX2, pixY2], ...]
            outs = np.array([thetax, thetay]).transpose()
        return outs

        
class Ui_MainWindow(object):
    # This generically sets up the main components of the GUI (was largely
    # produced using the pyqt developer)    
    def setupUi(self, MainWindow, nSats):
        # Set up the main window and its properties 
        MainWindow.setObjectName("MainWindow")
        wsize = 375
        MainWindow.resize(260+wsize*nSats, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        
        # Set up one to three plotting windows depending on nSats -----------------------|
        ypos = int((600-wsize)/2)-75
        self.graphWidget1 = pg.PlotWidget(self.centralwidget)
        self.graphWidget1.setGeometry(QtCore.QRect(220, ypos, wsize, wsize))
        self.labelSat1 = QtWidgets.QLabel(self.centralwidget)
        self.labelSat1.setGeometry(QtCore.QRect(220, 20, 300, 16))
        if nSats >=2:
            self.graphWidget2 = pg.PlotWidget(self.centralwidget)
            self.graphWidget2.setGeometry(QtCore.QRect(230+wsize, ypos, wsize, wsize))
            self.labelSat2 = QtWidgets.QLabel(self.centralwidget)
            self.labelSat2.setGeometry(QtCore.QRect(230+wsize, 20, 300, 16))
        if nSats == 3:
            self.graphWidget3 = pg.PlotWidget(self.centralwidget)
            self.graphWidget3.setGeometry(QtCore.QRect(240+2*wsize, ypos, wsize, wsize))
            self.labelSat3 = QtWidgets.QLabel(self.centralwidget)
            self.labelSat3.setGeometry(QtCore.QRect(240+2*wsize, 20, 300, 16))

        # Set up the individual sliders, their text boxes, and their labels -------------|
        # GCS shell parameters here
        # Latitude
        self.sliderLat = QtWidgets.QSlider(self.centralwidget)
        self.sliderLat.setGeometry(QtCore.QRect(30, 130, 160, 22))
        self.sliderLat.setOrientation(QtCore.Qt.Horizontal)
        self.sliderLat.setMinimum(-90)
        self.sliderLat.setMaximum(90)        
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(30, 90, 59, 16))
        self.leLat = QtWidgets.QLineEdit(self.centralwidget)
        self.leLat.setGeometry(QtCore.QRect(30, 110, 113, 21))
        self.leLat.setText('0')         # set to a default value, may get replaced
        # Longitude 
        self.sliderLon = QtWidgets.QSlider(self.centralwidget)
        self.sliderLon.setGeometry(QtCore.QRect(30, 60, 160, 22))
        self.sliderLon.setOrientation(QtCore.Qt.Horizontal)
        self.sliderLon.setMinimum(-180)
        self.sliderLon.setMaximum(180)  
        self.sliderLon.setValue(0) 
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(30, 20, 71, 16))
        self.leLon = QtWidgets.QLineEdit(self.centralwidget)
        self.leLon.setGeometry(QtCore.QRect(30, 40, 113, 21))
        self.leLon.setText('0')        
        # Longitude coord sys button set up
        self.LonButGroup = QtWidgets.QButtonGroup()
        self.stonyBut = QtWidgets.QRadioButton(self.centralwidget)
        self.stonyBut.setChecked(True)
        self.LonButGroup.addButton(self.stonyBut)
        self.stonyBut.setGeometry(QtCore.QRect(100, 20, 15, 15))
        self.stonyButLab = QtWidgets.QLabel(self.centralwidget)
        self.stonyButLab.setGeometry(QtCore.QRect(120, 20, 35, 16))
        self.carrBut = QtWidgets.QRadioButton(self.centralwidget)
        self.LonButGroup.addButton(self.carrBut)
        self.carrBut.setGeometry(QtCore.QRect(150, 20, 15, 15))
        self.carrButLab = QtWidgets.QLabel(self.centralwidget)
        self.carrButLab.setGeometry(QtCore.QRect(170, 20, 35, 16))
        # Tilt
        self.sliderTilt = QtWidgets.QSlider(self.centralwidget)
        self.sliderTilt.setGeometry(QtCore.QRect(30, 200, 160, 22))
        self.sliderTilt.setOrientation(QtCore.Qt.Horizontal)
        self.sliderTilt.setMinimum(-90)
        self.sliderTilt.setMaximum(90) 
        self.sliderTilt.setValue(0)                       
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(30, 160, 71, 16))
        self.leTilt = QtWidgets.QLineEdit(self.centralwidget)
        self.leTilt.setGeometry(QtCore.QRect(30, 180, 113, 21))
        self.leTilt.setText('0')        
        # Height
        self.sliderHeight = QtWidgets.QSlider(self.centralwidget)
        self.sliderHeight.setGeometry(QtCore.QRect(30, 270, 160, 22))
        self.sliderHeight.setOrientation(QtCore.Qt.Horizontal)
        self.sliderHeight.setMinimum(11)
        self.sliderHeight.setMaximum(250)
        self.sliderHeight.setValue(50)                        
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(30, 230, 71, 16))
        self.leHeight = QtWidgets.QLineEdit(self.centralwidget)
        self.leHeight.setGeometry(QtCore.QRect(30, 250, 113, 21))
        self.leHeight.setText('5.0')        
        # Angular Width
        self.sliderAW = QtWidgets.QSlider(self.centralwidget)
        self.sliderAW.setGeometry(QtCore.QRect(30, 340, 160, 22))
        self.sliderAW.setOrientation(QtCore.Qt.Horizontal)
        self.sliderAW.setMinimum(5)
        self.sliderAW.setMaximum(90)
        self.sliderAW.setValue(30)                        
        self.leAW = QtWidgets.QLineEdit(self.centralwidget)
        self.leAW.setGeometry(QtCore.QRect(30, 320, 113, 21))
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(30, 300, 91, 16))
        self.leAW.setText('30')       
        # Kappa/ratio
        self.sliderK = QtWidgets.QSlider(self.centralwidget)
        self.sliderK.setGeometry(QtCore.QRect(30, 420, 160, 22))
        self.sliderK.setOrientation(QtCore.Qt.Horizontal)
        self.sliderK.setMinimum(5)
        self.sliderK.setMaximum(90)
        self.sliderK.setValue(20)                        
        self.leK = QtWidgets.QLineEdit(self.centralwidget)
        self.leK.setGeometry(QtCore.QRect(30, 400, 113, 21))
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(30, 380, 71, 16))
        self.leK.setText('0.20')
        
        # Button to save the GCS parameters ---------------------------------------------|
        self.savePButton = QtWidgets.QPushButton(self.centralwidget)
        self.savePButton.setGeometry(QtCore.QRect(30, 515, 131, 32))

        # Button to save the figure(s) --------------------------------------------------|
        self.saveFButton = QtWidgets.QPushButton(self.centralwidget)
        self.saveFButton.setGeometry(QtCore.QRect(30, 545, 131, 32))
        
        
        # Button for turining the wireframe on or off -----------------------------------|
        self.wireButton = QtWidgets.QPushButton(self.centralwidget)
        self.wireButton.setGeometry(QtCore.QRect(30, 485, 131, 32))
        
        # Sliders and drop menu for scaling parameters ----------------------------------|
        # Matches number of plot windows
        # Drop menu
        self.menuScale = QtWidgets.QComboBox(self.centralwidget)
        self.menuScale.setGeometry(QtCore.QRect(30, 450, 90, 22))
        self.menuScale.addItems(["Linear", "Log", "Sqrt"])
        # Sat1 minimum brightness
        self.slSat1low = QtWidgets.QSlider(self.centralwidget)
        self.slSat1low.setGeometry(QtCore.QRect(220, ypos+wsize+50, 160, 22))
        self.slSat1low.setOrientation(QtCore.Qt.Horizontal)
        self.leSat1low = QtWidgets.QLineEdit(self.centralwidget)
        self.leSat1low.setGeometry(QtCore.QRect(220, ypos+wsize+30, 160, 22))
        self.labelSat1low = QtWidgets.QLabel(self.centralwidget)
        self.labelSat1low.setGeometry(QtCore.QRect(220, ypos+wsize+10, 160, 22))
        # Sat1 maximum brightness        
        self.slSat1hi = QtWidgets.QSlider(self.centralwidget)
        self.slSat1hi.setGeometry(QtCore.QRect(220, ypos+wsize+120, 160, 22))
        self.slSat1hi.setOrientation(QtCore.Qt.Horizontal)
        self.leSat1hi = QtWidgets.QLineEdit(self.centralwidget)
        self.leSat1hi.setGeometry(QtCore.QRect(220, ypos+wsize+100, 160, 22))
        self.labelSat1hi = QtWidgets.QLabel(self.centralwidget)
        self.labelSat1hi.setGeometry(QtCore.QRect(220, ypos+wsize+80, 160, 22))
        if nSats > 1:
            # Sat2 minimum brightness 
            self.slSat2low = QtWidgets.QSlider(self.centralwidget)
            self.slSat2low.setGeometry(QtCore.QRect(220+wsize+10, ypos+wsize+50, 160, 22))
            self.slSat2low.setOrientation(QtCore.Qt.Horizontal)
            self.leSat2low = QtWidgets.QLineEdit(self.centralwidget)
            self.leSat2low.setGeometry(QtCore.QRect(220+wsize+10, ypos+wsize+30, 160, 22))
            self.labelSat2low = QtWidgets.QLabel(self.centralwidget)
            self.labelSat2low.setGeometry(QtCore.QRect(220+wsize+10, ypos+wsize+10, 160, 22))
            # Sat2 maximum brightness 
            self.slSat2hi = QtWidgets.QSlider(self.centralwidget)
            self.slSat2hi.setGeometry(QtCore.QRect(220+wsize+10, ypos+wsize+120, 160, 22))
            self.slSat2hi.setOrientation(QtCore.Qt.Horizontal)
            self.leSat2hi = QtWidgets.QLineEdit(self.centralwidget)
            self.leSat2hi.setGeometry(QtCore.QRect(220+wsize+10, ypos+wsize+100, 160, 22))
            self.labelSat2hi = QtWidgets.QLabel(self.centralwidget)
            self.labelSat2hi.setGeometry(QtCore.QRect(220+wsize+10, ypos+wsize+80, 160, 22))        
        if nSats == 3:
            # Sat3 minimum brightness 
            self.slSat3low = QtWidgets.QSlider(self.centralwidget)
            self.slSat3low.setGeometry(QtCore.QRect(220+2*wsize+20, ypos+wsize+50, 160, 22))
            self.slSat3low.setOrientation(QtCore.Qt.Horizontal)
            self.leSat3low = QtWidgets.QLineEdit(self.centralwidget)
            self.leSat3low.setGeometry(QtCore.QRect(220+2*wsize+20, ypos+wsize+30, 160, 22))
            self.labelSat3low = QtWidgets.QLabel(self.centralwidget)
            self.labelSat3low.setGeometry(QtCore.QRect(220+2*wsize+20, ypos+wsize+10, 160, 22))
            # Sat3 maximum brightness 
            self.slSat3hi = QtWidgets.QSlider(self.centralwidget)
            self.slSat3hi.setGeometry(QtCore.QRect(220+2*wsize+20, ypos+wsize+120, 160, 22))
            self.slSat3hi.setOrientation(QtCore.Qt.Horizontal)
            self.leSat3hi = QtWidgets.QLineEdit(self.centralwidget)
            self.leSat3hi.setGeometry(QtCore.QRect(220+2*wsize+20, ypos+wsize+100, 160, 22))
            self.labelSat3hi = QtWidgets.QLabel(self.centralwidget)
            self.labelSat3hi.setGeometry(QtCore.QRect(220+2*wsize+20, ypos+wsize+80, 160, 22))
                        
        MainWindow.setCentralWidget(self.centralwidget)
        self.retranslateUi(MainWindow, nSats)
  
    def retranslateUi(self, MainWindow, nSats):
        # This takes the generic widgets and renames them what we want -----------------|
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "pyGCS"))
        self.label.setText(_translate("MainWindow", "Latitude"))
        self.label_2.setText(_translate("MainWindow", "Longitude"))
        self.stonyButLab.setText(_translate("MainWindow", "Sto"))
        self.carrButLab.setText(_translate("MainWindow", "Car"))
        self.label_3.setText(_translate("MainWindow", "Tilt"))
        self.label_4.setText(_translate("MainWindow", "Height"))
        self.label_5.setText(_translate("MainWindow", "Half Angle"))
        self.label_6.setText(_translate("MainWindow", "Ratio"))
        #self.saveLabel.setText(_translate("MainWindow", "Right click to save"))
        self.savePButton.setText(_translate("MainWindow", "Save GCS values"))
        self.saveFButton.setText(_translate("MainWindow", "Save Figure(s)"))
        self.wireButton.setText(_translate("MainWindow", "Wireframe On/Off"))
        self.labelSat1.setText(_translate("MainWindow", "Sat 1"))
        self.labelSat1low.setText(_translate("MainWindow", "Min Brightness"))
        self.labelSat1hi.setText(_translate("MainWindow", "Max Brightness"))
        if nSats > 1:
            self.labelSat2.setText(_translate("MainWindow", "Sat 2"))
            self.labelSat2low.setText(_translate("MainWindow", "Min Brightness"))
            self.labelSat2hi.setText(_translate("MainWindow", "Max Brightness"))
        if nSats == 3:
            self.labelSat3.setText(_translate("MainWindow", "Sat 3"))
            self.labelSat3low.setText(_translate("MainWindow", "Sat 3 Min Brightness"))
            self.labelSat3hi.setText(_translate("MainWindow", "Sat 3 Max Brightness"))


class mywindow(QtWidgets.QMainWindow):
    # This takes the generic but properly labeled window and adapts it to ---------------| 
    # our specific needs
    def __init__(self, diffMapsIn, nsIn=[5,20,30]):
        # Set up globals for the number of sats, plotranges, original images
        # the actual images displayed in the GUI, and the wireframe point density
        global nSats, plotranges, imgOrig, imgOut, ns, diffMaps
        ns = nsIn
        diffMaps = diffMapsIn
        #plotranges = plotrangesIn
        imgOrig = []
        imgOut  = []
        for i in range(nSats):
            if diffMaps[i].meta['obsrvtry'] == 'Parker Solar Probe':
                # Not sure about these values, just copied as place holder
                scl2ints = 100/np.median(np.abs(diffMaps[i].data[~np.isnan(diffMaps[i].data)]))
                thisIm = diffMaps[i].data * scl2ints
                thisIm[np.where(thisIm == np.nan)] = 0
            elif 'HI' in diffMaps[i].meta['detector']:
                thisIm = 1000*(diffMaps[i].data -0.5)
            else:
                scl2ints = 100/np.median(np.abs(diffMaps[i].data))
                thisIm = diffMaps[i].data * scl2ints
            imgOrig.append(np.transpose(thisIm))
            imgOut.append(np.transpose(thisIm))
        
        # -------------------------------------------------------------------------------|            
        # ESSENTIAL GUI SETUP THAT I TOTALLY UNDERSTAND! --------------------------------|           
        super(mywindow, self).__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)
        
        # -------------------------------------------------------------------------------|            
        # Get a generic MainWindow then add our labels ----------------------------------|
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self,nSats)
        
        # -------------------------------------------------------------------------------|            
        # Give nice titles (sat+instr) to each plot -------------------------------------|
        for i in range(nSats):
            myhdr = diffMaps[i].meta
            if 'telescop' in myhdr.keys():
                myScope = myhdr['telescop']
            else:
                myScope = myhdr['obsrvtry']
            print (myScope)
            if myScope in ['STEREO', 'STEREO_A', 'STEREO_B']:
                if '_' in myScope:
                    diffMaps[i].meta['telescop'] = 'STEREO'
                mySat = myhdr['obsrvtry'] + ' ' + myhdr['instrume'] + ' ' + myhdr['detector']
                myDate = myhdr['date-avg']
            if myScope == 'SOHO':
                mySat = myhdr['telescop']  + ' '+ myhdr['instrume'] + ' ' + myhdr['detector']
                myDate = myhdr['date-obs'].replace('/','-')+'T'+myhdr['time-obs']
            if myScope == 'Parker Solar Probe':
                mySat   = myhdr['obsrvtry'] + '_' + myhdr['instrume'] + '_HI' + myhdr['detector']
                myDate = myhdr['date-avg']
                
            
            #myDate = myhdr['date']
            # take off decimal secs
            dotIdx = myDate.find('.')
            myDate = myDate[:dotIdx]
            myName = mySat + ' ' + myDate
            if i == 0:
                self.ui.labelSat1.setText(myName)
            if i == 1:
                self.ui.labelSat2.setText(myName)              
            if i == 2:
                self.ui.labelSat3.setText(myName)
                                
        # -------------------------------------------------------------------------------|            
        # Make a mask for the occulter and outside circular FOV -------------------------|
        # Again being lazy and making the same Sun centered approx
        # !!!!! Probly want to add this back in, here or elsewhere
        global masks, innerR, scaleNshifts, wcss
        # Occulter distance for each satellite
        masks = []
        innerR = []
        cents = []
        scaleNshifts = []
        wcss = []
        for idx in range(nSats):  
            diffMap = diffMaps[idx]
            myhdr   = diffMap.meta
            if myhdr['obsrvtry'] != 'Parker Solar Probe':
                myTag   = myhdr['telescop'] + '_' + myhdr['instrume'] + '_' + myhdr['detector']
            else:
                myTag   = myhdr['obsrvtry'] + '_' + myhdr['instrume'] + '_HI' + myhdr['detector']
            print (myTag)
            obsR    = diffMap.observer_coordinate.radius.m / 7e8 # in Rs
            # Going to assume xy scales are equal for now
            obsScl  = diffMap.scale[0].to_value() # arcsec / pix
            
            mySnS = np.zeros([4,2])
            # Reference pixel, which should be img center now
            cx,cy = int(myhdr['crpix1'])-1, int(myhdr['crpix2'])-1
            obsLon = diffMap.observer_coordinate.lon.degree
            obsLat = diffMap.observer_coordinate.lat.degree
            obsR = diffMap.observer_coordinate.radius.m

            # Old version with skycoord
            #skyPt = SkyCoord(obsLon*u.deg, obsLat*u.deg, 1*u.solRad,frame="heliographic_stonyhurst", obstime=diffMap.date)
            #centS = diffMap.wcs.world_to_pixel(skyPt)

            # This is exact match for COR, HI diffs slightly from map version but match to IDL
            myWCS = fitshead2wcs(myhdr)
            wcss.append(myWCS)
            centS = wcs_get_pixel(myWCS, [0.,0.])
            sx, sy = centS[0], centS[1]
                        
            mySnS[0,0], mySnS[0,1] = cx, cy
            mySnS[2,0], mySnS[2,1] = sx, sy
            # 1 Rs in pix
            if 'rsun' in diffMap.meta:
                myRs = diffMap.meta['rsun'] # in arcsec
            else:
                myDist = diffMap.observer_coordinate.radius.m / 7e8
                myRs   = np.arctan2(1, myDist) * 206265
            mySnS[1,0] = myRs/diffMap.scale[0].to_value()
            mySnS[1,1] = myRs/diffMap.scale[1].to_value()

            mask = np.zeros(diffMap.data.shape)
            
            if ('HI' not in diffMaps[idx].meta['detector']) & (diffMaps[idx].meta['obsrvtry'] != 'Parker Solar Probe'):   
                myOccR  = occultDict[myTag][0] # radius of the occulter in Rs
                occRpix = int(myOccR * mySnS[1,0])
                mySnS[3,0] = mySnS[1,0] * occultDict[myTag][0]
                mySnS[3,1] = mySnS[1,1] * occultDict[myTag][0]
                # Fill in a circle around the occulter center
                for i in range(occRpix):
                    j = int(np.sqrt(occRpix**2 - i**2))
                    lowY = np.max([0,cy-j])
                    hiY  = np.min([diffMap.meta['naxis2']-2, cy+j])
                    #print (cx+i, lowY,hiY+1)
                    if cx+i <= diffMap.meta['naxis2']-1:
                        mask[cx+i, lowY:hiY+1] = 1
                    if cx-i >=0:
                        mask[cx-i, lowY:hiY+1] = 1    
            
                # Fill in outside FoV
                outRpix = int(occultDict[myTag][1] * mySnS[1,0]) 
                for i in range(diffMap.meta['naxis1']):
                    myHdist = np.abs(cx-i)
                    if myHdist >= outRpix:
                        mask[i,:] = 1
                    else:
                        possY = int(np.sqrt(outRpix**2 - myHdist**2))
                        lowY = np.max([0,cy - possY])
                        hiY  = np.min([diffMap.meta['naxis2'],cy + possY])
                        mask[i,:lowY+1] = 1
                        mask[i,hiY:] = 1
            
            masks.append(mask)
            scaleNshifts.append(mySnS)
            
        # -------------------------------------------------------------------------------|            
        # Set up the image spots in the GUI and make an array ---------------------------|
        # holding the graphWidgets so we can access elsewhere
        global images 
        images=[]
        image1 = pg.ImageItem()
        self.ui.graphWidget1.addItem(image1)
        self.ui.graphWidget1.setRange(xRange=(0,diffMaps[0].data.shape[0]), yRange=(0,diffMaps[0].data.shape[1]), padding=0)
        self.ui.graphWidget1.hideAxis('bottom')
        self.ui.graphWidget1.hideAxis('left')
        images.append(image1)        
        if nSats > 1:
            image2 = pg.ImageItem()
            self.ui.graphWidget2.addItem(image2)
            self.ui.graphWidget2.setRange(xRange=(0,diffMaps[1].data.shape[0]), yRange=(0,diffMaps[1].data.shape[1]), padding=0)
            self.ui.graphWidget2.hideAxis('bottom')
            self.ui.graphWidget2.hideAxis('left')
            images.append(image2)        
        if nSats ==3:
            image3 = pg.ImageItem()
            self.ui.graphWidget3.addItem(image3)
            self.ui.graphWidget3.setRange(xRange=(0,diffMaps[2].data.shape[0]), yRange=(0,diffMaps[2].data.shape[1]), padding=0)
            self.ui.graphWidget3.hideAxis('bottom')
            self.ui.graphWidget3.hideAxis('left')
            images.append(image3)
            
            
        # -------------------------------------------------------------------------------|            
        # Get the Carrington longitude of the Earth so can convert lon ------------------|
        # assume all times are similar and use only the first panel time ----------------|
        global Elon
        myStony = diffMaps[0].observer_coordinate.lon.deg
        myCarr  = diffMaps[0].carrington_longitude.deg
        Elon = myCarr - myStony
        if Elon > 180:
            Elon -= 360.
        if Elon < -180:
            Elon += 360.
        Elon = float('{:.2f}'.format(Elon))
                
        # -------------------------------------------------------------------------------|            
        # Check if there is an file with previous values to load ------------------------|
        minmaxesIn = [[-9999, -9999],[-9999, -9999], [-9999, -9999]]
        #if (os.path.exists(fname)): 
        #minmaxesIn =self.initShellValues()
        
        # -------------------------------------------------------------------------------|            
        # Initialize CME variables to the defaults or whatever was ----------------------|
        # pulled in from the input file
        global CMElat, CMElon, CMEtilt, height, k, ang, satpo
        CMElon = float(self.ui.leLon.text())
        CMElat = float(self.ui.leLat.text())
        CMEtilt = float(self.ui.leTilt.text())
        height = float(self.ui.leHeight.text())
        ang = float(self.ui.leAW.text())
        k = float(self.ui.leK.text())
                       
        # -------------------------------------------------------------------------------|            
        # Connect the sliders and textboxes to their actions (and each other!) ----------|
        self.ui.sliderLon.valueChanged[int].connect(self.slLon)
        self.ui.leLon.returnPressed.connect(self.allGCSText)
        self.ui.sliderLat.valueChanged[int].connect(self.slLat)
        self.ui.leLat.returnPressed.connect(self.allGCSText)
        self.ui.sliderTilt.valueChanged[int].connect(self.slTilt)
        self.ui.leTilt.returnPressed.connect(self.allGCSText)
        self.ui.sliderHeight.valueChanged[int].connect(self.slHeight)
        self.ui.leHeight.returnPressed.connect(self.allGCSText)
        self.ui.sliderAW.valueChanged[int].connect(self.slAW)
        self.ui.leAW.returnPressed.connect(self.allGCSText)
        self.ui.sliderK.valueChanged[int].connect(self.slK)
        self.ui.leK.returnPressed.connect(self.allGCSText)
        self.ui.slSat1low.valueChanged[int].connect(self.slBmin1)
        self.ui.leSat1low.returnPressed.connect(self.allImText)
        self.ui.slSat1hi.valueChanged[int].connect(self.slBmax1)
        self.ui.leSat1hi.returnPressed.connect(self.allImText)
        if nSats > 1:
                self.ui.slSat2low.valueChanged[int].connect(self.slBmin2)
                self.ui.leSat2low.returnPressed.connect(self.allImText)
                self.ui.slSat2hi.valueChanged[int].connect(self.slBmax2)
                self.ui.leSat2hi.returnPressed.connect(self.allImText)
        if nSats == 3:
                self.ui.slSat3low.valueChanged[int].connect(self.slBmin3)
                self.ui.leSat3low.returnPressed.connect(self.allImText)
                self.ui.slSat3hi.valueChanged[int].connect(self.slBmax3)
                self.ui.leSat3hi.returnPressed.connect(self.allImText)
        
        # -------------------------------------------------------------------------------|            
        # Connect the menu Box ----------------------------------------------------------|
        self.ui.menuScale.currentIndexChanged.connect(self.selectionchange)    
       
        # -------------------------------------------------------------------------------|            
        # Connect the wireframe button --------------------------------------------------|
        global wireShow
        wireShow = True
        self.ui.wireButton.clicked.connect(self.wireOO)

        # -------------------------------------------------------------------------------|            
        # Connect the lon coords toggle -------------------------------------------------|
        global isStony
        isStony = True
        self.ui.stonyBut.toggled.connect(self.swapLonSys)
        
        
        
        # -------------------------------------------------------------------------------|            
        # Connect the save buttons ------------------------------------------------------|
        self.ui.savePButton.clicked.connect(self.saveParams)
        self.ui.saveFButton.clicked.connect(self.saveFigs)
                    
        # -------------------------------------------------------------------------------|            
        # Set up a global for the minmax brightness -------------------------------------|
        # This can either come from the input file or be
        # calculated by initBrange
        global minmaxes 
        minmaxes = np.array([[-1,1], [-1,1], [-1,1]])
        if minmaxesIn[0][0] == -9999:
            for i in range(nSats):
                minmaxes[i] = self.initBrange(imgOut[i],i)
            if 'HI' in diffMaps[i].meta['detector']:   
                minmaxes[i] = [-500,500]                
        else:
            minmaxes = minmaxesIn        
        self.resetBrights(minmaxes)
            
        # -------------------------------------------------------------------------------|            
        # Add the mask to the images, then the images to the ----------------------------|
        # widgets (with appropriate levels)
        for i in range(nSats):   
            # !!!! add back in mask  
            if 'HI' not in diffMaps[i].meta['detector']:            
                imgOut[i][np.where(masks[i] == 1)] = np.min(imgOut[i])
            images[i].setImage(imgOut[i], levels=minmaxes[i])
 
        # -------------------------------------------------------------------------------|            
        # Set up the scatter items that will hold the GCS shells ------------------------|        
        global scatters
        scatters = []
        scatter1 = pg.ScatterPlotItem(pen=pg.mkPen(width=1, color='g'), symbol='o', size=2)
        self.ui.graphWidget1.addItem(scatter1)
        scatters.append(scatter1)
        if nSats > 1:
            scatter2 = pg.ScatterPlotItem(pen=pg.mkPen(width=1, color='g'), symbol='o', size=2)
            self.ui.graphWidget2.addItem(scatter2)
            scatters.append(scatter2)
        if nSats ==3:
            scatter3 = pg.ScatterPlotItem(pen=pg.mkPen(width=1, color='g'), symbol='o', size=2)
            self.ui.graphWidget3.addItem(scatter3)
            scatters.append(scatter3)
            
        # -------------------------------------------------------------------------------|            
        # Draw a white circle showing the outline of the Sun ----------------------------|
        thetas = np.linspace(0, 2*3.14159)
        mySnS = scaleNshifts[0]
        self.ui.graphWidget1.plot(mySnS[1,0]*np.cos(thetas)+mySnS[2,0], mySnS[1,1]*np.sin(thetas)+mySnS[2,1])
        if nSats > 1:
            mySnS = scaleNshifts[1]            
            self.ui.graphWidget2.plot(mySnS[1,0]*np.cos(thetas)+mySnS[2,0], mySnS[1,1]*np.sin(thetas)+mySnS[2,1])
        if nSats==3:
            mySnS = scaleNshifts[2]
            self.ui.graphWidget3.plot(mySnS[1,0]*np.cos(thetas)+mySnS[2,0], mySnS[1,1]*np.sin(thetas)+mySnS[2,1])
        
        # -------------------------------------------------------------------------------|
        # Calculate the GCS shell using pyGCS and add to the figures --------------------|
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])       
        #data = getGCS(40, 10, 0, 10, 0.8, 60)  
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])
    
    
    
    # -----------------------------------------------------------------------------------|            
    # -----------------------------------------------------------------------------------|            
    # -----------------------------------------------------------------------------------|            
    # -----------------------------------------------------------------------------------|            
    
    # Various functions that init will call ---------------------------------------------|
    # Most actions for sliders and text boxes but a few actual things
    
    def swapLonSys(self):
        global isStony, CMElon, Elon
        val = float(self.ui.leLon.text())
        isStony = self.ui.stonyBut.isChecked()  
        if not isStony: # convert from stony to carr
            newval = val + Elon
            #self.ui.sliderLon.setValue(int(newval))
        else: # convert from carr to stony
            newval = val - Elon
        if newval > 180:
            newval -= 360
        elif newval < -180:
            newval += 360
        self.ui.sliderLon.setValue(int(newval))
        self.ui.leLon.setText('{:.2f}'.format(newval))
        if isStony:
            CMElon = val - Elon
        else:
            CMElon = val 
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])
          

    # -----------------------------------------------------------------------------------|            
    def resetBrights(self, minmaxes): #--------------------------------------------------|
        # Take the minmax values and rescale the images
        self.ui.slSat1low.setValue(int(minmaxes[0][0]))
        self.ui.leSat1low.setText(str(int(minmaxes[0][0])))
        self.ui.slSat1hi.setValue(int(minmaxes[0][1]))
        self.ui.leSat1hi.setText(str(int(minmaxes[0][1])))
        if nSats > 1:
            self.ui.slSat2low.setValue(int(minmaxes[1][0]))
            self.ui.leSat2low.setText(str(int(minmaxes[1][0])))
            self.ui.slSat2hi.setValue(int(minmaxes[1][1]))
            self.ui.leSat2hi.setText(str(int(minmaxes[1][1])))
        if nSats == 3:
            self.ui.slSat3low.setValue(int(minmaxes[2][0]))
            self.ui.leSat3low.setText(str(int(minmaxes[2][0])))
            self.ui.slSat3hi.setValue(int(minmaxes[2][1]))
            self.ui.leSat3hi.setText(str(int(minmaxes[2][1])))
    
    # -----------------------------------------------------------------------------------|            
    def initBrange(self, imIn, idx): # --------------------------------------------------|
        # Make a guess at a good range for each plot based on
        # the current scaling method and rescale to that value
        goodIm = imIn[np.isfinite(imIn)]
        absMed = (np.median(np.abs(goodIm[goodIm != 0])))
        if self.ui.menuScale.currentText() == 'Linear':
            Bmin, Bmax = -10*absMed, 10*absMed
            slLow, slHigh = -30*absMed, 30*absMed
        if self.ui.menuScale.currentText() in ['Sqrt', 'Log']:
            Bmin, Bmax = int(absMed), int(1.25*absMed)
            slLow, slHigh = 0, 3*absMed
        if 'HI' in diffMaps[idx].meta['detector']:
            Bmin, Bmax = -250, 250
            slLow, slHigh = -500,500
        # Figure out which slider and box        
        if idx==0:
            sls = [self.ui.slSat1low, self.ui.slSat1hi]    
            les = [self.ui.leSat1low, self.ui.leSat1hi]    
        if idx==1:
            sls = [self.ui.slSat2low, self.ui.slSat2hi]    
            les = [self.ui.leSat2low, self.ui.leSat2hi]    
        if idx==2:
            sls = [self.ui.slSat3low, self.ui.slSat3hi]    
            les = [self.ui.leSat3low, self.ui.leSat3hi]    
        # Reset things                
        sls[0].setMinimum(int(slLow))
        sls[0].setMaximum(int(slHigh))
        sls[1].setMinimum(int(slLow))
        sls[1].setMaximum(int(slHigh))
        sls[0].setValue(int(Bmin)) 
        sls[1].setValue(int(Bmax)) 
        les[0].setText(str(int(Bmin)))
        les[1].setText(str(int(Bmax)))
        return Bmin, Bmax
        
    # -----------------------------------------------------------------------------------|            
    def selectionchange(self): # --------------------------------------------------------#
        # When the menu box changes reprocess imgOut from imgOrig 
        # however is needed.  Also calls initBrange to attempt to
        # make images pretty again
        global imgOut, imgOrig, images, masks
        for i in range(nSats):
            imgOut[i] = np.copy(imgOrig[i])
            # for Linear only need to copy to orig so nothing else...
            if self.ui.menuScale.currentText() == 'Sqrt':
                medval =(np.median(np.abs(imgOrig[i][np.isreal(imgOrig[i])])))
                imgOut[i][np.isreal(imgOrig[i])==False] = medval
                imgOut[i] += int(30.*medval)
                imgOut[i][np.where(imgOut[i]<0)] = 0
                # spread the sqrt out for more integers for slider
                imgOut[i] = 10*np.sqrt(imgOut[i])
            if self.ui.menuScale.currentText() == 'Log':  
                medval =(np.median(np.abs(imgOrig[i][np.isreal(imgOrig[i])])))
                imgOut[i][np.isreal(imgOrig[i])==False] = medval
                imgOut[i] += int(30.*medval)
                imgOut[i][np.where(imgOut[i]<0)] = medval
                # spread out again
                imgOut[i] = 10*np.log(imgOut[i])
                        
            Bmin, Bmax = self.initBrange(imgOut[i],i) 
            # !!! ADD masks back in
            #imgOut[i][np.where(masks[i] == 1)] = np.min(imgOut[i])
            images[i].updateImage(image=imgOut[i], levels=(Bmin, Bmax))
                        
    # -----------------------------------------------------------------------------------|            
    def wireOO(self): #------------------------------------------------------------------|
        # Turn the wireframe on or off
        global wireShow
        if wireShow:
            for i in range(nSats):
                scatters[i].setData() 
            wireShow = False
        else:
            data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])                
            for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])
            wireShow = True 

    # -----------------------------------------------------------------------------------|            
    def saveParams(self): #--------------------------------------------------------------|
        print('Saving output in GCSvals.txt')
        f1 = open('GCSvals.txt', 'w')
        if isStony:
            f1.write('Stony Lon:     '+self.ui.leLon.text()+'\n')
            carLon = float(self.ui.leLon.text()) + Elon
            if carLon > 180:
                carLon -= 360.
            if carLon < -180:
                carLon += 360
            f1.write('Carr Lon:     '+'{:.2f}'.format(carLon) +'\n')
        else:    
            stoLon = float(self.ui.leLon.text()) - Elon
            if stoLon > 180:
                stoLon -= 360.
            if stoLon < -180:
                stoLon += 360
            f1.write('Stony Lon:     '+'{:.2f}'.format(stoLon) +'\n')
            f1.write('Carr Lon:     '+self.ui.leLon.text() +'\n')
        f1.write('Lat:     '+self.ui.leLat.text()+'\n')
        f1.write('Tilt:    '+self.ui.leTilt.text()+'\n')
        f1.write('Height:  '+self.ui.leHeight.text()+'\n')
        f1.write('HalfAng: '+self.ui.leAW.text()+'\n')
        f1.write('Ratio:   '+self.ui.leK.text()+'\n')
        scDict = {'Linear':'0', 'Log':'1', 'Sqrt':'2'}
        f1.write('Scaling: '+scDict[self.ui.menuScale.currentText()]+'\n')
        f1.write('Sat1min: '+self.ui.leSat1low.text()+'\n')
        f1.write('Sat1max: '+self.ui.leSat1hi.text()+'\n')
        if nSats > 1:
            f1.write('Sat2min: '+self.ui.leSat2low.text()+'\n')
            f1.write('Sat2max: '+self.ui.leSat2hi.text()+'\n')
        if nSats == 3:
            f1.write('Sat3min: '+self.ui.leSat3low.text()+'\n')
            f1.write('Sat3max: '+self.ui.leSat3hi.text()+'\n')
        f1.close()

    # -----------------------------------------------------------------------------------|            
    def saveFigs(self): #--------------------------------------------------------------|
        pixmap = self.ui.graphWidget1.grab()
        pixmap.save('GCSframe1.png')
        print ('Saving panel 1 as GCSframe1.png')
        if nSats > 1:
            pixmap2 = self.ui.graphWidget2.grab()
            pixmap2.save('GCSframe2.png')
            print ('Saving panel 2 as GCSframe2.png')
        if nSats == 3:
            pixmap = self.ui.graphWidget3.grab()
            pixmap.save('GCSframe3.png') 
            print ('Saving panel 3 as GCSframe3.png')
                
    # -----------------------------------------------------------------------------------|            
    def plotGCSscatter(self, scatter, data, diffMap, mySnS, mywcs): #----------------------------------|
        # Take the data from pyGCS, shift and scale to match coronagraph
        # range and put it in a happy pyqt format.  If statement will hide
        # any values behind the Sun to try and help projection confusion
        # Can switch full for loop to single pos = [] line to turn off
        pos = []
        
        pixCloud = []
        obs = [diffMap.observer_coordinate.lat.deg, diffMap.observer_coordinate.lon.deg, diffMap.observer_coordinate.radius.m]
        obsScl = [diffMap.scale[0].to_value(), diffMap.scale[1].to_value()]  # arcsec/pix for COR
        if 'HI' in diffMap.meta['detector']: # then assuming deg/pix, true for STEREO
            obsScl = [diffMap.scale[0].to_value()*3600, diffMap.scale[1].to_value()*3600]
        cent   = mySnS[2,:]
        occultR = mySnS[3,0] * diffMap.scale[0].to_value()
        # Can probably unpack this, built pts2proj for arrays
        counter = 0
        for pt in data:
            counter += 1
            # Old version, new fast matches within ~ 1 pix
            # Keepin in here to validate against as needed
            #skyPt = SkyCoord(x=pt[0], y=pt[1], z=pt[2], unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
            # This is what is slowing everything down.
            #myPt2 = diffMap.world_to_pixel(skyPt)
            #if counter < 20:
            #    print ('a', myPt2.x.to_value(), myPt2.y.to_value())
            #pos.append({'pos': [myPt2.x.to_value(), myPt2.y.to_value()]})
           
            r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)
            lat = np.arcsin(pt[2]/r) * 180/np.pi
            lon = np.arctan2(pt[1],pt[0]) * 180 / np.pi
            pt = [lat, lon, r*7e8] 

            # Old approx version, new is more accurate
            '''myPt = pts2projOLD(pt, obs, obsScl, center=cent, occultR=occultR)
            if len(myPt) > 0:          
                pos.append({'pos': [myPt[0][0], myPt[0][1]]})
                if counter < 20:
                    print ('b',myPt, myPt2)'''


            myPt = pts2proj(pt, obs, obsScl, mywcs, center=cent, occultR=occultR)
            if len(myPt) > 0:          
                pos.append({'pos': [myPt[0][0], myPt[0][1]]})
                #if counter < 20:
                #    print ('c',myPt, myPt2)

        scatter.setData(pos)
        
            
        
    # -----------------------------------------------------------------------------------|            
    # -----------------------------------------------------------------------------------|            
    # -----------------------------------------------------------------------------------|            
    # All the slider things -------------------------------------------------------------|
    # brightness minmax for each spacecraft            
    def slBmin1(self, value):
        global minmaxes, images, imgOut
        minmaxes[0][0] = value
        self.ui.leSat1low.setText(str(value))
        images[0].updateImage(image=imgOut[0], levels=minmaxes[0])           
    def slBmax1(self, value):
        global minmaxes, images, imgOut
        minmaxes[0][1] = value
        self.ui.leSat1hi.setText(str(value))
        images[0].updateImage(image=imgOut[0], levels=minmaxes[0])           
    def slBmin2(self, value):
        global minmaxes, images, imgOut
        minmaxes[1][0] = value
        self.ui.leSat2low.setText(str(value))
        images[1].updateImage(image=imgOut[1], levels=minmaxes[1])           
    def slBmax2(self, value):
        global minmaxes, images, imgOut
        minmaxes[1][1] = value
        self.ui.leSat2hi.setText(str(value))
        images[1].updateImage(image=imgOut[1], levels=minmaxes[1])           
    def slBmin3(self, value):
        global minmaxes, images, imgOut
        minmaxes[2][0] = value
        self.ui.leSat3low.setText(str(value))
        images[2].updateImage(image=imgOut[2], levels=minmaxes[2])           
    def slBmax3(self, value):
        global minmaxes, images, imgOut
        minmaxes[2][1] = value
        self.ui.leSat3hi.setText(str(value))
        images[2].updateImage(image=imgOut[2], levels=minmaxes[2])           
        
    # Wireframe values    
    def slLon(self, value):
        global CMElon, isStony
        if isStony:
            CMElon = value
            self.ui.leLon.setText(str(int(CMElon)))
        else:
            CMElon = value - Elon
            if CMElon > 180:
                CMElon -= 360
            elif CMElon < -180:
                CMElon += 360
            self.ui.leLon.setText('{:.2f}'.format(value))
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])      
    def slLat(self, value):
        global CMElat
        CMElat = value
        self.ui.leLat.setText(str(CMElat))
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])          
    def slTilt(self, value):
        global CMEtilt
        CMEtilt = value
        self.ui.leTilt.setText(str(CMEtilt))
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])         
    def slHeight(self, value):
        global height
        # scale the height to a larger range to make the slider
        # have better resolution
        height = 0.1*value
        self.ui.leHeight.setText('{:6.2f}'.format(height))
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])           
    def slAW(self, value):
        global ang
        ang = value
        self.ui.leAW.setText(str(ang))
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])  
    def slK(self, value):
        global k
        # scale the ratio to a larger range to make the slider
        # have better resolution
        k = 0.01*value
        self.ui.leK.setText('{:3.2f}'.format(k))
        data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i], scaleNshifts[i], wcss[i])          
        
    # All the text things ---------------------------------------------------------------|
    def allImText(self):
        global minmaxes, images, imgOut, CMElon, CMElat, CMEtilt, height, ang, k
        # Pull all the brightness values
        minmaxes[0][0] = float(self.ui.leSat1low.text())
        minmaxes[0][1] = float(self.ui.leSat1hi.text())
        minmaxes[1][0] = float(self.ui.leSat2low.text())
        minmaxes[1][1] = float(self.ui.leSat2hi.text())
        minmaxes[2][0] = float(self.ui.leSat3low.text())
        minmaxes[2][1] = float(self.ui.leSat3hi.text())
        # Update brighteness sliders
        self.ui.slSat1low.setValue(minmaxes[0][0])
        self.ui.slSat1hi.setValue(minmaxes[0][1])
        self.ui.slSat2low.setValue(minmaxes[1][0])
        self.ui.slSat2hi.setValue(minmaxes[1][1])
        self.ui.slSat3low.setValue(minmaxes[2][0])
        self.ui.slSat3hi.setValue(minmaxes[2][1])
        # Update images
        images[0].updateImage(image=imgOut[0], levels=minmaxes[0])        
        images[1].updateImage(image=imgOut[1], levels=minmaxes[1])        
        images[2].updateImage(image=imgOut[2], levels=minmaxes[2])        
                
    def allGCSText(self):
        # Pull GCS values
        CMElon = float(self.ui.leLon.text())
        CMElat = float(self.ui.leLat.text())
        CMEtilt = float(self.ui.leTilt.text())
        height = float(self.ui.leHeight.text())
        ang = float(self.ui.leAW.text())
        k = float(self.ui.leK.text())
        # Reset Sliders
        self.ui.sliderLon.setValue(int(CMElon))
        self.ui.sliderLat.setValue(int(CMElat))
        self.ui.sliderTilt.setValue(int(CMEtilt))
        self.ui.sliderHeight.setValue(int(height*10))
        self.ui.sliderAW.setValue(int(ang))
        self.ui.sliderK.setValue(int(k*100))
        # Make new wirerframes and plot
        #data = getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=ns[0], ncirc=ns[1], ncross=ns[2])
        #for i in range(nSats):  self.plotGCSscatter(scatters[i], data, diffMaps[i])

                        
# Simple code to set up and run the GUI -------------------------------------------------|        
def runGCSgui(diffMaps, ns=[3,10,31]):
    # Make an application
    app = QtWidgets.QApplication([])
    # Make a widget
    global nSats
    nSats = len(diffMaps)
    application = mywindow(diffMaps, nsIn=ns)
    # Run it
    application.show()
    # Exit nicely
    sys.exit(app.exec())
    
             
           