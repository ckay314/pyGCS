import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QGridLayout, QTabWidget, QSlider, QComboBox, QLineEdit, QPushButton
from PyQt5 import QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import wombatWF as wf
import pyqtgraph as pg

from astropy.coordinates import SkyCoord


from wcs_funs import fitshead2wcs, wcs_get_pixel
#from GCSgui import pts2proj

import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('QGridLayout')
slogger.setLevel(logging.ERROR)

# Maintain a limited number of global variables to make passing things easier
# (had to move into releaseTheWombats after function-ifying)
# main window = the parameter window
# pws = array of plot windows
# wfs = array of wireframes in theoryland coords
# bmodes = array of integer background scaling modes, defaults to linear = 0
#global mainwindow, pws, nSats, wfs, nwfs, bmodes


global occultDict
# Nominal radii (in Rs) for the occulters for each instrument. Pulled from google so 
# generally correct (hopefully) but not the most precise
occultDict = {'STEREO_SECCHI_COR2':[3,14], 'STEREO_SECCHI_COR1':[1.1,4], 'SOHO_LASCO_C1':[1.1,3], 'SOHO_LASCO_C2':[2,6], 'SOHO_LASCO_C3':[3.7,32], 'STEREO_SECCHI_HI1':[15,80], 'STEREO_SECCHI_HI2':[80,215]} 

class ParamWindow(QMainWindow):
    def __init__(self, nTabs):
        super().__init__()
        if nTabs > 10:
            sys.exit('Do you really need to fit >10 wireframes at once? If so figure out where the upper limit of 10 is hardcoded in the ParamWindow class in wombatGUI.py')
            
        self.nTabs = nTabs
        self.setWindowTitle('Wombat Parameters')
        self.setGeometry(100, 100, 300, 900)

        # Create a QTabWidget
        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)
        
        # Create individual tabs (pages)
        self.tabs =[]
        
        # Create holder for the WF types
        self.WFtypes = np.zeros(nTabs)
        self.WFnum2type = ['None', 'GCS', 'Torus', 'Sphere', 'Half Sphere', 'Ellipsoid', 'Half Ellipsoid', 'Slab']
        self.WFshort = {'GCS':'GCS', 'Torus':'Tor', 'Sphere':'Sph', 'Half Sphere':'HSph', 'Ellipsoid':'Ell', 'Half Ellipsoid':'HEll', 'Slab':'Slab'}
        
        # Create holder for the WF and the params
        #self.WFs = theWFs #np.array([None for i in range(nTabs)])
        self.WFparams = np.array([np.zeros(10) for i in range(nTabs)])
        
        # holders for the param widgets so we can rm them
        self.WFLays = []
        self.layouts = []

        for i in range(nTabs):
            aTab = QWidget()
            layout, WFlay = self.paramLayout(i)
            aTab.setLayout(layout)
            self.tabs.append(aTab)
            self.tab_widget.addTab(aTab, "Tab" + str(i))
            self.layouts.append(layout)
            self.WFLays.append(WFlay)
        

    #|------------------------------| 
    #|----------- Layout -----------|
    #|------------------------------| 

    def paramLayout(self, i):
        # Start with fake massive layout to setup grid size
        layout = QGridLayout()
        label = QLabel('')
        layout.addWidget(label, 0,0,40,11,alignment=QtCore.Qt.AlignCenter)
        
        # Time slider/label
        Tlabel = QLabel('Time selection ')
        layout.addWidget(Tlabel,0,0,1,11,alignment=QtCore.Qt.AlignCenter)
        
        slider1 = QSlider()
        #slider1.setGeometry(QtCore.QRect(0, 100, 160, 25))
        slider1.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(slider1, 1,1,1,11)

        label = QLabel('Wireframe Type')
        layout.addWidget(label, 2,0,1,11,alignment=QtCore.Qt.AlignCenter)
        
        # Drop down box
        cbox = self.wfComboBox(i)
        layout.addWidget(cbox,3,0,1,11,alignment=QtCore.Qt.AlignCenter)
        
        label = QLabel('Wireframe Parameters')
        layout.addWidget(label, 5,0, 1, 11, alignment=QtCore.Qt.AlignCenter)
        
        # Dummy empty label to space it nicely before making wireframe specific
        #label = QLabel('')
        #layout.addWidget(label, 6,1,20,1)
        
        # Put a layout within the layout for the slider friends
        # It's like inception but without Elliot Page explaining everything
        WFLay = QGridLayout()
        layout.addLayout(WFLay, 7,1,25,11)
        
        
        # Background mode drop down box
        label = QLabel('Background Scaling Type')
        layout.addWidget(label, 40,0,1,11,alignment=QtCore.Qt.AlignCenter)
        cbox = self.bgComboBox()
        layout.addWidget(cbox,41,0,1,11,alignment=QtCore.Qt.AlignCenter)
        
        label = QLabel('')
        layout.addWidget(label, 40,0,2,11,alignment=QtCore.Qt.AlignCenter)
        

        saveBut = QPushButton('Save')
        saveBut.released.connect(self.SBclicked)
        layout.addWidget(saveBut, 43, 0, 1,3,alignment=QtCore.Qt.AlignCenter)

        massBut = QPushButton('Mass')
        massBut.released.connect(self.MBclicked)
        layout.addWidget(massBut, 43, 4, 1,3,alignment=QtCore.Qt.AlignCenter)

        # Add things at the bottom
        exitBut = QPushButton('Exit')
        exitBut.released.connect(self.EBclicked)
        exitBut.setStyleSheet("background-color: red")
        layout.addWidget(exitBut, 43, 8, 1,3,alignment=QtCore.Qt.AlignCenter)
        
        return layout, WFLay
    
    def wfComboBox(self,i):
        cbox = QComboBox()
        cbox.addItem('|----None/Select One----|')
        cbox.addItem('GCS')
        cbox.addItem('Torus')
        cbox.addItem('Sphere')
        cbox.addItem('Half Sphere')
        cbox.addItem('Ellipse')
        cbox.addItem('Half Ellipse')
        cbox.addItem('Slab')
        cbox.currentIndexChanged.connect(lambda x: self.cb_index_changed(x,i))
        return cbox
        
    def bgComboBox(self):
        cbox = QComboBox()
        cbox.addItem('Linear')
        cbox.addItem('Log')
        cbox.addItem('SQRT')
        cbox.currentIndexChanged.connect(self.back_changed)
        return cbox
        
    def WFparamLayout(self, myWF):
        WFLay = QGridLayout()
        widges = [[], []]
        i2f = []
        nSliders = 201
        for i in range(9):
            if i < len(myWF.labels):
                # Get the conversion factor between the slider integers and float vals
                myRng = myWF.ranges[i]
                i2f.append((myRng[1] - myRng[0]) / (nSliders - 1))
                # Label
                myDef = myWF.params[i]
                label = QLabel(myWF.labels[i]) 
                WFLay.addWidget(label, 3*i,1,1,3)   
                # Text box
                wBox = QLineEdit()
                #wBox.setText(str(myDef)) 
                WFLay.addWidget(wBox, 3*i,7,1,3)  
                widges[0].append(wBox)
                
                # Slider
                #myRng = myWF.ranges[i]
                slider = QSlider()
                slider.setOrientation(QtCore.Qt.Horizontal)
                slider.setRange(0,nSliders)
                #slider.setValue(int(myDef))
                WFLay.addWidget(slider, 3*i+1,1,1,9)  
                widges[1].append(slider)
                
            else:
                # Throw in the same size things but hidden to make it happy
                label = QLabel('')
                WFLay.addWidget(label, 3*i,1,1,3)   
                label = QLabel('')
                WFLay.addWidget(label, 3*i,1,1,3)   
                # Need to add something same height as box that doesnt complain
                #wBox = QLineEdit()
                #wBox.setFixedSize(1, 30)
                #WFLay.addWidget(wBox, 3*i,7,0,1)
                slider = QSlider()
                slider.setOrientation(QtCore.Qt.Horizontal)
                # Setting these to zero makes it disappear
                slider.setMinimum(0)
                slider.setMaximum(0)
                WFLay.addWidget(slider, 3*i+1,1,1,9)  
        
        # Need to do this explicit for each one for some reason other wise gets
        # upset about the looped index variable
        # Params 1-4 always happen
        widges[1][0].valueChanged.connect(lambda x: self.s2b(x, widges[0][0], i2f[0], myWF.ranges[0][0], myWF, widges))  
        widges[0][0].returnPressed.connect(lambda: self.b2s(widges[1][0], widges[0][0], i2f[0], myWF.ranges[0][0],nSliders, myWF, widges))     
        widges[1][1].valueChanged.connect(lambda x: self.s2b(x, widges[0][1], i2f[1], myWF.ranges[1][0], myWF, widges))  
        widges[0][1].returnPressed.connect(lambda: self.b2s(widges[1][1], widges[0][1], i2f[1], myWF.ranges[1][0],nSliders, myWF, widges))
        widges[1][2].valueChanged.connect(lambda x: self.s2b(x, widges[0][2], i2f[2], myWF.ranges[2][0], myWF, widges))  
        widges[0][2].returnPressed.connect(lambda: self.b2s(widges[1][2], widges[0][2], i2f[2], myWF.ranges[2][0],nSliders, myWF, widges))
        widges[1][3].valueChanged.connect(lambda x: self.s2b(x, widges[0][3], i2f[3], myWF.ranges[3][0], myWF, widges))  
        widges[0][3].returnPressed.connect(lambda: self.b2s(widges[1][3], widges[0][3], i2f[3], myWF.ranges[3][0],nSliders, myWF, widges))
        # Have to check remaining
        myNP = len(myWF.labels)
        # At least 5 params
        if myNP > 4:
            widges[1][4].valueChanged.connect(lambda x: self.s2b(x, widges[0][4], i2f[4], myWF.ranges[4][0], myWF, widges))  
            widges[0][4].returnPressed.connect(lambda: self.b2s(widges[1][4], widges[0][4], i2f[4], myWF.ranges[4][0],nSliders, myWF, widges))
        # At least 6 params    
        if myNP > 5:
            widges[1][5].valueChanged.connect(lambda x: self.s2b(x, widges[0][5], i2f[5], myWF.ranges[5][0], myWF, widges))  
            widges[0][5].returnPressed.connect(lambda: self.b2s(widges[1][5], widges[0][5], i2f[5], myWF.ranges[5][0],nSliders, myWF, widges))
        # At least 7 params    
        if myNP > 6:
            widges[1][6].valueChanged.connect(lambda x: self.s2b(x, widges[0][6], i2f[6], myWF.ranges[6][0], myWF, widges))  
            widges[0][6].returnPressed.connect(lambda: self.b2s(widges[1][6], widges[0][6], i2f[6], myWF.ranges[6][0],nSliders, myWF, widges))
        # At least 8 params           
        if myNP > 7:
            widges[1][7].valueChanged.connect(lambda x: self.s2b(x, widges[0][7], i2f[7], myWF.ranges[7][0], myWF, widges))  
            widges[0][7].returnPressed.connect(lambda: self.b2s(widges[1][7], widges[0][7], i2f[7], myWF.ranges[7][0],nSliders, myWF, widges))
        # At least 9 params    
        if myNP > 8:
            widges[1][8].valueChanged.connect(lambda x: self.s2b(x, widges[0][8], i2f[8], myWF.ranges[8][0], myWF, widges))  
            widges[0][8].returnPressed.connect(lambda: self.b2s(widges[1][8], widges[0][8], i2f[8], myWF.ranges[8][0],nSliders, myWF, widges))
            
        # Set things to the values the WF has
        for i in range(myNP):
            myVal = myWF.params[i]
            slidx = int((myVal - myWF.ranges[i][0])/ i2f[i])
            if slidx > nSliders -1:
                slidx = nSliders -1
                myVal = myWF.ranges[i][1]
            elif slidx < 0:
                slidx = 0
                myVal = myWF.ranges[i][0]
            widges[1][i].setValue(slidx)
            widges[0][i].setText(str(myVal))
        
        
        
         
        # This isn't exact same across all WF but close
        for i in range(WFLay.rowCount()):
                WFLay.setRowStretch(i, 1)
        return WFLay, widges        
        
   
    #|------------------------------| 
    #|----------- Events -----------|
    #|------------------------------| 
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Q: 
            self.close()
        elif event.key() == QtCore.Qt.Key_Escape:
            sys.exit()
            
    def s2b(self, x=None, b=None, dx=None, x0=None, myWF=None, widges=None):
        myVal = x0+dx*x
        if dx < 0.005:
            myStr = '{:.3f}'.format(myVal) 
        else:
            myStr = '{:.2f}'.format(myVal) 
        b.setText(myStr)
        self.updateWFpoints(myWF, widges)

    #def b2s(self, x=None, s=None):
    #    s.setValue(int(x))
    def b2s(self,s,b, dx=None, x0=None, nSli=None, myWF=None, widges=None):
        temp = b.text()
        slidx = int((float(b.text()) - x0)/dx)
        if slidx > nSli -1:
            slidx = nSli -1
            temp = str(x0 + (nSli-1) * dx)
            print('Value too high, maximum value allowed is ', temp)
        elif slidx < 0:
            slidx = 0
            temp = str(x0)
            print('Value too low, minimum value allowed is ', x0)
        s.setValue(slidx)
        # Reset it to what we actual wanted instead of slider rounded val
        # since this will trigger s2b as it hits valueChanged
        b.setText(temp)
        self.updateWFpoints(myWF, widges)
        
    def cb_index_changed(self, a='None',idx=-10):
        self.WFtypes[idx] = a

        # Check if making a new wf
        if type(wfs[idx].WFtype) == type(None):
            myType = self.WFnum2type[a]
            wfs[idx] = wf.wireframe(myType)
            self.tab_widget.setTabText(idx,self.WFshort[myType])
            
            WFLay, widges = self.WFparamLayout(wfs[idx])
                        
            self.layouts[idx].addLayout(WFLay, 7,1,30,11)
            self.WFLays[idx] = WFLay

            
        elif a == 0:
            # Set back to none if didn't select a WF
            wfs[idx] = wf.wireframe(None)
            self.tab_widget.setTabText(idx,'None')
            thisLay = self.cleanLayout(self.WFLays[idx])
            for aPW in pws:
                aPW.scatters[idx].setData([])
        else:
            # Create a new wf object but pass it any matching
            # parameters from the previous version
            ogLabs = wfs[idx].labels
            ogParams = wfs[idx].params
            myType = self.WFnum2type[a]
            newWF = wf.wireframe(self.WFnum2type[a])
            newLabs = newWF.labels
            for iii in range(len(ogLabs)):
                aLab = ogLabs[iii]
                if aLab in newLabs:
                    pidx = np.where(newLabs == aLab)[0]
                    newWF.params[pidx] = ogParams[iii]
                    
            # Change the tab text        
            self.tab_widget.setTabText(idx,self.WFshort[self.WFnum2type[a]])
            
            # Update the slider layout
            thisLay = self.cleanLayout(self.WFLays[idx])
            WFLay, widges = self.WFparamLayout(newWF)
            self.layouts[idx].addLayout(WFLay, 7,1,30,11)
            self.WFLays[idx] = WFLay
           
            # Give the structure the new wf
            wfs[idx] = newWF
            
        for aPW in pws:
            aPW.plotBackground()
            aPW.plotWFs()
            
    def back_changed(self,text):
        for aPW in pws:
             aPW.cbox.setCurrentIndex(text)         
            
    def EBclicked(self):
        sys.exit()

    def SBclicked(self):
        print('Save not coded yet')

    def MBclicked(self):
        print('Mass not coded yet')
            
    #|------------------------------| 
    #|----------- Others -----------|
    #|------------------------------| 
    def cleanLayout(self,lay):
        for i in reversed(range(lay.count())): 
            item = lay.takeAt(i)
            widget = item.widget()
            widget.deleteLater()
        return lay
    
    def updateWFpoints(self, aWF, widges):
        # Got to check if all the points are set or this
        # will blow up on the first run through before panel is built
        flagIt = False
        for i in range(len(widges[0])):
            if widges[0][i].text() != '':
                aWF.params[i] = float(widges[0][i].text())
            else:
                flagIt = True
        if not flagIt:
            aWF.getPoints()
            for ipw in range(nSats):
                pws[ipw].plotWFs()
    
        
class FigWindow(QWidget):
    def __init__(self, satName, myObs, myScls, satStuff, myNum=0):
        super().__init__()
        # Obs are [[ims], [hdrs]] as passed to the GUI (prob MSB)
        # Scls are [[lin, log, sqrt]] on -100,100 range
        
        self.setGeometry(550*(myNum+1), 350, 350, 450)
        
        self.satStuff = satStuff
        self.satName = satStuff[0]['OBS'] +' '+ satStuff[0]['INST']
        self.setWindowTitle(self.satName)
        
        self.OGims = myObs[0]
        self.hdrs = myObs[1]
        self.myScls2 = myScls
        self.tidx = 0
        self.sclidx = 0

        layoutP =  QGridLayout()
        
        self.pWindow = pg.PlotWidget()
        self.pWindow.setMinimumSize(400, 400)
        layoutP.addWidget(self.pWindow,0,0,11,11,alignment=QtCore.Qt.AlignCenter)
        
        # make an image item
        self.image = pg.ImageItem()
        self.pWindow.addItem(self.image)
        self.pWindow.setRange(xRange=(0,myObs[0][0].data.shape[0]), yRange=(0,myObs[0][0].data.shape[1]), padding=0)
        
        # Hide the axes
        self.pWindow.hideAxis('bottom')
        self.pWindow.hideAxis('left')
        
        # Set up scatters for the WF points so can adjust without clearing
        self.scatters = []
        for i in range(nwfs):
            aScat = pg.ScatterPlotItem(pen=pg.mkPen(width=1, color='g'), brush=pg.mkBrush(color='g'),symbol='o', size=2.5)
            self.scatters.append(aScat)
            self.pWindow.addItem(aScat)
        
        # Background mode drop down box
        label = QLabel('Background Scaling Type')
        layoutP.addWidget(label, 12,0,1,5,alignment=QtCore.Qt.AlignCenter)
        self.cbox = self.bgComboBox()
        layoutP.addWidget(self.cbox,12,5,1,5,alignment=QtCore.Qt.AlignCenter)
        
        # Min/max brightness sliders
        minL = QLabel('Min Value:     ')
        layoutP.addWidget(minL, 13,0,1,9)
        self.MinSlider = QSlider()
        self.MinSlider.setOrientation(QtCore.Qt.Horizontal)
        self.MinSlider.setMinimum(-50)
        self.MinSlider.setMaximum(150)
        self.MinSlider.setValue(0)
        layoutP.addWidget(self.MinSlider, 13,3,1,9)
        self.MinSlider.valueChanged.connect(lambda x: self.s2l(x, minL, 'Min Value: '))  
        
        maxL = QLabel('Max Value:     ')
        layoutP.addWidget(maxL, 15,0,1,9)
        self.MaxSlider = QSlider()
        self.MaxSlider.setOrientation(QtCore.Qt.Horizontal)
        self.MaxSlider.setMinimum(-50)
        self.MaxSlider.setMaximum(150)
        self.MaxSlider.setValue(100)
        layoutP.addWidget(self.MaxSlider, 15,3,1,9)
        self.MaxSlider.valueChanged.connect(lambda x: self.s2l(x, maxL, 'Max Value: '))  
        
        saveBut = QPushButton('Save')
        saveBut.released.connect(self.SBclicked)
        layoutP.addWidget(saveBut, 17, 0, 1,3,alignment=QtCore.Qt.AlignCenter)

        massBut = QPushButton('Mass')
        massBut.released.connect(self.MBclicked)
        layoutP.addWidget(massBut, 17, 4, 1,3,alignment=QtCore.Qt.AlignCenter)

        # Add things at the bottom
        exitBut = QPushButton('Exit')
        exitBut.released.connect(self.EBclicked)
        exitBut.setStyleSheet("background-color: red")
        layoutP.addWidget(exitBut, 17, 8, 1,3,alignment=QtCore.Qt.AlignCenter)
        
        self.setLayout(layoutP)
        
        self.plotBackground()
        
        
    #|------------------------------| 
    #|----------- Layout -----------|
    #|------------------------------| 

    def bgComboBox(self):
        cbox = QComboBox()
        cbox.addItem('Linear')
        cbox.addItem('Log')
        cbox.addItem('SQRT')
        cbox.currentIndexChanged.connect(self.back_changed)
        return cbox
        
    #|------------------------------| 
    #|----------- Events -----------|
    #|------------------------------| 
    
    def s2l(self, x=None, l=None, pref=None):
        l.setText(pref + str(x))
        self.plotBackground()

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Q: 
            self.close()
        elif event.key() == QtCore.Qt.Key_Escape:
            sys.exit()
    

    def back_changed(self,text):
        self.sclidx = text     
        self.plotBackground()
            
       
    def EBclicked(self):
        sys.exit()

    def SBclicked(self):
        print('Save not coded yet')

    def MBclicked(self):
        print('Mass not coded yet')
        
    #|------------------------------| 
    #|----------- Others -----------|
    #|------------------------------| 
    def plotWFs(self, justN=0):
        # This is the slow version but keeping syntax around if need to check anything
        #skyPt = SkyCoord(x=0, y=0, z=1, unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
        #myPt2 = self.OGims[self.tidx].world_to_pixel(skyPt)
        
        for i in range(nwfs):
            if type(wfs[i].WFtype) != type(None):                
                # Set up the input parameters for pts2proj
                pos = []
                obs = self.satStuff[self.tidx]['POS']
                obsScl = [self.satStuff[self.tidx]['SCALE'], self.satStuff[self.tidx]['SCALE']]
                if 'HI' in self.satStuff[self.tidx]['INST']:
                    obsScl = [self.satStuff[self.tidx]['SCALE'] * 3600, self.satStuff[self.tidx]['SCALE'] * 3600]
                cent = self.satStuff[self.tidx]['SUNPIX']
                occultR = self.satStuff[self.tidx]['OCCRARC']
                mywcs  = self.satStuff[self.tidx]['WCS']
                myColor =wfs[i].WFcolor
                for jj in range(len(wfs[i].points[:,0])):
                    pt = wfs[i].points[jj,:]
                    r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)
                    lat = np.arcsin(pt[2]/r) * 180/np.pi
                    lon = np.arctan2(pt[1],pt[0]) * 180 / np.pi
                    pt = [lat, lon, r*7e8]
                    myPt = pts2proj(pt, obs, obsScl, mywcs, center=cent, occultR=occultR)
                    if len(myPt) > 0:          
                        pos.append({'pos': [myPt[0][0], myPt[0][1]], 'pen':{'color':myColor, 'width':1}, 'pen':{'color':myColor}})
                        
                    #pos.append({'pos': [ys[jj], zs[jj]], 'pen':{'color':myColor, 'width':2}})
                self.scatters[i].setData(pos)
 
    def plotBackground(self):
        #self.ui.graphWidget1.addItem(image1)
        myIm = self.myScls2[self.tidx][self.sclidx]
        slMin = self.MinSlider.value()
        slMax = self.MaxSlider.value()
        self.image.updateImage(image=myIm, levels=(slMin, slMax))
        self.pWindow.plot(self.satStuff[self.tidx]['SUNCIRC'][0], self.satStuff[self.tidx]['SUNCIRC'][1])
        self.pWindow.plot(self.satStuff[self.tidx]['SUNNORTH'][0], self.satStuff[self.tidx]['SUNNORTH'][1], symbolSize=3, symbolBrush='w', pen=pg.mkPen(color='w', width=1))

def makeNiceMMs(imIn, hdr):
    # Clean out any nans
    im = imIn.data
    im[np.where(im == np.nan)] = 0
    # Transpose it 
    im = np.transpose(im)
    
    medval = np.median(np.abs(im))
    
    # Linear Image
    sclLin = 1 / medval
    linIm  = im*sclLin
    
    
    # Log Im
    sclMed = np.median(np.abs(linIm))
    cpLin = np.copy(linIm)
    minVal = np.min(np.abs(linIm))
    cpLin[np.where(cpLin < minVal)] = minVal
    logIm = np.log(np.copy(cpLin))
    perc95 = np.percentile(logIm, 95)
    logIm = 100 * logIm / perc95   
   
    # SQRT im
    # mostly the same as log prep
    sqrtIm = np.sqrt(cpLin)
    percX = np.percentile(sqrtIm, 90)
    sqrtIm = 100 * sqrtIm / percX
    
    # Scale the lin img down here so others can use before
    percX = np.percentile(linIm, 90)
    linIm = 100 * linIm / percX
    
    
    sclIms = [linIm, logIm,  sqrtIm]
    return sclIms
    
def getSatStuff(imMap):
    # Set up satDict as a micro header that we will package the mask in
    # Keys are OBS, INST, POS, SCALE, CRPIX, WCS, SUNPIX, ONERSUN, MASK
    # OCCRPIX, OCCRARC, SUNCIRC, SUNNORTH
    satDict = {}
    
    # Get the name 
    myhdr   = imMap.meta
    satDict['OBS'] =  myhdr['obsrvtry']
    
    # PSP format
    if myhdr['obsrvtry'] == 'Parker Solar Probe':
        satDict['OBS'] =  myhdr['obsrvtry']
        satDict['INST'] =  myhdr['instrume'] + '_HI' + myhdr['detector']
        myTag   = myhdr['obsrvtry'] + '_' + myhdr['instrume'] + '_HI' + myhdr['detector']
    # SolO format
    elif myhdr['obsrvtry'] == 'Solar Orbiter':
        satDict['OBS'] =  myhdr['obsrvtry']
        satDict['INST'] = myhdr['instrume'] 
        myTag   = myhdr['obsrvtry'] + '_' + myhdr['instrume']
    elif myhdr['telescop'] == 'STEREO':
        satDict['OBS'] =  myhdr['obsrvtry'] 
        satDict['INST'] = myhdr['instrume'] + '_' + myhdr['detector']
        myTag   = myhdr['telescop'] + '_' + myhdr['instrume'] + '_' + myhdr['detector']
    # Other less picky sats
    else:
        satDict['OBS'] =  myhdr['telescop']
        satDict['INST'] = myhdr['instrume'] + '_' + myhdr['detector']
        myTag   = myhdr['telescop'] + '_' + myhdr['instrume'] + '_' + myhdr['detector']

    # Get satellite info    
    obsLon = imMap.observer_coordinate.lon.degree
    obsLat = imMap.observer_coordinate.lat.degree
    obsR = imMap.observer_coordinate.radius.m
    satDict['POS'] = [obsLat, obsLon,  obsR]
    
    # Plate scale in arcsec/pix
    # Check to make sure same in x/y since we will assume as much
    if (imMap.scale[0].to_value() != imMap.scale[1].to_value()):
        sys.exit('xy scales not equilent. Not set up to handle this')    
    obsScl  = imMap.scale[0].to_value()
    satDict['SCALE'] = obsScl
    
    # Reference pixel
    cx,cy = int(myhdr['crpix1'])-1, int(myhdr['crpix2'])-1
    satDict['CRPIX'] = [cx, cy]
    
    # Get WCS header and pix location of Sun
    myWCS = fitshead2wcs(myhdr)
    satDict['WCS'] = myWCS
    centS = wcs_get_pixel(myWCS, [0.,0.])
    sx, sy = centS[0], centS[1]
    satDict['SUNPIX'] = [sx, sy]
    
    # Get size of 1 Rs in pix
    if 'rsun' in imMap.meta:
        myRs = imMap.meta['rsun'] # in arcsec
    else:
        myDist = imMap.observer_coordinate.radius.m / 7e8
        myRs   = np.arctan2(1, myDist) * 206265
    oners = myRs/imMap.scale[0].to_value()
    satDict['ONERSUN'] = oners
    
    # Actually make the mask
    mask = np.zeros(imMap.data.shape)
    # Check that not SolO/PSP or STEREO HI
    if ('HI' not in imMap.meta['detector']) & (imMap.meta['obsrvtry'] not in ['Parker Solar Probe', 'Solar Orbiter']):   
        myOccR  = occultDict[myTag][0] # radius of the occulter in Rs
        occRpix = int(myOccR * oners)
        satDict['OCCRPIX'] = myOccR * oners
        satDict['OCCRARC'] = myOccR * oners * imMap.scale[0].to_value()
        
         # Fill in a circle around the occulter center
        for i in range(occRpix):
            j = int(np.sqrt(occRpix**2 - i**2))
            lowY = np.max([0,cy-j])
            hiY  = np.min([imMap.meta['naxis2']-2, cy+j])
            #print (cx+i, lowY,hiY+1)
            if cx+i <= imMap.meta['naxis2']-1:
                mask[cx+i, lowY:hiY+1] = 1
            if cx-i >=0:
                mask[cx-i, lowY:hiY+1] = 1    
    
        # Fill in outside FoV
        outRpix = int(occultDict[myTag][1] * oners) 
        for i in range(imMap.meta['naxis1']):
            myHdist = np.abs(cx-i)
            if myHdist >= outRpix:
                mask[i,:] = 1
            else:
                possY = int(np.sqrt(outRpix**2 - myHdist**2))
                lowY = np.max([0,cy - possY])
                hiY  = np.min([imMap.meta['naxis2'],cy + possY])
                mask[i,:lowY+1] = 1
                mask[i,hiY:] = 1
        satDict['MASK'] = mask
        
        # Get the sun outline
        thetas = np.linspace(0, 2.1*3.14159,100)
        xs = oners * np.cos(thetas) + sx
        ys = oners * np.sin(thetas) + sy
        satDict['SUNCIRC'] = [xs, ys]
        
        # Get the north line
        skyPt = SkyCoord(x=0, y=0, z=1, unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
        myPt2 = imMap.world_to_pixel(skyPt)
        satDict['SUNNORTH'] = [[sx, myPt2[0].to_value()], [sy, myPt2[1].to_value()]]

    return satDict
    
    
def pts2proj(pts_in, obs, scale, mywcs, center=[0,0], occultR=None):
    #  
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
    if mywcs['cunit'][0].lower() == 'arcsec':
        rad2unit = rad2arcsec
    elif mywcs['cunit'][0].lower() == 'deg':
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

def releaseTheWombat(obsFiles, nWFs=1):
    
    global mainwindow, pws, nSats, wfs, nwfs, bmodes
    
    # obsFiles should have 1 array for each satellite
    nSats = len(obsFiles)
    # each sat array then [[ims], [hdrs]]
    nwfs = nWFs
    wfs = [wf.wireframe(None) for i in range(nWFs)]
    
    # Find the min/max for each type of plot range
    sclIms = []    
    satStuff = []
    for i in range(nSats):
        satScls = []
        someStuff = []
        for j in range(len(obsFiles[i][0])):
            sclIm = makeNiceMMs(obsFiles[i][0][j], obsFiles[i][1][j])
            mySatStuff = getSatStuff(obsFiles[i][0][j])
            # Check if it made a mask and just use it now if so
            if 'MASK' in mySatStuff:
                midx = np.where(mySatStuff['MASK'] == 1)
                for k in range(len(sclIm)):
                    # black out all occulted
                    sclIm[k][midx] = -100. # might need to change if adjust plot ranges
                    # assigning sun pixels instead of just drawing on it makes disjointed so do later
                    #sx, sy = mySatStuff['SUNCIRC'][0].astype(int), mySatStuff['SUNCIRC'][1].astype(int)
                    #sclIm[k][sy,sx] = 200
            satScls.append(sclIm)
            someStuff.append(mySatStuff)
            
        sclIms.append(satScls)
        satStuff.append(someStuff)
        
    
    # Start the application
    app = QApplication(sys.argv)
    
    # Launch obs windows
    pws = []
    for i in range(nSats):
        pw = FigWindow('Sat1', obsFiles[i], sclIms[i], satStuff[i], myNum=i)
        pw.show()
        pws.append(pw) 
    
    
    # Launch the parameter panel    
    mainwindow = ParamWindow(nWFs)
    mainwindow.show()

    sys.exit(app.exec_())
    


if __name__ == "__main__":
    app = QApplication(sys.argv)
    global mainwindow, pws, nSats, wfs, nwfs, bmodes
    
    '''screen = app.primaryScreen()
    size = screen.size()
    width, height = size.width(), size.height()
    print (width, height)'''
    
    nwfs = 3

    wfs = [wf.wireframe(None) for i in range(nwfs)]
    
    pws = []
    pw = FigWindow('Sat1')
    pw.show()
    pws.append(pw)
    pw = FigWindow('Sat2', myNum=1)
    pw.show()
    pws.append(pw)
    nSats = len(pws)
    

    mainwindow = ParamWindow(3)
    mainwindow.show()


    sys.exit(app.exec_())