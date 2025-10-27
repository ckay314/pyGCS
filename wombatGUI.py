import sys, os
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QGridLayout, QTabWidget, QSlider, QComboBox, QLineEdit, QPushButton
from PyQt5 import QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import wombatWF as wf
import pyqtgraph as pg
from sunpy.visualization.colormaps import color_tables


from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u


from wombatLoadCTs import *

from wcs_funs import fitshead2wcs, wcs_get_pixel, wcs_get_coord
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


global occultDict, WFname2id
# Nominal radii (in Rs) for the occulters for each instrument. Pulled from google so 
# generally correct (hopefully) but not the most precise
occultDict = {'STEREO_SECCHI_COR2':[3,14], 'STEREO_SECCHI_COR1':[1.1,4], 'SOHO_LASCO_C1':[1.1,3], 'SOHO_LASCO_C2':[2,6], 'SOHO_LASCO_C3':[3.7,32], 'STEREO_SECCHI_HI1':[15,80], 'STEREO_SECCHI_HI2':[80,215], 'STEREO_SECCHI_EUVI':[0,1.7],'SDO_AIA':[0,1.3]} 
WFname2id = {'GCS':1, 'Torus':2, 'Sphere':3, 'Half Sphere':4, 'Ellipse':5, 'Half Ellipse':6, 'Slab':7}

        
class ParamWindow(QMainWindow):
    def __init__(self, nTabs):
        super().__init__()
        if nTabs > 10:
            sys.exit('Do you really need to fit >10 wireframes at once? If so figure out where the upper limit of 10 is hardcoded in the ParamWindow class in wombatGUI.py')
            
        self.nTabs = nTabs
        self.setWindowTitle('Wombat Parameters')
        self.setGeometry(100, 100, 300, 950)
        self.setFixedSize(300, 950) 

        # Create a QTabWidget
        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)
        
        # Create individual tabs (pages)
        self.tabs =[]
        
        # Create holder for the WF types
        self.WFtypes = np.zeros(nTabs)
        self.WFnum2type = ['None', 'GCS', 'Torus', 'Sphere', 'Half Sphere', 'Ellipse', 'Half Ellipse', 'Slab']
        self.WFshort = {'GCS':'GCS', 'Torus':'Tor', 'Sphere':'Sph', 'Half Sphere':'HSph', 'Ellipse':'Ell', 'Half Ellipse':'HEll', 'Slab':'Slab'}
        
        # Create holder for the WF and the params
        #self.WFs = theWFs #np.array([None for i in range(nTabs)])
        self.WFparams = np.array([np.zeros(10) for i in range(nTabs)])
        
        # holders for the param widgets so we can rm them
        self.WFLays = []
        self.widges = [None for i in range(nTabs)]
        self.layouts = []
        self.cbs = []
        
        # Number of points in the sliders
        self.nSliders = 201

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
        layout.addWidget(Tlabel,0,0,1,10,alignment=QtCore.Qt.AlignCenter)
        
        slider1 = QSlider()
        #slider1.setGeometry(QtCore.QRect(0, 100, 160, 25))
        slider1.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(slider1, 1,0,1,11)

        label = QLabel('Wireframe Type')
        layout.addWidget(label, 2,0,1,11,alignment=QtCore.Qt.AlignCenter)
        
        # Drop down box
        cbox = self.wfComboBox(i)
        self.cbs.append(cbox)
        layout.addWidget(cbox,3,0,1,11,alignment=QtCore.Qt.AlignCenter)
        
        label = QLabel('Parameters')
        layout.addWidget(label, 5,0, 1, 5,alignment=QtCore.Qt.AlignCenter)
        
        hideBut = QPushButton('Show/Hide WF')
        hideBut.released.connect(lambda: self.HBclicked(i))
        layout.addWidget(hideBut, 5, 6, 1,5)
        
        
        # Put a layout within the layout for the slider friends
        # It's like inception but without Elliot Page explaining everything
        WFLay = QGridLayout()
        layout.addLayout(WFLay, 7,0,25,11)
        
        
        # Background mode drop down box
        label = QLabel('Background Scaling')
        layout.addWidget(label, 41,0,1,6,alignment=QtCore.Qt.AlignCenter)
        cbox = self.bgComboBox()
        layout.addWidget(cbox,41,5,1,6,alignment=QtCore.Qt.AlignCenter)
        
        label = QLabel('')
        layout.addWidget(label, 42,0,2,11,alignment=QtCore.Qt.AlignCenter)
        

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
        nSliders = self.nSliders
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
            wfs[idx] = wf.wireframe(myType, WFidx=idx+1)
            self.tab_widget.setTabText(idx,self.WFshort[myType])
            
            WFLay, widges = self.WFparamLayout(wfs[idx])
                        
            self.layouts[idx].addLayout(WFLay, 7,0,30,11)
            self.WFLays[idx] = WFLay
            self.widges[idx] = widges
            
        elif a == 0:
            # Set back to none if didn't select a WF
            wfs[idx] = wf.wireframe(None)
            self.tab_widget.setTabText(idx,'None')
            thisLay = self.cleanLayout(self.WFLays[idx])
            for aPW in pws:
                aPW.scatters[idx].setData([])
            ovw.arrows[idx].setStyle(angle=0, headWidth=0, headLen=0, tailLen=0, tailWidth=0, pxMode=False,  pen={'color': color, 'width': 0}, brush=color)    
        else:
            # Create a new wf object but pass it any matching
            # parameters from the previous version
            ogLabs = wfs[idx].labels
            ogParams = wfs[idx].params
            myType = self.WFnum2type[a]
            newWF = wf.wireframe(self.WFnum2type[a], WFidx=idx+1)
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
            self.layouts[idx].addLayout(WFLay, 7,0,30,11)
            self.WFLays[idx] = WFLay
            self.widges[idx] = widges
           
            # Give the structure the new wf
            wfs[idx] = newWF
            
        for aPW in pws:
            aPW.plotBackground()
            aPW.plotWFs(justN=idx)
            
    def back_changed(self,text):
        for aPW in pws:
             aPW.cbox.setCurrentIndex(text)         
            
    def EBclicked(self):
        sys.exit()

    def HBclicked(self, i):
        if wfs[i].showMe: 
            wfs[i].showMe = False
            for aPW in pws:
                aPW.scatters[i].setData([])
        else:
            wfs[i].showMe = True
            for aPW in pws:
                aPW.plotWFs(justN=i)

    def SBclicked(self, singleSat=None):
        # Single sat is integer corresponding to a plot window/satellite number
        fileName = 'wombatSummaryFile.txt'
        # Get the filename 
        
        outFile = open('wboutputs/'+fileName, 'w')
        print ('Saving results in wboutputs/'+fileName)
        # Save the wireframe points
        for j in range(nwfs):
            aWF = wfs[j]
            if aWF.WFtype:
                outFile.write('WFtype'+str(j+1)+': ' + str(aWF.WFtype)+'\n')
                for i in range(len(aWF.labels)):
                    thisLab = aWF.labels[i]
                    spIdx = thisLab.find(' ')
                    if spIdx > 0:
                        outStr = thisLab[:spIdx]+str(j+1)+': ' + str(aWF.params[i])
                    else:
                        outStr = thisLab+str(j+1)+': ' + str(aWF.params[i])
                    outFile.write(outStr+'\n')
                    
        # Save the background plot info
        if type(singleSat) != type(None):
            toDo = [singleSat]
        else:
            toDo = range(nSats)
        for j in toDo:
            aPW = pws[j]
            tidx = aPW.tidx
            outStr = 'ObsTime'+str(j+1)+': ' + aPW.satStuff[tidx]['DATEOBS']
            outFile.write(outStr+'\n')
            if 'DATEOBS0' in aPW.satStuff[tidx]:
                outStr = 'ObsTimeZero'+str(j+1)+': ' + aPW.satStuff[tidx]['DATEOBS0']
                outFile.write(outStr+'\n')
            
            outStr = 'ObsType'+str(j+1)+': ' + aPW.satStuff[tidx]['MYTAG']
            outFile.write(outStr+'\n')
            
            outStr = 'Scaling'+str(j+1)+': ' +str(aPW.sclidx)
            outFile.write(outStr+'\n')
            outStr = 'MinVal'+str(j+1)+': ' +str(aPW.MinSlider.value())
            outFile.write(outStr+'\n')
            outStr = 'MaxVal'+str(j+1)+': ' +str(aPW.MaxSlider.value())
            outFile.write(outStr+'\n')
        
        outFile.close()

        # Save figures            
        for j in toDo:
            aPW = pws[j]
            tidx = aPW.tidx
            figName = 'wombat_'+ aPW.satStuff[tidx]['DATEOBS'] + '_' +  aPW.satStuff[tidx]['MYTAG'] +'.png'
            figGrab = aPW.pWindow.grab()
            figGrab.save('wboutputs/'+figName)
            print ('Saving figure in wboutputs/'+figName )
        if ovw:
            figName = 'wombat_'+ pws[0].satStuff[0]['DATEOBS'] + '_overview.png'
            figGrab = ovw.pWindow.grab()
            figGrab.save('wboutputs/'+figName)
            print ('Saving figure in wboutputs/'+figName )
            
        
        # Save Files
        for j in toDo:
            aPW = pws[j]
            tidx = aPW.tidx
            fitsName = 'wombat_'+ aPW.satStuff[tidx]['DATEOBS'] + '_' +  aPW.satStuff[tidx]['MYTAG'] +'.fits'
            if not os.path.exists('wbfits/'+fitsName):
                fitsdata = aPW.OGims[tidx].data               
                fitshdr  = aPW.hdrs[tidx]
                
                # Don't know why this gets tripped for some LASCO cases but making
                # it a string makes it happy
                if 'OBT_TIME' in fitshdr:
                    fitshdr['OBT_TIME'] = str(fitshdr['OBT_TIME'])
                print ('Saving fits file as wbfits/'+fitsName)
                hdu = fits.PrimaryHDU(fitsdata, header=fitshdr)
                hdu.writeto('wbfits/'+fitsName, overwrite=True)
        

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
                pws[ipw].plotWFs(justN=aWF.WFidx-1)
            if ovw:
                ovw.updateArrow(aWF.WFidx-1,color=aWF.WFcolor)
    
        
class FigWindow(QWidget):
    def __init__(self, satName, myObs, myScls, satStuff, myNum=0, labelPW=True):
        super().__init__()
        # Obs are [[ims], [hdrs]] as passed to the GUI (prob MSB)
        # Scls are [[lin, log, sqrt]] on -100,100 range
        
        self.setGeometry(550*(myNum+1), 350, 350, 450)
        self.winidx = myNum
        self.labelIt = labelPW
        
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
        self.pWindow.scene().sigMouseClicked.connect(self.mouse_clicked)
        layoutP.addWidget(self.pWindow,0,0,11,11,alignment=QtCore.Qt.AlignCenter)
        
        # make an image item
        self.image = pg.ImageItem()

        hasCT = check4CT(satStuff[0])
        if type(hasCT) != type(None):
            self.image.setLookupTable(hasCT)
        
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
        label = QLabel('Background Scaling')
        layoutP.addWidget(label, 12,0,1,5,alignment=QtCore.Qt.AlignCenter)
        self.cbox = self.bgComboBox()
        layoutP.addWidget(self.cbox,12,5,1,5,alignment=QtCore.Qt.AlignCenter)
        
        # Min/max brightness sliders
        minL = QLabel('Min Value:     ')
        layoutP.addWidget(minL, 13,0,1,9)
        self.MinSlider = QSlider()
        self.MinSlider.setOrientation(QtCore.Qt.Horizontal)
        self.MinSlider.setMinimum(0)
        self.MinSlider.setMaximum(255)
        self.MinSlider.setValue(satStuff[0]['SLIVALS'][0][0])
        layoutP.addWidget(self.MinSlider, 13,3,1,9)
        self.MinSlider.valueChanged.connect(lambda x: self.s2l(x, minL, 'Min Value: '))  
        
        maxL = QLabel('Max Value:     ')
        layoutP.addWidget(maxL, 15,0,1,9)
        self.MaxSlider = QSlider()
        self.MaxSlider.setOrientation(QtCore.Qt.Horizontal)
        self.MaxSlider.setMinimum(0)
        self.MaxSlider.setMaximum(255)
        self.MaxSlider.setValue(satStuff[0]['SLIVALS'][1][0])
        layoutP.addWidget(self.MaxSlider, 15,3,1,9)
        self.MaxSlider.valueChanged.connect(lambda x: self.s2l(x, maxL, 'Max Value: '))  
        
        # If EUV switch to log at the start. Have to do after
        # weve addded the sliders since will adjust them
        if self.satStuff[0]['OBSTYPE'] == 'EUV':
            self.cbox.setCurrentIndex(1)
        
        
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
        self.MinSlider.setValue(self.satStuff[0]['SLIVALS'][0][int(text)])  
        self.MaxSlider.setValue(self.satStuff[0]['SLIVALS'][1][int(text)])  
        self.plotBackground()
            
       
    def EBclicked(self):
        sys.exit()

    def SBclicked(self):
        mainwindow.SBclicked(singleSat=self.winidx)
        

    def MBclicked(self):
        print('Mass not coded yet')
    
    def mouse_clicked(self,event):
        scene_pos = event.scenePos()
        
        view_pos = self.pWindow.plotItem.vb.mapSceneToView(scene_pos)
        pix = [view_pos.x(), view_pos.y()]
        prefA = self.satStuff[self.tidx]['MYTAG'] + ' pix:'
        print (prefA, str(int(pix[0])).rjust(8), str(int(pix[1])).rjust(8))
        
        # Convert to ra/dec
        skyres = self.OGims[self.tidx].pixel_to_world(pix[0]*u.pixel, pix[1]*u.pixel)
        Tx, Ty = skyres.Tx.to_value(), skyres.Ty.to_value()
        print ('Tx, Ty (arcsec):'.rjust(len(prefA)), str(int(Tx)).rjust(8), str(int(Ty)).rjust(8))
        
        # Convert to proj Rsun/PA
        Rarc = np.sqrt(Tx**2 + Ty**2)
        Rpix = Rarc / self.satStuff[self.tidx]['SCALE']
        if self.satStuff[self.tidx]['OBSTYPE'] == 'HI':
            Rpix = Rpix / 3600
        RRSun = Rpix /  self.satStuff[self.tidx]['ONERSUN']
        # PA define w/ N as 0 and E (left) as 90
        PA = (np.arctan2(-Tx,Ty) * 180 / np.pi) % 360.
        print ('Proj R (Rs), PA (deg):'.rjust(len(prefA)), '{:8.2f}'.format(RRSun), '{:8.1f}'.format(PA))
        
    #|------------------------------| 
    #|----------- Others -----------|
    #|------------------------------| 
    def plotWFs(self, justN=0):
        # This is the slow version but keeping syntax around if need to check anything
        #skyPt = SkyCoord(x=0, y=0, z=1, unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
        #myPt2 = self.OGims[self.tidx].world_to_pixel(skyPt)
        
        toDo = range(nwfs)
        if justN:
            toDo = [justN]
            
        for i in toDo:
            if wfs[i].showMe & (type(wfs[i].WFtype) != type(None)):    
            #if type(wfs[i].WFtype) != type(None):                
                # Set up the input parameters for pts2proj
                pos = []
                obs = self.satStuff[self.tidx]['POS']
                obsScl = [self.satStuff[self.tidx]['SCALE'], self.satStuff[self.tidx]['SCALE']]
                if self.satStuff[self.tidx]['OBSTYPE'] == 'HI':
                    obsScl = [self.satStuff[self.tidx]['SCALE'] * 3600, self.satStuff[self.tidx]['SCALE'] * 3600]
                cent = self.satStuff[self.tidx]['SUNPIX']
                if 'OCCRARC' in self.satStuff[self.tidx]:
                    occultR = self.satStuff[self.tidx]['OCCRARC']
                else:
                    occultR = None
                mywcs  = self.satStuff[self.tidx]['WCS']
                myColor =wfs[i].WFcolor
                # change pen wid if HI
                penwid =1
                if self.satStuff[self.tidx]['OBSTYPE'] == 'HI':
                    penwid = 4
                
                
                # For the EUV panels, check if the WF is much higher
                # than the FOV and just project it onto the surface
                # instead if it is
                flatEUV = False
                if self.satStuff[self.tidx]['OBSTYPE'] == 'EUV':
                    pts = wfs[i].points
                    rs = np.sqrt(pts[:,0]**2 + pts[:,1]**2 + pts[:,2]**2)
                    if np.mean(rs) > 1.5*self.satStuff[self.tidx]['FOV']:
                        flatEUV = True
                # Downselect to fewer points for proj EUV
                toShow = range(len(wfs[i].points[:,0]))
                if flatEUV:
                    toShow = toShow[::2]
                    myColor = '#C81CDE'
                    occultR = 1. * self.satStuff[self.tidx]['ONERSUN']
                
                for jj in toShow:
                    pt = wfs[i].points[jj,:]
                    r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)
                    lat = np.arcsin(pt[2]/r) * 180/np.pi
                    lon = np.arctan2(pt[1],pt[0]) * 180 / np.pi
                    pt = [lat, lon, r*7e8]
                    if flatEUV:
                        pt = [lat, lon, 7e8]
                    myPt = pts2proj(pt, obs, obsScl, mywcs, center=cent, occultR=occultR)
                   
                    if len(myPt) > 0:          
                        pos.append({'pos': [myPt[0][0], myPt[0][1]], 'pen':{'color':myColor, 'width':penwid}, 'brush':pg.mkBrush(myColor)})
                        
                self.scatters[i].setData(pos)
                

 
    def plotBackground(self):
        #self.ui.graphWidget1.addItem(image1)
        myIm = self.myScls2[self.tidx][self.sclidx]
        slMin = self.MinSlider.value()
        slMax = self.MaxSlider.value()
        self.image.updateImage(image=myIm, levels=(slMin, slMax))
        if self.satStuff[self.tidx]['OBSTYPE'] != 'EUV':
            if 'SUNCIRC' in self.satStuff[self.tidx]:
                self.pWindow.plot(self.satStuff[self.tidx]['SUNCIRC'][0], self.satStuff[self.tidx]['SUNCIRC'][1])
            if 'SUNNORTH' in self.satStuff[self.tidx]:
                self.pWindow.plot(self.satStuff[self.tidx]['SUNNORTH'][0], self.satStuff[self.tidx]['SUNNORTH'][1], symbolSize=3, symbolBrush='w', pen=pg.mkPen(color='w', width=1))
        if self.labelIt:
            geom = self.pWindow.visibleRange()
            wid = geom.width()
            text_item1 = pg.TextItem(self.satStuff[self.tidx]['MYTAG'], anchor=(0, 1), fill='k')
            text_item1.setPos(0.001*wid, 0.001*wid)
            self.pWindow.addItem(text_item1)
            text_item2 = pg.TextItem(self.satStuff[self.tidx]['DATEOBS'], anchor=(1, 1), fill='k')
            text_item2.setPos(0.999*wid, 0.001*wid)
            self.pWindow.addItem(text_item2)
            



class OverviewWindow(QWidget):
    def __init__(self, satStuff):
        super().__init__()
        # Shows the position of the satellites and the slider longitude
        self.setFixedSize(400, 400) 
        self.setWindowTitle('Polar View')

        layoutOV =  QGridLayout()
        
        self.pWindow = pg.PlotWidget()
        self.pWindow.setMinimumSize(350, 350)
        layoutOV.addWidget(self.pWindow,0,0,11,11,alignment=QtCore.Qt.AlignCenter)
        
        # make an image item
        #self.image = pg.ImageItem()
        #self.pWindow.addItem(self.image)
        self.pWindow.setRange(xRange=(-1.2, 1.2), yRange=(-1.2,1.2), padding=0)
        
        #self.pWindow.plot([-1,1], [-1,1])
        twopi = np.linspace(0, 2.01*np.pi, 200)
        x_data = np.cos(twopi)
        y_data = np.sin(twopi)
        self.pWindow.plot(x_data, y_data, pen=pg.mkPen('w', width=3))
        self.pWindow.plot([0], [0], symbol='o', symbolSize=20, symbolBrush=pg.mkBrush(color='y'))
        self.pWindow.plot([0], [-1], symbol='o', symbolSize=15, symbolBrush=pg.mkBrush(color='blue'))
        
        # Hide the axes
        self.pWindow.hideAxis('bottom')
        self.pWindow.hideAxis('left')
        
        
        # Set up scatters for the sat points so can adjust without clearing
        self.scatters = []
        self.satLabs = []
        for i in range(nSats):
            # would be good to update this to adjust with time
            myPos = satStuff[i][0]['POS']
            myName = satStuff[i][0]['OBS']
            myR = myPos[2] / 1.496e+11 
            myLon = myPos[1] * np.pi / 180.
            y = - myR * np.cos(myLon)
            x = myR * np.sin(myLon)
            aScat = pg.ScatterPlotItem(pen=pg.mkPen(width=1, color='w'), brush=pg.mkBrush(color='w'),symbol='o', size=8)
            
            self.scatters.append(aScat)
            pos = []
            pos.append({'pos': [x,y]})
            self.scatters[i].setData(pos)
            self.pWindow.addItem(aScat)
            
            myName = satStuff[i][0]['SHORTNAME']
            text_item = pg.TextItem(myName, anchor=(0.5, 0.5))
            ysat = - 0.85*myR * np.cos(myLon)
            xsat = 0.85*myR * np.sin(myLon)
            text_item.setPos(xsat, ysat)
            
            self.satLabs.append(text_item)
            self.pWindow.addItem(text_item)
        
        # Set up arrow for the WF lons
        self.arrows = []
        for i in range(nwfs):
            arrow = pg.ArrowItem(angle=-45, tipAngle=0, headLen=0, tailLen=0, tailWidth=0, pen={'color': 'w', 'width': 2}, brush='b')
            arrow.setPos(0, 0)
            self.pWindow.addItem(arrow)
            self.arrows.append(arrow)
        self.setLayout(layoutOV)
    
    def updateArrow(self, i, color='w'):
        mywf = wfs[i]
        lon  = mywf.params[1]
        rlon = lon * np.pi /180.
        hL, tL = 0.1, 0.3
        aL = hL+tL
        xh = aL * np.sin(rlon)
        yh = -aL * np.cos(rlon)
        ang = -np.arctan2(yh, xh) * 180 / np.pi
        self.arrows[i].setStyle(angle=lon-270, headWidth=0.05, headLen=hL, tailLen=tL, tailWidth=0.03, pxMode=False,  pen={'color': color, 'width': 2}, brush=color)
        tail_len = self.arrows[i].opts['tailLen']
        self.arrows[i].setPos(xh, yh)
        
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Q: 
            self.close()
        elif event.key() == QtCore.Qt.Key_Escape:
            sys.exit()

def makeNiceMMs(imIn, hdr, satStuff):
    # Transpose and get perc values ignoring the NaNs
    im = np.transpose(imIn.data)
    imNonNaN = im[~np.isnan(im)]
    medval  = np.median(np.abs(imNonNaN))
    
    myInst = satStuff['INST']
    # mins/maxs on percentiles by instrument [[lower], [upper]] with [lin, log, sqrt]
    pMMs = {'AIA':[[0.001,10,1], [99,99,99]], 'SECCHI_EUVI':[[0.001,10,1], [99,99,99]], 'LASCO_C2':[[15,1,15], [97,99,97]], 'LASCO_C3':[[40,1,10], [99,99,90]], 'SECCHI_COR1':[[20,1,10], [90,99,90]], 'SECCHI_COR2':[[20,1,10], [90,99,90]], 'WISPR_HI1':[[1,40,1], [99.9,80,99.9]], 'WISPR_HI2':[[1,40,1], [99.9,80,99.9]], 'SoloHI':[[1,40,1], [99.5,80,99.5]] }
    
    sliVals = {'AIA':[[0,0,0], [191,191,191]], 'SECCHI_EUVI':[[0,32,0], [191,191,191]], 'LASCO_C2':[[0,0,21],[191,191,191]], 'LASCO_C3':[[37,0,37],[191,191,191]], 'SECCHI_COR1':[[0,0,21],[128,191,191]], 'SECCHI_COR2':[[0,0,21],[128,191,191]], 'WISPR_HI1':[[0,0,21],[128,191,191]], 'WISPR_HI2':[[0,0,21],[128,191,191]], 'SoloHI':[[0,0,21],[128,191,191]]}
    
    #pMMs = {'EUV':[[0.001,10,1], [99,99,99]], 'COR':[[20,1,10], [90,99,90]], 'HI':[[1,1,1], [99.9,99,99]]}
    #sliVals = {'EUV':[[0,0,0], [191,191,191]], 'COR':[[0,0,21],[128,191,191]], 'HI':[[1,1,1], [99.9,99,99]]}
    myMM = pMMs[myInst]
    mySliVals = sliVals[myInst]
    satStuff['SLIVALS'] = mySliVals
    
    # remove the nans
    im[np.isnan(im)] = 0
    
    # Linear Image
    linMin, linMax = np.percentile(imNonNaN, myMM[0][0]), np.percentile(imNonNaN, myMM[1][0])  
    rng = linMax- linMin
    linIm = (im - linMin) * 255 / rng
    
    # Log im
    tempIm = im / medval
    minVal = np.percentile(np.abs(tempIm),myMM[0][1])
    pidx = np.where(tempIm > minVal)
    nidx = np.where(tempIm < -minVal)
    logIm = np.zeros(tempIm.shape)
    logIm[np.where(np.abs(tempIm) < minVal)] = 1
    logIm[pidx] = np.log(tempIm[pidx] - minVal + 1)  
    logIm[nidx] = -np.log(-tempIm[nidx] - minVal + 1)  
    #minLog, maxLog = np.percent(logIm, myMM[1][0]), np.percent(logIm, myMM[1][1])
    perc95 = np.percentile(logIm, myMM[1][1])
    logIm = 191 * logIm / perc95
    
    # SQRT image
    tempIm = im / medval
    minVal = np.percentile(tempIm,myMM[0][2])
    tempIm = tempIm - minVal # set minVal at zero
    tempIm[np.where(tempIm < 0)] = 0
    sqrtIm = np.sqrt(tempIm)
    
    percX = np.percentile(sqrtIm, myMM[1][2])
    sqrtIm = 191 * sqrtIm / percX
      
    
    sclIms = [linIm, logIm,  sqrtIm]
    return sclIms, satStuff
    
def getSatStuff(imMap, diffDate=None):
    # Set up satDict as a micro header that we will package the mask in
    # Keys are OBS, INST, MYTAG, POS, SCALE, CRPIX, WCS, SUNPIX, ONERSUN, MASK
    # OCCRPIX, OCCRARC, SUNCIRC, SUNNORTH, OBSTYPE, FOV, DATEOBS, DATEOBS0
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
    elif myhdr['obsrvtry'] == 'SDO':
        satDict['OBS'] =  myhdr['obsrvtry'] 
        satDict['INST'] = myhdr['detector']
        myTag   = myhdr['obsrvtry'] + '_' + myhdr['detector']
    # Other less picky sats
    else:
        satDict['OBS'] =  myhdr['telescop']
        satDict['INST'] = myhdr['instrume'] + '_' + myhdr['detector']
        myTag   = myhdr['telescop'] + '_' + myhdr['instrume'] + '_' + myhdr['detector']
    satDict['MYTAG'] = myTag
        
    # Obs type - flag between HI, COR, EUV
    if satDict['OBS'] in ['Parker Solar Probe', 'Solar Orbiter']:
        satDict['OBSTYPE'] = 'HI'
    elif satDict['OBS'] in ['STEREO_A', 'STEREO_B']:
        if myhdr['detector'] in ['COR1', 'COR2']:
            satDict['OBSTYPE'] = 'COR'
        elif myhdr['instrume'] in ['HI1', 'HI2']:
            satDict['OBSTYPE'] = 'HI'
        else:
            satDict['OBSTYPE'] = 'EUV'
    elif satDict['OBS'] == 'SOHO':
        if myhdr['detector'] in ['C2', 'C3']:
            satDict['OBSTYPE'] = 'COR'
        else:
            satDict['OBSTYPE'] = 'EUV'
    elif satDict['OBS'] == 'SDO':
         satDict['OBSTYPE'] = 'EUV'

    # Add the wavelength if EUV
    if satDict['OBSTYPE'] == 'EUV':
        satDict['WAVE'] = str(myhdr['WAVELNTH'])
        satDict['MYTAG'] = satDict['MYTAG'] + '_' + satDict['WAVE']
            
    shortNames = {'Parker Solar Probe':'PSP', 'Solar Orbiter':'SolO', 'STEREO_A':'STA', 'STEREO_B':'STB', 'SOHO':'SOHO', 'SDO':'SDO'}
    satDict['SHORTNAME'] = shortNames[satDict['OBS']]
    
    # get a color table if we can
    if satDict['OBS'] == 'SDO':
        ct = color_tables.aia_color_table(myhdr['WAVELNTH']*u.angstrom)
        satDict['MYCT'] = ct
    
    # Get obs date/time
    if len(myhdr['date-obs']) > 13:
        satDict['DATEOBS'] = myhdr['date-obs']
    else:
        satDict['DATEOBS'] = myhdr['date-obs'] + 'T' + myhdr['time-obs']
    satDict['DATEOBS'] = satDict['DATEOBS'].replace('/','-')
    if '.' in satDict['DATEOBS']:
        dotidx = satDict['DATEOBS'].find('.')
        satDict['DATEOBS'] = satDict['DATEOBS'][:dotidx]

    if diffDate:
        satDict['DATEOBS0'] = diffDate
        satDict['DATEOBS0'] = satDict['DATEOBS0'].replace('/','-')
        if '.' in satDict['DATEOBS0']:
            dotidx = satDict['DATEOBS0'].find('.')
            satDict['DATEOBS0'] = satDict['DATEOBS0'][:dotidx]
     

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
    if imMap.scale[0].unit == 'deg / pix':
        oners = oners / 3600
    satDict['ONERSUN'] = oners
    
    # Get location of edges
    myFOV = 0
    for i in [0,imMap.data.shape[0]-1]:
        for j in [0,imMap.data.shape[1]-1]:
            coord = wcs_get_coord(myWCS, pixels = np.array([i,j]))
            edgeR = np.sqrt(coord[0]**2 + coord[1]**2)
            thisFOV = edgeR / obsScl / oners
            if thisFOV > myFOV: myFOV = thisFOV
    satDict['FOV'] = myFOV
    
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


def reloadIt(rD):
    # Set the wf params
    for i in range(nwfs):
        ii = str(i+1)
        if 'WFtype'+ii in rD:
            WFid = WFname2id[rD['WFtype'+ii]]
            wfs[i]  = wf.wireframe(rD['WFtype'+ii])
            mainwindow.cbs[i].setCurrentIndex(WFid)
            # get the expected labels
            for j in range(len(wfs[i].labels)):
                thisLab = wfs[i].labels[j]
                spIdx = thisLab.find(' ')
                shortStr = thisLab[:spIdx]+ii
                if shortStr in rD:
                    wfs[i].params[j] = float(rD[shortStr])
                    mainwindow.widges[i][0][j].setText(str(wfs[i].params[j]))
                myRng = wfs[i].ranges[j]
                dx = (myRng[1] - myRng[0]) / (mainwindow.nSliders - 1)
                x0 = myRng[0]
                slidx = int((float(wfs[i].params[j]) - x0)/dx)
                mainwindow.widges[i][1][j].setValue(slidx)
                #mainwindow.b2s(mainwindow.widges[i][1][0], mainwindow.widges[i][0][0], mainwindow.i2fs[i][0], wfs[i].ranges[0][0],mainwindow.nSliders, wfs[i], mainwindow.widges[i])
            mainwindow.updateWFpoints(wfs[i], mainwindow.widges[i])

    # Set the fig params
    for i in range(nSats):
        ii = str(i+1)
        myscl = int(rD['Scaling'+ii])
        myMin = int(rD['MinVal'+ii])
        myMax = int(rD['MaxVal'+ii])
        pws[i].cbox.setCurrentIndex(myscl)
        pws[i].MinSlider.setValue(myMin)
        pws[i].MaxSlider.setValue(myMax)
    
def releaseTheWombat(obsFiles, nWFs=1, overviewPlot=False, labelPW=True, reloadDict=None):
    
    global mainwindow, pws, nSats, wfs, nwfs, bmodes, ovw
    
    # obsFiles should have 1 array for each satellite
    nSats = len(obsFiles)
    # each sat array then [[ims], [hdrs]]
    if type(reloadDict) != type(None):
        nwfs = reloadDict['nWFs']
    else:
        nwfs = nWFs
    wfs = [wf.wireframe(None) for i in range(nwfs)]
    
    # Find the min/max for each type of plot range
    sclIms = []    
    satStuff = []
    for i in range(nSats):
        satScls = []
        someStuff = []
        for j in range(len(obsFiles[i][0])):
            mySatStuff = getSatStuff(obsFiles[i][0][j], diffDate=obsFiles[i][1][j]['DATE-OBS0'])
            sclIm, mySatStuff = makeNiceMMs(obsFiles[i][0][j], obsFiles[i][1][j], mySatStuff)
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
    
    # Get the max FoV from the satStuff so we can adjust height slider appropriately
    maxFoV = 0
    for stuff in satStuff:
        if stuff[0]['FOV'] > maxFoV: maxFoV = stuff[0]['FOV']
    # pad it a bit then round to a nice number
    maxFoV = int((1.1 * maxFoV) / 5) * 5
    # Edit the dictionaries in wombatWF
    wf.rngDict['Height (Rs)'] = [1,maxFoV]
    if maxFoV > 20:
        wf.defDict['Height (Rs)'] = 25
    elif maxFoV < 5:
        wf.defDict['Height (Rs)'] = 1.5
            
    
    # Start the application
    app = QApplication(sys.argv)
    
    # Launch obs windows
    pws = []
    for i in range(nSats):
        pw = FigWindow('Sat1', obsFiles[i], sclIms[i], satStuff[i], myNum=i, labelPW=labelPW)
        pw.show()
        pws.append(pw) 
    
    # Launch the overview panel (if turned on)
    if overviewPlot:
        ovw = OverviewWindow(satStuff)
        ovw.show()
    else:
        ovw = None
    
    
    # Launch the parameter panel    
    mainwindow = ParamWindow(nwfs)
    mainwindow.show()
    
    # Set to values from reload file if passed
    if type(reloadDict) != type(None):
        reloadIt(reloadDict)

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