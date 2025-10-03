import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QGridLayout, QTabWidget, QSlider, QComboBox, QLineEdit, QPushButton
from PyQt5 import QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import wombatWF as wf
import pyqtgraph as pg

import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('QGridLayout')
slogger.setLevel(logging.ERROR)

# Maintain a limited number of global variables to make passing things easier
# main window = the parameter window
# pws = array of plot windows
# wfs = array of wireframes in theoryland coords
# bmodes = array of integer background scaling modes, defaults to linear = 0
global mainwindow, pws, nwinds, wfs, nwfs, bmodes


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
            
        for ipw in range(nwinds):
            pws[ipw].plotWFs()
            
    def back_changed(self,text):
        print (text)
            
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
            for ipw in range(nwinds):
                pws[ipw].plotWFs()
    
        
class FigWindow(QWidget):
    def __init__(self, satName, myNum=0):
        super().__init__()
        self.setWindowTitle(satName)
        self.setGeometry(550*(myNum+1), 350, 350, 450)
        self.satName = satName

        layout =  QGridLayout()
        
        self.pWindow = pg.PlotWidget()
        self.pWindow.setMinimumSize(400, 400)
        layout.addWidget(self.pWindow,0,0,11,11,alignment=QtCore.Qt.AlignCenter)
        
        # attempt to get a nice size window
        self.pWindow.hideAxis('bottom')
        self.pWindow.hideAxis('left')
        self.pWindow.setRange(xRange=(-25,25), yRange=(-25,25), padding=0)
        
        # Background mode drop down box
        label = QLabel('Background Scaling Type')
        layout.addWidget(label, 12,0,1,5,alignment=QtCore.Qt.AlignCenter)
        cbox = self.bgComboBox()
        layout.addWidget(cbox,12,5,1,5,alignment=QtCore.Qt.AlignCenter)
        
        # Min/max brightness sliders
        minL = QLabel('Min Value:     ')
        layout.addWidget(minL, 13,0,1,9)
        slider = QSlider()
        slider.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(slider, 13,3,1,9)
        slider.valueChanged.connect(lambda x: self.s2l(x, minL, 'Min Value: '))  
        
        maxL = QLabel('Max Value:     ')
        layout.addWidget(maxL, 15,0,1,9)
        slider = QSlider()
        slider.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(slider, 15,3,1,9)
        slider.valueChanged.connect(lambda x: self.s2l(x, maxL, 'Max Value: '))  
        
        saveBut = QPushButton('Save')
        saveBut.released.connect(self.SBclicked)
        layout.addWidget(saveBut, 17, 0, 1,3,alignment=QtCore.Qt.AlignCenter)

        massBut = QPushButton('Mass')
        massBut.released.connect(self.MBclicked)
        layout.addWidget(massBut, 17, 4, 1,3,alignment=QtCore.Qt.AlignCenter)

        # Add things at the bottom
        exitBut = QPushButton('Exit')
        exitBut.released.connect(self.EBclicked)
        exitBut.setStyleSheet("background-color: red")
        layout.addWidget(exitBut, 17, 8, 1,3,alignment=QtCore.Qt.AlignCenter)
        
        self.setLayout(layout)
        #self.plotWFs(wfs)
        
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

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Q: 
            self.close()
        elif event.key() == QtCore.Qt.Key_Escape:
            sys.exit()
    

    def back_changed(self,text):
        print (text)
        
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
        #self.pWindow.plot(50*np.array(range(10)),50*np.array(range(10)), pen=None, symbol='o', symbolSize=350, symbolBrush='g', symbolPen=None)
        # No project yet, just plot yz
        self.pWindow.clear() # might have to change to avoid clearing background every time if slow
        for i in range(nwfs):
            if type(wfs[i].WFtype) != type(None):
                if self.satName == 'Sat1':
                    ys = wfs[i].points[:,1]
                    zs = wfs[i].points[:,2]
                elif self.satName == 'Sat2':
                    ys = wfs[i].points[:,0]
                    zs = wfs[i].points[:,2]
                myColor =wfs[i].WFcolor
                self.pWindow.plot(ys,zs, pen=None, symbol='o', symbolSize=3, symbolBrush=myColor, symbolPen=None)
    

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
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
    nwinds = len(pws)
    

    mainwindow = ParamWindow(3)
    mainwindow.show()


    sys.exit(app.exec_())