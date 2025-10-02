import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QGridLayout, QTabWidget, QSlider, QComboBox, QLineEdit
from PyQt5 import QtCore
import wombatWF as wf

import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('QGridLayout')
slogger.setLevel(logging.ERROR)

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
        self.WFs = np.array([None for i in range(nTabs)])
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
        layout.addWidget(label, 0,0,40,11)
        
        # Time slider/label
        Tlabel = QLabel('Time selection ')
        layout.addWidget(Tlabel,0,1,1,9,alignment=QtCore.Qt.AlignCenter)
        
        slider1 = QSlider()
        #slider1.setGeometry(QtCore.QRect(0, 100, 160, 25))
        slider1.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(slider1, 1,1,1,1)

        label = QLabel('Wireframe Type')
        layout.addWidget(label, 2,1,1,1,alignment=QtCore.Qt.AlignCenter)
        
        # Drop down box
        cbox = self.wfComboBox(i)
        layout.addWidget(cbox,3,0,1,5,alignment=QtCore.Qt.AlignCenter)
        
        label = QLabel('Wireframe Parameters')
        layout.addWidget(label, 5,1,alignment=QtCore.Qt.AlignCenter)
        
        # Dummy empty label to space it nicely before making wireframe specific
        #label = QLabel('')
        #layout.addWidget(label, 6,1,20,1)
        
        # Put a layout within the layout for the slider friends
        # It's like inception but without Elliot Page explaining everything
        WFLay = QGridLayout()
        layout.addLayout(WFLay, 7,1,25,11)
        
        '''slider1 = QSlider()
        slider1.setGeometry(QtCore.QRect(30, 130, 160, 22))
        slider1.setOrientation(QtCore.Qt.Horizontal)
        slider1.setMinimum(-90*i)
        slider1.setMaximum(90*i)
        layout.addWidget(slider1)'''
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
        widges[1][0].valueChanged.connect(lambda x: self.s2b(x, widges[0][0], i2f[0], myWF.ranges[0][0], myWF))  
        widges[0][0].returnPressed.connect(lambda: self.b2s(widges[1][0], widges[0][0], i2f[0], myWF.ranges[0][0],nSliders, myWF))     
        widges[1][1].valueChanged.connect(lambda x: self.s2b(x, widges[0][1], i2f[1], myWF.ranges[1][0], myWF))  
        widges[0][1].returnPressed.connect(lambda: self.b2s(widges[1][1], widges[0][1], i2f[1], myWF.ranges[1][0],nSliders, myWF))
        widges[1][2].valueChanged.connect(lambda x: self.s2b(x, widges[0][2], i2f[2], myWF.ranges[2][0], myWF))  
        widges[0][2].returnPressed.connect(lambda: self.b2s(widges[1][2], widges[0][2], i2f[2], myWF.ranges[2][0],nSliders, myWF))
        widges[1][3].valueChanged.connect(lambda x: self.s2b(x, widges[0][3], i2f[3], myWF.ranges[3][0], myWF))  
        widges[0][3].returnPressed.connect(lambda: self.b2s(widges[1][3], widges[0][3], i2f[3], myWF.ranges[3][0],nSliders, myWF))
        # Have to check remaining
        myNP = len(myWF.labels)
        # At least 5 params
        if myNP > 4:
            widges[1][4].valueChanged.connect(lambda x: self.s2b(x, widges[0][4], i2f[4], myWF.ranges[4][0], myWF))  
            widges[0][4].returnPressed.connect(lambda: self.b2s(widges[1][4], widges[0][4], i2f[4], myWF.ranges[4][0],nSliders, myWF))
        # At least 6 params    
        if myNP > 5:
            widges[1][5].valueChanged.connect(lambda x: self.s2b(x, widges[0][5], i2f[5], myWF.ranges[5][0], myWF))  
            widges[0][5].returnPressed.connect(lambda: self.b2s(widges[1][5], widges[0][5], i2f[5], myWF.ranges[5][0],nSliders, myWF))
        # At least 7 params    
        if myNP > 6:
            widges[1][6].valueChanged.connect(lambda x: self.s2b(x, widges[0][6], i2f[6], myWF.ranges[6][0], myWF))  
            widges[0][6].returnPressed.connect(lambda: self.b2s(widges[1][6], widges[0][6], i2f[6], myWF.ranges[6][0],nSliders, myWF))
        # At least 8 params           
        if myNP > 7:
            widges[1][7].valueChanged.connect(lambda x: self.s2b(x, widges[0][7], i2f[7], myWF.ranges[7][0], myWF))  
            widges[0][7].returnPressed.connect(lambda: self.b2s(widges[1][7], widges[0][7], i2f[7], myWF.ranges[7][0],nSliders, myWF))
        # At least 9 params    
        if myNP > 8:
            widges[1][8].valueChanged.connect(lambda x: self.s2b(x, widges[0][8], i2f[8], myWF.ranges[8][0], myWF))  
            widges[0][8].returnPressed.connect(lambda: self.b2s(widges[1][8], widges[0][8], i2f[8], myWF.ranges[8][0],nSliders, myWF))
            
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
        if event.key() in [QtCore.Qt.Key_Escape, QtCore.Qt.Key_Q]:
            self.close()
            
    def s2b(self, x=None, b=None, dx=None, x0=None, myWF=None):
        myVal = x0+dx*x
        if dx < 0.005:
            myStr = '{:.3f}'.format(myVal) 
        else:
            myStr = '{:.2f}'.format(myVal) 
        b.setText(myStr)

    #def b2s(self, x=None, s=None):
    #    s.setValue(int(x))
    def b2s(self,s,b, dx=None, x0=None, nSli=None, myWF=None):
        temp = b.text()
        slidx = int((float(b.text()) - x0)/dx)
        if slidx > nSli -1:
            slidx = nSli -1
            temp = x0 + nSli * dx
        elif slidx < 0:
            slidx = 0
            temp = x0
        s.setValue(slidx)
        # Reset it to what we actual wanted instead of slider rounded val
        # since this will trigger s2b as it hits valueChanged
        b.setText(temp)
        
    def cb_index_changed(self, a='None',idx=-10):
        self.WFtypes[idx] = a
        if type(self.WFs[idx]) == type(None):
            myType = self.WFnum2type[a]
            self.WFs[idx] = wf.wireframe(myType)
            self.tab_widget.setTabText(idx,self.WFshort[myType])
            
            WFLay, widges = self.WFparamLayout(self.WFs[idx])
                        
            self.layouts[idx].addLayout(WFLay, 7,1,30,11)
            self.WFLays[idx] = WFLay

            
        elif a == 0:
            # Set back to none if didn't select a WF
            self.WFs[idx] = None
            self.tab_widget.setTabText(idx,'None')
            thisLay = self.cleanLayout(self.WFLays[idx])
        else:
            # Create a new wf object but pass it any matching
            # parameters from the previous version
            ogLabs = self.WFs[idx].labels
            ogParams = self.WFs[idx].params
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
            self.WFs[idx] = newWF
            
            
    #|------------------------------| 
    #|----------- Others -----------|
    #|------------------------------| 
    def cleanLayout(self,lay):
        for i in reversed(range(lay.count())): 
            item = lay.takeAt(i)
            widget = item.widget()
            widget.deleteLater()
        return lay
    
        
class FigWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Secondary Window")
        self.setGeometry(550, 100, 300, 150)

        layout =  QVBoxLayout()
        label = QLabel("This is a secondary window.")
        layout.addWidget(label)
        self.setLayout(layout)
    

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    '''screen = app.primaryScreen()
    size = screen.size()
    width, height = size.width(), size.height()
    print (width, height)'''

    main_window = ParamWindow(3)
    main_window.show()

    #another_window = FigWindow()
    #another_window.show()

    sys.exit(app.exec_())