#from PyQt5 import QtWidgets, QtGui, QtCore
from astropy.io import fits
import astropy.units as u
from pyGCS import *
from GCSgui import *
from scipy import ndimage
import sunpy.map
import sys, os

sys.path.append('/Users/kaycd1/STEREO_Mass/IDLport') 
from secchi_prep import secchi_prep
from wispr_prep import wispr_prep
from lasco_prep import c2_prep, c3_prep
from solohi_prep import solohi_fits2grid
from aia_prep import aia_prep

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import matplotlib.image as mpimg

from wombatGUI import releaseTheWombat

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)


def processReload(fileIn):
    if os.path.exists(fileIn):
        print ('Loading configuration from '+fileIn)
        data = np.genfromtxt(fileIn,dtype=str)
        # Get the number of wfs and sats while making Dict
        nWFs, nSats = 0, 0
        reloadDict = {}
        for i in range(len(data[:,0])):
            mykey = data[i,0]
            if 'WFtype' in mykey: nWFs += 1
            if 'ObsType' in mykey: nSats += 1
            reloadDict[str(mykey).replace(':','')] = str(data[i,1])
        reloadDict['nWFs'] = nWFs
        reloadDict['nSats'] = nSats
        
        allFH = []
        for i in range(nSats):
            fitsname = 'wombat_'+reloadDict['ObsTime'+str(i+1)]+'_' + reloadDict['ObsType'+str(i+1)]+ '.fits'
            if os.path.exists('wbfits/'+fitsname):
                #with fits.open('wbfits/'+fitsname) as hdulist:
                #    im  = np.asarray(hdulist[0].data)
                #    hdr = hdulist[0].header
                hdul = fits.open('wbfits/'+fitsname)
                im = np.asarray(hdul[0].data)
                hdr = hdul[0].header
                hdul.close()
                diff = sunpy.map.Map(im, hdr)
                allFH.append([[diff], [hdr]])    # double pack for single time
            else:
                sys.exit('Cannot find wbfits/'+fitsname+' Make sure it is in the correct folder.')
        return allFH, reloadDict
        
    else:
        sys.exit('Reload file not found, cannot launch')


def pullProcFiles(theFile, diffMode='RD', diffEUV=False):
    # allFH is all the fits and headers packaged into an array
    # allFH = [[instr1], [instr2], ...]
    # where each instr is [[data_t1, data_t2...], [hdr_t1, hdr_t2]]
    # and the data and hdrs are assume to be sorted in time order
    # (but will not actually fail if they are not)
    
    # diffMode can be RD or BD for running diff or base diff
    
    # Set up empty holders
    allFH0 = []
    names = []
    counter = -1
    
    pFiles = np.genfromtxt(theFile, dtype=str)
    for apF in pFiles:
        if '.fits' not in apF:
            print ('Reading in ', apF, ' files')
            allFH0.append([[], []])
            names.append(apF)
            counter += 1
        else:
            print(' ', apF)
            with fits.open(apF) as hdulist:
                im  = hdulist[0].data
                hdr = hdulist[0].header
                
            allFH0[counter][0].append(im)
            allFH0[counter][1].append(hdr)
    
    # Convert the fits to maps and take differences
    allFH = []
    for i in range(len(allFH0)):
        # Make sure it has more than one file 
        if (len(allFH0[i][0]) > 1) or (('EUVI' in names[i]) & (not diffEUV)) or (('AIA' in names[i]) & (not diffEUV)):
            allFH.append([[], []])
            # EUV no diff cases
            if (('EUVI' in names[i]) & (not diffEUV)) or (('AIA' in names[i]) & (not diffEUV)):            
                for j in range(len(allFH0[i][0])):
                    myData = allFH0[i][0][j]
                    myHdr  = allFH0[i][1][j]
                    myMap  = sunpy.map.Map(myData, myHdr)
                    # Add date obs0 into header 
                    myHdr['date-obs0'] = None
                    
                    allFH[i][0].append(myMap)
                    allFH[i][1].append(myHdr)
            # Any form of difference
            else:
                for j in range(len(allFH0[i][0])-1):
                    # Get the data for a time step
                    myData = allFH0[i][0][j+1]
                    myHdr  = allFH0[i][1][j+1]
                    # Get the 
                    if diffMode == 'BD':
                        myBase = allFH0[i][0][0]
                        hdr0   = allFH0[i][1][0]
                    else:
                        myBase = allFH0[i][0][j]
                        hdr0   = allFH0[i][1][j]
                        
                    if (myData.shape == myBase.shape):
                        diffData = myData - myBase
                        diffMap = sunpy.map.Map(diffData, myHdr)
                        
                        if len(hdr0['date-obs']) > 13:
                            myHdr['date-obs0'] = hdr0['date-obs']
                        else:
                             myHdr['date-obs0'] = hdr0['date-obs'] + 'T' + hdr0['time-obs']
                        
                        allFH[i][0].append(diffMap)
                        allFH[i][1].append(myHdr)
                    else:
                        print ('Skipping file -- mismatch in size for ' + names[i] + allFH0[i][1][j+1]['DATE-OBS'])
            
        else:
            print ('Cannot make diff for '+ names[i]+ ' from only a single file')
            

    return allFH
    
# Check if we were passed a file
# Could be either a reload file or a list of pre processed obs
if len(sys.argv) > 1:
    theFile = sys.argv[1]
    
    doReload = False
    preProc  = False
    
    if 'reload' in theFile:
        doReload = True
    elif 'obslist' in theFile:
        preProc = True
    else:
        print ('Need to determine if reload/obs, tbd')
        
    if doReload:    
        allFH, reloadDict = processReload(theFile)
    else:
        reloadDict = None
        
    if preProc:
        allFH = pullProcFiles(theFile)
    

    releaseTheWombat(allFH,overviewPlot=False, reloadDict=reloadDict)

    print (sd)

else:
    # STEREO
    fnameB1 = 'fits/20120712_162400_d4c2A.fts'
    fnameB2 = 'fits/20120712_172400_d4c2A.fts' 

    #fnameB1 = 'fits/20120712_162400_d4c2B.fts'
    #fnameB2 = 'fits/20120712_172400_d4c2B.fts' 

    #fnameB1 = 'fits/20120713_004901_s4h1A.fts'
    #fnameB2 = 'fits/20120713_072901_s4h1A.fts'

    # LASCO
    #fnameB1 = '/Users/kaycd1/wombat/fits/C2_22800178.fts'
    #fnameB2 = '/Users/kaycd1/wombat/fits/C2_22800186.fts'

    #fnameA1 = '/Users/kaycd1/wombat/fits/C3_32048310.fts'
    #fnameA2 = '/Users/kaycd1/wombat/fits/C3_32048311.fts'

    # PSP     
    fnameW1 = 'fits/psp_L2_wispr_20210121T055621_V1_1211.fits'
    fnameW2 = 'fits/psp_L2_wispr_20210121T062821_V1_1211.fits'
    #fnameW1 = 'fits/testing/psp_L2_wispr_20250322T010204_V0_2222.fits'
    #fnameW2 = 'fits/testing/psp_L2_wispr_20250322T010702_V0_2222.fits'


    # SoloHI
    file1A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T053609_V02.fits'
    file2A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-2ft_20220329T053809_V02.fits'
    file3A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-3fg_20220329T054157_V02.fits'
    file4A = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T055357_V02.fits'

    file1B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-1ft_20220329T052409_V02.fits'
    file2B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-2ft_20220329T052609_V02.fits'
    file3B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-3fg_20220329T051757_V02.fits'
    file4B = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T052957_V02.fits'
    
    fnameC1 = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T052957_V02.fits'
    fnameC2 = '/Users/kaycd1/wombat/fits/solo_L2_solohi-4fg_20220329T055357_V02.fits'

    '''filenames12 = ['/Users/kaycd1/wombat/obsFiles/SECCHI/COR1_20230304_120600_n4c1a.fts', '/Users/kaycd1/wombat/obsFiles/SECCHI/COR1_20230304_120618_n4c1a.fts', '/Users/kaycd1/wombat/obsFiles/SECCHI/COR1_20230304_120636_n4c1a.fts']

    ims12, hdrs12 = secchi_prep(filenames12, polarizeOn=True)
    
    filenames13 = ['/Users/kaycd1/wombat/obsFiles/SECCHI/COR1_20230304_130600_n4c1a.fts', '/Users/kaycd1/wombat/obsFiles/SECCHI/COR1_20230304_130618_n4c1a.fts', '/Users/kaycd1/wombat/obsFiles/SECCHI/COR1_20230304_130636_n4c1a.fts']

    ims13, hdrs13 = secchi_prep(filenames13, polarizeOn=True)
    diff = ims13[0] - ims12[0]
    diff = sunpy.map.Map(ims3[0], hdrs3[0])
    hdrs = [hdrs12[0], hdrs13[0]]'''
    
    # STEREO EUVI
    '''aEUV = ['/Users/kaycd1/wombat/obsFiles/SECCHI/EUVI_304a_20230304_144545_n4eua.fts']
    ims, hdrs = secchi_prep(aEUV)
    #diff = ims[1] - ims[0]
    diff = sunpy.map.Map(ims[0], hdrs[0])
    # hack to make it run with single file
    hdrs = [hdrs[0], hdrs[0]]'''

    '''bPair = [fnameB1, fnameB2]
    ims2, hdrs2 = secchi_prep([bPair[0], bPair[1]])
    diff2 = ims2[1] - ims2[0]
    diff2 = sunpy.map.Map(diff2, hdrs2[1])'''
    
    # SDO AIA 
    fname = ['/Users/kaycd1/wombat/obsFiles/AIA/aia_lev1_193a_2023_03_04t14_20_04_84z_image_lev1.fits']
    ims = aia_prep(fname)
    ims[0].meta['OBSRVTRY'] = 'SDO'
    ims[0].meta['DETECTOR'] = 'AIA'
    ims[0].meta['SC_ROLL'] = 0. # is rotated in aia_prep
    hdrs = [ims[0].meta, ims[0].meta]
    diff = ims[0]
    
    # LASCO 
    '''aPair = [fnameA1, fnameA2]
    ims, hdrs = c3_prep([aPair[0], aPair[1]])
    if hdrs[1]['DETECTOR'] not in ['HI1', 'HI2']:
         diff = ims[1] - ims[0]
         hdrs[1]['OBSRVTRY'] = 'SOHO'
         diff = sunpy.map.Map(diff, hdrs[1])
    
    bPair = [fnameB1, fnameB2]
    ims2, hdrs2 = c2_prep([bPair[0], bPair[1]])
    if hdrs2[1]['DETECTOR'] not in ['HI1', 'HI2']:
         diff2 = ims2[1] - ims2[0]
         hdrs2[1]['OBSRVTRY'] = 'SOHO'
         diff2 = sunpy.map.Map(diff2, hdrs2[1])'''
    '''if True:
         norm = mpl.colors.Normalize(vmin=0, vmax=255)
         scalar_mappable = cm.ScalarMappable(norm=norm, cmap=cm.gist_heat)
         r = ''
         b = ''
         g = ''
         for i in range(256):
             #print (scalar_mappable.to_rgba(i)[0]*255, scalar_mappable.to_rgba(i)[1]*255,)
             r = r + str(int(255*scalar_mappable.to_rgba(i)[0])) + ' '
             g = g + str(int(255*scalar_mappable.to_rgba(i)[1])) + ' '
             b = b + str(int(255*scalar_mappable.to_rgba(i)[2])) + ' '
            # print (int(255*scalar_mappable.to_rgba(i)[0]), int(255*scalar_mappable.to_rgba(i)[1]), int(255*scalar_mappable.to_rgba(i)[2]))
         
         print (g)'''

    # PSP
    '''ims2, hdrs2 = wispr_prep([fnameW1, fnameW2], straylightOff=True)
    diff2 = ims2[1] - ims2[0]
    # detector as int not string so fix that
    for hdr in hdrs2:
        hdr['detector'] = str(hdr['detector'])
    diff2 = sunpy.map.Map(diff2, hdrs2[1])'''

    # SolO Quad
    '''aPair = [[file1B, file2B, file3B, file4B], [file1A, file2A, file3A, file4A]]
    im0, hdr0 = solohi_fits2grid(aPair[0])
    im1, hdr1 = solohi_fits2grid(aPair[1])
    diff2 = im1 - im0
    hdrs2 = [hdr0, hdr1] 
    diff2 = sunpy.map.Map(diff2, hdrs2[1])'''

    # SolO single (or anything fully prepped)
    aPair = [fnameC1, fnameC2]
    with fits.open(aPair[0]) as hdulist:
        im1  = hdulist[0].data
        hdr1 = hdulist[0].header
    with fits.open(aPair[1]) as hdulist:
        im2  = hdulist[0].data
        hdr2 = hdulist[0].header
    diff2 = im2 - im1
    hdrs2 = [hdr1, hdr2]
    diff2 = sunpy.map.Map(diff2, hdrs2[1])

    if len(hdrs[0]['date-obs']) > 13:
        hdrs[1]['date-obs0'] = hdrs[0]['date-obs']
    else:
         hdrs[1]['date-obs0'] = hdrs[0]['date-obs'] + 'T' + hdrs[0]['time-obs']
    if len(hdrs2[0]['date-obs']) > 13:
        hdrs2[1]['date-obs0'] = hdrs2[0]['date-obs']
    else:
         hdrs2[1]['date-obs0'] = hdrs2[0]['date-obs'] + 'T'  + hdrs2[0]['time-obs']

    releaseTheWombat([[[diff],[hdrs[1]]], [[diff2],[hdrs2[1]]]], nWFs=2, overviewPlot=False)