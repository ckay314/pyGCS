import sys, os
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.time import TimeDelta
from sunpy.time import parse_time


sys.path.append('/Users/kaycd1/STEREO_Mass/IDLport') 
from secchi_prep import secchi_prep
from wispr_prep import wispr_prep
from lasco_prep import c2_prep, c3_prep
from solohi_prep import solohi_fits2grid

# This routine takes already downloaded files and performs all the
# calibration so that they are ready for the wombat gui

# We assume that the files are stored separated into folders in the 
# same format as wombatPullObs uses. The indiv folders are AIA, LASCO
# SECCHI, SoloHI, and WISPR


# |-------------------------------|
# |------- Setup everything ------|
# |-------------------------------|

# Base folder for where the files are stored 
obsFiles = '/Users/kaycd1/wombat/obsFiles/'

# Start and end time of the period of interest
#startT = '2023/03/02T12:00'
#endT   = '2023/03/04T16:00'

startT = '2023/12/31T22:00'
endT   = '2024/01/01T02:00'

# Difference mode - option to save as base or running diff or just keep as is
diffMode = 'AsIs' # Select from 'AsIs', 'Base', 'Run'
noDiffEUV = True  # keep the EUV images as is instead of base/running diffs 

# AIA setup
doAIA = False
AIAwav = [171, 193, 304] # Select from 94, 131, 171*, 193*, 211, 304*, 335, 1600, 1700 (* most common)

# LASCO setup
doLASCO = False
whichLASCO = ['C2', 'C3'] # Select from 'C2' and 'C3'

# SECCHI setup
doSECCHI = False
whichSECCHI = ['EUVI', 'COR1', 'COR2',  'HI1', 'HI2'] # Select from 'EUVI', 'COR1', 'COR2',  'HI1', 'HI2'
EUVIwav     = [195, 171] # Select from 171, 195, 284, 304

# SoloHI setup
doSoloHI = True
whichSoloHI = ['Quad'] # options of 'Quad' or 'Single'

# WISPR setup
doWISPR = False
whichWISPR = ['In', 'Out'] # options of 'In' and 'Out


# |-------------------------------|
# |------- Setup Time Stuff ------|
# |-------------------------------|
# Convert strings to astropy time objects
startAPT = parse_time(startT)
endAPT = parse_time(endT)


# Get day of year
stDOY = int(startAPT.yday.split(':')[1])
# End easy if same year
if startAPT.datetime.year == endAPT.datetime.year:
    enDOY = int(endAPT.yday.split(':')[1])
# Split year case
else:
    nye = parse_time(str(startAPT.datetime.year) + '/12/31')
    nye_doy = int(nye.yday.split(':')[1]) # get last doy (good for leap years)
    enDOY = int(endAPT.yday.split(':')[1]) + nye_doy
    

# Check if same and doing single day 
if stDOY == enDOY:
    singleDay = True
    ymds = [str(startAPT)[:10].replace('-','')]
# Otherwise pull ymd strings for each date
else:
    singleDay = False
    nDays = enDOY - stDOY + 1
    ymds = []
    for i in range(nDays):
        nowDay = startAPT + TimeDelta(i * u.day)
        ymds.append(str(nowDay)[:10].replace('-',''))
        
# Get time on starting/ending days
hm0 = str(startAPT)[11:16].replace(':','')
hmf = str(endAPT)[11:16].replace(':','')


# |-------------------------------|
# |--------- Process AIA ---------|
# |-------------------------------|
if doAIA:
    # Get all the AIA files
    AIAfiles = os.listdir(obsFiles+'AIA/')
    
    # Format the date string as expected in AIA file names
    AIAdatestrs = [ymd[:4]+'_'+ymd[4:6]+'_'+ymd[6:] for ymd in ymds]
    
    # Check all the files to see if they have a matching date str
    # If date matches then check the hour/min on first/last date
    # Add into separate arrays for each wavelength
    nWav = len(AIAwav)
    goodFiles = [[] for i in range(nWav)]
    for aF in AIAfiles:
        if singleDay:
            if AIAdatestrs[0] in aF:
                hm = aF[25:30]
                if (hm >= hm0) & (hm <= hmf):
                    for i in range(nWav):
                        if '_'+str(AIAwav[i])+'a_' in aF:
                            goodFiles[i].append(aF)
        else:
            for j in range(nDays):
                if AIAdatestrs[j] in aF:
                    hm = aF[25:30]
                    addIt = True
                    if (j == 0) & (hm < hm0):
                        addIt = False
                    elif (j == nDays-1) & (hm > hmf):
                        addIt = False
                    if addIt:
                        for i in range(nWav):
                            if '_'+str(AIAwav[i])+'a_' in aF:
                                goodFiles[i].append(aF)
    # Make an array and make sure sorted               
    for i in range(nWav):
        goodFiles[i] = np.sort(np.array(goodFiles[i]))
        
    # Add in the actual processing, saving, and output to runFile  
 

# |-------------------------------|
# |------- Process LASCO ---------|
# |-------------------------------|
if doLASCO:
    # Get all the LASCO files
    LASCOfiles = os.listdir(obsFiles+'LASCO/')
    
    # date strs already formatted how we want bc we added
    # them onto the lasco names to start with 
    
    # Check all the files to see if they have a matching date str
    # If date matches then check the hour/min on first/last date
    # Sort all by C2/C3, will only process what we want later
    goodFiles =[[], []]
    for aF in LASCOfiles:
        if singleDay:
            if ymds[0] in aF: 
                hm = aF[9:13]
                if (hm >= hm0) & (hm <= hmf):
                    if '_C2_' in aF:
                        goodFiles[0].append(obsFiles+'LASCO/'+aF)
                    elif '_C3_' in aF:
                        goodFiles[1].append(obsFiles+'LASCO/'+aF)
        else:
            for j in range(nDays):
                if ymds[j] in aF: 
                    hm = aF[9:13]
                    addIt = True
                    if (j == 0) & (hm < hm0):
                        addIt = False
                    elif (j == nDays-1) & (hm > hmf):
                        addIt = False
                    if addIt:
                        if '_C2_' in aF:
                            goodFiles[0].append(obsFiles+'LASCO/'+aF)
                        elif '_C3_' in aF:
                            goodFiles[1].append(obsFiles+'LASCO/'+aF)
                            
    # Sort and array-ify
    goodFiles[0] = np.sort(np.array(goodFiles[0]))
    goodFiles[1] = np.sort(np.array(goodFiles[1]))
                        
    # Add in the actual processing, saving, and output to runFile  
    if 'C2' in whichLASCO:
        print ('|---- Processing LASCO C2 ----|')
        imsc2, hdrsc2 = c2_prep(goodFiles[0])    
    if 'C3' in whichLASCO:
        print ('|---- Processing LASCO C3 ----|')
        imsc3, hdrsc3 = c3_prep(goodFiles[1])    
                               
            
# |-------------------------------|
# |------- Process SECCHI --------|
# |-------------------------------|

if doSECCHI:
    # Get all the SECCHI files
    SECCHIfiles = os.listdir(obsFiles+'SECCHI/')
    
    # Sort everything first and just process what we want
    goodFiles = {'a':[[] for i in range(5)], 'b':[[] for i in range(5)]}
    
    # Dictionary to convert prefix to index for good f
    pre2ind = {'EUVI':0, 'COR1':1, 'COR2':2, 'HI1_':3, 'HI2_':4}
        
    for aF in SECCHIfiles:
        if singleDay:
            if ymds[0] in aF:
                # Check if EUVI which has two _ before date
                if aF[0] == 'E':
                    stidx = 9
                # Otherwise COR/HI have single _
                else:
                    stidx = aF.find('_')
                hm = aF[stidx+10:stidx+14]
                if (hm >= hm0) & (hm <= hmf):
                    whichInst = pre2ind[aF[:4]]
                    goodFiles[aF[-5]][whichInst].append(obsFiles+'SECCHI/'+aF)
        else:
            for j in range(nDays):
                if ymds[j] in aF:
                    # Check if EUVI which has two _ before date
                    if aF[0] == 'E':
                        stidx = 9
                    # Otherwise COR/HI have single _
                    else:
                        stidx = aF.find('_')
                    hm = aF[stidx+10:stidx+14]
                    addIt = True
                    if (j == 0) & (hm < hm0):
                        addIt = False
                    elif (j == nDays-1) & (hm > hmf):
                        addIt = False
                    if addIt:
                        whichInst = pre2ind[aF[:4]]
                        goodFiles[aF[-5]][whichInst].append(obsFiles+'SECCHI/'+aF)
    
    for i in range(len(goodFiles['a'])):
        goodFiles['a'][i] = np.sort(np.array(goodFiles['a'][i]))
    for i in range(len(goodFiles['b'])):
        goodFiles['b'][i] = np.sort(np.array(goodFiles['b'][i]))
                
    # Add in the actual processing, saving, and output to runFile  
    #if 'EUVI' in whichSECCHI:
    if 'COR1' in whichSECCHI:   
        if len(goodFiles['a'][1]) > 0:
            ims1a, hdrs1a = secchi_prep(goodFiles['a'][1]) 
        if len(goodFiles['b'][1]) > 0:
            ims1b, hdrs1b = secchi_prep(goodFiles['b'][1])
    if 'COR2' in whichSECCHI:   
        if len(goodFiles['a'][2]) > 0:
            print ('|--- Processing STEREO COR2A ---|')
            imsc2a, hdrsc2a = secchi_prep(goodFiles['a'][2]) 
        if len(goodFiles['b'][2]) > 0:
            print ('|--- Processing STEREO COR2B ---|')
            imsc2b, hdrsc2b = secchi_prep(goodFiles['b'][2])
    if 'HI1' in whichSECCHI:   
        if len(goodFiles['a'][3]) > 0:
            print ('|--- Processing STEREO HI1A ---|')
            imsh1a, hdrsh1a = secchi_prep(goodFiles['a'][3]) 
        if len(goodFiles['b'][3]) > 0:
            print ('|--- Processing STEREO HI1B ---|')
            imsh1b, hdrsh1b = secchi_prep(goodFiles['b'][3]) 
    if 'HI2' in whichSECCHI:   
        if len(goodFiles['a'][4]) > 0:
            print ('|--- Processing STEREO HI2A ---|')
            imsh2a, hdrsh2a = secchi_prep(goodFiles['a'][4]) 
        if len(goodFiles['b'][4]) > 0:
            print ('|--- Processing STEREO HI2B ---|')
            imsh2b, hdrsh2b = secchi_prep(goodFiles['b'][4]) 
        
# |-------------------------------|
# |------- Process SoloHI --------|
# |-------------------------------|
if doSoloHI:
    # Get all the SoloHI files
    SoloHIfiles = os.listdir(obsFiles+'SoloHI/')
    # Check all the files to see if they have a matching date str
    # If date matches then check the hour/min on first/last date
    goodFiles =[[], [], [], []]
    for aF in SoloHIfiles:
        if singleDay:
            if ymds[0] in aF:
                hm = aF[28:32]
                if (hm >= hm0) & (hm <= hmf):
                    if 'hi-1' in aF:
                        goodFiles[0].append(aF)
                    elif 'hi-2' in aF:
                        goodFiles[1].append(aF)
                    elif 'hi-3' in aF:
                        goodFiles[2].append(aF)
                    elif 'hi-4' in aF:
                        goodFiles[3].append(aF)
        else:
            for j in range(nDays):
                if ymds[j] in aF:
                    hm = aF[28:32]
                    addIt = True
                    if (j == 0) & (hm < hm0):
                        addIt = False
                    elif (j == nDays-1) & (hm > hmf):
                        addIt = False
                    if addIt:
                        if 'hi-1' in aF:
                            goodFiles[0].append(aF)
                        elif 'hi-2' in aF:
                            goodFiles[1].append(aF)
                        elif 'hi-3' in aF:
                            goodFiles[2].append(aF)
                        elif 'hi-4' in aF:
                            goodFiles[3].append(aF)
                        
    # Make an array and make sure sorted               
    for i in range(4):
        goodFiles[i] = np.sort(np.array(goodFiles[i]))
    
    # Need to match images to quad for this one with diff time cadences..
    if 'Quad' in whichSoloHI:
        shTimes = [[], [], [], []]
        prefixs = ['solo_L2_solohi-1ft_', 'solo_L2_solohi-2ft_', 'solo_L2_solohi-3fg_', 'solo_L2_solohi-4fg_']
        for i in range(4):
            for j in range(len(goodFiles[i])):
                thisT = parse_time(goodFiles[i][j].replace(prefixs[i],'').replace('_V01.fits',''))
                shTimes[i].append(thisT)
    # Figure out which has the most obs
    nEach =[len(shTimes[i]) for i in range(4)]
    nMost = np.where(nEach == np.max(nEach))[0][0]
    nQuads = nEach[nMost] # total number of quad images
    quadFiles = [[] for i in range(4)]
    # Find the closest time match for each panel
    # This will duplicate panels as needed if nEach < nQuads
    for i in range(4):
        if i == nMost:
            quadFiles[i] = goodFiles[i]
        else:
            quadFiles[i] = np.copy(goodFiles[nMost])
            for j in range(nQuads):
                maxDiff = 9999
                for k in range(nEach[i]):
                    nowDiff = np.abs(shTimes[i][k] - shTimes[nMost][j])
                    if nowDiff < maxDiff:
                        quadFiles[i][j] = goodFiles[i][k]
                        maxDiff = nowDiff
            quadFiles[i] = np.array(quadFiles[i])
    quadFiles = np.transpose(quadFiles)
    
    # Add in the full path for each file
    # Doesn't like just updating quadFiles so make new array
    fullQuadFiles = [[] for i in range(nQuads)]
    for j in range(nQuads):
        fullLine = []
        for i in range(4):    
            fullLine.append(obsFiles+'SoloHI/'+ quadFiles[j][i])
        fullQuadFiles[j] = fullLine
        
    # Run the processing
    ims = []
    hdrs = []
    for i in range(nQuads):
        im, hdr = solohi_fits2grid(fullQuadFiles[i])
        ims.append(im)
        hdrs.append(hdr)

    # Add in the actual processing, saving, and output to runFile  
                        
        
# |-------------------------------|
# |-------- Process WISPR --------|
# |-------------------------------|
if doWISPR:    
    # Get all the WISPR files
    WISPRfiles = os.listdir(obsFiles+'WISPR/')
    # Check all the files to see if they have a matching date str
    # If date matches then check the hour/min on first/last date
    # Sort all by in/out, will only process what we want later
    goodFiles =[[], []]
    for aF in WISPRfiles:
        if singleDay:
            if ymds[0] in aF:
                hm = aF[22:26]
                if (hm >= hm0) & (hm <= hmf):
                    if '_V1_1' in aF:
                        goodFiles[0].append(obsFiles+'WISPR/'+aF)
                    elif '_V1_2' in aF:
                        goodFiles[1].append(obsFiles+'WISPR/'+aF)
        else:
            for j in range(nDays):
                if ymds[j] in aF:
                    hm = aF[22:26]
                    addIt = True
                    if (j == 0) & (hm < hm0):
                        addIt = False
                    elif (j == nDays-1) & (hm > hmf):
                        addIt = False
                    if addIt:
                        if '_V1_1' in aF:
                            goodFiles[0].append(obsFiles+'WISPR/'+aF)
                        elif '_V1_2' in aF:
                            goodFiles[1].append(obsFiles+'WISPR/'+aF)
                                    
    # Make an array and make sure sorted               
    for i in range(2):
        goodFiles[i] = np.sort(np.array(goodFiles[i]))
                        
    # Process, save, and add names to runFile  
    if 'In' in whichWISPR:
        print ('|--- Processing WISPR Inner ---|')
        imsI, hdrsI = wispr_prep(goodFiles[0], straylightOff=True)
    if 'Out' in whichWISPR:
        print ('|--- Processing WISPR Outer ---|')
        imsO, hdrsO = wispr_prep(goodFiles[0], straylightOff=True)
        
 