from sunpy.net import Fido, attrs as a
import astropy.units as u
import numpy as np
import wget
from sunpy.time import parse_time
from astropy.io import fits
import os

obsFiles = '/Users/kaycd1/wombat/obsFiles/'

startT = '2012/07/12T16:00'
endT   = '2012/07/12T22:00'

# AIA setup
doAIA = False
AIAwav = [171] # Select from 94, 131, 171*, 193*, 211, 304*, 335, 1600, 1700 (* most common)
AIAtime = 20 # resolution to download in minutes 

# LASCO setup
doLASCO = False
whichLASCO = ['C2', 'C3'] # Select from 'C2' or 'C3'
LASCOtime = 20 # resolution to download in minutes 

# SECCHI setup
doSECCHI = True
whichSECCHI = ['COR2'] # Select from 'EUVI', 'COR1', 'COR2',  'HI1', 'HI2'
EUVIwav     = [171] # Select from 171, 195, 284, 304
EUVItime    = 20 # resolution to download in minutes 
CORtime     = 20 # resolution to download in minutes 

# SoloHI setup
doSoloHI = False

# WISPR setup
doWISPR = False
whichWISPR = ['In', 'Out'] # options of 'In' and 'Out
WISPRtime = 30

# |---------------------------------|
# |----- Set up the time range -----|
# |---------------------------------|
# Get time range in minutes, which we may use later
# since the Sample option in Fido is some times problematic
# make astropy time objects
startAPT = parse_time(startT)
endAPT = parse_time(endT)
timeRange = (endAPT - startAPT).to(u.min).to_value()




# |---------------------------|
# |-- Pull the SDO AIA data --|
# |-------  (level 1) --------|
# |---------------------------|

if doAIA:
    for wv in AIAwav:
        result = Fido.search(a.Time(startT, endT), a.Instrument.aia, a.Wavelength(wv*u.angstrom), a.Sample(AIAtime*u.minute))
        # fileids have format 'aia__lev1:1700:1457020818:1457020818'
        # seems like all level 1 so don't need to filter
        if len(result) > 0:
            print('Downloading AIA files...')
            downloaded_files = Fido.fetch(result[0,:], path=obsFiles+'AIA/{file}') 
        

# |------------------------------|
# |-- Pull the SOHO LASCO data --|
# |-------- (level 0.5) ---------|
# |------------------------------|

if doLASCO:    
    result = Fido.search(a.Time(startT, endT), a.Instrument.lasco, a.Sample(LASCOtime*u.minute))
    # fileids have format /archive/soho/private/data/processed/lasco/level_05/230304/c3/32733517.fts  
    whichC = [[], []]
    ymdts = [[], []]
    
    for i in range(len(result[0]['fileid'])):
        # Don't pull the 512x512 images, have to check file
        if ('/c2/' in result[0]['fileid'][i]) & (result[0]['Size'][i].to_value() > 1.1):
            whichC[0].append(i)
            ymdts[0].append(str(result[0]['Start Time'][i])[0:10].replace('-','') + 'T' + str(result[0]['Start Time'][i])[11:16].replace(':',''))
        elif ('/c3/' in result[0]['fileid'][i]) & (result[0]['Size'][i].to_value() > 1.1):
            whichC[1].append(i)
            ymdts[1].append(str(result[0]['Start Time'][i])[0:10].replace('-','') + 'T' + str(result[0]['Start Time'][i])[11:16].replace(':',''))
            
    if 'C2' in whichLASCO:
        if len(whichC[0]):
            print('Downloading LASCO C2 files...')
            for i in range(len(whichC[0])):
                downloaded_files = Fido.fetch(result[0,whichC[0][i]], path=obsFiles + 'LASCO/' + ymdts[0][i] + '_C2_{file}') 
    if 'C3' in whichLASCO:
        if len(whichC[1]):
            print('Downloading LASCO C3 files...')
            for i in range(len(whichC[1])):
                downloaded_files = Fido.fetch(result[0,whichC[1][i]], path=obsFiles + 'LASCO/' + ymdts[1][i] + '_C3_{file}') 



# |----------------------------|
# |--- Pull the SECCHI data ---|
# |-------  (level ?) ---------|
# |----------------------------|

if doSECCHI:
    # Secchi doesn't like using Sample, randomly starts yeeting things so no files left
    result = Fido.search(a.Time(startT, endT), a.Instrument.secchi)

    # Start pulling EUVI if needed
    if 'EUVI' in whichSECCHI:
        wavidx = [[] for i in range(len(EUVIwav))]
        
        # Sometimes Fido decides to give no wavelengths for EUVI so have to work around
        # Test if all entries have two values, the EUVI none flags will trip this
        waveIssue = False
        try:
            temp = np.where(result[0]['Wavelength'][:,0] == EUVIwav[0]*u.AA)[0]
        except:
            waveIssue = True
            print("Fido has decided not to provide EUVI wavelengths so have to download them all")
        if not waveIssue:
            for i in range(len(EUVIwav)):
                wavidx[i]=np.where(result[0]['Wavelength'][:,0] == EUVIwav[i]*u.AA)[0]
                # Manually downselect 
                fullN = len(wavidx[i])
                expN  = timeRange / EUVItime
                if fullN > expN:
                    downN = int(fullN / expN)
                    wavidx[i] = wavidx[i][::downN]
        else:
            wavidx = []
            for i in range(len(result[0]['fileid'])):
                if ('/euvi/' in result[0]['fileid'][i]) & ('_n7eu' not in result[0]['fileid'][i]):
                    wavidx.append(i)       
                    
            # Downselect because Fido will def lag with too many    
            # There's no good way of picking which ones we keep    
            nCrit = 30
            if len(wavidx) > nCrit:  
                downSel = int(len(wavidx) / nCrit)
                wavidx = wavidx[::downSel]
    
    whichC = [[], [], [], []]    
    for i in range(len(result[0]['fileid'])):
        if '/cor1/' in result[0]['fileid'][i]:
            whichC[0].append(i)
        elif '/cor2/' in result[0]['fileid'][i]:
            whichC[1].append(i)
        elif '/hi_1/' in result[0]['fileid'][i]:
            whichC[2].append(i)
        elif '/hi_2/' in result[0]['fileid'][i]:
            whichC[3].append(i)

    # Need to filter COR1 but keep full pB series
    # the series will have s4c1 tags early that 
    # at some point in switch to n4c1 
    whichC[0] = whichC[0][2:]
    nTimes = len(whichC[0]) / 3
    if nTimes - int(nTimes) != 0:
        if (whichC[0][1]-whichC[0][0]) != 1:
            whichC[0] = whichC[0][1:]
        elif (whichC[0][2]-whichC[0][1]) != 1:
            whichC[0] = whichC[0][2:]
        nTimes = len(whichC[0]) / 3
    expN  = timeRange / CORtime
    if nTimes > expN:
       downN = int(nTimes / expN) * 3
       partial = whichC[0][::downN]    
       full = np.zeros(len(partial)*3, dtype=int)
       for i in range(len(partial)):
           idx = int(partial[i])
           full[3*i: 3*(i+1)]= [idx, idx+1, idx+2]
       whichC[0] = full
    
    # Filter out pB images based on file tags
    # d4c tag is totB, n4c is pB and skipping those (for now at least)
    justTot = []
    for idx in whichC[1]:
        # d4c tag is totB, n4c is pB and skipping those (for now at least)
        # Also only want to pull the ~8 mb images
        if ('d4c' in result[0]['fileid'][idx]) & (result[0]['Size'][idx].to_value() > 4.):
            if len(justTot) == 0:
                justTot.append(idx)
            else:
                if (result[0]['Start Time'][idx] - result[0]['Start Time'][justTot[-1]]).to(u.min).to_value() >= CORtime:
                     justTot.append(idx)
    whichC[1] = justTot    
    
    # Do the DLs        
    if 'EUVI' in whichSECCHI:
        if not waveIssue:
            for i in range(len(wavidx)):
                idxs = wavidx[i]
                if len(idxs) > 0:
                    print('Downloading STEREO EUVI ' + str(EUVIwav[i]) + ' files...')
                    downloaded_files = Fido.fetch(result[0,idxs], path=obsFiles + 'SECCHI/EUVI_' + str(EUVIwav[i]) + 'a_{file}')
        else:
            if len(wavidx)>0:
                print('Downloading Mystery STEREO EUVI files...')
                downloaded_files = Fido.fetch(result[0,wavidx], path=obsFiles + 'SECCHI/EUVI_UNKa_{file}')
                for aF in downloaded_files:
                    with fits.open(aF) as hdulist:
                        im  = hdulist[0].data
                        hdr = hdulist[0].header
                        newName = aF.replace('UNK',str(hdr['WAVELNTH']))
                        print(aF, 'has wavelength ', str(hdr['WAVELNTH']))
                        os.replace(aF,newName)
                
    if 'COR1' in whichSECCHI:
        #print(result[0,whichC[0]])
        if len(whichC[0]) > 0:
            print('Downloading STEREO COR1 files...')
            downloaded_files = Fido.fetch(result[0,whichC[0]], path=obsFiles+'SECCHI/COR1_{file}') 
    if 'COR2' in whichSECCHI:
        #print(result[0,whichC[1]])
        if len(whichC[0]) > 0:
            print('Downloading STEREO COR2 files...')
            downloaded_files = Fido.fetch(result[0,whichC[1]], path=obsFiles+'SECCHI/COR2_{file}') 
    if 'HI1' in whichSECCHI:
        #print(result[0,whichC[2]])
        if len(whichC[0]) > 0:
            print('Downloading STEREO HI1 files...')
            downloaded_files = Fido.fetch(result[0,whichC[2]], path=obsFiles+'SECCHI/HI1_{file}') 
    if 'HI2' in whichSECCHI:
        #print(result[0,whichC[3]])
        if len(whichC[0]) > 0:
            print('Downloading STEREO HI2 files...')
            downloaded_files = Fido.fetch(result[0,whichC[3]], path=obsFiles+'SECCHI/HI2_{file}') 
                
# |----------------------------|
# |--- Pull the SoloHI data ---|
# |-------  (level 2) ---------|
# |----------------------------|
if doSoloHI:
    result = Fido.search(a.Time(startT, endT), a.Instrument.solohi)
    goodIdx = []
    for i in range(len(result[result.keys()[0]]['fileid'])):
        # fileid form data/so/solohi/L2/2023/03/04/solo_L2_solohi-4fg_20230304T155426_V02.fits
        # only take L2
        if '/L2/' in result[result.keys()[0]]['fileid'][i]:
            goodIdx.append(i)
    # CK has given up on the vso downloading bc its almost always a 404 or 403 error
    # but we can wget these files
    # base path = 'https://solohi.nrl.navy.mil/so_data/L2/'
    # followed by YYYYMMDD/samefilename.fits
    basePath = 'https://solohi.nrl.navy.mil/so_data/L2/'
    if len(goodIdx) > 0:
        print('Downloading SoloHI files...')
        for i in goodIdx:
            ogname = result[result.keys()[0]]['fileid'][i]
            fname = ogname[ogname.rfind('/')+1:]
            ymd = fname[19:27]
            myPath = basePath + ymd + '/' + fname
            temp = wget.download(myPath, out=obsFiles+'SoloHI/')
 

# |----------------------------|
# |---- Pull the WISPR data ---|
# |-------  (level 2) ---------|
# |----------------------------|
# The NRL repo does have some V2 files where the L2 data has been reprocessed for 
# straylight (?) (effective 2024/03). Isn't consistent when and where they show up
# so only pulling what Fido/VSO suggests
if doWISPR:       
    result = Fido.search(a.Time(startT, endT), a.Instrument.wispr)
    goodIdx = []
    for i in range(len(result[0]['fileid'])):
        if '/L2/' in result[0]['fileid'][i]:
            goodIdx.append(i)
            #print (result[0]['fileid'][i])
    
    # CK has given up on the vso downloading bc its almost always a 404 or 403 error
    # but we can wget these files
    # base path = 'https://wispr.nrl.navy.mil/data/rel/fits/L2/'
    # followed by YYYYMMDD/samefilename.fits
    basePath  = 'https://wispr.nrl.navy.mil/data/rel/fits/L2/'
    if len(goodIdx) > 0:
        print('')
        print('Downloading PSP WISPR files...')
        #print (len(goodIdx), ' files')
        
        if type(WISPRtime) == type(None):
            WISPRtime = 1 # set at min of one minute res, should grab all
        lastTime = [[],[]]
        
        for i in goodIdx:
            ogname = result[result.keys()[0]]['fileid'][i]
            fname = ogname[ogname.rfind('/')+1:]
            ymd = fname[13:21]
            myPath = basePath + ymd + '/' + fname
            if ('V1_1' in fname) & ('In' in whichWISPR):
                if len(lastTime[0]) == 0:
                    if os.path.exists(obsFiles+'WISPR/'+fname):
                        os.remove(obsFiles+'WISPR/'+fname)
                    temp = wget.download(myPath, out=obsFiles+'WISPR/')
                    lastTime[0].append(i)
                else:
                    if (result[0]['Start Time'][i] - result[0]['Start Time'][lastTime[0][-1]]).to(u.min).to_value() >= WISPRtime:
                        if os.path.exists(obsFiles+'WISPR/'+fname):
                            os.remove(obsFiles+'WISPR/'+fname)
                        temp = wget.download(myPath, out=obsFiles+'WISPR/'+fname)
                        print ('')
                        lastTime[0].append(i)
                    else:
                        print ('Skipping ', result[0]['Start Time'][i])
            elif ('V1_2' in fname) & ('Out' in whichWISPR):
                if len(lastTime[1]) == 0:
                    if os.path.exists(obsFiles+'WISPR/'+fname):
                        os.remove(obsFiles+'WISPR/'+fname)
                    temp = wget.download(myPath, out=obsFiles+'WISPR/')
                    lastTime[1].append(i)
                else:
                    if (result[0]['Start Time'][i] - result[0]['Start Time'][lastTime[1][-1]]).to(u.min).to_value() >= WISPRtime:
                        if os.path.exists(obsFiles+'WISPR/'+fname):
                            os.remove(obsFiles+'WISPR/'+fname)
                        temp = wget.download(myPath, out=obsFiles+'WISPR/'+fname)
                        print ('')
                        lastTime[1].append(i)
                    else:
                        print ('Skipping ', result[0]['Start Time'][i])
            