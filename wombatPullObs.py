#!/usr/bin/env python

"""
Module to pull observations for the wombat GUI using Sunpy Fido/wget

Note that this only works as well as Fido works and there seems to be a lot
of issues in terms of communication with JSOC, especially with VPNs and firewalls
and whatnot. We already pull PSP/SolO using wget as they are available online
and we can take the Fido find results and avoid the download issues. Fixing Fido
is outside of wombat scope.


Inputs:
    times: an array with [startTime, endTime] where both
           times are of the format YYYY-MM-DDTHH:MM 
           (or any format that sunpy parse_time likes)
    
    inst: an array with the tags for which data to pull
        Available tags:
            AIAnum  = SDO AIA where num represents a wavelength from [94, 131, 171*, 
                      193*, 211, 304*, 335, 1600, 1700] with * most common
            C2      = LASCO C2
            C3      = LASCO C3
            COR1    = STEREO COR1
            COR2    = STEREO COR2    
            EUVInum = STEREO EUVI where num is a wavelength from [171, 195, 284, 304]
            HI1     = STEREO HI1
            HI2     = STEREO H2
            SolOHI  = All quadrants from Solar Orbiter HI
            WISPR   = Both inner and outer from PSP WISPR
            WISPRI  = Inner only from PSP WISPR
            WISRPO  = Outer only from PSP WISPR
            * all STEREO values will pull A and B (as available)
    
Optional Inputs:
    EUVtime: time resolution in minutes to use for EUV images
    
    CORtime: time resolution in minutes to use for coronagraph images

    HItime:  time resolution in minutes to use for HI images
    
    ---> default values of 10, 20, 30. Will downselect only if avail at higher res

Outputs:
    The code will set up a folder structure at outFolder with nested satellite
    and instument folders (if it doesn't already exist). It will then dump the 
    unprocessed files into the appropriate folders

Usage:
    from wombatPullObs import pullObs
    times = ['YYYY/MM/DDTHH:MM', 'YYYY/MM/DDTHH:MM']
    sats  = ['AIA', 'COR2', 'WISPR']
    waves = [171]
    outFolder = '/Users/kaycd1/wombat/obsFiles/'
    pullObs(times, sats, waves, outFolder=outFolder)

"""


from sunpy.net import Fido, attrs as a
import astropy.units as u
import numpy as np
import wget
from sunpy.time import parse_time
from astropy.io import fits
import os, sys


if False:
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


# |------------------------------------------------------------|
# |------------ Check/Make Output Folder Structure ------------|
# |------------------------------------------------------------|
def setupFolderStructure(topName):
    """
    Helper function to set up the folder structure for the wombat observations

    Inputs:
        topName: the name of the top level directory

    Outputs:
        None
    
    """
    # |---------------------------------------|
    # |---- Check the top level directory ----|
    # |---------------------------------------|
    if not os.path.exists(topName):
        # Make if if it doesn't exist
        os.mkdir(topName)
    
    # |---------------------------------------|
    # |------------ Check subdirs ------------|
    # |---------------------------------------|
    # Satellite directory names 
    satDirs = ['PSP', 'SDO', 'SOHO', 'SolO', 'STEREO']
    # Instrument/wavelength directory names
    instDirs = [['WISPR','WISPR/Inner', 'WISPR/Outer'], 
                ['AIA', 'AIA/94', 'AIA/131', 'AIA/171', 'AIA/193', 'AIA/211', 'AIA/304', 'AIA/335', 'AIA/1600', 'AIA/1700'], 
                ['LASCO', 'LASCO/C2', 'LASCO/C3'],
                ['SoloHI'], 
                ['EUVI', 'EUVI/171', 'EUVI/195', 'EUVI/284', 'EUVI/304', 'COR1', 'COR2', 'HI1', 'HI2']]    
    # Check each dir/subdir
    for i in range(len(satDirs)):
        if not os.path.exists(topName+'/'+satDirs[i]):
            os.mkdir(topName+'/'+satDirs[i])
        for subDir in instDirs[i]:
            if not os.path.exists(topName+'/'+satDirs[i]+'/'+subDir):
                os.mkdir(topName+'/'+satDirs[i]+'/'+subDir)


# |------------------------------------------------------------|
# |------------------ Grab AIA Observations -------------------|
# |------------------------------------------------------------|
def pullAIA(times, wavs, EUVtime=10, outFolder='obsFiles/'):
    """
    Function to pull level 1 AIA observations using FIDO

    Inputs:
        times: an array with [startTime, endTime] where both
        
        wavs:  an array of wavelength strings
    
    Optional Inputs:
        EUVtime: time resolution in minutes, will set min spacing between images
        
        outFolder: top folder, results will be saved in outFolder/SDO/AIA/wav/

    Outputs:
        The downloaded fits files will be placed in the appropriate folders within outFolder

    """
    
    for wv in wavs:
        # Search Fido for a wavelength, sample time works here
        result = Fido.search(a.Time(times[0], times[1]), a.Instrument.aia, a.Wavelength(int(wv)*u.angstrom), a.Sample(EUVtime*u.minute))
        # If we find things then download
        if len(result) > 0:
            print('Downloading AIA files...')
            downloaded_files = Fido.fetch(result[0,:], path=outFolder+'SDO/AIA/'+wv+'{file}') 
        else:
            print ('Cannot find any files for AIA '+wv)
    

# |------------------------------------------------------------|
# |----------------- Grab LASCO Observations ------------------|
# |------------------------------------------------------------|
def pullLASCO(times, insts, CORtime=20, outFolder='obsFiles/'):
    """
    Function to pull level 1 AIA observations using FIDO

    Inputs:
        times: an array with [startTime, endTime] where both
        
        insts: strings for selected instruments (C2 or C3)
    
    Optional Inputs:
        CORtime: time resolution in minutes, will set min spacing between images
        
        outFolder: top folder, results will be saved in outFolder/LASCO/C#/

    Outputs:
        The downloaded fits files will be placed in the appropriate folders within outFolder

    """
    # |-----------------------------|
    # |--------- Searching ---------|
    # |-----------------------------|
    # Search Fido, sample works here
    result = Fido.search(a.Time(times[0], times[1]), a.Instrument.lasco, a.Sample(CORtime*u.minute))

    
    # |-----------------------------|
    # |---------- Sorting ----------|
    # |-----------------------------|
    # arrays to hold file names and string datetimes
    whichC = [[], []]
    ymdts = [[], []]


    # Loop through what we found and store it if we like it
    for i in range(len(result[0]['fileid'])):
        # Don't pull the 512x512 images, have to check file size
        # Check if c2 via the file path
        if ('/c2/' in result[0]['fileid'][i]) & (result[0]['Size'][i].to_value() > 1.1):
            whichC[0].append(i)
            ymdts[0].append(str(result[0]['Start Time'][i])[0:10].replace('-','') + 'T' + str(result[0]['Start Time'][i])[11:16].replace(':',''))
            
        # Check if c3 via the file path
        elif ('/c3/' in result[0]['fileid'][i]) & (result[0]['Size'][i].to_value() > 1.1):
            whichC[1].append(i)
            ymdts[1].append(str(result[0]['Start Time'][i])[0:10].replace('-','') + 'T' + str(result[0]['Start Time'][i])[11:16].replace(':',''))
    
    
    # |-----------------------------|
    # |-------- Downloading --------|
    # |-----------------------------|
    # LASCO C2   
    if 'C2' in insts:
        if len(whichC[0]):
            print('Downloading LASCO C2 files...')
            c2path = outFolder + 'LASCO/C2/'
            for i in range(len(whichC[0])):
                downloaded_files = Fido.fetch(result[0,whichC[0][i]], path= c2path + ymdts[0][i] + '_C2_{file}') 
        else:
            print('Cannot find any files for LASCO C2')
    # LASCO C3         
    if 'C3' in insts:
        if len(whichC[1]):
            print('Downloading LASCO C3 files...')
            c3path = outFolder + 'LASCO/C3/'
            for i in range(len(whichC[1])):
                downloaded_files = Fido.fetch(result[0,whichC[1][i]], path= c3path + ymdts[1][i] + '_C3_{file}') 
        else:
            print('Cannot find any files for LASCO C3')


# |------------------------------------------------------------|
# |----------------- Grab WISPR Observations ------------------|
# |------------------------------------------------------------|
def pullWISPR(times, insts, HItime=20, outFolder='obsFiles/'):
    """
    Function to pull level 2 WISPR observations using a combo of Fido search
    and wget from the NRL site directly
    
    The NRL repo does have some V2 files where the L2 data has been reprocessed for 
    straylight (?) (effective 2024/03). Isn't consistent when and where they show up
    so only pulling what Fido/VSO suggests

    Inputs:
        times: an array with [startTime, endTime] where both
        
        insts: strings for selected instruments (WISPRI, WISPRO)
    
    Optional Inputs:
        HItime: time resolution in minutes, will set min spacing between images
        
        outFolder: top folder, results will be saved in outFolder/PSP/WISPR/

    Outputs:
        The downloaded fits files will be placed in the appropriate folders within outFolder

    """
    # Fido usually ok searching at least
    result = Fido.search(a.Time(times[0], times[1]), a.Instrument.wispr)

    
    # Pull the L2 data
    goodIdx = []
    if len(result) > 0:
        for i in range(len(result[0]['fileid'])):
            if '/L2/' in result[0]['fileid'][i]:
                goodIdx.append(i)
    else:
        print ('Cannot find any files for LASCO C2')
 
    # CK has given up on the vso downloading bc its almost always a 404 or 403 error
    # but we can wget these files
    # base path = 'https://wispr.nrl.navy.mil/data/rel/fits/L2/'
    # followed by YYYYMMDD/samefilename.fits
    basePath  = 'https://wispr.nrl.navy.mil/data/rel/fits/L2/'
    if len(goodIdx) > 0:
        print('')
        print('Downloading PSP WISPR files...')
           
        # Track the last time we pulled so we can get the res we want
        lastTime = [[],[]]
    
        for i in goodIdx:
            # Full path name from Fido
            ogname = result[result.keys()[0]]['fileid'][i]
            # Just the file name
            fname = ogname[ogname.rfind('/')+1:]
            # Just the date
            ymd = fname[13:21]
            # Build the NRL paths
            myPath = basePath + ymd + '/' + fname
            # Local Paths
            inPath = outFolder+'PSP/WISPR/Inner/'
            outPath = outFolder+'PSP/WISPR/Outer/'
            
            # Inner cases
            if ('V1_1' in fname) & ('WISPRI' in insts):
                # Save the first one, always
                if len(lastTime[0]) == 0:
                    # Make sure we don't already have it
                    print ('Downloading ', fname)
                    if not os.path.exists(inPath+fname):
                        temp = wget.download(myPath, out=inPath)
                    lastTime[0].append(i)
                # Otherwise check how long its been since last
                else:
                    if (result[0]['Start Time'][i] - result[0]['Start Time'][lastTime[0][-1]]).to(u.min).to_value() >= HItime:
                        print ('Downloading ', fname)
                        if not os.path.exists(inPath+fname):
                           temp = wget.download(myPath, out=inPath)
                        lastTime[0].append(i)
                   
            # Outer cases
            elif ('V1_2' in fname) & ('WISPRO' in insts):
                if len(lastTime[1]) == 0:
                    print ('Downloading ', fname)
                    if not os.path.exists(outPath+fname):
                        temp = wget.download(myPath, out=outPath)
                    lastTime[1].append(i)
                else:
                    if (result[0]['Start Time'][i] - result[0]['Start Time'][lastTime[1][-1]]).to(u.min).to_value() >= HItime:
                        print ('Downloading ', fname)
                        if not os.path.exists(outPath+fname):
                            temp = wget.download(myPath, out=outPath)
                        lastTime[1].append(i)
    



# |------------------------------------------------------------|
# |------------ Main function to pull observations ------------|
# |------------------------------------------------------------|
def pullObs(times, insts, outFolder='obsFiles/', EUVtime=10, CORtime=20, HItime=30):
    """
    Main wrapper function for pulling observations

    Inputs:
        times: an array with [startTime, endTime] where both
            times are of the format YYYY-MM-DDTHH:MM 
            (or any format that sunpy parse_time likes)
    
        inst: an array with the tags for which data to pull
            Available tags:
                AIAnum  = SDO AIA where num represents a wavelength from [94, 131, 171*, 
                          193*, 211, 304*, 335, 1600, 1700] with * most common
                C2      = LASCO C2
                C3      = LASCO C3
                COR1    = STEREO COR1
                COR2    = STEREO COR2    
                EUVInum = STEREO EUVI where num is a wavelength from [171, 195, 284, 304]
                HI1     = STEREO HI1
                HI2     = STEREO H2
                SolOHI  = All quadrants from Solar Orbiter HI
                WISPR   = Both inner and outer from PSP WISPR
                WISPRI  = Inner only from PSP WISPR
                WISRPO  = Outer only from PSP WISPR
                * all STEREO values will pull A and B (as available)
    
    Optional Inputs:
        EUVtime: time resolution in minutes to use for EUV images
    
        CORtime: time resolution in minutes to use for coronagraph images

        HItime:  time resolution in minutes to use for HI images
    
        ---> default values of 10, 20, 30. Will downselect only if avail at higher res
    

    Outputs:
        The downloaded fits files will be placed in the appropriate folders within outFolder

    """

    # |---------------------------------------| 
    # |---- Check the top level directory ----|
    # |---------------------------------------| 
    setupFolderStructure(outFolder)
    
    
    # |---------------------------------------| 
    # |---- Check all inst keys are valid ----|
    # |---------------------------------------| 
    goods = ['AIA94', 'AIA131', 'AIA171','AIA193','AIA211','AIA304','AIA335','AIA1600','AIA1700', 'C2', 'C3', 'COR1', 'COR2', 'EUVI171', 'EUVI195', 'EUVI284', 'EUVI304', 'HI1', 'HI2', 'SoloHI', 'WISPR', 'WISPRI', 'WISPRO']
    quitIt = False
    for inst in insts:
        if inst not in goods:
            print ('Unknown instrument tag', inst)
            quitIt = True
    if quitIt:
        sys.exit('Quitting pullObs since passed invalid instrument tag')
        
        
    # |---------------------------------------| 
    # |--------- Set up time things ----------|
    # |---------------------------------------|    
    # Get time range in minutes, which we may use later
    # since the Sample option in Fido is some times problematic
    # make astropy time objects
    startAPT = parse_time(times[0])
    endAPT = parse_time(times[1])
    timeRange = (endAPT - startAPT).to(u.min).to_value()
    
    
    # |-------------------------------------|    
    # |------- Loop through each sat -------|   
    # |-------------------------------------|   
        
    # |------------- AIA -------------|
    doAIA = []
    for inst in insts:
        if 'AIA' in inst:
            # If found just save the wavelength
            doAIA.append(inst.replace('AIA',''))
    if len(doAIA) > 0:
        pullAIA(times, doAIA, EUVtime=EUVtime)
            
            
    
    # |------------ SOHO -------------|
    doLASCO = []
    if 'C2' in insts: doLASCO.append('C2')
    if 'C3' in insts: doLASCO.append('C3')
    if len(doLASCO) > 0:
        pullLASCO(times, doLASCO, CORtime=CORtime)

    # |------------ STEREO -----------|


    # |------------ WISPR ------------|
    doWISPR = []
    if 'WISPR' in insts: doWISPR = ['WISRPI', 'WISPRO']
    elif 'WISPRI' in insts: doWISPR.append('WISPRI')
    elif 'WISPRO' in insts: doWISPR.append('WISPRO')
    if len(doWISPR) > 0:
        pullWISPR(times, doWISPR, HItime=HItime)
    

    # |------------- SolOHI -------------|


if __name__ == '__main__':
    startT = '2023/09/24T16:00'
    endT   = '2023/09/24T20:00'
    times = [startT, endT]
    sats  = ['WISPR', 'COR2']
    pullObs(times, sats)
            