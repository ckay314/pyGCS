from sunpy.net import Fido, attrs as a
import astropy.units as u
import numpy as np
import wget

obsFiles = '/Users/kaycd1/wombat/obsFiles/'

startT = '2023/3/4T12:00'
endT   = '2023/3/4T16:00'

# AIA setup
doAIA = False
AIAwav = [171, 191, 304] # Select from 94, 131, 171*, 193*, 211, 304*, 335, 1600, 1700 (* most common)
AIAtime = 20 # resolution to download in minutes 

# LASCO setup
doLASCO = False
whichLASCO = ['C2', 'C3'] # Select from 'C2' or 'C3'
LASCOtime = 20 # resolution to download in minutes 

# SECCHI setup
doSECCHI = True
whichSECCHI = ['EUVI', 'COR1', 'COR2',  'HI1', 'HI2'] # Select from 'EUVI', 'COR1', 'COR2',  'HI1', 'HI2'
EUVIwav     = [195, 171] # Select from 171, 195, 284, 304

# SoloHI setup
doSoloHI = False

# WISPR setup
doWISPR = False


# |---------------------------|
# |-- Pull the SDO AIA data --|
# |-------  (level 1) --------|
# |---------------------------|

if doAIA:
    for wv in AIAwav:
        result = Fido.search(a.Time(startT, endT), a.Instrument.aia, a.Wavelength(wv*u.angstrom), a.Sample(AIAtime*u.minute))
        # fileids have format 'aia__lev1:1700:1457020818:1457020818'
        # seems like all level 1 so don't need to filter
        #downloaded_files = Fido.fetch(result[0,:], path=FidoFiles+'AIA/{file}') 
        

# |------------------------------|
# |-- Pull the SOHO LASCO data --|
# |-------- (level 0.5) ---------|
# |------------------------------|

if doLASCO:    
    result = Fido.search(a.Time(startT, endT), a.Instrument.lasco, a.Sample(LASCOtime*u.minute))
    # fileids have format /archive/soho/private/data/processed/lasco/level_05/230304/c3/32733517.fts    
    whichC = [[], []]
    for i in range(len(result[result.keys()[0]]['fileid'])):
        if '/c2/' in result[result.keys()[0]]['fileid'][i]:
            whichC[0].append(i)
        elif '/c3/' in result[result.keys()[0]]['fileid'][i]:
            whichC[1].append(i)
    #if 'C2' in whichLASCO:
    #    downloaded_files = Fido.fetch(result[0,whichC[0]], path=FidoFiles+'LASCO/{file}') 
    #if 'C3' in whichLASCO:
    #    downloaded_files = Fido.fetch(result[0,whichC[1]], path=FidoFiles+'LASCO/{file}') 


# |----------------------------|
# |--- Pull the SECCHI data ---|
# |-------  (level ?) ---------|
# |----------------------------|

if doSECCHI:
    # Secchi doesn't like using Sample, randomly starts yeeting things so no files left
    # 
    result = Fido.search(a.Time(startT, endT), a.Instrument.secchi)

    # Start pulling EUVI if needed
    if 'EUVI' in whichSECCHI:
        wavidx = [[] for i in range(len(EUVIwav))]
        for i in range(len(EUVIwav)):
            wavidx[i]=np.where(resultEUV[0]['Wavelength'][:,0] == EUVIwav[i]*u.AA)[0]
    
    result = Fido.search(a.Time(startT, endT), a.Instrument.secchi)
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
     
    # Do the DLs        
    if 'EUVI' in whichSECCHI:
        for idxs in wavidx:
            print (resultEUV[0,idxs])
            #downloaded_files = Fido.fetch(result[0,idxs[:3]], path=FidoFiles+'SECCHI/{file}')          
    if 'COR1' in whichSECCHI:
        print(result[0,whichC[0]])
        #downloaded_files = Fido.fetch(result[0,whichC[0][:3]], path=FidoFiles+'SECCHI/{file}') 
    if 'COR2' in whichSECCHI:
        print(result[0,whichC[1]])
        #downloaded_files = Fido.fetch(result[0,whichC[1][:3]], path=FidoFiles+'SECCHI/{file}') 
    if 'HI1' in whichSECCHI:
        print(result[0,whichC[2]])
        #downloaded_files = Fido.fetch(result[0,whichC[2][:3]], path=FidoFiles+'SECCHI/{file}') 
    if 'HI2' in whichSECCHI:
        print(result[0,whichC[3]])
        #downloaded_files = Fido.fetch(result[0,whichC[3][:3]], path=FidoFiles+'SECCHI/{file}') 
                
# |----------------------------|
# |--- Pull the SoloHI data ---|
# |-------  (level 2) ---------|
# |----------------------------|
if doSoloHI:
    result = Fido.search(a.Time(startT, endT), a.Instrument.solohi)
    if len(result.keys()) > 1:
        print ('Warning -- more than one source found, running with first one ', result.keys()[0])
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
    for i in goodIdx:
        ogname = result[result.keys()[0]]['fileid'][i]
        fname = ogname[ogname.rfind('/')+1:]
        ymd = fname[19:27]
        myPath = basePath + ymd + '/' + fname
        temp = wget.download(myPath, out=FidoFiles+'SoloHI')
 

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
    for i in goodIdx:
        ogname = result[result.keys()[0]]['fileid'][i]
        fname = ogname[ogname.rfind('/')+1:]
        ymd = fname[13:21]
        myPath = basePath + ymd + '/' + fname
        temp = wget.download(myPath, out=FidoFiles+'WISPR')
            