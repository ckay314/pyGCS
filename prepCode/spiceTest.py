import spiceypy as spice
import os, sys
import numpy as np

global topDir, LSKfile
topDir = '/Users/kaycd1/ssw/psp/gen/data/spice'
otherDir = '/Users/kaycd1/ssw/packages/sunspice/data'
soloDir = '/Users/kaycd1/ssw/so/gen/data/sunspice'
stereoDir = '/Users/kaycd1/ssw/stereo/gen/data/spice'

def setupPSPkernels():
    # Get any loaded files
    num_kernels = spice.ktotal('ALL')
    loadKerns = []
    for i in range(num_kernels):
        filename, kind, source, handle,  = spice.kdata(i, 'ALL')
        loadKerns.append(filename)
    loadKerns = np.array(loadKerns)
    #print (loadKerns)
    
    if os.path.isdir(topDir):
        files = os.listdir(topDir)
        
        # Leap second kernel
        LSKfile = otherDir + '/naif0012.tls'
        if LSKfile not in loadKerns:
            spice.furnsh(LSKfile)
        
        # Load the file that makes ckgp happy
        ckgpFile = otherDir + '/pck00011_n0066.tpc'
        if ckgpFile not in loadKerns:
            spice.furnsh(ckgpFile)
        
        # Load the file that makes ckgp happy
        deFile = otherDir + '/de421.bsp'
        if deFile not in loadKerns:
            spice.furnsh(deFile)
        
        # Gen folder
        if 'gen' in files:
            genF = os.listdir(topDir+'/'+'gen')
            for aFile in genF:
                if aFile[-3:] in ['.tf', '.ti']:
                    thisKern = topDir+'/'+'gen/'+aFile
                    if thisKern not in loadKerns:
                        spice.furnsh(thisKern)
        else:
            print ('Missing gen files')

        # Operations
        if 'operations_sclk_kernel' in files:
            oskF = os.listdir(topDir+'/'+'operations_sclk_kernel')
            allFs = []
            for aFile in oskF:
                if aFile[-4:] == '.tsc':
                    allFs.append(aFile)
            # Just keep the last one
            aFile = np.sort(np.array(allFs))[-1]
            thisKern = topDir+'/'+'operations_sclk_kernel/'+aFile
            if thisKern not in loadKerns:
                spice.furnsh(thisKern)
        else:
            print ('Missing operatiosn_sclk_kernel')
            
        
        # Orbit
        if 'orbit' in files:
            orbF = os.listdir(topDir+'/'+'orbit')
            if len(orbF) == 1:
                thisKern = topDir+'/'+'orbit/'+orbF[0]
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            else:
                print('Multiple orbit files, need to sort out')
                
        # Ephemeris Predict
        if 'ephemeris_predict' in files:
            epF = os.listdir(topDir+'/'+'ephemeris_predict')
            if len(epF) == 1:
                thisKern = topDir+'/'+'ephemeris_predict/'+epF[0]
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            else:
                print('Multiple ephemeris_predict files, need to sort out')
        
        # Recon ephemeris
        if 'reconstructed_ephemeris' in files:
            rephF = os.listdir(topDir+'/'+'reconstructed_ephemeris')
            if len(rephF) > 0:
                for yr in rephF:
                    yrF = os.listdir(topDir+'/'+'reconstructed_ephemeris/'+yr)
                    for aF in yrF:
                        if aF[-4:] == '.bsp':
                            thisKern = topDir+'/'+'reconstructed_ephemeris/'+yr+'/'+aF
                            if thisKern not in loadKerns:
                                spice.furnsh(thisKern)
                            
        # Attitude long term
        if 'attitude_long_term_predict' in files:
            altF = os.listdir(topDir+'/'+'attitude_long_term_predict')
            # loading more than IDL but run with it for now
            for aF in altF:
                thisKern = topDir+'/'+'attitude_long_term_predict/'+aF
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
        
        # Attitude short term
        if 'attitude_short_term_predict' in files:
            astF = os.listdir(topDir+'/'+'attitude_short_term_predict')
            for aF in astF:
                thisKern = topDir+'/'+'attitude_short_term_predict/'+aF
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
        
        # Bonus things
        moreFiles = ['de421.bsp', 'pck00011_n0066.tpc', 'heliospheric.tf',  'sdo_body_name.tf']
        for aFile in moreFiles:
            thisKern = otherDir + '/' + aFile
            if thisKern not in loadKerns:
                spice.furnsh(thisKern)
                
                
    # Solar Orbiter things
    # (also needs the naif0012, de421, pck00011_n0066, heliospheric, sdo_body_name)        
    # Most of this needs to be cleaned up to read in newest version not hardcoded
    if os.path.isdir(soloDir):
        files = os.listdir(soloDir)
        
        if 'gen' in files:
            genF = os.listdir(soloDir+'/'+'gen')
            if 'solo_rtn.tf' in genF:
                thisKern = soloDir+'/'+'gen/solo_rtn.tf'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            if 'solo_ANC_soc-sc-fk_V09.tf' in genF:
                thisKern = soloDir+'/'+'gen/solo_ANC_soc-sc-fk_V09.tf'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            if 'solo_ANC_soc-ops-fk_V02.tf' in genF:
                thisKern = soloDir+'/'+'gen/solo_ANC_soc-ops-fk_V02.tf'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            if 'solo_ANC_soc-sci-fk_V08.tf' in genF:
                thisKern = soloDir+'/'+'gen/solo_ANC_soc-sci-fk_V08.tf'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
        if 'orbit' in files:
            orbF = os.listdir(soloDir+'/'+'orbit')
            if 'solo_ANC_soc-orbit-stp_20200210-20301120_379_V1_00476_V01.bsp' in orbF:
                thisKern = soloDir+'/'+'orbit/solo_ANC_soc-orbit-stp_20200210-20301120_379_V1_00476_V01.bsp'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
        
        if 'att' in files:
            attF = os.listdir(soloDir+'/'+'att')
            for aF in attF:
                if '.bc' in aF:
                    thisKern = soloDir+'/'+'att/'+aF
                    if thisKern not in loadKerns:
                        spice.furnsh(thisKern)
                        
        if 'sclk' in files:
            sclkF = os.listdir(soloDir+'/'+'sclk')
            if 'solo_ANC_soc-sclk-fict_20000101_V01.tsc' in sclkF:
                thisKern = soloDir+'/'+'sclk/'+'solo_ANC_soc-sclk-fict_20000101_V01.tsc'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
        
        if 'ik' in files:
            ikF = os.listdir(soloDir+'/'+'ik')
            for aF in ikF:
                if '.ti' in aF:
                    thisKern = soloDir+'/'+'ik/'+aF
                    if thisKern not in loadKerns:
                        spice.furnsh(thisKern)
      
    # General STEREO things
    if os.path.isdir(stereoDir):
        files = os.listdir(stereoDir)   
        if 'stereo_rtn.tf' in files:     
            thisKern = stereoDir+'/'+'stereo_rtn.tf'
            if thisKern not in loadKerns:
                spice.furnsh(thisKern)
                
        if 'sclk' in files:   
            # Ahead
            moreFs = os.listdir(stereoDir+'/sclk/ahead/')
            asFs = []
            for aF in moreFs:
                if 'ahead_science' in aF:
                    asFs.append(aF)
            asFs = np.sort(np.array(asFs))
            if len(asFs) > 0:
                thisKern = stereoDir+'/sclk/ahead/' + asFs[-1] # take the newest
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            # Behind 
            moreFs = os.listdir(stereoDir+'/sclk/behind/')
            asFs = []
            for aF in moreFs:
                if 'ahead_science' in aF:
                    asFs.append(aF)
            asFs = np.sort(np.array(asFs))
            if len(asFs) > 0:
                thisKern = stereoDir+'/sclk/behind/' + asFs[-1] # take the newest
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
                         
        if 'epm' in files:   
            moreFs = os.listdir(stereoDir+'/epm/ahead/')
            if 'ahead_2017_061_5295day_predict.epm.bsp' in moreFs:
                thisKern = stereoDir+'/epm/ahead/ahead_2017_061_5295day_predict.epm.bsp'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            if 'ahead_2024_226_01.epm.bsp' in moreFs:
                thisKern = stereoDir+'/epm/ahead/ahead_2024_226_01.epm.bsp'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            moreFs = os.listdir(stereoDir+'/epm/behind/')
            if 'behind_2009_049_definitive_predict.epm.bsp' in moreFs:
                thisKern = stereoDir+'/epm/behind/behind_2009_049_definitive_predict.epm.bsp'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            if 'behind_2020_301_baseline_1460day_01.epm.bsp' in moreFs:
                thisKern = stereoDir+'/epm/behind/behind_2020_301_baseline_1460day_01.epm.bsp'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            if 'behind_2024_026_01.epm.bsp' in moreFs:
                thisKern = stereoDir+'/epm/behind/behind_2024_026_01.epm.bsp'
                if thisKern not in loadKerns:
                    spice.furnsh(thisKern)
            
        if 'depm' in files:   
            moreFs = os.listdir(stereoDir+'/depm/ahead/')
            for aF in moreFs:
                if '.depm.bsp' in aF:
                    thisKern = stereoDir+'/depm/ahead/'+aF
                    if thisKern not in loadKerns:
                        spice.furnsh(thisKern)
                moreFs = os.listdir(stereoDir+'/depm/behind/')
                for aF in moreFs:
                    if '.depm.bsp' in aF:
                        thisKern = stereoDir+'/depm/behind/'+aF
                        if thisKern not in loadKerns:
                            spice.furnsh(thisKern)
        
        
def loadSomeSTEREO(yd):
    # STEREO things
    if os.path.isdir(stereoDir):
        files = os.listdir(stereoDir)
        if 'ah' in files:
            moreFs = os.listdir(stereoDir+'/ah/')
            if 'ahead' in moreFs:
                aKs = os.listdir(stereoDir+'/ah/ahead/')
                for aF in aKs:
                    if yd in aF:
                        thisKern = stereoDir+'/ah/ahead/'+aF
                        spice.furnsh(thisKern)

if False:
    setupPSPkernels()
    date = '2025-06-10T00:01:34.971'
    frame = 'PSPHGRTN'
    et = spice.str2et(date)
    sc = -96
    sclkdp = spice.sce2c(sc, et)
    tol = spice.sctiks(int(sc), "1:000")
    cmat, clkout = spice.ckgp(int(sc)*1000, sclkdp, tol, frame)



# Have included
'''
/Users/kaycd1/ssw/packages/sunspice/data/naif0012.tls
/Users/kaycd1/ssw/packages/sunspice/data/de421.bsp
/Users/kaycd1/ssw/packages/sunspice/data/pck00011_n0066.tpc
/Users/kaycd1/ssw/packages/sunspice/data/heliospheric.tf
/Users/kaycd1/ssw/packages/sunspice/data/sdo_body_name.tf
/Users/kaycd1/ssw/psp/gen/data/spice/gen/psp_rtn.tf
/Users/kaycd1/ssw/psp/gen/data/spice/gen/spp_dyn_v201.tf
/Users/kaycd1/ssw/psp/gen/data/spice/gen/spp_v300.tf
/Users/kaycd1/ssw/psp/gen/data/spice/gen/spp_wispr_v002.ti
/Users/kaycd1/ssw/psp/gen/data/spice/gen/spp_sweap_v100.ti
/Users/kaycd1/ssw/psp/gen/data/spice/gen/epihi_v100.ti
/Users/kaycd1/ssw/psp/gen/data/spice/gen/spp_epilo_v100.ti
/Users/kaycd1/ssw/psp/gen/data/spice/gen/spp_fields_v100.ti
/Users/kaycd1/ssw/psp/gen/data/spice/operations_sclk_kernel/spp_sclk_1499.tsc
/Users/kaycd1/ssw/psp/gen/data/spice/orbit/spp_nom_20180812_20250831_v040_RO7.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/ephemeris_predict/spp_pred_20240205_20240625_od203_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2018/spp_recon_20180812_20181008_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2018/spp_recon_20181008_20190120_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2019/spp_recon_20190120_20190416_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2019/spp_recon_20190416_20190914_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2019/spp_recon_20190914_20200101_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2020/spp_recon_20200101_20200301_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2020/spp_recon_20200301_20200505_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2020/spp_recon_20200505_20200705_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2020/spp_recon_20200705_20200802_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2020/spp_recon_20200802_20201016_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2020/spp_recon_20201016_20210101_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20210101_20210226_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20210226_20210325_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20210325_20210525_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20210524_20210723_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20210723_20210904_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20210904_20211104_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20211104_20211217_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2021/spp_recon_20211217_20220329_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2022/spp_recon_20220329_20220620_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2022/spp_recon_20220620_20220725_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2022/spp_recon_20220725_20220923_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2022/spp_recon_20220923_20221030_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2022/spp_recon_20221030_20230124_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2023/spp_recon_20230124_20230402_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2023/spp_recon_20230402_20230522_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2023/spp_recon_20230522_20230723_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2023/spp_recon_20230723_20231008_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2023/spp_recon_20231008_20240201_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/reconstructed_ephemeris/2024/spp_recon_20240201_20240414_v001.bsp
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_long_term_predict/spp_nom_20180812_20250831_v035_RO2_20190129_20250710.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_long_term_predict/spp_2023_184_2024_043_05_RO6_contacts_only.asp.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_long_term_predict/spp_2023_317_2024_141_00_CaseV04_15d_outside_trks_040_to_057_epiLoProtect_start_058.alp.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_short_term_predict/spp_2024_155_2024_176_00.asp.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_short_term_predict/spp_2024_176_2024_197_00.asp.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_short_term_predict/spp_2024_197_2024_218_00.asp.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_short_term_predict/spp_2024_218_2024_232_00.asp.bc
/Users/kaycd1/ssw/psp/gen/data/spice/attitude_short_term_predict/spp_2024_232_2024_246_00.asp.bc
'''