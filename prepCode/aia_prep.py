import sunpy.map
import astropy.units as u
from aiapy.calibrate import register, update_pointing, respike, correct_degradation
from aiapy.calibrate.util import get_pointing_table
from aiapy.psf import deconvolve
                
def aia_prep(filesIn, downSize=1024):
    # No reason to port idl, we have aiapy and wont be doing any WL
    # mass calc so assuming a close enough match
    
    
    
    # |------------------------------------------------------|
    # |------------- Loop to process each image -------------|
    # |------------------------------------------------------|
    num = len(filesIn)
    maps_out = []

    # Assume were working from files and not something loaded from read_sdo
    for i in range(num):
        #print ('Processing AIA image '+str(i+1) + ' out of '+str(num))
        # aia files are compressed/different from secchi so doensn't work with
        # straight up fits read, but using sunpy map equiv to read_sdo
        aia_map = sunpy.map.Map(filesIn[i])
                
        # Make range wide enough to get closest 3-hour pointing
        pointing_table = get_pointing_table("JSOC", time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h))
        aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)
        
        # Respike
        # This fails bc of a JSOC TLS/SSL certificate error
        #aia_map_respike = respike(aia_map_updated_pointing)
        
        # Deconvolve
        # This either runs inifinitely slow or catches in loop so not running
        #aia_map_decon = deconvolve(aia_map_updated_pointing)
        
        # Registration
        aia_map_registered = register(aia_map_updated_pointing)
        
        if aia_map_registered.dimensions.x.to_value() > downSize:
            new_dimensions = [downSize, downSize] * u.pixel
            aia_map_registered = aia_map_registered.resample(new_dimensions)
        
        # Normalize by exposure time
        aia_map_normed = aia_map_registered / aia_map_registered.exposure_time
        
        
        maps_out.append(aia_map_normed)
    return maps_out 
    
#fname = ['/Users/kaycd1/wombat/obsFiles/AIA/aia_lev1_304a_2023_03_04t14_00_05_15z_image_lev1.fits']
#aiaMap = aia_prep(fname)
