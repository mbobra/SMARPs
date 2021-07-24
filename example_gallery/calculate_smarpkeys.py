"""
Purpose:   To calculate the following spaceweather parameters SoHO/MDI line-of-sight magnetic field data:
           USFLUXL Total unsigned flux in Maxwells
           MEANGBL Mean value of the line-of-sight field gradient, in Gauss/Mm
           CMASKL  Number of pixels used in the USFLUXL and MEANGBL calculation
           R_VALUE Flux along gradient-weighted neutral-line length in Maxwells

Inputs:    All SDO/HMI data is stored in a pSQL database; the web interface is here:
           http://jsoc.stanford.edu/ajax/lookdata.html.

           The data used for this code is available in the DRMS series mdi.smarp_cea_96m,
           which is documented extensively in Bobra et al., ApJ, 2021, an open-access publication:
           [insert DOI here].

           We use the following segments:
           
           [example filename]                 --> [description]
           mdi.smarp_cea_*.bitmap.fits        --> bits indicate result of automatic detection algorithm
           mdi.smarp_cea_*.magnetogram.fits   --> line-of-sight component of the magnetic field

Examples:  
           > python calculate_smarpkeys.py --file_bitmap=files/mdi.smarp_cea_96m.13520.20100714_124800_TAI.bitmap.fits --file_los=files/mdi.smarp_cea_96m.13520.20100714_124800_TAI.magnetogram.fits
"""

# import some modules
import sunpy.map
import scipy.ndimage
import numpy as np
import sys
import math
import argparse
from skimage.measure import block_reduce

# define some constants
radsindeg = np.pi/180.
munaught  = 0.0000012566370614

#===========================================

def main():

    file_bitmap = ''
    file_los     = ''
    
    parser = argparse.ArgumentParser(description='calculate spaceweather keywords from vector magnetic field data')
    parser.add_argument('-a', '--file_bitmap', type=str, help='FITS file with bits identifying the active region', required=True)
    parser.add_argument('-b', '--file_los', type=str, help='FITS file containing line-of-sight component of magnetic field', required=True)
    parser._optionals.title = "flag arguments"
    args = parser.parse_args()

    file_bitmap        = args.file_bitmap
    file_los           = args.file_los

    print('')
    print('These are the files:')
    print('file_bitmap is', file_bitmap)
    print('file_los is', file_los)
    print('')
    
    # get the data
    print('Getting the data.')    
    bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, los, los_err = get_data(file_bitmap, file_los)

    print('These are the keyword values:')
    # compute the gradient-weighted neutral line length
    Rparam, Rparam_err = computeR(los, los_err, nx, ny, cdelt1_arcsec)
    print('R_VALUE ', Rparam,'Mx')
    print('The error in R_VALUE is', Rparam_err)

    # compute mean line-of-sight field gradient
    mean_derivative_blos, mean_derivative_blos_err = computeLOSderivative(los, los_err, nx, ny, bitmap, rsun_ref, rsun_obs, cdelt1_arcsec)
    print('MEANGBL ', mean_derivative_blos,'G * Mm^(-1)')
    print('The error in MEANGBL is', mean_derivative_blos_err)

    # compute the total unsigned flux and associated errors
    mean_vf, mean_vf_err, count_mask  = compute_abs_flux_los(los, los_err, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec)
    print('USFLUXL ',mean_vf,'Mx')
    print('The error in USFLUXL is', mean_vf_err)
    print('CMASKL ', count_mask,'pixels')

    print('Note that the calculation for R_VALUE uses a slightly different method than applied for the mdi.smarp*_96m series. The results, however, should be identical or within a log(R) value of 0.1. ')
    print('All the other keyword calculations use an identical method, and the results are identical. ')

def get_data(file_bitmap, file_los):

    """function: get_data

    This function reads the appropriate data and metadata.
    """

    try:
        bitmap_map = sunpy.map.Map(file_bitmap)
    except:
        print("Could not open the bitmap fits file")
        sys.exit(1)

    try:
        los_map = sunpy.map.Map(file_los)
    except:
        print("Could not open the LoS fits file")
        sys.exit(1)
    
    # get array data
    bitmap            = bitmap_map.data
    los               = los_map.data

    # get metadata
    header = los_map.meta
    
    # get fits header key information
    rsun_ref = header['rsun_ref']
    dsun_obs = header['dsun_obs']
    rsun_obs = header['rsun_obs']
    cdelt1   = header['cdelt1']

    # Note that the value of CDELT1 in mdi.smarp_cea_720s is in units of degrees per pixel.
    # The following calculation converts CDELT1 into arcseconds.
    # Therefore the variable cdelt1_arcseconds is in units of arseconds per pixel.
    cdelt1_arcsec = (math.atan((rsun_ref*cdelt1*radsindeg)/(dsun_obs)))*(1/radsindeg)*(3600.)
    print("The native CDELT1 coordinates for this series is in units of degrees per pixel", cdelt1)
    print("Here is the conversion of CDELT1 to units of arseconds per pixel", cdelt1_arcsec)

    # get dimensions
    nx     = los.shape[1]
    ny     = los.shape[0]

    # Create an error array. Liu et al. (2012) [DOI: 10.1007/s11207-012-9976-x] determined
    # the median noise in the MDI one-minute full-disk magnetic field maps is 26.4 Mx cm^(−2) 
    # (see Figure 9). We will assume this noise value is homogeneous throughout the disk to
    # estimate the error in the keyword quantities. Here, 1 Gauss = 1 Mx cm^(−2).
    los_err = np.ndarray(shape=(ny,nx), dtype=float)
    los_err.fill(26.4)

    return [bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, los, los_err] 

#===========================================

def computeR(los, los_err, nx, ny, cdelt1_arcsec):
    """
    function: computeR

    This function computes R, or the log of the gradient-weighted neutral line length. 
    So the output is unitless.

    This function also computes the error in R. The general formula for an error of a
    function is ERR(f(x)) = d/dx f(x) * ERR(x). Thus
                ERR(R) = d/dx (log_10 R) * ERR(x)
                       = [1/ln(10)] * [1/x] * ERR(x)
                       = ERR(x) / ln(10)*x
    """

    sum   = 0.0
    err   = 0.0
    sigma = 10.0/2.3548
    scale = int(round(2.0/cdelt1_arcsec))

    # =============== [STEP 1] =============== 
    # bin the line-of-sight magnetogram down by a factor of scale
    rim = block_reduce(los, block_size=(scale,scale), func=np.mean)

    # =============== [STEP 2] =============== 
    # identify positive and negative pixels greater than +/- 150 gauss
    # and label those pixels with a 1.0 in arrays p1p0 and p1n0

    nx1  = rim.shape[1]
    ny1  = rim.shape[0]
    p1p0 = np.zeros([ny1,nx1])
    p1n0 = np.zeros([ny1,nx1])

    for j in range(ny1):
        for i in range (nx1):
            if (rim[j,i] > 150):
                p1p0[j,i]=1.0
            else:
                p1p0[j,i]=0.0
            if (rim[j,i] < -150):
                p1n0[j,i]=1.0
            else:
                p1n0[j,i]=0.0

    # =============== [STEP 3] =============== 
    # smooth each of the negative and positive pixel bitmaps by convolving with a boxcar     

    # set up the convolution kernel
    boxcar_kernel = np.zeros([ny1,nx1])
    midpoint_ny1  = int(round(ny1/2))
    midpoint_nx1  = int(round(nx1/2))

    for j in range(midpoint_ny1,midpoint_ny1+3):
        for i in range(midpoint_nx1,midpoint_nx1+3):
            boxcar_kernel[j,i]=0.1111

    p1p = scipy.ndimage.convolve(p1p0,boxcar_kernel)
    p1n = scipy.ndimage.convolve(p1n0,boxcar_kernel)

    # =============== [STEP 4] =============== 
    # find the pixels for which p1p and p1n are both equal to 1. 
    # this defines the polarity inversion line

    p1 = np.zeros([ny1,nx1])
    for j in range(ny1):
        for i in range (nx1):
            if ((p1p[j,i] > 0.0) and (p1n[j,i] > 0.0)):
                p1[j,i]=1.0
            else:
                p1[j,i]=0.0
                
    # =============== [STEP 5] =============== 
    # convolve the polarity inversion line map with a gaussian
    # to identify the region near the plarity inversion line
    # the resultant array is called pmap

    pmap = scipy.ndimage.gaussian_filter(p1,sigma,order=0)

    # =============== [STEP 6] =============== 
    # the R parameter is calculated

    for j in range(ny1):
        for i in range (nx1):
            if np.isnan(pmap[j,i]):
                continue
            if np.isnan(rim[j,i]):
                continue
            sum += pmap[j,i]*abs(rim[j,i])
            err += pmap[j,i]*abs(los_err[j,i])

    if (sum < 1.0):
        Rparam = 0.0
        Rparam_err = 0.0
    else:
        Rparam = math.log10(sum)
        Rparam_err = err / (math.log(10) * sum) # note that math.log is a natural log by default

    return [Rparam, Rparam_err]
    
#===========================================

def computeLOSderivative(los, los_err, nx, ny, bitmap, rsun_ref, rsun_obs, cdelt1_arcsec):

    """function: computeLOSderivative

    This function computes the derivative of the line-of-sight field, or sqrt[(dB/dx)^2 + (dB/dy)^2].
    The native units in the series mdi.smarp_96m and mdi.smarp_cea_96m are in Gauss/pixel.

    Here are the steps to convert from Gauss/pixel to Gauss/Mm:
    The units of the magnetic field, or dB, are in Gauss.
    The units of length, i.e. dx or dy, are in pixels.
    Therefore, the units of dB/dx or dB/dy = (Gauss/pix)(pix/arcsec)(arsec/meter)(meter/Mm), or
                                           = (Gauss/pix)(1/cdelt1_arcsec)(RSUN_OBS/RSUN_REF)(1000000)
                                           = Gauss/Mm

    In other words, multiply MEANGBL by a factor of (1/cdelt1_arcsec)*(RSUN_OBS/RSUN_REF)*(1000000).
    
    Note that cdelt1_arcsec is defined in the get_data function above.
    """

    count_mask = 0
    sum        = 0.0
    err        = 0.0

    derx_blos  = np.zeros([ny,nx])
    dery_blos  = np.zeros([ny,nx])
    err_term1  = np.zeros([ny,nx])
    err_term2  = np.zeros([ny,nx])

    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1,nx-1):
        for j in range(0,ny):
           derx_blos[j,i]   = (los[j,i+1] - los[j,i-1])*0.5
           err_term1[j,i] = ( ((los[j,i+1] - los[j,i-1])*(los[j,i+1]-los[j,i-1])) * (los_err[j,i+1]*los_err[j,i+1] + los_err[j,i-1]*los_err[j,i-1]) )

    #brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0,nx):
        for j in range(1,ny-1):
           dery_blos[j,i]   = (los[j+1,i] - los[j-1,i])*0.5
           err_term2[j,i] = ( ((los[j+1,i]-los[j-1,i])*(los[j+1,i]-los[j-1,i])) * (los_err[j+1,i]*los_err[j+1,i] + los_err[j-1,i]*los_err[j-1,i]) )
    
    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero. 
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the conf_disambig and bitmap arrays.

    i=0
    for j in range(ny):
        derx_blos[j,i] = ( (-3*los[j,i]) + (4*los[j,i+1]) - (los[j,i+2]) )*0.5
        
    i=nx-1
    for j in range(ny):
        derx_blos[j,i] = ( (3*los[j,i]) + (-4*los[j,i-1]) - (-los[j,i-2]) )*0.5
    
    j=0
    for i in range(nx):
        dery_blos[j,i] = ( (-3*los[j,i]) + (4*los[j+1,i]) - (los[(j+2),i]) )*0.5
    
    j=ny-1
    for i in range(nx):
        dery_blos[j,i] = ( (3*los[j,i]) + (-4*los[j-1,i]) - (-los[j-2,i]) )*0.5

    # Calculate the sum only
    for j in range(1,ny-1):
        for i in range (1,nx-1):
            if ( bitmap[j,i] < 36 ):
                continue
            if np.isnan(los[j,i]):
                continue
            if np.isnan(los[j+1,i]):
                continue
            if np.isnan(los[j-1,i]):
                continue
            if np.isnan(los[j,i-1]):
                continue
            if np.isnan(los[j,i+1]):
                continue
            if np.isnan(derx_blos[j,i]):
                continue
            if np.isnan(dery_blos[j,i]):
                continue
            sum += np.sqrt( derx_blos[j,i]*derx_blos[j,i]  + dery_blos[j,i]*dery_blos[j,i]  )
            denominator_1 = 16.0*( derx_blos[j,i]*derx_blos[j,i] + dery_blos[j,i]*dery_blos[j,i])
            denominator_2 = 16.0*( derx_blos[j,i]*derx_blos[j,i] + dery_blos[j,i]*dery_blos[j,i])
            if np.isnan(denominator_1):
                continue
            if np.isnan(denominator_2):
                continue
            if denominator_1 == 0:
                continue
            if denominator_2 == 0:
                continue       
            err += (err_term2[j,i] / denominator_1) + (err_term1[j,i] / denominator_2)
            count_mask += 1
            
    if count_mask == 0:
        mean_derivative_blos = 0.0
        mean_derivative_blos_err = 0.0
    else:
        mean_derivative_blos = (sum)/(count_mask)
        mean_derivative_blos_err = (np.sqrt(err))/(count_mask)

    return [mean_derivative_blos, mean_derivative_blos_err]

#===========================================

def compute_abs_flux_los(los, los_err, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec):

    """function: compute_abs_flux_los

    This function computes the total unsigned flux, on the line-of-sight field, in units of G/cm^2.
    It also returns the number of pixels used in this calculation in the keyword CMASK.
    
    To compute the unsigned flux, we simply calculate
       flux = surface integral [(vector Blos) dot (normal vector)],
            = surface integral [(magnitude Blos)*(magnitude normal)*(cos theta)].

    However, since the field is radial, we will assume cos theta = 1.
    Therefore, the pixels only need to be corrected for the projection.

    To convert G to G*cm^2, simply multiply by the number of square centimeters per pixel: 
       (Gauss/pix^2)(CDELT1)^2(RSUN_REF/RSUN_OBS)^2(100.cm/m)^2
       =Gauss*cm^2
    """

    count_mask = 0
    sum        = 0.0
    err        = 0.0
    
    for j in range(ny):
        for i in range(nx):
            if ( bitmap[j,i] < 36 ):
                continue
            if np.isnan(los[j,i]):
                continue
            sum += abs(los[j,i])
            err += los_err[j,i]*los_err[j,i]
            count_mask += 1

    mean_vf     = sum*cdelt1_arcsec*cdelt1_arcsec*(rsun_ref/rsun_obs)*(rsun_ref/rsun_obs)*100.0*100.0
    mean_vf_err = (np.sqrt(err))*abs(cdelt1_arcsec*cdelt1_arcsec*(rsun_ref/rsun_obs)*(rsun_ref/rsun_obs)*100.0*100.0)

    return [mean_vf, mean_vf_err, count_mask]

#===========================================

if __name__ == "__main__":
    main()
    
__author__ = 'Monica Bobra'
