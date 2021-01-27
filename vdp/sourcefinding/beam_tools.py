"""Functions for computing the primary beam correction
and expected source counts.

"""
import os
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import  EarthLocation, AltAz
from astropy.time import Time
from astropy import wcs
from datetime import datetime

TWOPI = 2*np.pi
RAD2DEG = 180.0/np.pi
DEG2RAD = np.pi/180.0

locVLA=EarthLocation.of_site('vla')
VLALON = locVLA.lon.deg
VLALAT = locVLA.lat.deg


# Define beam offset & dictionary for each primary observing frequency.
# Beams are offset from phase center by amount "offset" in given direction.
# Angles are in deg, "E" of zenith direction (counter-clkwise looking into dish)
# *Note P band primary observing images do NOT need offset correction
offset = 6.5/60  # deg
offset_dict = {'1.5': -174.1,
               '3': 11.6,
               '6': 75.2,
               '10': 113.7,
               '15': -42.4,
               '22': -64.1,
               '33': -106.9,
               '45': -85.5}

#primary beam class
class pribeam(object):
    def __init__(self):
        self.power = None
        self.beamstep = None
        self.phistep = None
        self.nbeam = None
        self.beamsteprad = None
        self.phisteprad = None
        self.error = None

def EQtoPIX(RA, DEC, w):
    """Convert from RA, DEC to pixels
    """
    pixcrd=w.wcs_world2pix([[RA,DEC]],1)
    X = pixcrd[0][0]
    Y = pixcrd[0][1]
    return X,Y

def PIXtoEQ(X, Y, w):
    """Convert from pixels to RA, DEC
    """
    pixcrd=w.wcs_pix2world([[X,Y]],1)
    ra = pixcrd[0][0]
    dec = pixcrd[0][1]
    return ra,dec


def Read_Fitted_Beam(priband):
    """Reads the fitted beam file which contains
    the correction factor as a function of distance
    from beam center at every arcsecond.

    Parameters
    ----------
    priband : str
        Primary frequency of the observations in GHz.
    
    Returns
    -------
    pbobj : object
        Object containing the primary beam
        correction factor ('power') as a function
        of distance from the image center and 
        angle phi from zenith.
        Uncertainty in fitted beam set to 3% for PBmap v3 beams
    """
    pb = pribeam()
    #for PBmap v3 beams:
    beamstep= 1. #arcsec
    nbeam = 16201 #4.5 deg, including end points
    phistep = 5. #deg
    nphi  = 19 #includes end points
    #set file name
    if float(priband) > 14.9:
        fname='FITBEAM2D_15+22+33+45GHz_0.10.txt'
    elif float(priband) > 5.9:
        fname='FITBEAM2D_6+10GHz_0.10.txt'
    elif float(priband) > 2.9:
        fname='FITBEAM2D_3GHz_0.15.txt'
    else:
        fname='FITBEAM2D_1.5GHz_0.10.txt'
    print(fname)
    beamfile=path.join('/home/vpipe/VLITE/vdp/resources',fname)
    f = open(beamfile,'r')
    b = f.readlines()
    f.close() 
    n = len(b)
    if n != nphi:
        print('ERROR! %d lines in %s, expected %d!' % (n,beamfile,nphi))
        sys.exit(1)
    FitBeam = []
    for i in range(nphi): #loop over phi lines
        d = b[i].split()# split line
        for j in d:
            FitBeam.append(float(j))
    #now fill in phi = 90-180 deg
    for i in range(nphi-2,-1,-1):
        d = b[i].split()# split line
        for j in d:
            FitBeam.append(float(j))
    #phi = 180-270
    for i in range(1,nphi):
        d = b[i].split()# split line
        for j in d:
            FitBeam.append(float(j))
    #phi = 270-360
    for i in range(nphi-2,-1,-1):
        d = b[i].split()# split line
        for j in d:
            FitBeam.append(float(j))
    pb.power = np.array(FitBeam)
    pb.beamstep = beamstep #arcsec
    pb.phistep = phistep #deg
    pb.beamsteprad = beamstep*DEG2RAD/3600 #rad
    pb.phisteprad = phistep*DEG2RAD #rad
    pb.nbeam = nbeam
    pb.error = 0.03
    return pb


def Get_Perley_Beam(freq=344.): #freq in MHz
    """Returns Perley P band primary beam from EVLA Memo 195 
    as numpy array. This beam is 1D (symmetric in phi)
    Currently only returns 344 MHz beam
    Returns
    -------
    pbobj : object
        Object containing the primary beam
        correction factor ('power') as a function
        of distance from the image center 
        Uncertainty in fitted beam set to 3%
    """
    pb = pribeam()
    beamstep=1.0 #arcsec
    nbeam = 16201 #4.5 deg, including end points
    freqGHz=freq*1e-3 #GHz
    arcsec=np.arange(nbeam)
    theta=arcsec/3600. #sample pts, angle in deg. Perley fits only to <~ 3deg
    capr=theta*60. #arcmin
    r= freqGHz*capr #Perley's "normalized" angle units
    FitBeam = []
    ##coeffs at 344 MHz (Table 2 in memo)
    a0 =  1.0
    a2 = -0.974e-3
    a4 =  4.09e-7
    a6 = -0.76e-10
    a8 =  0.53e-14
    for i in range(nbeam):
        power  = a0+(a2*np.power(r[i],2))+(a4*np.power(r[i],4))+(a6*np.power(r[i],6))+(a8*np.power(r[i],8))
        FitBeam.append(power)
    pb.power = np.array(FitBeam)
    pb.beamstep = beamstep
    pb.phistep = None #Beam is 1D symmetric for dedicated P band observing
    pb.beamsteprad = beamstep*DEG2RAD/3600
    pb.phisteprad = None
    pb.nbeam = nbeam
    pb.error = 0.03
    return pb


def Get_Primary_Beams():
    """Returns dictionary with primary beam for each
    primary observing band
    """
    pribands=['0.3','1.5','3','6','10','15','22','33','45']
    pbdic={}
    for priband in pribands:
        if priband=='0.3':
            pbdic[priband] = Get_Perley_Beam()
        else:
            pbdic[priband] = Read_Fitted_Beam(priband)
    return pbdic



def Find_Beam_Center(imobj, parang, dxpix=0, dypix=0):
    """Return (RA,Dec) as SkyCoord object & X,Y [pixels] of beam center
    for image at time of parallactic angle parang.
    * dxpix is optional offset added to account for RA drift (i.e. VCSS snapshots)
    * dypix is optional offset added to account for RA drift (i.e. VCSS snapshots)
    """
    if imobj.priband == '0.3': #300 MHz needs no offset correction
        xbeam = imobj.xref
        ybeam = imobj.yref
        beam_center = SkyCoord(imobj.obs_ra,imobj.obs_dec,unit='deg')
        return beam_center,xbeam,ybeam
    #find coords of beam center
    x0 = imobj.xref
    y0 = imobj.yref
    theta = parang - offset_dict[priband]
    if theta > 180: theta -= 360
    if theta < -180: theta += 360
    ##beam center is 'offset' deg in direction theta
    offset_pix = offset*3600/imobj.pixel_scale # [pixels]
    delx = -1*offset_pix*np.sin(theta*DEG2RAD)
    dely = offset_pix*np.cos(theta*DEG2RAD)
    xbeam = x0 + delx + dxpix
    ybeam = y0 + dely + dypix
    rabeam,decbeam = PIXtoEQ(xbeam,ybeam,imobj.wcsobj)
    beam_center = SkyCoord(rabeam,decbeam,unit='deg')
    return beam_center,xbeam,ybeam



def Calc_Beampix_One(cenx, ceny, x, y, pbdic, imobj, parang):
    """Calculates primary beam value at all pixels 
    at one time given by the parallactic angle

    Parameters
    ----------
    cenx, ceny : float
        Pixel coords of beam center (theta=0)
    x, y : np.ndarray
        x, y values of pixels
    pbdic : dictionary
        Primary beam dictionary
    imobj : object
        Image object
    parang : float
        Parallactic angle [rad]

    Returns
    -------
    theta : np.ndarray
        Primary beam values at all pixels
    """
    pb = pbdic[imobj.priband]
    thetapix = np.sqrt(np.power(x-cenx,2)+np.power(y-ceny,2))
    theta    = ((thetapix*imobj.pixel_scale)+(0.5*pb.beamstep))/pb.beamstep #arcsec
    if imobj.priband=='0.3':
        theta = pb.power[np.asarray(theta,dtype=int)]
    else:
        phi = np.arctan2(cenx-x,y-ceny) - parang #radians
        phi = np.mod(phi,TWOPI)
        phi = (phi+(0.5*pb.phisteprad))/pb.phisteprad
        index = np.asarray(phi,dtype=int)*pb.nbeam + np.asarray(theta,dtype=int)
        theta = pb.power[np.asarray(index,dtype=int)] #1:12
    return theta


def zeroimg(x, y):
    """Sets image pixels to 0
    """
    val = 0.*x + 0.*y
    return val


def Calc_Parang(mjd, ra, dec):
    """Calculates parallactic angle of (RA,Dec) at time mjd
    """
    t = Time(mjd,format='mjd')
    lst = t.sidereal_time('apparent',VLALON).hour #hrs
    hrang = (15*lst) - ra #deg
    if hrang > 180: hrang -= 360
    if hrang < -180: hrang += 360
    tmp1 = np.sin(hrang*DEG2RAD)
    tmp2 = np.tan(VLALAT*DEG2RAD)*np.cos(dec*DEG2RAD) - np.sin(dec*DEG2RAD)*np.cos(hrang*DEG2RAD)
    parang = np.arctan2(tmp1,tmp2)*RAD2DEG # [deg]
    return parang 


def Calc_Beam_Image(imobj, pbdic, smeartime):
    """Calculates primary beam image

    Parameters
    ----------
    imobj : object
        Image object
    pbdic : dictionary
        Primary beam dictionary
    smeartime : float
        Max time step for smearing beam [s]

    Returns
    -------
    bmimg : np.ndarray
        Primary beam image
    """
    y,x = np.indices((imobj.naxis2,imobj.naxis1))
    ntime    = int(math.ceil(imobj.duration/smeartime)) + 1
    deltat = duration/(ntime-1) #[s]
    deltatday = deltat/86400

    #initialize beam image to 0
    bmimg=zeroimg(x,y)
    #loop over times
    for i in range(ntime):
        #Set mjd time
        mjdtime0 = imobj.mjdtime + (i*deltatday)
        #Calc parang
        parang = Calc_Parang(mjdtime0,imobj.obs_ra,imobj.obs_dec)
        #Calc beam center
        beam_center,xbeam,ybeam = Find_Beam_Center(imobj,parang)
        #Calc beam at this time, add to beam image
        bmimg += Calc_Beampix_One(xbeam,ybeam,x,y,pbdic,imobj,parang*DEG2RAD)
    #normalize beam image
    bmimg /= ntime    
    return bmimg


'''
def find_nearest_pbcorr(angle, pri_freq):
    """Finds the primary beam correction factor
    corresponding to the angle value closest to 
    the given distance from the image center.

    Parameters
    ----------
    angle : float
        Distance from the image center in degrees.
    pri_freq : float
        Primary frequency of the observations in GHz.

    Returns
    -------
    pbdata['power'][idx] : float
        Primary beam correction factor.
    err : float
        Uncertainity of the primary beam correction. 
        Add in quadrature to flux uncertainty
    """
    pbdata,err = read_fitted_beam(pri_freq)
    idx = (np.abs(np.array(pbdata['angle']) - angle)).argmin()

    return pbdata['power'][idx],err
'''

def expected_nsrc(pri_freq, rms, max_angle=1.5, sigma=5.):
    """Calculates expected number of sources that would
    have been detected within a maximum distance from the
    image center (max_angle) based on the detection threshold
    (sigma) in the image from logN-logS fit of WENSS catalog 
    scaled to 341 MHz from 325 MHz with alpha= -0.7, rms at 
    center of image (mJy), and fitted beam.

    Parameters
    ----------
    rms : float
        Estimate of the noise in the center of the image
        in mJy/beam.
    max_angle : float, optional
        Maximum angular distance in degrees from image center
        within which to count/estimate the number of sources.
        Default value is 1.5 degrees.
    sigma : float, optional
        The detection threshold used when source finding.
        Default value is 5.

    Returns
    -------
    nexp : int
        The expected number of sources within 'max_angle'
        given the detection threshold, 'sigma', and the
        image noise.
    """
    wenss_n0 = 10.637 # number per square degree
    wenss_S0 = 247.569 # mJy
    wenss_alpha1 = 1.846
    wenss_alpha2 = 0.332
    wenss_n = 0.497
    wenss_n2 = -1.0 / wenss_n
    # Read fitted beam file
    pbdata,err = read_fitted_beam(pri_freq)
    nexp = 0
    for i in range(len(pbdata['angle'])):
        if pbdata['angle'][i] > max_angle:
            break
        flux = ((sigma * rms) / pbdata['power'][i]) / wenss_S0
        nsrc = pbdata['domega'][i] * wenss_n0 * (
            (flux**(wenss_alpha1*wenss_n) \
             + flux**(wenss_alpha2*wenss_n))**wenss_n2)
        nexp += nsrc

    return nexp
