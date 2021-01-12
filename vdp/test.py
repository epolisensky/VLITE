import math
from math import sqrt
from math import exp
import os.path  as path
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft
from astropy.convolution import convolve
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


print(astropy.__version__)


filename='/roadraid/vpipe/processed/vcss/Mosaics/T09t01/Images/T09t01.J000003-030.IMSC.fits'
f1=filename
hdu_list=fits.open(f1)
hdr=fits.getheader(filename, 0)

#print(hdu_list[0].header['DATE-OBS'])
print(hdr['DATE-OBS'])
try:
	print(hdr['STARTIME'])
except:
	print('pooop')
