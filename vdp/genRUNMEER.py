import numpy as np
import astropy.io
from astropy.io import fits
from scipy.stats import sigmaclip
import sys,os


fout=open('run_meerkat.sh','w')

fields = ['G000+11.5','G000-14.5','G000-2.5','G000+8.5','G003-2.5','G357+2.5','G357-5.5','G000-11.5','G000-8.5','G003+5.5','G357-2.5','G000-17.5','G000+2.5','G000-5.5','G003+2.5','G003-5.5','G357+5.5']

#polarization?
allpol = ['P','I','V','U','Q']

#bands?
allbands = ['12','4','2','FULL']

for field in fields:
    for pol in allpol:
        for band in allbands:
            if band=='FULL': 
                filename = 'meerkat_FULLBAND_'+pol+'.yaml'
            else:
                filename = 'meerkat_'+band+'bands_'+pol+'.yaml'
            fout.write('python3 vdp.py config_files/meerkat/%s/%s --ignore_prompt\n' % (field,filename))
