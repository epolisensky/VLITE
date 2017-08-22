"""radio_cat_match.py matches a catalog or database of radio
sources to other large radio survey catalogs. Matched sources
can be output as text files and/or inserted into a database.

Adapted from Emil's VSLOW.py.

Survey catalogs for matching (listed in the order they are checked):
1. TGSS - GMRT 150 MHz, 25" res.
2. NVSS - VLA 1.4 GHz, 45" res.
3. FIRST - VLA 1.4 GHz, 5" res.
4. SUMSS - MOST 843 MHz, 45" res.
5. WENSS - WSRT 325 MHz, 54" res."""


import os
import sys
import sqlite3
import catalogio
import crossmatch


def dbconnect(dbname):
    """Creates an sqlite3 cursor object and specifies
    row_factory so that fetch commands return the rows
    as dictionaries."""
    conn = sqlite3.connect(dbname)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    return cur


def search_range(image):
    """Returns the minimum & maximum RA, Dec of an image (in degrees).
    Used to limit the number of catalog sources extracted for 
    cross-matching."""
    try:
        # image must be a dictionary
        imsize = float(image['imsize'][1:].split(',')[0])
        fov = (image['pixel_scale'] * imsize) / 3600. # deg
        ramin = image['obs_ra'] - fov
        ramax = image['obs_ra'] + fov
        decmin = image['obs_dec'] - fov
        decmax = image['obs_dec'] + fov
    except TypeError:
        # If imsize, pixel_scale, obs_ra, or obs_dec is 'None'
        ramin = None
        ramax = None
        decmin = None
        decmax = None
    return ramin, ramax, decmin, decmax


def catalog_extract(cursor, catalogs, *p):
    """Returns a dictionary containing catalogs and their sources
    which lie in the specified range of RA, Dec.""" 
    ramin, ramax, decmin, decmax = p
    catdict = {}
    for catalog in catalogs:
        query = 'SELECT ra, e_ra, dec, e_dec FROM %s WHERE ra BETWEEN '\
                '? AND ? AND dec BETWEEN ? AND ?' % catalog
        cursor.execute(query, (ramin, ramax, decmin, decmax))
        catdict[catalog] = cursor.fetchall()
    return catdict


def match(image, sources, catalogs):              
    """Performs the steps necessary to cross-match a 
    provided list of sources from the specified image
    to a given list of catalogs. The list of catalog
    names specifies which tables of the SkyCatalogs
    database to look in and in which order. The image
    is used to define a search range when extracting
    sources from the SkyCatalog database tables. The
    sources are cross-matched using the deRuiter
    radius criterion. A list of matched sources,
    non-matched sources, and the matched catalog
    sources is returned."""
    # Connect the sky catalog database
    skycat = os.path.join(catalogio.catalogdir, 'SkyCatalogs.sqlite')
    catcur = dbconnect(skycat)

    # Extract catalog sources
    catdict = catalog_extract(catcur, catalogs, search_range(image))

    # Define beam major axis FWHM:
    bmaj = image['bmaj'] / 3600. # deg

    # Cross-match each source to the catalog sources
    matches = []
    non_matches = []
    cat_matches = []
    num_match = 0
    for src in sources:
        idx = 0
        match, min_der, cat_match = crossmatch.deRuitermatch(
            src, catdict[catalogs[idx]], bmaj)
        while match != 1:
            idx += 1
            if idx < len(catalogs):
                match, min_der, cat_match = crossmatch.deRuitermatch(
                    src, catdict[catalogs[idx]], bmaj)
            else:
                # Ran out of catalogs to search -- no match found
                non_matches.append(src)
                break
        else:
            # Found a match!
            matches.append(src)
            cat_matches.append(cat_match)
            num_match += 1

    return matches, non_matches, cat_matches
    print('\nMatched {}/{} sources.'.format(num_match, len(sources)))


def main(catalogs=['NVSS'], **kwargs):
    """Prepares all possible inputs for catalog cross-matching.
    Input options are either a database of images and their sources
    (like the one written by database.py) or an ImageTable object
    and list of DetectedSource objects which are used to populate
    the database tables. If the input is the former, the database
    keyword must be set to True and the full path to the database 
    along with a list of image names for which sources are 
    to be cross-matched must be provided (i.e. database=True,
    file='/path/to/database.sqlite', images=['']). If the input is
    the latter, then the database keyword must be set to False and
    the ImageTable and DetectedSource objects must be provided
    (i.e. database=False, objects=(imgtbl, sources)). In both
    cases, the list of catalogs to check can be specified through
    the catalogs keyword with the default being only NVSS."""  
    # main(catalogs=catalogs, database=False, objects=(imgtbl, sources))
    # main(catalogs=catalogs, database=True, images=imglist,
    #      file='/data3/erichards/vlite/allvlite.sqlite')
    if kwargs['database'] is True:
        try:
            dbname = kwargs['file']
            if not os.path.exists(dbname):
                print('ERROR: File or path {} does not exist!'.format(dbname))
                sys.exit(0)
            else: pass
            srccur = dbconnect(dbname)
        except KeyError:
            print('If database is True, a file must be provided.')
        try:
            imglist = kwargs['images']
            # Select sources from one image at a time
            for img in imglist:
                srccur.execute('''SELECT id, bmaj, obs_ra, obs_dec, imsize, 
                    pixel_scale FROM Image WHERE name = ?''', (img, ))
                image = srccur.fetchone()
                srccur.execute('''SELECT ra, e_ra, dec, e_dec FROM Source 
                    WHERE image_id = ?''', (image['id'], ))
                sources = srccur.fetchall()
                matches, non_matches, cat_matches = match(image, sources,
                                                          catalogs)
        except KeyError:
            print('If database is True, a list of images must be provided.')
    else:
        try:
            image, sources = kwargs['objects']
            # Make sure we're dealing with dictionaries for consistency
            # across data input types
            try:
                image['bmaj']
                dimage = image
            except TypeError:
                dimage = image.__dict__
            dsources = []
            for src in sources:
                try:
                    src['ra']
                    dsources.append(src)
                except TypeError:
                    dsources.append(src.__dict__)
            matches, non_matches, cat_matches = match(dimage, dsources,
                                                      catalogs)
        except KeyError:
            print('If database is False, an ImageTable object and list of '
                  'DetectedSource objects must be provided.')


if __name__ == '__main__':
         main()


'''
#this program reads *.IPln1.fits and *.pybdsm.srl files
#INPUTS: 
# main pipeline processed dir (home dir with YYYY-MM subdirs)
# 1, 2, 3: YYYY, MM, DD
#

# for now, procdir is local directory containing output from runPyBDSM.py
procdir='/data3/erichards/testing/'
#procdir='/nfsshare/vpipe/processed/'
yr=int(sys.argv[1])
mon=int(sys.argv[2])
day=int(sys.argv[3])

#procdir='/extraid/vpipe/processed_test/'
#yr=2014
#mon=12
#day=3



#define startdir
if procdir[-1]=='/':
    procdir=procdir[:-1]
startdir='%s/%04d-%02d/%02d' % (procdir,yr,mon,day)
print startdir
if os.path.exists(startdir)==False:
    print 'ERROR! directory %s does not exist!' % startdir
    sys.exit(1)

#check to make sure the Images directory isn't empty
#(carryover from runPyBDSF not copying the images)
vimdir ='%s/Images' % startdir
vcatdir='%s/Catalogs' % startdir
#print vimdir
#print vcatdir
files=os.listdir(vimdir)
if not files:
    print 'Empty Images directory. Exiting.\n'
    sys.exit(1)

#make vslow_pybdsm subdir for output files
vslowdir=startdir+'/vslow_pybdsf'
print vslowdir
if os.path.isdir(vslowdir)==False:
    os.system('mkdir '+vslowdir)


#output files
logfile     ='%d-%02d-%02d_log.txt'             % (yr,mon,day)
sourcesfname='%d-%02d-%02d_sources_summary.txt' % (yr,mon,day)
imagesfname ='%d-%02d-%02d_images_summary.txt'  % (yr,mon,day)
matchesfname='%d-%02d-%02d_matches_summary.txt' % (yr,mon,day)
regfname    ='%d-%02d-%02d_matches.reg'         % (yr,mon,day)

globals.logf =open(logfile,'w') 
srcf =open(sourcesfname,'w')
imgf =open(imagesfname,'w')
matf =open(matchesfname,'w')
regf =open(regfname,'w')

line='#VLITE SLOW PyBDSF VERSION %s\n' % VERSION
srcf.write(line)
imgf.write(line)
matf.write(line)
regf.write(line)
globals.logf.write(line)

line='#Sources summary for %d-%02d-%02d\n' % (yr,mon,day)
srcf.write(line)
line='#srcID, islID, RA [deg], Dec [deg], src grade, image name,  src_sigma,  deltaRA [deg], deltaDec[deg],  angle from center [deg], int flux, int flux ERR [mJy], peak flux, peak flux ERR [mJy/bm], bmaj [arcsec], bmin [arcsec], bpa [deg]  bmaj/CLEANBMJ, bmin/CLEANBMN, image ACTNOISE [mJy/bm], image TAU_TIME [s], image RADIUS [deg], min deRuiter, name of catsrc min deRuiter, name of closest problem src, angle to closest probsrc [deg]\n'
srcf.write(line)

line='#Images summary for %d-%02d-%02d\n' % (yr,mon,day)
imgf.write(line)
line='#Name, RA [deg], Dec [deg],  NAXIS2, CDELT2 [deg/pix], Image Radius [deg],  BMAJ [arcsec], BMIN [arcsec], BPA [deg],  TAU_TIME [s], ACTNOISE [mJy/bm],  total # srcs, # center sources, # srcs in cutangle, expected # srcs in cutangle,  NVIS,  quality grade,  DURATION [s], MJDSTART, LST_START, ALT [deg], AZ [deg]\n'
imgf.write(line)

line='#Matched sources summary for %d-%02d-%02d\n' % (yr,mon,day)
matf.write(line)
line='#Source_id, Isl_id, RA [deg], E_RA, DEC [deg], E_DEC, Grade, Image, src_sigma, Total_flux [mJy], E_Total_flux, Peak_flux [mJy/bm], E_Peak_flux, RA_max, E_RA_max, DEC_max, E_DEC_max, Maj [arcsec], E_Maj, Min, E_Min, PA [deg], E_PA, Isl_Total_flux [mJy], E_Isl_Total_flux, Isl_rms [mJy/bm], Isl_mean, Resid_Isl_rms, Resid_Isl_mean, min_deRuiter, deRuiter_catsrc\n'
matf.write(line)

line='#Matched sources ds9 regions file for %d-%02d-%02d\n' % (yr,mon,day)
regf.write(line)
line='global color=blue font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n'
regf.write(line)
line='fk5\n'
regf.write(line)

#read sky catalogs
print 'reading tgss...'
tgss=iofuncs.readTGSS()
print 'reading nvss...'
#nvss=iofuncs.readNVSSerrs_SHORT()
nvss=iofuncs.readNVSSerrs()
print 'reading first...'
first=iofuncs.readFIRST()
print 'reading sumss...'
sumss=iofuncs.readSUMSS()
print 'reading wenss...'
wenss=iofuncs.readWENScomplete()
print 'reading gleam...'
gleam=iofuncs.readGLEAM()
#print 'reading VLITE catalog...'
#vlitecat=iofuncs.readVLITEcatalog()

#read fitted beam file
fitbeam=beamtools.ReadFittedBeamFile()

#define problem sources
psrc=iofuncs.defineproblemsources()


#loop over images in list "files"
ntmp=0
for jj in files:
    if jj[-10:] == 'IPln1.fits':
        print 'starting %s...' % jj
        nimages+=1
        fname=jj[:-5]
        vcatname=fname+'.pybdsm.srl'
        vcatfile=path.join(vcatdir,vcatname)

        #read VLITE CATALOG -- vltable 
        print 'reading source catalog...'
        VSOU=iofuncs.readPybdsfCatFile(vcatfile) #fluxes in mJy, angles in deg

        #read VLITE IMAGE
        print 'reading image...'
        VIMG=iofuncs.readVLITEimage(jj,vimdir)

        #how many VLITE sources total?
        VIMG.NSRCTOT=len(VSOU)
        print '  nsource= %d' % VIMG.NSRCTOT

        #how many VLITE sources w/in ANGLE_EXPECTED?
        #VIMG.NSRCCUT=transtools.countVLITEsources(VIMG,VSOU)
        for i in VSOU:
            i.alpha=(transtools.angdist(VIMG.CRVAL1,VIMG.CRVAL2,i.RA,i.Dec))
            if i.alpha < ANGLE_EXPECTED:
                VIMG.NSRCCUT+=1

        #calc expected number of srcs w/in ANGLE_EXPECTED
        VIMG.NSRCCUTEXPECT=beamtools.CalcNsrc(fitbeam,VIMG.NOISE,ANGLE_EXPECTED)

        #check image quality
        VIMG.QUALITY=transtools.checkVLITEimageQuality(VIMG,SRCMETRIC_CUT)


        
        #loop over sources:
        #  1) check for matches to sky catalogs & assign letter code to each source
        #  2) check for problem sources
        #  3) calc alpha angles & count number of center sources (should never be > 1, hopefully)
        #  4) calc phi angles of sources --> NOT ANYMORE; NO X,Y PIXEL VALUES IN PYBDSM CATALOGS
        #  5) write to files
        cutangle=0.95*VIMG.RADIUS   
        for i in VSOU:
            if VIMG.BMAJ>0.0:
                #count number of srcs at center
                if i.alpha <= (centercutangle*VIMG.BMAJ): #alpha in deg, VIMG.BMAJ in deg
                    VIMG.NSRCCENTER+=1
            #calc phi
#            tmpphi=math.atan2(i.X-VIMG.CRPIX1,i.Y-VIMG.CRPIX2) #rad
#            i.phi = (tmpphi*globals.RAD2DEG) #deg

            #is this near a problem source?
            probID,probangle,probflag=iofuncs.closestproblem(i.RA,i.Dec,psrc) #probflag=1 if near a problem src

            #check for matches to sky catalogs
            minder=999999.9
            minderid='dummy'
            #start w/ TGSS
            match,minder,minderid,matRA,matDec,matBMAJ,matBMIN,matPA=transtools.checkCatalogMatch(i,VIMG,tgss,DERUITER,minder,minderid)
            #if unmatched, check NVSS
            if match==0:
                match,minder,minderid,matRA,matDec,matBMAJ,matBMIN,matPA=transtools.checkCatalogMatch(i,VIMG,nvss,DERUITER,minder,minderid)
            #if unmatched, check FIRST
            if match==0:
                match,minder,minderid,matRA,matDec,matBMAJ,matBMIN,matPA=transtools.checkCatalogMatch(i,VIMG,first,DERUITER,minder,minderid)
            #if unmatched, check SUMSS
            if match==0:
                match,minder,minderid,matRA,matDec,matBMAJ,matBMIN,matPA=transtools.checkCatalogMatch(i,VIMG,sumss,DERUITER,minder,minderid)
            #if unmatched, check WENSS
            if match==0:
                match,minder,minderid,matRA,matDec,matBMAJ,matBMIN,matPA=transtools.checkCatalogMatch(i,VIMG,wenss,DERUITER,minder,minderid)
            #if unmatched, check VLITE catalog
            #if match==0:
                #match,minder,minderid,matRA,matDec,matBMAJ,matBMIN,matPA=transtools.checkCatalogMatch(i,VIMG,vlitecat,DERUITER,minder,minderid)
            #if unmatched, check GLEAM but do not categorized as "matched" (because GLEAM resolution poor)
            gleammatch=0
            if match==0:
                for j in gleam:
                    #first do a quickcheck:
                    if transtools.quickcheck(i.Dec,j.Dec,0.25)==1:
                        #now do deruiter calculation
                        der=transtools.deruiter(i.RA,i.Dec,i.dRA,i.dDec,j.RA,j.Dec,j.dRA,j.dDec)
                        if der<minder:
                            minder   = der
                            minderid = j.ID
                        if der<DERUITER:
                            #check angular distance too
                            if transtools.angdist(i.RA,i.Dec,j.RA,j.Dec)<(0.5*VIMG.BMAJ):
                                gleammatch=1
                                break   

            #assign src category
            srccode=VIMG.QUALITY
            if match>0:                   # matched
                srccode='M'+srccode
            elif i.alpha < cutangle:# unmatched & not near edge
                if minder>10.0:      
                    srccode='A'+srccode   #   not near a sky cat source
                else:
                    srccode='B'+srccode   #   near a sky cat source
            else:                         # unmatched & near edge
                if minder>10.0:
                    srccode='C'+srccode   #   not near a sky cat source
                else:
                    srccode='D'+srccode   #   near a sky cat source
            if probflag==1:
                srccode=srccode+'P'
            if gleammatch==1:
                srccode=srccode+'G'

            #write src info to file
            srcf.write('%s %s %f %f  %s  %s  %.1f  %.3e %.3e  %.3f  %.2f %.2f %.2f %.2f  %.3f %.3f %.3f %.2f %.2f  %.2f %.2f %.2f  %.2f %s  %s %.1f\n' % (i.srcID,i.islID,i.RA,i.Dec,srccode,VIMG.ID,i.sig,i.dRA,i.dDec,i.alpha,i.Total,i.dTotal,i.Peak,i.dPeak,i.BMAJ*3600.0,i.BMIN*3600.0,i.PA,i.BMAJ/VIMG.BMAJ,i.BMIN/VIMG.BMIN,VIMG.NOISE,VIMG.TINT,VIMG.RADIUS,minder,minderid,probID,probangle))

            #write matched src info to file
            if match>0:
                matf.write('%s %s %f %.3e %f %.3e %s %s %.1f %.2f %.2f %.2f %.2f %f %.3e %f %.3e %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %s\n' % (i.srcID,i.islID,i.RA,i.dRA,i.Dec,i.dDec,srccode,VIMG.ID,i.sig,i.Total,i.dTotal,i.Peak,i.dPeak,i.RAmax,i.dRAmax,i.Decmax,i.dDecmax,i.BMAJ*3600.0,i.dBMAJ*3600.0,i.BMIN*3600.0,i.dBMIN*3600.0,i.PA,i.dPA,i.islTotal,i.dislTotal,i.islRMS,i.islMEAN,i.RRMS,i.RMEAN,minder,minderid))
                regf.write('ellipse(%f,%f,%.2f",%.2f",%.1f) # text={%s}\n' % (matRA,matDec,matBMAJ*3600.0,matBMIN*3600.0,matPA+90.0,minderid))

        #END OF LOOP OVER SOURCES
        #
        #write image info to file
        imgf.write('%s %f %f  %d %e %.2f  %.3f %.3f %.2f  %.2f %.3f  %d %d %d %.3f  %d  %s  %.3f %f %s %.3f %.3f\n' % (VIMG.ID,VIMG.RA,VIMG.DEC,VIMG.NAXIS2,VIMG.CDELT2,VIMG.RADIUS,VIMG.BMAJ*3600.0,VIMG.BMIN*3600.0,VIMG.BPA,VIMG.TINT,VIMG.NOISE,VIMG.NSRCTOT,VIMG.NSRCCENTER,VIMG.NSRCCUT,VIMG.NSRCCUTEXPECT,VIMG.NVIS,VIMG.QUALITY,VIMG.DURATION,VIMG.MJDIMAGE,VIMG.LSTSTART,VIMG.ALT,VIMG.AZ))



t2=time.clock()
print 'time: %f [s]' % (round(t2-t1,3))
globals.logf.write('run time= %f [hrs]' % (round(t2-t1,3)/3600.0)) 

#close files
srcf.close()
matf.close()
regf.close()
imgf.close()
globals.logf.close()
#move files
#  sort sources file before moving
os.system('(head -n 3 '+sourcesfname+' && tail -n +4 '+sourcesfname+' | sort -k +5) > tmp.txt')
os.system('mv tmp.txt '+sourcesfname)
os.system('mv '+sourcesfname+' '+vslowdir)
os.system('mv '+matchesfname+' '+vslowdir)
os.system('mv '+regfname+' '+vslowdir)
os.system('mv '+imagesfname+' '+vslowdir)
os.system('mv '+logfile+' '+vslowdir)
sys.exit(0)
'''
