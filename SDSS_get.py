#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2017"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "v1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"

#################################################
## Importing Important Python Libraries
#################################################
import sys
import os
import subprocess
import numpy as np
import pyfits
from optparse import OptionParser
import sqlcl    # tools to querry SDSS DR12
#################################################


def arg_parser():
    parser = OptionParser(usage="""\
\n
 This code downloads SDSS images, based on the coordinates of an object
 and the given size of the field of view. User inputs the coordinates of an object,
 and the desired size of the image, and/or the SDSS bands. The final image has a
 square shape.
 
 * Dependencies:
 
 1) Python packages:
   - pyfits
   - numpy
   - optparse
   - sqlcl  (provided)
 
 2) Montage software (An astronomical Image Mosaic Engine)
   Download: http://montage.ipac.caltech.edu/
   - mImgtbl
   - mProjExec
   - mImgtbl
   - mAdd
   - mConvert
 
 3) Linux tools
   - bzip2
   - wget
   - gzip
  
 * How it works:
   - Based on the given size and coordinates, it queries SDSS DR12 data-base for the
   parameters of the SDSS footprint images. This query returns all relevant run, camcol, and field parameters. Then it downloads all images from  
   SDSS DR12 repository: 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'
  
   - Using the software Montage to stack all the retrieved images and to finalize the desired image
 
   - The magnitude zero point of the final image is given in the header of the output fits files, where ZP = 16.40006562. This means that the physical unit of pixel values is mJy/pix.
  
   - All the output images are normalized by Exposure time, i.e. EXPTIME=1.
 
 * How to run:
 
           1$ %prog -r [RA] -d [DEC] -w [FOV-width]
   example  $ %prog -r 20.5604 -d 31.2547 -w 0.02
       -w 0.02   --> The final image would be a square, size: 0.02 degree on each side
     
      where RA and DEC are in degree. FOV-width is the
      width of the final image in degree. The final image is square shape.
    
           2$ %prog -r [RA] -d [DEC] -w [FOV-width] -b [SDSS-bands] -o [output_dir] -n [Object_name]
   example  $ %prog -r 20.5604 -d 31.2547 -w 0.02 -b ugr -n pgc5057732 -o /home/user/

       -b ugr    --> To only download u,g,r-bands
       -b r      --> Only downloads r-band
       not given --> Downloads all bands: ugriz

       -n <name> --> Sets the name of final files and folders based on the given name
       -o /home/user/ --> The final image(s) would be stored in the folder: /home/user/pgc5057732

  
 * Copyright 2017
 * Author: Ehsan Kourkchi <ehsan@ifa.hawaii.edu>
 * Feel free to distribute and modify for a better performance
  
""")
    
    parser.add_option('-r', '--ra',
                      type='float', action='store',
                      help="""Right Ascension of the object center [degree]""")
    parser.add_option('-d', '--dec',
                      type='float', action='store',
                      help="""Declination of the object center [degree]""")    
    parser.add_option("-w", "--width",
                      type='float', action='store',
                      help="""The width of fied of view [degree]""") 
    parser.add_option("-n", "--name",
                      type='string', action='store',
                      help="""object name (optional)""")   
    parser.add_option("-o", "--outfolder",
                      type='string', action='store',
                      help="""output folder (optional)""")    
    parser.add_option("-b", "--bands",
                      type='string', action='store',
                      help="""SDSS bands (optional)""")   
   
    (opts, args) = parser.parse_args()
    return opts, args

#################################################

def xcmd(cmd,verbose):

  if verbose: print '\n'+cmd

  tmp=os.popen(cmd)
  output=''
  for x in tmp: output+=x
  if 'abort' in output:
    failure=True
  else:
    failure=tmp.close()
  if False:
    print 'execution of %s failed' % cmd
    print 'error is as follows',output
    sys.exit()
  else:
    return output

#################################################

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        
#################################################
def getSDSSfields(ra, dec, size):  # all in degree
  
  
  delta = 0.5*size+0.13
  ra_max  =  ra+1.5*delta
  ra_min  =  ra-1.5*delta
  dec_max =  dec+delta
  dec_min =  dec-delta
  
  querry = """
  
 SELECT
 fieldID,
 run, 
 camCol, 
 field,
 ra, 
 dec,
 run,
 rerun 
 FROM Field   
  """
  
  querry += "WHERE ra BETWEEN "+str(ra_min)+" and "+str(ra_max)+" and dec BETWEEN "+str(dec_min)+" and "+str(dec_max)
  
  lines = sqlcl.query(querry).readlines()
  N = len(lines)

  
  field_lst = []
  for i in np.arange(2,N):
      line = lines[i]
      line = line.split(',')
      run    = line[1]
      camcol = line[2]
      field  = line[3]
      ra_    = line[4]
      dec_   = line[5]
      field_lst.append([run, camcol, field])
  
  return field_lst

#################################################
## fieldInfo = [run, camcol, field]
## band = 'u', 'g', 'r', 'i', 'z'
## folde: <where to store files>

def getSDSSfiles(fieldInfo, band, folder):   
  
  run    = fieldInfo[0]
  camcol = fieldInfo[1]
  field  = fieldInfo[2]
  
  fileName = 'frame-'+band + '-''{0:06d}'.format(int(run))+'-'+camcol+'-'+'{0:04d}'.format(int(field))+'.fits.bz2'
  filename = 'frame-'+band + '-''{0:06d}'.format(int(run))+'-'+camcol+'-'+'{0:04d}'.format(int(field))+'.fits'
  http = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'
  
  http += run + '/'
  http += camcol + '/'
  http += fileName
  
  with cd(folder):
    xcmd("wget " + http, True)
    xcmd("bzip2 -d " + fileName, True)
    fits_prep(filename, edit_zp=True)

#################################################
# removing the siffix of a filename
def removefix(filename):
  
  
  name_list = filename.split('.')
  N = len(name_list)
  if name_list[N-1] == 'fits':
    
    name =''
    for i in range(N-2):
      name = name + name_list[i] + '.'
    name += name_list[N-2]

  return name
#################################################
# preparing fits files before stacking
# Editing the zero points 
def fits_prep(filename, edit_zp=True):
  
  name = removefix(filename)
  xcmd("cp " + filename + ' ' + name+'.tmp.fits', True) 
  
  hdulist = pyfits.open(filename)
  prihdr = hdulist[1].header    
      
  if edit_zp:
     
     hdulist = pyfits.open(name+'.tmp.fits', mode='update')
     prihdr = hdulist[0].header
     prihdr['BSCALE'] = 1.
     prihdr['BZERO'] = 0.

     expTime = float(prihdr['EXPTIME'])
     nMgy = float(prihdr['NMGY'])

     new_zp = 16.40006562  # mJy / pix
     
     ### References from SDSS official site
     ### https://data.sdss.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
     ### http://www.sdss.org/dr12/imaging/images/
     zeroPoint = 22.5
     alpha = 10**(-0.4*(zeroPoint-new_zp))
     
     #hdulist[0].data = data = hdulist[0].data * (3631E-6 )
     hdulist[0].data = data = hdulist[0].data * (alpha )

     prihdr['EXPTIME'] = 1.0
     prihdr['ZP'] = new_zp

     hdulist.flush() 

     cmd = 'mv '+name+'.tmp.fits ' + filename
     xcmd(cmd, True)    
     
#################################################
# making mosaics
def make_sdss(objName, ra, dec, size, band):
    
    band_root ='.'
    cmd = 'mkdir '+band_root+'/tmp'
    xcmd(cmd, True)
    
    cmd = 'mv '+band_root+'/frame-'+band+'*.fits '+band_root+'/tmp/.'
    xcmd(cmd, True)
    
#################################################
    with cd(band_root):
        
      # How to get Montage to make the mosaics
      # preparing csh file 
      fout=open('montage.csh','w')
      montage = '''
      rm *tbl
      rm -rf projected
      mkdir  projected
      mImgtbl tmp rimages.tbl
      mProjExec -p tmp rimages.tbl tmp.hdr projected stats.tbl
      mImgtbl projected pimages.tbl
      mAdd -p projected pimages.tbl tmp.hdr tmp.fits
      mConvert -b -32 tmp.fits tmp.fits
      rm *area*fits
      rm *tbl
      rm -rf projected
      '''
      fout.write(montage+'\n')
      fout.close()
      
      # The header of the template fits file
      naxis = int(size*3600 / 0.396)
      crpix = naxis/2
      hdr = ''
      hdr += 'SIMPLE = T' + '\n'
      hdr += 'BITPIX = -32' + '\n'
      hdr += 'NAXIS = 2' + '\n'
      hdr += 'NAXIS1 = '+str(naxis) + '\n'
      hdr += 'NAXIS2 = '+str(naxis) + '\n'
      hdr += 'CTYPE1 = \'RA---TAN\'' + '\n'
      hdr += 'CTYPE2 = \'DEC--TAN\'' + '\n'
      hdr += 'CRVAL1 = ' +str(ra)+ '\n'
      hdr += 'CRVAL2 = ' +str(dec)+ '\n'
      hdr += 'CRPIX1 = '+str(crpix) + '\n'
      hdr += 'CRPIX2 = '+str(crpix) + '\n'
      hdr += 'CDELT1 = -1.1E-04' + '\n'
      hdr += 'CDELT2 = 1.1E-04' + '\n'
      hdr += 'CROTA2 = 0.000000' + '\n'
      hdr += 'EQUINOX = 2000.0' + '\n'
      hdr += 'BSCALE = 1' + '\n'
      hdr += 'BZERO = 0' + '\n'
      hdr += 'EXPTIME = 1.0' + '\n'
      hdr += 'ZP = 16.40006562 ' + '\n'
      hdr += 'OBJECT = ' + objName + '\n'
      hdr += 'HISTORY = \'email: <ehsan@ifa.hawaii.edu>\'' + '\n'  
      hdr += 'HISTORY = \'by: SDSS_get.py\'' + '\n'  
      hdr += 'HISTORY = \'count unit = [mJy/pixel]\'' + '\n'
      hdr += 'END' + '\n'

      fout=open('tmp.hdr','w')
      fout.write(hdr)
      fout.close()
      
      # executing all montage commands
      xcmd('csh montage.csh', True)
      
      # dealing with files and folders
      xcmd('mv tmp.fits '+objName+'_'+band+'.fits', True)
      xcmd('rm montage.csh', True)
      xcmd('rm tmp.hdr', True)
      xcmd('rm -rf tmp', True)
      
#################################################
# Downloading SDSS images and then making mosaics
def sdssget(objName, ra, dec, size, repository): 
    
    
  fields = getSDSSfields(ra, dec, size)
  
  folder = repository+'tmp'
  if os.path.isdir(folder):
      xcmd('rm -rf '+folder, True)
  xcmd('mkdir '+folder, True)
  
  if not os.path.isdir(repository+objName):
      xcmd('mkdir '+repository+objName, True)
  
  for band in ['u','g','r','i','z']:
      
      for i in range(len(fields)):
        print ' '  
        print ' # '+str(i+1)+'/'+str(len(fields))+'  ... obj: '+objName+'  band: '+band
        
        getSDSSfiles(fields[i], band, folder)
      
      with cd(folder):
        make_sdss(objName, ra, dec, size, band)       
      
      xcmd('mv '+folder+'/*fits '+repository+objName, True)
      xcmd('gzip '+repository+objName+'/*fits', True)

  xcmd('rm -rf '+folder, True)
#################################################

  
if __name__ == '__main__':

#################################################
  warning = False
  
  if (len(sys.argv) < 2): 
    print "\nNot enough input arguments ..."
    print >> sys.stderr, "Use -h option for help ... \n"
    exit(1)
  
  
  opts, args =  arg_parser()
  print "\n------------------------------------"
  print "Input Arguments (provided by User)"
  print "------------------------------------"
  print "R.A.       [deg]:", opts.ra
  print "Dec.       [deg]:", opts.dec
  print "FOV-width  [deg]:", opts.width
  print "Object Name  (optional):", opts.name
  print "Output Foldr (optional):", opts.outfolder
  print "SDSS Bands   (optional):", opts.bands
  print "------------------------------------"
  
  ra  = opts.ra 
  dec = opts.dec
  width = opts.width
  
  if ra == None or dec == None or width == None:
    print "\nWrong input parameters ..."
    print >> sys.stderr, "Use -h option for more information ... \n"
    exit(1)      


  bands = []

  if opts.bands != None:
      for b in opts.bands:
        if b in ['u', 'g', 'r', 'i', 'z'] and not b in bands:
            bands.append(b)
  if len(bands) == 0:
     print "\n[Warning] No output SDSS band is given ..."
     print "[Warning] Downloading all bands (ugriz) ..."
     bands = ['u', 'g', 'r', 'i', 'z']
     warning = True
  
  
  if  opts.outfolder == None:
      print "\n[Warning] No output folder were provided ..."
      print "[Warning] Using the current directory to store results ..."
      warning = True
      repository = '.'
  else:
      repository = opts.outfolder
      if not os.path.isdir(repository):
        print "\n[Error] no such directory: " + repository
        print >> sys.stderr, "Use -h option for help ... \n"
        exit(1)
         
      
      
  if  opts.name == None:
      print "\n[Warning] No object name were given ..."
      warning = True
      ra_st = '{:.4f}'.format(ra)
      dec_st = '{:.4f}'.format(dec)
      if dec < 0:
         objName = 'Obj'+ra_st+'-'+dec_st
      else: 
         objName = 'Obj'+ra_st+'+'+dec_st
      print "Using this name for the object: "+objName
  else:
      objName = opts.name  
   
  
  if warning:
      
    print "\n------------------------------------"
    print "[Warning] Use the following Parameters ?"
    print "------------------------------------"
    print "R.A.       [deg]:", ra
    print "Dec.       [deg]:", dec
    print "FOV-width  [deg]:", width
    print "Object Name     :", objName
    if repository == '.':
       print "Output Foldr    : <current folder>"
    else: print "Output Foldr    :", repository
    print "SDSS Bands    :", bands
    print "------------------------------------"  
      
    input_var = raw_input("Continue (yes,no)?: ")
    if input_var.upper() == 'YES': 
       warning = False
  
  if not warning:  
      sdssget(objName, ra, dec, width, repository+'/')
  else: 
      print "\nTry again ..."
      print >> sys.stderr, "Use the '-h' option for more information ... \n"


