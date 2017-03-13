# SDSS_get
A program to get SDSS images for given coordinates and size.

This code downloads SDSS images, based on the coordinates of an object
 and the given size of field of view. User inputs the coordinates of an objectsds
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
   - bzip2   --> decompressing SDSS data
   - wget    --> downloading SDSS data
   - gzip    --> compressing the final images
   - csh     --> running the middle scripts
  
 * How it works:
   - Based on the given size and coordinates, it queries SDSS DR12 data-base for the
   parameters of the SDSS footprint images. This query returns all relevant run, camcol, and field parameters. Then it downloads all images from  
   SDSS DR12 repository: 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'
  
   - Using the software Montage to stack all the retrieved images and to finalize the desired image
 
   - The magnitude zero point of the final image is given in the header of the output fits files, where ZP = 16.40006562. This means that the physical unit of pixel values is mJy/pix.
  
   - All the output images are normalized by Exposure time, i.e. EXPTIME=1.
 
 * How to run:
 
           1$ python SDSS_get.py -r [RA] -d [DEC] -w [FOV-width]
   example  $ python SDSS_get.py -r 20.5604 -d 31.2547 -w 0.02
       -w 0.02   --> The final image would be a square, size: 0.02 degree on each side
     
      where RA and DEC are in degree. FOV-width is the
      width of the final image in degree. The final image is square shape.
    
           2$ python SDSS_get.py -r [RA] -d [DEC] -w [FOV-width] -b [SDSS-bands] -o [output_dir] -n [Object_name]
   example  $ python SDSS_get.py -r 20.5604 -d 31.2547 -w 0.02 -b ugr -n pgc5057732 -o /home/user/

       -b ugr    --> To only download u,g,r-bands
       -b r      --> Only downloads r-band
       not given --> Downloads all bands: ugriz

       -n <name> --> Sets the name of final files and folders based on the given name
       -o /home/user/ --> The final image(s) would be stored in the folder: /home/user/pgc5057732

  
 * Copyright 2017
 * Author: Ehsan Kourkchi <ehsan@ifa.hawaii.edu>sd
 * Feel free to distribute and modify for a better performance
