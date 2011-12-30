# gdal_copy_gt_srs.py - copy geotransform and spatial reference from a source
# raster to a target raster... target raster will be modified!

# Chris Toney (christoney at fs.fed.us)

###############################################################################

version = '0.1'

###############################################################################

import sys

try:
	from osgeo import gdal
	from osgeo.gdalconst import *
except ImportError:
	import gdal

import numpy

# =============================================================================
def Usage():
	print 'Usage: gdal_copy_gt_srs.py src_raster trg_raster'
	print
	print 'WARNING: target raster will be modified!'
	print
	sys.exit(1)

# =============================================================================

if __name__ == '__main__':

	if len(sys.argv) <> 3:
		Usage()

	src_file = sys.argv[1]
	trg_file = sys.argv[2]

	gdal.AllRegister()
    
	# get reference info from source
	gd_src = gdal.Open(src_file, GA_ReadOnly)
	if gd_src is None:
		print 'Could not open the source raster'
		sys.exit(1)
	#refXSize = gd_src.RasterXSize
	#refYSize = gd_src.RasterYSize
	refProj = gd_src.GetProjection()
	refGeotransform = gd_src.GetGeoTransform()
	if refGeotransform is not None:
		refOriginX = float(refGeotransform[0])
		refOriginY = float(refGeotransform[3])
		refPixelXSize = abs(float(refGeotransform[1]))
		refPixelYSize = abs(float(refGeotransform[5]))
	else:
		print 'Error getting geotransform'
		gd_src = None
		sys.exit(1)
        
	# open target raster
	gd_trg = gdal.Open(trg_file, GA_Update)
        
	gd_trg.SetGeoTransform(refGeotransform)
	gd_trg.SetProjection(refProj)

	gd_trg = None
	gd_src = None

	print 'Done'

