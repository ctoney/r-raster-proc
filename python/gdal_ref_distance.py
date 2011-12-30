# gdal_ref_distance.py - for each pixel in a grid of reference plot ids,
# calculates the distance to the actual reference plot location

# Chris Toney (christoney at fs.fed.us)

###############################################################################

version = '0.1'

###############################################################################

import sys
import csv
import math

try:
	from osgeo import gdal
	from osgeo.gdalconst import *
except ImportError:
	import gdal

import numpy

# =============================================================================
def Usage():
	print 'Usage: gdal_ref_distance.py -o out_filename -of out_format [-co NAME=VALUE]*'
	print '                     -n_src nodata_value_source -n_dst nodata_dest'
	print '                     -csv plot_idxy_csv -src plot_id_raster'
	print
	sys.exit(1)

# =============================================================================

# =============================================================================

def get_distance(gridid, x, y):
	# assumes refidxy dictionary is populated
	if gridid == n_src: # nodata value?
		return n_dst

	elif not (refidxy.has_key(gridid)):
		return n_dst
    
	else:
		refx = refidxy[gridid][0]
		refy = refidxy[gridid][1]
		dist = math.sqrt( (x-refx)*(x-refx) + (y-refy)*(y-refy) )

        return dist

# =============================================================================

if __name__ == '__main__':

	src_file = None
	format = None

	n_src = -9999
	n_dst = -9999
	create_options = []
	band_type = gdal.GDT_Int32

	src = None
	csv_file = None

	gdal.AllRegister()
	argv = gdal.GeneralCmdLineProcessor(sys.argv)
	if argv is None:
		sys.exit(0)

	# Parse command line arguments.
	i = 1
	while i < len(argv):
		arg = argv[i]

		if arg == '-o':
			i = i + 1
			out_file = argv[i]

		elif arg == '-n_src':
			i = i + 1
			n_src = float(argv[i])
            
		elif arg == '-n_dst':
			i = i + 1
			n_dst = float(argv[i])  

		elif arg == '-of':
			i = i + 1
			format = argv[i]

		elif arg == '-co':
			i = i + 1
			create_options.append(argv[i])

		elif arg == '-src':
			i = i + 1
			src = argv[i]

		elif arg == '-csv':
			i = i + 1
			csv_file = argv[i]
            
		elif arg[:1] == '-':
			print 'Unrecognised command option: ', arg
			Usage()
			sys.exit(1)

		else:
			print 'Unrecognised command option: ', arg
			Usage()
			sys.exit(1)
            
		i = i + 1


	if src is None:
		print 'Source raster file has not been specified.'
		Usage()
		sys.exit(1)

	if csv_file is None:
		print 'Plot idxy file has not been specified.'
		Usage()
		sys.exit(1)
        
	Driver = gdal.GetDriverByName(format)
	if Driver is None:
		print 'Format driver %s not found, pick a supported driver.' % format
		sys.exit(1)

	DriverMD = Driver.GetMetadata()
	if not DriverMD.has_key('DCAP_CREATE'):
		print 'Format driver %s does not support creation and piecewise writing.\nPlease select a format that does, such as GTiff (the default) or HFA (Erdas Imagine).' % format
		sys.exit(1)

	# open the plot idxy file (id,x,y) and populate a dictionary
	refidxy = {}
	pfile = open(csv_file, 'r')
	reader = csv.reader(pfile)
	for line in reader:
		refidxy[int(line[0])] = [float(line[1]),float(line[2])]
    
	# get reference info from source raster
	gd = gdal.Open(src, gdal.GA_ReadOnly)
	if gd is None:
		print 'Could not open the source raster'
		sys.exit(1)
	refXSize = gd.RasterXSize
	refYSize = gd.RasterYSize
	refProj = gd.GetProjection()
	refGeotransform = gd.GetGeoTransform()
	if refGeotransform is not None:
		refOriginX = float(refGeotransform[0])
		refOriginY = float(refGeotransform[3])
		refPixelXSize = abs(float(refGeotransform[1]))
		refPixelYSize = abs(float(refGeotransform[5]))
	else:
		print 'Error getting geotransform'
		gd = None
		sys.exit(1)

	# Create the output layer
	print 'Creating the output raster...'
	outdataset = Driver.Create(out_file, refXSize, refYSize, 1, band_type, create_options)
	if outdataset is None:
		print 'Creation of output raster failed, terminating gdal_mask.'
		sys.exit(1)
	outdataset.SetGeoTransform(refGeotransform)
	outdataset.SetProjection(refProj)

	# Line by line, read the id grid and write out distance to the reference plot location
	print 'Writing distances to the output raster...'
	src_band = gd.GetRasterBand(1)
	out_band = outdataset.GetRasterBand(1)
	# make a line of cell center x coordinates...
	x_line = []
	for n in range(refXSize):
		x_line.append( (refOriginX+(refPixelXSize/2)) + n*refPixelXSize )
	x_array = numpy.array(x_line, dtype=numpy.float32)
	x_array.shape = (1, refXSize)
	# vectorize the distance function...
	get_distance_vect = numpy.vectorize(get_distance)
	for y in xrange(refYSize):
		id_line = src_band.ReadAsArray(0, y, refXSize, 1)
		# make a line of cell center y coordinates...
		y_line = [( (refOriginY-(refPixelYSize/2)) - refPixelYSize*y )] * refXSize
		y_array = numpy.array(y_line, dtype=numpy.float32)
		y_array.shape = (1, refXSize)
		outline = get_distance_vect(id_line, x_array, y_array)
		out_band.WriteArray(outline, 0, y)

	print 'Raster output written to: %s' % out_file

	print 'Done'
