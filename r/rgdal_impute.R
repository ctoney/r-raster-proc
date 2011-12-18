# impute reference plot ids into target pixels using yaImpute
# Chris Toney (christoney at fs.fed.us)

# imputation of reference plots can be automatically constrained to specific
# landscape strata

# raster IO using rgdal
# parallel using snowfall

#******************************************************************************
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#******************************************************************************

# uses global variables, functions assume certain variables are present

run_parallel = FALSE
ncpus = 2

# format of training data is observation id, x1, x2, ...
# header has column names that match the variable names in raster lut
# id should be a 16-bit int:
train_data_fn <- "/home/ctoney/work/rf/test/VModelMapData_nntest.csv"

# if use_xy = TRUE then train_data_fn must include columns x and y
# x,y of pixel centers will be calulated automatically
use_xy = TRUE

# if slp_asp_transform = TRUE, aspect will be transformed to cartesian coordinates
# train_data_fn must contain columns asp (degrees from north) and slp (slope percent)
slp_asp_transform = FALSE

# method for yai object in package yaImpute
yai_method = "mahalanobis"

# format of raster lut (no header): raster file path, var name, band num:
raster_lut_fn <- "/home/ctoney/work/rf/test/VModelMapData_LUT.csv"
out_raster_fn <- "/home/ctoney/work/rf/test/nn_test.img"
out_raster_fmt <- "HFA"
out_raster_dt <- "Int16"
nodata_value <- -9999

library(snowfall)

sfInit(parallel=run_parallel, cpus=ncpus)

sfLibrary(yaImpute)
library(rgdal)

# read in training data and generate yai object(s) .. 1 yai for now
print("generating yai object(s)...")
df_tr <- read.csv(train_data_fn)
row.names(df_tr) <- df[,1]
if (slp_asp_transform) {
	# convert slope/aspect to cartesian coordinates.. Stage (1976) transformation
	# following the example in yaIpmute doc...
	polar <- data.frame( c(df_tr$slp*.01, df_tr$asp*(pi/180)) )
	cartesian <- t(apply(polar,1,function (x) {return (c(x[1]*cos(x[2]),x[1]*sin(x[2]))) }))
	df_tr$slp_asp_x <- cartesian[,1]
	df_tr$slp_asp_y <- cartesian[,2]
	df_tr$asp <- NULL
}
yai_list <- list()
yai_list[[1]] <- yai(x=df_tr[,-1], noTrgs=TRUE, noRefs=TRUE, method=yai_method)
print(yai_list[[1]])

raster_lut <- read.table(file=raster_lut_fn,sep=",",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)

# assuming all input rasters have same extent and cell size
# TO DO: check extents of input rasters

# open first raster and read extents
sp.rast <- open.SpatialGDAL(raster_lut[1,1], silent=TRUE)
ncols <- sp.rast@grid@cells.dim[1]
nrows <- sp.rast@grid@cells.dim[2]
#print(sp.rast@bbox)
xmin <- sp.rast@bbox[1,1]
ymax <- sp.rast@bbox[2,1]
cellsize <- sp.rast@grid@cellsize[1]
close(sp.rast)
if (use_xy) {
	cell_centers_x <- seq(from=xmin+(cellsize/2), by=cellsize, length.out=ncols)
}

# build a list of GDAL datasets for the rasters
gd_list <- list()
for (r in 1:length(raster_lut[,1])) {
	print(raster_lut[r,1])
	gd_list[raster_lut[r,2]] <- GDAL.open(raster_lut[r,1], read.only=TRUE, silent=TRUE)
	
}

# not using - this code sets up gd_list on the cluster for parallel input
#gd_list <- list()
#r <- 1
#sfExportAll() #("raster_lut", "gd_list")

#for (n in 1:length(raster_lut[,1])) {
#	#print(raster_lut[r,1])
#	sfClusterEval( gd_list[raster_lut[r,2]] <- GDAL.open(raster_lut[r,1], read.only=TRUE, silent=TRUE) )
#	sfClusterEval( r <- r + 1 )
#}


# a function to read a row from the input raster stack
read_input_row <- function(scanline) {
	# scanline is a scanline (row) number to process
	# returns a dataframe with the input rows in named colums

	ncols <- dim(gd_list[[1]])[2]

	df <- data.frame(matrix(-9999, ncols, length(raster_lut[,2])))
	names(df) <- raster_lut[,2]
	# need unique row names that are different from refs, for yai newtargets...
	row_names <- rep("t", ncols)
	row.names(df) <- make.unique(row_names, sep="")

	for (r in 1:length(gd_list)) {
		src <- gd_list[[r]]
		b <- raster_lut[r,3]
		df[raster_lut[r,2]] <- getRasterData(src, band=b, offset=c(scanline,0), region.dim=c(1,ncols))
	}

	if (use_xy) {
		df$x <- cell_centers_x
		df$y <- rep( (ymax-(cellsize/2) - (cellsize*scanline)), ncols)
	}

	if (slp_asp_transform) {
		# convert slope/aspect to cartesian coordinates.. Stage (1976) transformation
		# following the example in yaIpmute doc...
		polar <- data.frame( c(df$slp*.01, df$asp*(pi/180)) )
		cartesian <- t(apply(polar,1,function (x) {return (c(x[1]*cos(x[2]),x[1]*sin(x[2]))) }))
		df$slp_asp_x <- cartesian[,1]
		df$slp_asp_y <- cartesian[,2]
		df$asp <- NULL
	}

	return(df)
}

# a wraper function for a predict method
predict.wrapper <- function(df) {
# assumes the yai object list has been exported to the cluster
	yai.nnref <- newtargets(yai_list[[1]], df)
	return(yai.nnref$neiIdsTrgs)
}

# a function to calculate row values
# needs yai object available on the cluster
process_row <- function(scanline) {
	# impute across the input vectors

	df <- read_input_row(scanline)

	if (sfParallel()) {

		# chunk each row:
		# parallelize in n chunks, n = number of cpus
		# the chunking adds overhead, is this faster in parallel?...
		df$chunk <- rep(1:ncpus, length.out=ncols)
		chunks <- split(df, df$chunk)
		chunks_out <- sfClusterApplyLB(chunks, predict.wrapper)
		df_out <- unsplit(chunks_out, df$chunk)
		outline <- as.numeric(df_out[,1])

	}
	else {
		# sequential execution:
		yai.nnref <- newtargets(yai_list[[1]], df)
		outline <- as.numeric(yai.nnref$neiIdsTrgs[,1])
	}

	return(outline)
}


sfExportAll()
print("processing...")

# process by rows:
t <- system.time( output <- lapply(0:(nrows-1), process_row) )
print(t)

print("cluster processing done...")

# close the GDAL datasets
lapply(gd_list, GDAL.close)

sfStop()

# **output has a 1-based index** (index-1 to get raster row numbers)

# write the ouput
print("writing output raster...")
dst <- new('GDALTransientDataset', driver=new('GDALDriver', out_raster_fmt), rows=nrows, cols=ncols, bands=1, type=out_raster_dt)
for (y in 0:(nrows-1)) {
	a <- array(unlist(output[[(y+1)]]), dim=c(ncols, 1))
	a[is.na(a)] <- nodata_value
	putRasterData(dst, a, band=1, offset=c(y,0))
}
saveDataset(dst, out_raster_fn)
GDAL.close(dst)

rm(output)

print("done.")

