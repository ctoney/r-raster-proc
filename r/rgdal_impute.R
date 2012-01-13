# Impute reference plot ids into target pixels using yaImpute
# Chris Toney (christoney at fs.fed.us)

# Imputation of subsets of reference plots can be automatically constrained to 
# specified landscape strata including nodata regions for improved performance.
# i.e., groups of reference plots can be targeted to specific pixels in the 
# output raster.

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

# Uses global variables... functions assume certain objects are present.

# BEGIN configuration section

ncpus = 2
if (ncpus > 1) {run_parallel = TRUE} else {run_parallel = FALSE}

# if use_strata = TRUE, second column of training data should have strata ids
# strata ids should be 16-bit integers
# the first input raster in the raster lut should be the strata id grid
use_strata = TRUE
nlevels = 1

# format of training data is observation id, [strata_id], x1, x2, ...
# header has column names that match the variable names in raster lut
# id should be a 16-bit int:
#train_data_fn <- "/home/ctoney/work/rf/test/VModelMapData_nntest.csv"
train_data_fn <- "/home/ctoney/work/rf/test/z19_imp_ref_plots_v1.csv"

# if use_xy = TRUE then train_data_fn must include columns x and y
# x,y of pixel centers will be calculated automatically
use_xy = TRUE
write_xy_grids = FALSE # write out xy grids if use_xy = TRUE
xy_grid_dt = "Int32" # round to nearest meter, or use Float32 instead

# if slp_asp_transform = TRUE, aspect will be transformed to Cartesian
# coordinates. train_data_fn must contain columns asp (degrees from north)
# and slpp (slope percent).
slp_asp_transform = TRUE

# parameters for yai object in package yaImpute
yai_method = "euclidean"
yai_ann = TRUE # approximate nearest neighbor search, FALSE for brute force

# format of raster lut (no header): raster file path, var name, band num:
#raster_lut_fn <- "/home/ctoney/work/rf/test/VModelMapData_LUT.csv"
raster_lut_fn <- "/home/ctoney/work/rf/test/z19_raster_lut.csv"
#out_raster_fn <- "/home/ctoney/work/rf/test/nn_test_seq.img"
out_raster_fn <- "/home/ctoney/work/rf/test/z19_piece_test_imp_2.img"
out_raster_fmt <- "HFA"
out_raster_dt <- "Int16"
nodata_value <- 0

# END configuration section

library(snowfall)

sfInit(parallel=run_parallel, cpus=ncpus)

sfLibrary(yaImpute)
library(rgdal)

# read in training data and generate yai object(s) .. 1 yai for now
print("generating yai object(s)...")
df_tr <- read.csv(train_data_fn)
row.names(df_tr) <- df_tr[,1]

#if (use_xy) { ... TO DO: check whether columns x and y are in ref data

transform.slp_asp <- function(slpp, asp) {
	# convert slope/aspect to Cartesian coordinates.. Stage (1976) transformation
	# following the example in yaImpute doc...
	polar <- data.frame( slpp*.01, asp*(pi/180) )
	cartesian <- t(apply(polar,1,function (x) {return (c(x[1]*cos(x[2]),x[1]*sin(x[2]))) }))
	return(cartesian)
}

if (slp_asp_transform) {
	print("transforming slp/asp to Cartesian...")
	cartesian <- transform.slp_asp(df_tr$slpp, df_tr$asp)
	df_tr$slp_asp_x <- cartesian[,1]
	df_tr$slp_asp_y <- cartesian[,2]
	df_tr$asp <- NULL
}

if (use_strata) {
	# columns 2 through <nlevels+1> of df_tr have the strata ids

	# a yai object for all ref plots
	yai.allrefs <- yai(x=df_tr[,-1:-(nlevels+1)], noTrgs=TRUE, noRefs=TRUE, method=yai_method, ann=yai_ann)
	#print(yai.allrefs)

	# yai objects for strata subsets
	yai.strata <- list()
	yai.strata.ids <- list()
	idref.strata <- list()
	idref.strata.ids <- list()
	for (n in 1:nlevels) {
		subsets.df_tr <- split(df_tr, df_tr[,n+1])
		strata_ids <- unique(df_tr[,n+1])
		yai.strata[[n]] <- list()
		yai.strata.ids[[n]] <- vector(mode="integer", 0)
		idref.strata[[n]] <- list()
		idref.strata.ids[[n]] <- vector(mode="integer", 0)
		for (id in strata_ids) {
			idx <- which( names(subsets.df_tr) == as.character(id) )
			if ( length(subsets.df_tr[[idx]][,1]) > 1 ) {
				# a yai object for the subset
				yai.strata[[n]][[id]] <- yai(x=subsets.df_tr[[idx]][,-1:-(nlevels+1)], noTrgs=TRUE, noRefs=TRUE, method=yai_method, ann=yai_ann)
				yai.strata.ids[[n]] <- append(yai.strata.ids[[n]], id)
			} else {
				# the reference plot id
				idref.strata[[n]][[id]] <- subsets.df_tr[[idx]][1,1]
				idref.strata.ids[[n]] <- append(idref.strata.ids[[n]], id)
			}
		}
	}
	#print(yai.strata)

} else {
	yai.allrefs <- yai(x=df_tr[,-1], noTrgs=TRUE, noRefs=TRUE, method=yai_method, ann=yai_ann)
	print(yai.allrefs)
}

# populate a table of raster names and corresponding variable names
raster_lut <- read.table(file=raster_lut_fn,sep=",",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)

# assuming all input rasters have same extent and cell size
# TO DO: check extents of input rasters

# open first raster and read extents
sp.rast <- open.SpatialGDAL(raster_lut[1,1], silent=TRUE)
ncols <- sp.rast@grid@cells.dim[1]
nrows <- sp.rast@grid@cells.dim[2]
#print(sp.rast@bbox)
xmin <- sp.rast@bbox[1,1]
ymax <- sp.rast@bbox[2,2]
cellsize <- sp.rast@grid@cellsize[1]
close(sp.rast)
if (use_xy) {
	print("using x,y...")
	cell_centers_x <- seq(from=xmin+(cellsize/2), by=cellsize, length.out=ncols)
}

# open a list of GDAL datasets for the rasters
gd_list <- list()
for (r in 1:length(raster_lut[,1])) {
	print(raster_lut[r,1])
	gd_list[raster_lut[r,2]] <- GDAL.open(raster_lut[r,1], read.only=TRUE, silent=TRUE)
	
}

# not using - this sets up gd_list on the cluster for parallel input
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
	# scanline is a row number to process
	# assumes gd_list is populated
	# returns a dataframe with the input rows in named columns

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
		df$y <- rep( ( ymax-(cellsize/2) - (cellsize*scanline) ), ncols)
	}

	return(df)
}

# a wrapper function for a predict method.. yaImpute in this case
predict.wrapper <- function(df) {
# assumes the yai object list has been exported to the cluster

	if (slp_asp_transform) {
		cartesian <- transform.slp_asp(df$slpp, df$asp)
		df$slp_asp_x <- cartesian[,1]
		df$slp_asp_y <- cartesian[,2]
		df$asp <- NULL
	}

	if (use_strata) {
		# values of the first <nlevels> rasters (df columns 1-nlevels) are the strata ids

		df$neiIdsTrgs <- as.vector(rep(NA, length(df[,1])), mode="character")

		# assign the nodata pixels
		df[df[,1]==nodata_value, "neiIdsTrgs"] <- nodata_value

		# strata level <nlevels> is most specific, strata level 1 is most general
		# try to impute from the most specific level if there are plots, then step 
		# back to the more general stratification levels...
		for (n in seq(nlevels, 1, -1)) {
			# strata ids having only one ref plot... assign that plot id
			ids <- unique(df[df[,n] %in% idref.strata.ids[[n]] == TRUE, n])
			for (id in ids) {
				df[df[,n]==id & is.na(df[, "neiIdsTrgs"]), "neiIdsTrgs"] <- idref.strata[[n]][[id]]
			}

			# strata ids having yai objects for nearest neighbor imputation
			ids <- unique(df[df[,n] %in% yai.strata.ids[[n]] == TRUE, n])
			for (id in ids) {
				try(df[df[,n]==id & is.na(df[, "neiIdsTrgs"]), "neiIdsTrgs"] <- newtargets(yai.strata[[n]][[id]], df[df[,n]==id & is.na(df[, "neiIdsTrgs"]), ])$neiIdsTrgs, silent=TRUE)
			}
		}

		# strata ids that are absent in the reference data... impute from all plots
		if ( any(is.na(df$neiIdsTrgs)) ) {
			df[is.na(df[, "neiIdsTrgs"]), "neiIdsTrgs"] <- newtargets(yai.allrefs, df[is.na(df[, "neiIdsTrgs"]), ])$neiIdsTrgs
		}

		return(df$neiIdsTrgs)

	} else {

		#yai.nnref <- newtargets(yai.allrefs, df)
		#return(yai.nnref$neiIdsTrgs)
		return(newtargets(yai.allrefs, df)$neiIdsTrgs)
	}
}

# a function to calculate row values
process_row <- function(scanline) {
# impute across the row vectors
# needs yai object available on the cluster
# assumes dst GDAL dataset is available for writing output

	df_in <- read_input_row(scanline)

	if (sfParallel()) {
		# split each row:
		# parallelize in n pieces, n = number of cpus
		# is this faster in parallel?...
		df_in$piece <- rep(1:ncpus, length.out=ncols)
		pieces <- split(df_in, df_in$piece)
		pieces_out <- sfClusterApplyLB(pieces, predict.wrapper)
		outline <- as.numeric(unsplit(pieces_out, df_in$piece))

	} else {
		# sequential execution:
		outline <- as.numeric(predict.wrapper(df_in))
	}

	# write a row of output to the transient dataset
	a <- array(unlist(outline), dim=c(ncols, 1))
	a[is.na(a)] <- nodata_value
	putRasterData(dst, a, band=1, offset=c(scanline,0))

	if (write_xy_grids) {
		# write xy rows
		a <- array(unlist(df_in$x), dim=c(ncols, 1))
		putRasterData(dst_x, a, band=1, offset=c(scanline,0))
		a <- array(unlist(df_in$y), dim=c(ncols, 1))
		putRasterData(dst_y, a, band=1, offset=c(scanline,0))
	}

	return()
}


sfExportAll()

# a transient dataset to hold the output
dst <- new('GDALTransientDataset', driver=new('GDALDriver', out_raster_fmt), rows=nrows, cols=ncols, bands=1, type=out_raster_dt)
if (write_xy_grids) {
	dst_x <- new('GDALTransientDataset', driver=new('GDALDriver', out_raster_fmt), rows=nrows, cols=ncols, bands=1, type=xy_grid_dt)
	dst_y <- new('GDALTransientDataset', driver=new('GDALDriver', out_raster_fmt), rows=nrows, cols=ncols, bands=1, type=xy_grid_dt)
}

# process by rows:
print("processing...")
t <- system.time( lapply(0:(nrows-1), process_row) )
print(t)

print("cluster processing done...")

# save the output dataset(s)
saveDataset(dst, out_raster_fn)
GDAL.close(dst)
if (write_xy_grids) {
	x_fn <- sub( basename(out_raster_fn), sub(".", "_x.", basename(out_raster_fn), fixed=TRUE), out_raster_fn, fixed=TRUE )
	y_fn <- sub( basename(out_raster_fn), sub(".", "_y.", basename(out_raster_fn), fixed=TRUE), out_raster_fn, fixed=TRUE )
	saveDataset(dst_x, x_fn)
	saveDataset(dst_y, y_fn)
	GDAL.close(dst_x)
	GDAL.close(dst_y)
}

# close the input datasets
lapply(gd_list, GDAL.close)

sfStop()

print("done.")

