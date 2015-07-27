
#' Save a feature set
#'
#' Save the components of a feature set for later retrieval.
#' @param obj.to.save Feature set to be saved
#' @param name Name of files to be used to save the feature set. Should not include file extension, since it will be used to save multiple files
#' @export
#' @examples
#' save.features(genlis, "Genlis")
#' 
#' 
save.features <- function(obj.to.save, name) {
    writeRaster(obj.to.save$features, filename = paste(name, "-features.grd", sep = ""), overwrite = T)
    write.csv(obj.to.save$feature.types, file = paste(name, "-feature-types.csv", sep = ""), row.names = F)
}


#' Load a feature set
#'
#' Load previously-saved components of a feature set and assign to an object.
#' @param name Name of files that was to save the feature set. Should not include file extension, since it will be used to refer to multiple files
#' @return A new feature set: a list containing a raster of numbered features, and a matrix assigning a type to each feature.
#' @export
#' @examples
#' genlis <- load.features("Genlis")
#' 
#' 
load.features <- function(name) {
    feature.types <- read.csv(file = paste(name, "-feature-types.csv", sep = ""))
    features <- raster(paste(name, "-features.grd", sep = ""))
    
    list(features = features, feature.types = feature.types)
}


#' Plot features & feature types
#'
#' Plot all features in set, coloured according to their assigned type.
#' @param feature.set The feature set to be plotted
#' @export
#' @examples
#' show.features(genlis)
#' 
#' 
show.features <- function(feature.set){
    
    cols <- c("grey", "black", "red3", "red", "gold", "lightseagreen", "purple")
    desc <- c("Unclassified", "Post-hole", "Scale", "N-S axis", "Annotation", "Large feature", "Changed")
    
    ind <- sort(unique(feature.set$feature.types[,2])) + 1
    l.ind <- ind[(ind != 3) & (ind != 4)]
    
    plot(reclassify(feature.set$features, rcl = feature.set$feature.types), col = cols[ind], asp = T, legend = F, frame.plot = F)
    legend("bottom", legend = desc[l.ind], col = cols[l.ind], pch = 20, cex = 0.6, bty = "n")  
}


#' Import map from jpeg
#'
#' Load a jpeg map image and convert to a black-and-white binary raster
#' @param jpg.file String containing the full name (including path, if not in current working directory) of the jpeg to import.
#' @param plot Logical indicating whether to draw a plot of the imported data or not. Defaults to T.
#' @param threshold Numeric between 0 and 1, indicating the threshold at which points are divided into white or black. A higher threshold will only convert darker objects to black pixels, while a lower threshold will convert fainter objects to black. defaults to 0.2.
#' @return A new feature set: a list containing a raster of numbered features, and a matrix assigning a type to each feature.
#' @export
#' @examples
#' genlis <- import.map("Genlis-cropped.jpg", plot = T, threshold = 0.2)
#' 
#' 
import.map <- function(jpg.file, plot = T, threshold = 0.2) {
    
    # load the JPEG image
    jpegimage <- readJPEG(jpg.file)
    data <- matrix(1-jpegimage, nrow = nrow(jpegimage), ncol = ncol(jpegimage))
    
    # convert to raster
    r <- raster(data, xmn = 0, ymn = 0, xmx = 1, ymx = nrow(data)/ncol(data))
    
    # convert to binaary raster: threshold defaults to 0.5
    # (equivalent to rounding to 0dp):
    # should be adequate, since dealing with black & white images
    r0 <- r >= threshold
    
    # identify clumps of points
    cc <- clump(r0, dir = 8)
    
    # create table to store clump 'type'
    clump.types <- cbind(unique(getValues(cc)), 0)
    clump.types <- clump.types[!is.na(clump.types[,1]),]
    
    if (plot) {
        org.par <- par()
        par(mar = c(0,0,0,0))
        plot(cc, col = "black", asp = T, legend = F)
        par(mar = org.par$mar)
    }
    list(features = cc, feature.types = clump.types)
}

