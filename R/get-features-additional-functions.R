
#' Assign categories to features based on size & shape
#'
#' For each feature, get the ratio of the longest to shortest bounding edge. Small features with a ratio close to 1 are set as post-holes; slightly larger as annotations; the largest as linear features.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param lower Cell count threshold below which a feature is classed as noise.
#' @param mid Cell count threshold below which a feature is classed as a post-hole.
#' @param upper Cell count threshold above which a feature is classed as a linear feature of the map.
#' @param ratio Ratio (of longest to shortest edge) above which
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types. Also plots the new feature types.
#' @export
#' @examples
#' assign.by.size(genlis)
assign.by.size <- function(site, lower = 3, mid = 50, upper = 100, ratio = 0.8, density = 0.2, plot = T) {
    
    # get sizes of all clumps identified
    cc <- site$features
    fr <- freq(cc)
    fr <- cbind(fr[!is.na(fr[,1]),], ratio = 0, sq.density = 0)
    
    xy <- cbind(xyFromCell(cc, cell = 1:ncell(cc)), z = getValues(cc))
    xy <- xy[!is.na(xy[,3]),]
    
    for (i in 1:nrow(fr)) {
        coords <- matrix(xy[xy[,3] == i, ], ncol = 3)
        max.l <- max(length(unique(coords[,1])),length(unique(coords[,2])))
        min.l <- min(length(unique(coords[,1])),length(unique(coords[,2])))
        
        # ratio of shape sides
        fr[i,3] <- round(min.l / max.l,3)
        
        # shape density: if shape was bounded by a square, how much is filled?
        fr[i,4] <- round(fr[i,2] / max.l^2,3)
    }        
    
    ct <- site$feature.types
    
    # post-holes: small, round/square in shape
    ct[(fr[,3] >= ratio) & (fr[,2] <= mid) & (ct[,2] == 0), 2] <- 1
    
    # large features
    ct[(ct[,2] == 0) & (fr[,2] >= upper), 2] <- 5
    
    # annotations/smaller linear features: long, thin shapes 
    ct[(ct[,2] <= 1) & (fr[,4] < density), 2] <- 4
    
    # noise: very small features
    ct[(ct[,2] == 0) & (fr[,2] < lower), 2] <- 4
    
    # text: smallish, non-round features
    ct[(ct[,2] == 0) & (fr[,2] < upper) & (fr[,2] > mid), 2] <- 4
    
    new.types <<- list(features = cc, feature.types = ct)
    if (plot) {show.features(new.types)
               title (main = "Features reclassified according to size & compactness")}
}


#' Identify vertical features
#'
#' Identify features that contain a vertical strip of black pixels
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param l Length of black pixels to be tested for (should be an odd number, but function will add 1 if an even number is supplied). Default is 9.
#' @param plot Boolean indicator: should the new set of features be displayed on a plot?
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
#' @export
#' @examples
#' get.vertical(genlis)
get.verticals <- function(site, l = 9, plot = T) {
    
    # if l even, make odd
    if(l/2 == round(l/2,0)) {
        print(paste("l must be odd: set to ", l+1, ".", sep = ""), quote = F)
        l = l+1
    }
    line.filter <- t(matrix(c(0,1,0), nrow = 3, ncol = l))
    
    ct <- site$feature.types
    cc <- site$features
    
    # recreate binary data from clumps
    r0 <- reclassify(cc, rcl = rbind(cbind(ct[,1], 1), c(NA, 0)))
    
    # only pick up the highest-scoring regions
    f <- focal(r0, w = line.filter, pad = TRUE, padValue = 0) == l
    
    ct[(ct[,1] %in% unique(cc[f])) & (ct[,2] <= 1), 2] <- 4 
    
    verticals <<- list(features = cc, feature.types = ct)
    if (plot) {show.features(verticals)
               title (main = "Features containing vertical lines reclassified as annotations")
    }
}


#' Identify horizontal features
#'
#' Identify unclassified/post-hole features that contain a horizontal strip of black pixels, and set as annotations.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param l Length of black pixels to be tested for (should be an odd number, but function will add 1 if an even number is supplied). Default is 9.
#' @param plot Boolean indicator: should the new set of features be displayed on a plot?
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
#' @export
#' @examples
#' get.horizontal(genlis)
get.horizontals <- function(site, l = 7, plot = T) {
    
    # if l even, make odd
    if(l/2 == round(l/2,0)) {
        print(paste("l must be odd: set to ", l+1, ".", sep = ""), quote = F)
        l = l+1
    }
    
    ct <- site$feature.types
    cc <- site$features
    
    # recreate binary data from clumps
    r0 <- reclassify(cc, rcl = rbind(cbind(ct[,1], 1), c(NA, 0)))
    
    line.filter <- (matrix(c(0,1,0), nrow = 3, ncol = l))
    f <- (focal(r0, w = line.filter, pad = TRUE, padValue = 0)) == l
    
    # recreate binary data from clumps
    r0 <- reclassify(cc, rcl = rbind(cbind(ct[,1], 1), c(NA, 0)))
    
    # only pick up the highest-scoring regions
    f <- focal(r0, w = line.filter, pad = TRUE, padValue = 0) == l
    
    ct[(ct[,1] %in% unique(cc[f])) & (ct[,2] <= 1), 2] <- 4 
    
    horizontals <<- list(features = cc, feature.types = ct)
    if (plot) {show.features(horizontals)
               title (main = "Features containing horizontal lines reclassified as annotations")
    }
}



#' Manually reclassify a single feature
#'
#' A last resort: manually select a single feature from the map and reclassify it as a post-hole, an annotation, or a linear feature of the map.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param replot Boolean: do you want to redraw the site before selecting the feature to reclassify?
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
#' @export
#' @examples
#' manually.set.feature.type(genlis)
manually.set.feature.type <- function(site, replot = F) {
    cc <- site$features
    ct <- site$feature.types
    
    if (replot) {show.features(site)}
    cat("Please click on the feature to be manually reclassified...")    
    scale.clump <- click(cc, n = 1, xy = T, id = F, show = F)
    
    # If not on a point, get nearest point 
    if (is.na(scale.clump$value)) {
        coords <- xyFromCell(cc, cell = 1:ncell(cc))
        coords <- cbind(getValues(cc), coords)
        coords <- coords[!(is.na(coords[,1])), ]
        
        nn <- knnx.index(coords[,2:3], data.frame(scale.clump$x, scale.clump$y), k = 1)
        target <- coords[as.numeric(nn),1]
    } else {
        target <- scale.clump$value
    }
    sl <- cbind(ct[,1], 0)
    sl[target, 2] <- 1
    plot(reclassify(cc, rcl = sl), col = c("Grey", "Red"), legend = F, axes = F)
    
    upd <- F
    conf <- readline(prompt = "Is the highlighted object the one you want to change? (y/n) > ")
    if (tolower(substr(conf, 1, 1)) != "y") {
        cat("No change made.")
    } else {
        cat("Select new feature type: \n   1: post-hole \n   2: annotation \n   3: large linear map feature \n   x: exit without changing")
        new.type <- readline(prompt = "New feature type: > ")
        
        if (tolower(substr(new.type, 1, 1)) == "x") {
            cat("No change made.")
        } else {
            if (new.type == 1) {ct[target, 2] <- 1; upd <- T}     # post-hole
            if (new.type == 2) {ct[target, 2] <- 1; upd <- T}     # annotation
            if (new.type == 3) {ct[target, 2] <- 1; upd <- T}     # large linear map feature
        }   
    }
    if (upd) {
        cat(" New feature set output to 'updated.features'. \n Remember to rename this before making any further manual changes, as it will be overwritten.")
        updated.features <<- list(features = cc, feature.types = ct)
    }
}


#' Remove isolated features
#'
#' For each post-hole feature identified, find the distance to the nearest neighbour. Extremely isolated points (outliers) are defined in the same way as in a boxplot: any points whose nearest-neighbour distance lies more than 1.5 x the IQR above the third quartile.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param plot Boolean indicator: should the new set of features be displayed on a plot?
#' @export
#' @examples
#' remove.isolated(genlis)
remove.isolated <- function(site, plot = T) {
    
    cc <- site$features
    ct <- site$feature.types
    
    # if no points set as postholes yet, set all unclassified points as postholes
    if (length(ct[ct[,2]==1, 2]) == 0) {
        ct[ct[,2]==0, 2] <- 1
    } 
    # get centrepoints of all post-hole features
    xy <- data.frame(cbind(id = getValues(cc), xyFromCell(cc, 1:ncell(cc))))
    mids <- ddply(xy[!is.na(xy$id),], .(id), summarise, xm = mean(x), ym = mean(y))
    pts <- mids[ct[,2] == 1,2:3]
    
    d <- knn.dist(pts, k = 1)
    l <- quantile(d, 0.75) + (1.5 * IQR(d))
    ct[c(as.numeric(rownames(pts[d > l,]))), 2] <- 0
    
    density.filtered <<- list(features = cc, feature.types = ct)
    
    if (plot) {
        plot(site$features, col = "grey", asp = T, legend = F)
        points(pts[d <= l, 1], pts[d <= l, 2], pch = 20, cex = 0.5, asp = T)
        points(pts[d > l, 1], pts[d > l, 2], col = "blue", asp = T)
    }
}



#' Identify a feature from the map
#'
#' Get the id number and x, y coordinates from a particular feature by clicking on it on the map. Can be used in conjunction with \code{feature.dims} to check the height, density and other characteristics of a particular feature.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param replot Boolean indicator: should the map be redrawn, or can the function use the currently drawn plot?
#' @export
#' @examples
#' select.feature(genlis, replot = F)
select.feature <- function(site, replot = F) {
    if (replot) {
        plot(reclassify(site$features, site$feature.types), col = "grey", asp = T, legend = F)
    }
    print("Click on the map to select a feature...", quote = F)
    
    s <- click(site$features, n = 1, xy = T, id = F, show = F)
    
    # If not on a point, get nearest point 
    if (is.na(s$value)) {
        coords <- cbind(getValues(site$features),xyFromCell(site$features, cell = 1:ncell(site$features)))
        coords <- coords[!(is.na(coords[,1])), ]
        nn <- knnx.index(coords[,2:3], data.frame(s$x, s$y), k = 1)
        target <- coords[as.numeric(nn),1]
    } else {
        target <- s$value
    }
    list(x = s$x, y = s$y, index = target)
}



#' Get feature dimensions
#'
#' Creates a matrix containing the dimensions of every feature in the map: feature id, midpoint, max and min x and y coordinates; width, height and size in terms of cells covered; and ratio of longest to shortest edge.
#' Can be used to investigate plausible thresholds for the other functions.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @return A matrix containing the dimensions of every feature identified in the map.
#' @export
#' @examples
#' feature.dims(genlis)
feature.dims <- function(site) {
    
    xy <- data.frame(cbind(id = getValues(site$features), xyFromCell(site$features, 1:ncell(site$features))))
    xy <- xy[!is.na(xy[,1]),]
    
    dims <- ddply(xy, .(id), summarise, xm = mean(x), ym = mean(y),
                  x.max = max(x), x.min = min(x),
                  y.max = max(y), y.min = min(y),
                  freq = length(x))
    
    # add width & height of shape, in cells
    dims <- cbind(dims, width = (dims$x.max - dims$x.min) / res(site$features)[1] + 1,
                  height = (dims$y.max - dims$y.min) / res(site$features)[2] + 1)
    
    # get shape density: if it were bounded by a square, what proportion of cells are coloured?
    d.max <- apply(cbind(dims$width, dims$height), 1, max)
    d.min <- apply(cbind(dims$width, dims$height), 1, min)
    
    dims$density <- dims$freq/(d.max^2)
    dims$ratio <- d.min/d.max
    dims
}



