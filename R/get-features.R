
#' Identify map scale
#'
#' Identify the scale marker on the image. Will either output a list containing the details of the scale identified, or a new feature set object with the raster rescaled as directed by the user. If automatic detection fails, the legend can be selected automatically by the user. (Can be run more than once to rescale map repeatedly if necessary)
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param l Length of horizontal strip of pixels to search for. If the image contains a lot of long horizontal objects, this may need to be increased (default is 17).
#' @return New feature set list containing a rescaled raster of numbered features and a matrix of corresponding feature types.
#' @details Uses a focal raster to find clusters containing large horizontal strips. Function will also rescale map to use life-size coordinate system, but doesn't have to; can be run again at any time to rescale further).
#' @export
#' @examples
#' get.scale(genlis)
#' # will output a new feature set object, with coordinates rescaled
get.scale <- function(site, l = 17) {
    
    ct <- site$feature.types
    cc <- site$features
    cs <- cbind(ct[,1],0)
        
    # recreate binary data from clumps
    r0 <- reclassify(cc, rcl = rbind(cbind(ct[,1], 1), c(NA, 0)))
    
    org.par <- par()
    par(mar = c(0,0,0,0))
    
    if (length(ct[ct[,2] == 2, 1]) == 1) {
        candidate <- ct[ct[,2] == 2, 1]
        cs[candidate, 2] <- 2
        
        # show it on the map
        plot(reclassify(cc, rcl = cs), col = c("Black", "Red"), asp = T, legend = F)
        par(mar = org.par$mar)
    } else {
        # filter to look for long horizontal lines of black cells
        scale.filter <- matrix(c(0,1,0), nrow = 3, ncol = l)
        
        # only picks up the highest-scoring regions
        cf <- table(cc[focal(r0, w = scale.filter, pad = TRUE, padValue = 0) == l])
        
        # find the largest of the highest-scoring regions, set to 2
        candidate <- as.numeric(names(cf)[which.max(cf)])
        reset <- ct[as.numeric(names(cf)[which.max(cf)]),2]
        ct[candidate,2] <- 2 
        cs[candidate,2] <- 2
        
        # show it on the map
        plot(reclassify(cc, rcl = cs), col = c("Black", "Red"), asp = T, legend = F)
        par(mar = org.par$mar)
        
        # get confirmation: is this the right scale item?
        conf <- readline(prompt = "Is the highlighted object the map scale? (y/n) > ")
        if (tolower(substr(conf, 1, 1)) != "y") {
            ct[candidate,2] <- reset
            cs[candidate,2] <- 0
            print("Please click on the scale marker...", quote = F)    
            scale.clump <- click(cc, n = 1, xy = T, id = F, show = F)
            
            # If not on a point, get nearest point 
            if (is.na(scale.clump$value)) {
                coords <- cbind(getValues(cc), xyFromCell(cc, cell = 1:ncell(cc)))
                coords <- coords[!(is.na(coords[,1])), ]
                nn <- knnx.index(coords[,2:3], data.frame(scale.clump$x, scale.clump$y), k = 1)
                candidate <- coords[as.numeric(nn),1]
            } else {
                candidate <- scale.clump$value
            }  
            ct[candidate,2] <- 2
            cs[candidate,2] <- 2
            plot(reclassify(cc, rcl = cs), col = c("Black", "Red"), asp = T, legend = F)
        }
    }           
    # get real length of scale
    m.scale <- "x"
    while (is.na(suppressWarnings(as.numeric(m.scale)))) {
        m.scale <- readline(prompt = "What is the length, in metres, represented by the map scale? > ")
    }
    m.scale <- as.numeric(m.scale)
    
    # update feature found as type 2
    ct <- ct
    
    # get length of scale on plan
    sc <- xyFromCell(cc, 1:ncell(cc))
    sc <- data.frame(sc[(!is.na(getValues(cc))) & (getValues(cc) == candidate),])
    scale.length <- max(sc$x) - min(sc$x)
    sc.check <- min(sc$x) + scale.length
    
    rescale <- readline(prompt = "Do you want to rescale the map coordinates now? > ")
    if (tolower(substr(rescale, 1, 1)) != "y") {
        list(true.length = m.scale, scale.length = scale.length, clump.id = candidate)
    } else {
        scale <- scale.length / m.scale
        cc.xy <- rasterFromXYZ(cbind(xyFromCell(cc, cell = 1:ncell(cc))/ scale,
                                     z = getValues(cc)))
        print("New object 'rescaled' created containing rescaled raster.", quote = F)
        rescaled <<- list(features = cc.xy, feature.types = ct)
    }
}


#' Identify N-S marker on map
#'
#' Identify the N-S marker on the image and calculate the angle of the N-S axis. Since the marker shape varies from plan to plan, it must be selected manually.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @return The angle (in radians) of the N-S axis.
#' @export
#' @examples
#' get.orientation(genlis)
get.orientation <- function(site) {
    
    cc <- site$features
    ct <- site$feature.types
    coords <- xyFromCell(cc, cell = 1:ncell(cc))
    coords <- cbind(getValues(cc), coords)
    coords <- coords[!(is.na(coords[,1])), ]
    
    if (length(ct[ct[,2] == 3, 1]) == 1) {
        # already identified - just need to extract angle
        target <- ct[ct[,2] == 3, 1]
    } else {
        # manually select object to use as marker
        require(FNN)        # needed for knnx.index (to find nearest neighbour to click)
        show.features(site)
        print("Please click on the N-S marker...", quote = F)    
        scale.clump <- click(site$features, n = 1, xy = T, id = F, show = F)
        
        # If not on a point, get nearest point 
        if (is.na(scale.clump$value)) {
            nn <- knnx.index(coords[,2:3], data.frame(scale.clump$x, scale.clump$y), k = 1)
            target <- coords[as.numeric(nn),1]
        } else {
            target <- scale.clump$value
        }
    }
    
    # extract points in selected shape
    coords <- coords[coords[,1] == target, 2:3]
    
    # idenfity longest axis
    L <- coords[which.min(coords[,1]),]
    R <- coords[which.max(coords[,1]),]
    B <- coords[which.min(coords[,2]),]
    U <- coords[which.max(coords[,2]),]
    
    h <- sqrt(sum((L-R)^2))
    v <- sqrt(sum((U-B)^2))
    if (v > h) {
        axis <- atan2((U-B)[1], (U-B)[2])
    } else {
        axis <- atan2((R-L)[1], (R-L)[2])
    }
    
    ct[target, 2] <- 3
    
    angle.found <- list(features = site$features, feature.types = ct)
    show.features(angle.found)
    print("New object 'NS.marked' created with updated feature types.", quote = F)
    NS.marked <<- angle.found
    axis
}


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
assign.by.size <- function(site, lower = 3, mid = 50, upper = 100, ratio = 0.7) {
    
    # get sizes of all clumps identified
    cc <- site$features
    fr <- freq(cc)
    fr <- cbind(fr[!is.na(fr[,1]),],0)
    
    xy <- cbind(xyFromCell(cc, cell = 1:ncell(cc)), z = getValues(cc))
    xy <- xy[!is.na(xy[,3]),]
    
    for (i in 1:nrow(fr)) {
        coords <- matrix(xy[xy[,3] == i, ], ncol = 3)
        max.l <- max(length(unique(coords[,1])),length(unique(coords[,2])))
        min.l <- min(length(unique(coords[,1])),length(unique(coords[,2])))
        
        fr[i,3] <- min.l / max.l
    }        
    
    ct <- site$feature.types
    
    # post-holes: small, round/square in shape
    ct[(fr[,3] >= ratio) & (fr[,2] <= mid) & (ct[,2] == 0), 2] <- 1
    
    # large features
    ct[(ct[,2] == 0) & (fr[,2] >= upper), 2] <- 5
    
    # noise: very small features
    ct[(ct[,2] == 0) & (fr[,2] < lower), 2] <- 4
    
    # text: smallish, non-round features
    ct[(ct[,2] == 0) & (fr[,2] < upper) & (fr[,2] > mid), 2] <- 4
    
    new.types <<- list(features = cc, feature.types = ct)
    show.features(new.types)
    title (main = "Features reclassified according to size & compactness")
}


#' Identify vertical features
#'
#' Identify features that contain a vertical strip of black pixels
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param l Length of black pixels to be tested for (should be an odd number, but function will add 1 if an even number is supplied). Default is 9.
#' @param plot Logical indicator: should the new set of features be displayed on a plot?
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
#' @export
#' @examples
#' get.vertical(genlis)
get.vertical <- function(site, l = 9, plot = T) {
    
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
#' @param plot Logical indicator: should the new set of features be displayed on a plot?
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
#' @export
#' @examples
#' get.horizontal(genlis)
get.horizontal <- function(site, l = 9, plot = T) {
    
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


#' Identify any features that lie directly between two annotations
#'
#' For each feature classed as an annotation, the function extends the feature's footprint directly left, right, up and down. Any features that appear twice in these extensions (ie. features that are directly horizontally or vertically between two annotation features) are also classed as annotations. Helps with picking up broken line boundaries, dots in text etc.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @param plot Logical indicator: should the new set of features be displayed on a plot?
#' @return New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
#' @export
#' @examples
#' extend.annotations(genlis)
extend.annotations <- function(site, plot = T) {
    
    ct <- site$feature.types
    cc <- site$features
    
    # only consider annotation marks
    ann <- ct[ct[,2] == 4,1]
    extensions <- list()
    
    for (i in 1:length(ann)) {
        # create a solid block of points
        xy <- xyFromCell(cc, Which(cc == ann[i], cells=TRUE))
        all.x <- unique(xy[,1])
        all.y <- unique(xy[,2])
        block <- cbind(x = rep(all.x, length(all.y)),
                       y = sort(rep(all.y, length(all.x))))
        
        offset.x <- c(rep(max(block[,1]) - min(block[,1]), nrow(block)), rep(0, nrow(block)))        
        offset.y <- c(rep(0, nrow(block)), rep(max(block[,2]) - min(block[,2]), nrow(block)))
        
        neighbours <- unique(c(extract(cc, block - offset.x), extract(cc, block + offset.x),
                               extract(cc, block - offset.y), extract(cc, block + offset.y)))
        neighbours <- neighbours[(!is.na(neighbours)) & (neighbours != ann[i])]
        extensions[[i]] <- neighbours
    }
    ext <- table(unlist(extensions))
    
    show.features(site)
    
    ct[(ct[,1] %in% names(ext[ext == 2])) & (ct[,2] <= 1),2] <- 4
    with.extensions <<- list(features = cc, feature.types = ct)
    if (plot) {show.features(with.extensions)
               title (main = "Features directly between two other annotations reclassified as annotations")
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


#' Extract centre-points of post-holes
#'
#' For each post-hole feature identified, get the x and y coordinates of the centre of the feature.
#' @param site Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.
#' @return Matrix of coordinates showing the midpoints of the features selected
#' @export
#' @examples
#' get.postholes(genlis)
get.postholes <- function(site, plot = T) {
    cc <- site$features
    ct <- site$feature.types
    
    # get mid-points of all clumps
    xy <- data.frame(cbind(id = getValues(cc), xyFromCell(cc, 1:ncell(cc))))
    mids <- ddply(xy[!is.na(xy$id),], .(id), summarise, xm = mean(x), ym = mean(y))
    
    # filter to include only centroids of post-holes
    # can use 'extract(cc, centres)' to get feature id numbers if necessary
    centres <<- mids[ct[,2] == 1,2:3]
    if (plot) {
        x.lim <- c(floor(min(xy[,2])/10) * 10, ceiling(max(xy[,2])/10) * 10)
        y.lim <- c(floor(min(xy[,3])/10) * 10, ceiling(max(xy[,3])/10) * 10)
        
        plot(cc, col = "grey", asp = T, legend = F, xlim = x.lim, ylim = y.lim)
        points(centres[,1], centres[,2], pch = 20, cex = 0.5, asp = T, xlim = x.lim, ylim = y.lim)
    }   
}


