# setwd('~/Development/gitter/')
# #Import packages
# library('EBImage')
# library('jpeg')
# library('logging')
# library('multicore')
#
# source('Peaks.R')
# source('Help.R')

.defaultFormat <- function(record) {
  text <- paste(paste(record$timestamp, record$levelname, record$logger, record$msg, sep=':'))
}

.methods = c('kmeans', 'tophat')

gitter.batch <- function(image.files, ref.image.file=NA, 
                         failed.file=file.path(getwd(), 'gitter_failed_plates.txt'), logging=F, ...){
  # Set logging
  logReset()
  if(logging) addHandler(writeToConsole, formatter=.defaultFormat)
  
  if(!file.create(failed.file)) 
    stop(sprintf('Invalid failed file: "%s". The failed file is the path of the file which will contain file names of failed plates (if any) after batch processing is completed.', failed.file))
  
  is.dir = file.info(image.files[1])$isdir
  if(is.dir){
    image.files = image.files[1]
    writeLines(sprintf('Reading images from directory: %s', image.files))
    image.files = list.files(image.files, pattern='jpe?g$', full.names=T, ignore.case=T)
    if(length(image.files) == 0) stop('No images with JPEG or JPG extension found. Images must be JPG format, please convert any non-JPG images to JPG')
  }
  
  z = sapply(image.files, file.exists)
  if(!all(z)) stop(sprintf('Files "%s" do not exist', paste0(image.files[!z], collapse=', ')))
  
  params = NA
  is.ref = !is.na(ref.image.file)
  if(is.ref){
    loginfo('Processing reference image: %s', ref.image.file)
    r = gitter(ref.image.file, ..., .is.ref=T)
    params = attr(r, 'params')
  }
  
  failed.plates = c()
  for(image.file in image.files){
    result = tryCatch({ gitter(image.file, ..., .params=params, .is.ref=F) }, 
                error = function(e) { 
                  logerror('Failed to process "%s", skipping', image.file)
                  e
                })
    
    # If we have an error
    if('error' %in% class(result)){
      failed.plates = c(failed.plates, basename(image.file))
    }
  }
  
  # Save failed plates
  if(length(failed.plates) > 0){
    failed.plates = c('# gitter failed plates', failed.plates)
    writeLines(failed.plates, failed.file)
  }
  #dats = lapply(image.files, gitter, ..., .params=params, .is.ref=F)
  
  #return(dats)
}


gitter <- function(image.file, plate.format=c(32,48), remove.noise=F, autorotate=F, inverse=F, 
                         method="kmeans", smooth.factor=5, logging=T, contrast=NA, fast=F, fast.width=1500,
                         plate.edges=F, plot=F, gridded.save.dir=getwd(), dat.save.dir=getwd(), 
                         .is.ref=F, .params=NA){
  
  
  if(!method %in% .methods) 
    stop(sprintf('Invalid method, possible methods are: %s', paste(.methods, collapse=', ')))
  if(!is.na(contrast) & contrast < 0) 
    stop('Contrast cannot be a negative value')
  if(grepl('^gridded', basename(image.file))) 
    warning('Detected gridded image as input')
  if(!is.na(gridded.save.dir)){
    if(!file.info(gridded.save.dir)$isdir) 
      stop(sprintf('Invalid gridded directory "%s"', gridded.save.dir))
  }
  if(!is.na(dat.save.dir)){
    if(!file.info(dat.save.dir)$isdir) 
      stop(sprintf('Invalid dat directory "%s"', gridded.save.dir))
  }
  if(fast.width < 1000 | fast.width > 4000) 
    stop('Resize width must be between 1000-4000px')
  
  expf = 1.5
  params = as.list(environment(), all=TRUE)
  nrow = plate.format[1]
  ncol = plate.format[2]
  # Are we using a reference screen?
  is.ref = all(is.na(.params))
  ptm <- proc.time()
  
  # Set logging
  logReset()
  if(logging) addHandler(writeToConsole, formatter=.defaultFormat)
  
  # Read image
  loginfo('Reading image from: %s', image.file)
  im = readJPEG(image.file)
  is.color = (length(dim(im)) == 3)
  im.grey = im
  
  if(fast){
    loginfo('Resizing image...')
    im = resize(im, w=1000)
  }
  
  if(autorotate){
    loginfo('Autorotating image...')
    im = .autoRotateImage(im)
  }
  
  # Extract greyscale
  if(is.color){
    loginfo('\tDetected color image, extracting greyscale')
    # Luminosity grey scale from GIMP
    im.grey = (im[,,1]*0.72) + (im[,,2]*0.21) + (im[,,3]*0.07)
  }else{
    loginfo('\tDetected greyscale image')
  }
  
  if(!is.na(contrast)){
    loginfo('\tIncreasing image contrast with factor %s', contrast)
    im.grey = .setContrast(im.grey, contrast)
  }
  
  if(inverse){
    loginfo('Inversing image')
    im.grey = 1 - im.grey 
  }
  
  if(!is.ref){
    loginfo('Non-reference plate, registering image to reference')
    # Fix any shifts 
    im.grey = .register2d(im.grey, .params$row.sums, .params$col.sums)
  }
  
  if(grepl('tophat', method)){
    si = round((nrow(im.grey) / nrow) * 1.0)
    loginfo('Opening image with kernel size %s', si)
    kern = makeBrush(si, 'box')
    im.grey = whiteTopHatGreyScale(im.grey, kern)
    #im.grey = setContrast(im.grey, -10)
    thresh = .findOptimalThreshold(im.grey)
    im.grey = (im.grey >= thresh)+0
    
  }
  
  if(plate.edges){
    loginfo('Eroding plate edges for peak detection...')
    kern = matrix(0,3,9)
    kern[2,] = 1
    z1 = openingGreyScale(im.grey, kern)
    z2 = openingGreyScale(im.grey, t(kern))
    sum.y = rowSums(z2)
    sum.x = colSums(z1)
  }else{
    # Sum up rows / columns
    sum.y = rowSums(im.grey)
    sum.x = colSums(im.grey)
  }
  
  
  if(is.ref){
    # Get peaks of sums
    if(plot) par(mfrow=c(2,1), bty='n', las=1)
    loginfo('Getting row peaks...')
    cp.y = .getColonyPeaks(sum.y, n=nrow, smooth.factor, plot)
    loginfo('Getting column peaks...')
    cp.x = .getColonyPeaks(sum.x, n=ncol, smooth.factor, plot)
    
    # Average window (though they should be the same)
    w = round(mean( c(cp.x$window, cp.y$window) ))
    
    # Get center coordinates
    coords = expand.grid(cp.x$peaks, cp.y$peaks)
    names(coords) = c('x', 'y')
    
    params$window = w
    params$coords = coords
    params$row.sums = sum.y
    params$col.sums = sum.x
  }else{
    w = .params$window
    coords = .params$coords
  }
  
  # Window expansion factor
  d = round(w * expf)
  
  #Pad image
  im.pad =  .padmatrix(im.grey, w, 1)
  
  coords[,c('x','y')] = coords[,c('x','y')]+w
  coords$xl = coords$x - d
  coords$xr = coords$x + d
  coords$yt = coords$y - d
  coords$yb = coords$y + d
  
  if(! grepl('tophat', method)){
    loginfo('Binarizing image... ')
    im.pad = .getKmeansImage(im.pad, coords)
  }
  
  if(remove.noise){
    # Remove the noise
    loginfo('Denoising image... ')
    kern = makeBrush(3, 'diamond')
    im.pad = dilateGreyScale(erodeGreyScale(im.pad, kern), kern)
  }
  
  loginfo('Fitting rectangles... ')
  coords = .fitRects(coords, im.pad, w)
  
  im.grey = .unpadmatrix(im.pad, w)
  
  coords[,c('x','y','xl','xr','yt','yb')] = coords[,c('x','y','xl','xr','yt','yb')]-w
  coords[coords<0] = 1
  coords[,c('xl', 'xr')][coords[,c('xl', 'xr')] > ncol(im.grey)] = ncol(im.grey)
  coords[,c('yt', 'yb')][coords[,c('yt', 'yb')] > nrow(im.grey)] = nrow(im.grey)
  
  # Generate rows and columns to be bound to coords below
  rc = expand.grid(1:nrow, 1:ncol)
  names(rc) = c('row', 'col')
  results = cbind(rc, coords)
  
  # Print elapsed time
  elapsed = signif( (proc.time() - ptm)[['elapsed']], 5)
  
  # Rearrange
  results = results[,c('row', 'col', 'size', 'circularity',
                       'x', 'y', 'xl','xr', 'yt', 'yb')]
  
  class(results) = c('gitter', 'data.frame')
  
  # Save important attributes
  attr(results, 'params') = params
  attr(results, 'elapsed') = elapsed
  attr(results, 'call') = match.call()
  
  # Save gridded image
  if(!is.na(gridded.save.dir) & !.is.ref){
    #imr = drawPeaks(peaks.c=cp.y$peaks, peaks.r=cp.x$peaks, imr)
    imr = .drawRect(coords[,3:6], im.grey)
    save = file.path(gridded.save.dir, paste0('gridded_',basename(image.file)))
    loginfo('Saved gridded image to: %s', save)
    writeJPEG(imr, save, 1)
  }
  
  # Save dat file
  if(!is.na(dat.save.dir) & !.is.ref){
    save = file.path(dat.save.dir, paste0(basename(image.file), '.dat'))
    write.table(results, file=save, quote=F, sep='\t', row.names=F)
    loginfo('Saved dat file to: %s', save)
  }
  
  # Garbage collector
  gc(reset=T, verbose=F)
  
  loginfo('Time elapsed: %s seconds', elapsed)
  return(results)
}

.getKmeansImage <- function(im, coords){
  coords = coords[,c('xl', 'xr', 'yt', 'yb')]
  cent = .getMeans(im)
  ret = im
  t = matrix(T, nrow(im), ncol(im))
  for(i in 1:nrow(coords)){
    z = as.numeric(coords[i,])
    spot = im[z[3]:z[4],z[1]:z[2]]
    im.spot.vec = c(cent, spot)
    #cent = c(min(im.spot.vec), max(im.spot.vec))
    cl = kmeans(im.spot.vec, cent)
    cls = cl$cluster[-(1:2)] - 1
    spot.bw = matrix(cls, nrow=nrow(spot))
    ret[z[3]:z[4],z[1]:z[2]] = spot.bw
    t[z[3]:z[4],z[1]:z[2]] = F
  }
  
  ret[t] = ( ret[t] > 0.7 )+0
  return(ret)
}


.xl <- function(z, w){ m = which(z == min(z)); t = length(z) - m[length(m)]; if(t < w) t = w; return(t) }
.xr <- function(z, w){ t = which.min(z); if(t < w) t = w; return(t) }
.yt <- .xl
.yb <- .xr

.fitRects <- function(coords, im.kmeans, d){
  
  # Minimum border for really small colonies
  minb = d/3
  
  ret = lapply(1:nrow(coords), function(i){
    # Center of spot
    x = coords$x[i]
    y = coords$y[i]
    
    cent.pixel = im.kmeans[y,x]
    
    # Define expanded rectangle 
    rect = c(coords$xl[i], coords$xr[i], coords$yt[i], coords$yb[i])
    
    spot.bw = im.kmeans[rect[3]:rect[4], rect[1]:rect[2]]
    
    if(cent.pixel == 0){
      z = rep(minb*2, 4)
    }else{
      # Sum rows and cols, split sums in half 
      sp.y = .splitHalf(rowSums(spot.bw))
      sp.x = .splitHalf(colSums(spot.bw))
      # Get minimums relative to the spot and set minimums if too small
      z = c(.xl(sp.x$left, minb), 
            .xr(sp.x$right, minb),
            .yt(sp.y$left, minb), 
            .yb(sp.y$right, minb))
    }
    
    # Crop bw spot and compute new size
    x.rel = x - rect[1]
    y.rel = y - rect[3]
    
    rect.rel = c(x.rel - z[1], x.rel + z[2], y.rel - z[3], y.rel + z[4])
    spot.bw = spot.bw[rect.rel[3]:rect.rel[4], rect.rel[1]:rect.rel[2]]
    # Compute new spot relative to plate
    rect = c(x - z[1], x + z[2], y - z[3], y + z[4])
    
    c(x, y, rect, sum(spot.bw), .circularity(spot.bw))
  })
  
  coords = matrix(simplify2array(ret), nrow(coords), byrow=T)
  coords = as.data.frame(coords)
  names(coords) = c('x', 'y', 'xl', 'xr', 'yt', 'yb', 'size', 'circularity')
  return(coords)
}


.register2d <- function(im, r2, c2, lag.max=100){
  r1 = rowSums(im)
  c1 = colSums(im)
  
  z1 = ccf(r1, r2, lag.max, plot=F)
  y_s = as.vector(z1$lag)[ which.max(as.vector(z1$acf)) ]
  z2 = ccf(c1, c2, lag.max, plot=F)
  x_s = as.vector(z2$lag)[ which.max(as.vector(z2$acf)) ]
  
  return(translate(im, c(y_s,x_s)))
}