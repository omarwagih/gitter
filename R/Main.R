# z = setwd('~/Development/gitter/R')
# #Import packages
# library('EBImage')
# library('jpeg')
# library('logging')
# library('parallel')
# library('PET')
# library('ggplot2')
# 
# source('Peaks.R')
# source('Help.R')
# 
# setwd(z)

GITTER_VERSION = '1.0.1'
.methods = c('kmeans', 'tophat')
.pf = list('1536'=c(32,48), '384'=c(16,24), '96'=c(8,12))  

gitter.example <- function(eg='single'){
  if(!eg %in% c('single', 'ref')) stop('Invalid example')
  
  if(eg == 'single'){
    dat = gitter(system.file("extdata", "sample.jpg", package="gitter"))
    p <- plot.gitter(dat, title=sprintf('gitter v%s single image example', GITTER_VERSION))
    print(p)
    browseURL(f)
    summary(dat)
  }
  if(eg == 'ref'){
    f = system.file("extdata", "sample_dead.jpg", package="gitter")
    f.ref = system.file("extdata", "sample.jpg", package="gitter")
    gitter.batch(f, f.ref)
    warning(sprintf('NOTE: Output files were saved to working directory at %s', getwd()))
  }
}

.defaultFormat <- function(record) {
  text <- paste(paste(record$timestamp, record$levelname, record$logger, record$msg, sep=':'))
}

gitter.batch <- function(image.files, ref.image.file=NA, logging=F, ...){
  f = 'gitter_failed_images'
  ff = list.files(pattern=f)
  if(length(ff) > 0){
    f = paste0(f, format(Sys.time(), "_%d-%m-%y_%H-%M-%S"))
  }
  failed.file = paste0(f, '.txt')
  
  # Set logging
  logReset()
  if(logging) addHandler(writeToConsole, formatter=.defaultFormat)
  
  is.dir = file.info(image.files[1])$isdir
  if(is.dir){
    image.files = image.files[1]
    loginfo('Reading images from directory: %s', image.files)
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
    failed.plates = c('# gitter failed images', failed.plates)
    writeLines(failed.plates, failed.file)
  }
  #dats = lapply(image.files, gitter, ..., .params=params, .is.ref=F)
  
  #return(dats)
}


gitter <- function(image.file=file.choose(), plate.format=c(32,48), remove.noise=F, autorotate=F, inverse=F,
                   method="kmeans", logging=T, contrast=NA, fast=F, fast.width=1500, plot=F, 
                   gridded.save.dir=getwd(), dat.save.dir=getwd(), 
                   .is.ref=F, .params=NA){
  
  # Check if we have one number plate formats
  if(length(plate.format) == 1){
    t = as.character(plate.format)
    if(t %in% names(.pf)){
      plate.format = .pf[[t]]
    }else{
      stop('Invalid plate density, please use 1536, 384 or 96. If the density of your plate is not listed, you can specifcy a vector of the number of rows and columns in your plate (e.g. c(32,48))')
    }
  }
  # Check for incorrect plate formats
  if(length(plate.format) != 2){
    stop('Invalid plate format, plate formats must be a vector of the number of rows and columns (e.g. c(32,48)) or a value indicating the density of the plate (e.g 1536, 384 or 96) possible')
  }
  
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
    im = resize(im, h=fast.width)
  }
  
  # Extract greyscale
  if(is.color){
    loginfo('\tDetected color image, extracting greyscale')
    # Luminosity grey scale from GIMP
    im.grey = (im[,,1]*0.72) + (im[,,2]*0.21) + (im[,,3]*0.07)
  }else{
    loginfo('\tDetected greyscale image')
  }
  
  
  if(autorotate){
    loginfo('Autorotating image...')
    im.grey = .autoRotateImage2(im.grey)
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
  
  sum.y = rowSums(im.grey)
  sum.x = colSums(im.grey)
  
#     loginfo('Eroding plate edges for peak detection...')
#     kern = matrix(0,3,9)
#     kern[2,] = 1
#     z1 = openingGreyScale(im.grey, kern)
#     z2 = openingGreyScale(im.grey, t(kern))
#     sum.y = rowSums(z2)
#     sum.x = colSums(z1)

  
  
  if(is.ref){
    # Get peaks of sums
    if(plot) par(mfrow=c(2,1), bty='n', las=1)
    loginfo('Getting row peaks...')
    #cp.y = .getColonyPeaks(sum.y, n=nrow, smooth.factor, plot)
    cp.y = .getColonyPeaks2(sum.y, n=nrow, plot)
    loginfo('Getting column peaks...')
    #cp.x = .getColonyPeaks(sum.x, n=ncol, smooth.factor, plot)
    cp.x = .getColonyPeaks2(sum.x, n=ncol, plot)
    
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
  rc = expand.grid(1:ncol, 1:nrow)[,c(2,1)]
  names(rc) = c('row', 'col')
  results = cbind(rc, coords)
  
  # Print elapsed time
  elapsed = signif( (proc.time() - ptm)[['elapsed']], 5)
  
  # Rearrange
  results = results[,c('row', 'col', 'size', 'circularity',
                       'x', 'y', 'xl','xr', 'yt', 'yb')]
  
  class(results) = c('gitter', 'data.frame')
  
  
  # Save gridded image
  if(!is.na(gridded.save.dir) & !.is.ref){
    #imr = drawPeaks(peaks.c=cp.y$peaks, peaks.r=cp.x$peaks, imr)
    imr = .drawRect(coords[,3:6], im.grey)
    save = file.path(gridded.save.dir, paste0('gridded_',basename(image.file)))
    loginfo('Saved gridded image to: %s', save)
    writeJPEG(imr, save)
  }
  
  # Save dat file
  if(!is.na(dat.save.dir) & !.is.ref){
    save = file.path(dat.save.dir, paste0(basename(image.file), '.dat'))
    
    results[[4]] = round(results[[4]], 4)
    results = results[,1:4]
    
    write.table(results, file=save, quote=F, sep='\t', row.names=F)
    loginfo('Saved dat file to: %s', save)
  }
  
  # Garbage collector
  gc(reset=T, verbose=F)
  
  loginfo('Time elapsed: %s seconds', elapsed)
  
  # Save important attributes
  attr(results, 'params') = params
  attr(results, 'elapsed') = elapsed
  attr(results, 'call') = match.call()
  attr(results, 'file') = image.file
  attr(results, 'format') = plate.format
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



plot.gitter <- function(dat, title='', type='heatmap', low='turquoise', mid='black', high='yellow', 
                        show.text=F, text.color='white', norm=T, 
                        show.circ=T, circ.cutoff=0.6, circ.color='white'){
  
  if(is.character(dat)) dat = read.table(dat, stringsAsFactors=F, header=T)
  if(!is.data.frame(dat) & ! 'gitter' %in% class(dat)) stop('Data must be a data path or data frame generated by gitter')
  if(!type %in% c('heatmap', 'bubble')) stop('Invalid plot type. Use "heatmap" or "bubble"')
  
  if(length(dat) > 4 | length(dat) < 3) stop('Invalid number of columns for dat file')
  
  names(dat) = c('r', 'c', 's')
  dat.cs = NULL
  
  
  
  r = max(dat$r)
  c = max(dat$c)
  
  t = r:1
  names(t) = 1:r
  dat$r = t[as.character(dat$r)]
  
  
  m = mean(dat$s, na.rm=T)
  if(norm){
    z = quantile(1:nrow(dat), c(0.4, 0.6))
    m = mean(dat$s[z[1]:z[2]], na.rm=T)
    dat$s = dat$s / m 
    dat$s[dat$s > 2] = 2
    m = mean(dat$s[z[1]:z[2]], na.rm=T)
  }
  
  if(length(dat) == 4){
    names(dat)[4] = 'circ'
    dat.cs = dat[dat$circ < circ.cutoff & !is.na(dat$circ),] 
  }
  
  if(type == 'heatmap'){
    p <- ggplot(dat, aes(x = c, y = r, fill = s)) + 
      geom_tile(color='black') +
      scale_fill_gradient2(midpoint=m, low=low, high=high, mid=mid) + 
      scale_x_discrete(expand = c(0, 0), limits = as.character(1:c)) +
      scale_y_discrete(expand = c(0, 0), limits = as.character(r:1)) + 
      coord_equal() +
      theme_bw() +
      labs(list(title = title, x = "Column", y = "Row", fill = "Size")) + 
      theme(legend.position="right",title=element_text(size=14,face="bold"))
  }
  if(type == 'bubble'){
    p = ggplot(dat, aes(x = c, y = r)) + 
      geom_point(aes(x = c, y = r, size = s, colour = s),shape=16, alpha=0.80) +
      scale_colour_gradient(low=low, high=high) +
      scale_x_discrete(limits = as.character(1:c)) +
      scale_y_discrete(limits = as.character(r:1)) +
      coord_equal() +
      labs(list(title = title, x = "Column", y = "Row", color = "Size", size=""))+
      theme_bw() +
      theme(title=element_text(size=14,face="bold"))
  }
  if(show.circ & !is.null(dat.cs) & nrow(dat.cs) != 0){
    p = p + geom_point(aes(x=c, y=r), data=dat.cs, color=circ.color, size=1)
  }
  
  if(show.text) p = p + geom_text(color = text.color, aes(label=s), size=3)
  return(p)
}

summary.gitter <- function(d){
  pf = attr(d, 'format')
  call = attr(d, 'call')
  writeLines(sprintf('###########################\n# gitter V%s data file #\n###########################', GITTER_VERSION))
  writeLines(sprintf('Function call: %s', deparse(call)))
  writeLines(sprintf('Elapsed time: %s secs', attr(d, 'elapsed')))
  writeLines(sprintf('Plate format: %s × %s (%s)', pf[1], pf[2], prod(pf)))
  writeLines('Colony size statistics:')
  print(summary(d[[3]]))
  writeLines('Dat file (showing 6 rows):')
  print(head(d))
  
}
