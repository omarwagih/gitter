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

.GITTER_VERSION = '1.0.2'
.methods = c('kmeans', 'tophat')
.pf = list('1536'=c(32,48),'768'=c(32,48),'384'=c(16,24),'96'=c(8,12))  

gitter.demo <- function(eg='single'){
  if(!eg %in% c('single', 'ref')) stop('Invalid example')
  
  if(eg == 'single'){
    f = system.file("extdata", "sample.jpg", package="gitter")
    dat = gitter(f)
    p <- gitter.plot(dat, title=sprintf('gitter v%s single image example', .GITTER_VERSION))
    print(p)
    browseURL(f)
    gitter.summary(dat)
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

gitter.batch <- function(image.files, ref.image.file=NULL, verbose='l', ...){
  f = 'gitter_failed_images'
  ff = list.files(pattern=f)
  if(length(ff) > 0){
    f = paste0(f, format(Sys.time(), "_%d-%m-%y_%H-%M-%S"))
  }
  failed.file = paste0(f, '.txt')
  
  # Set verbose
  logReset()
  if(verbose == 'l') addHandler(writeToConsole, formatter=.defaultFormat)
  
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
  is.ref = !is.null(ref.image.file)
  if(is.ref){
    loginfo('Processing reference image: %s', ref.image.file)
    r = gitter(ref.image.file, verbose=verbose, .is.ref=T, ..., )
    params = attr(r, 'params')
  }
  
  failed.plates = c()
  for(image.file in image.files){
    result = tryCatch({ gitter(image.file, .params=params, .is.ref=F, verbose=verbose,...) }, 
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


gitter <- function(image.file=file.choose(), plate.format=c(32,48), remove.noise=F, autorotate=F, 
                   inverse=F, verbose='l', contrast=NULL, fast=NULL, plot=F, grid.save=getwd(), 
                   dat.save=getwd(), .is.ref=F, .params=NULL){
  
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
  
  if(!verbose %in% c('l', 'p', 'n')) 
    stop('Invalid verbose parameter, use "l" to show the log, "p" to show a progress bar or "n" for no verbose')
  if(!is.null(contrast)){
    if(contrast <= 0) stop('Contrast value must be positive')
  } 
  if(grepl('^gridded', basename(image.file))) 
    warning('Detected gridded image as input')
  if(!is.null(grid.save)){
    if(!file.info(grid.save)$isdir) 
      stop(sprintf('Invalid gridded directory "%s"', grid.save))
  }
  if(!is.null(dat.save)){
    if(!file.info(dat.save)$isdir) 
      stop(sprintf('Invalid dat directory "%s"', grid.save))
  }
  is.fast = !is.null(fast)
  if(is.fast){
    if(fast < 1500 | fast > 4000) stop('Fast resize width must be between 1500-4000px')
  } 
   
  
  expf = 1.5
  params = as.list(environment(), all=TRUE)
  nrow = plate.format[1]
  ncol = plate.format[2]
  # Are we using a reference screen?
  is.ref = all(is.null(.params))
  ptm <- proc.time()
  
  # Set verbose
  logReset()
  if(verbose == 'l') addHandler(writeToConsole, formatter=.defaultFormat)
  
  prog = verbose == 'p' 
  if(!prog) pb <- NULL
  if(prog){
    cat(file.path(dirname(image.file), basename(image.file)), '\n')
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
  }
  # Read image
  loginfo('Reading image from: %s', image.file)
  if(prog) setTxtProgressBar(pb, 5)
  
  im = readJPEG(image.file)
  is.color = (length(dim(im)) == 3)
  
  if(is.fast){
    if(prog) setTxtProgressBar(pb, 7)
    loginfo('Resizing image...')
    im = resize(im, h=fast)
  }
  
  # Extract greyscale
  if(is.color){
    if(prog) setTxtProgressBar(pb, 9)
    loginfo('\tDetected color image, extracting greyscale')
    # Luminosity grey scale from GIMP
    im.grey = (im[,,1]*0.72) + (im[,,2]*0.21) + (im[,,3]*0.07)
  }else{
    loginfo('\tDetected greyscale image')
    im.grey = im
  }
  
  
  if(autorotate){
    if(prog) setTxtProgressBar(pb, 14)
    loginfo('Autorotating image...')
    im.grey = .autoRotateImage2(im.grey)
  }
  
  if(!is.null(contrast)){
    if(prog) setTxtProgressBar(pb, 15)
    loginfo('\tIncreasing image contrast with factor %s', contrast)
    im.grey = .setContrast(im.grey, contrast)
  }
  
  if(inverse){
    if(prog) setTxtProgressBar(pb, 17)
    loginfo('Inversing image')
    im.grey = 1 - im.grey 
  }
  
  if(!is.ref){
    if(prog) setTxtProgressBar(pb, 18)
    loginfo('Non-reference plate, registering image to reference')
    # Fix any shifts 
    im.grey = .register2d(im.grey, .params$row.sums, .params$col.sums)
  }
  
  if(prog) setTxtProgressBar(pb, 20)
  im.grey = .threshold(im.grey, nrow, ncol, fast=T, f=1000, pb)
  
  if(remove.noise){
    # Remove the noise
    if(prog) setTxtProgressBar(pb, 50)
    loginfo('Denoising image... ')
    kern = makeBrush(3, 'diamond')
    im.grey = openingGreyScale(im.grey, kern)
  }
  
  if(prog) setTxtProgressBar(pb, 55)
  #sum.y = rowSums(im.grey)
  loginfo('Computing row sums')
  sum.y = .rmRle(im.grey, p=0.2, 1)
  
  if(prog) setTxtProgressBar(pb, 60)
  loginfo('Computing column sums')
  sum.x = .rmRle(im.grey, p=0.2, 2)
  
  #sum.x = colSums(im.grey)

  if(is.ref){
    # Get peaks of sums
    if(plot) par(mfrow=c(2,1), bty='n', las=1)
    z = nrow*ncol
    loginfo('Getting row peaks...')
    if(prog) setTxtProgressBar(pb, 65)
    cp.y = .colonyPeaks(sum.y, n=nrow,z, plot)
    
    loginfo('Getting column peaks...')
    if(prog) setTxtProgressBar(pb, 70)
    cp.x = .colonyPeaks(sum.x, n=ncol,z, plot)
    
    if(prog) setTxtProgressBar(pb, 73)
    # Average window (though they should be the same)
    w = round(mean( c(cp.x$window, cp.y$window) ))
      
    # Get center coordinates
    coords = expand.grid(cp.x$peaks, cp.y$peaks)
    names(coords) = c('x', 'y')
    
    params$window = w
    params$coords = coords
    params$row.sums = rowSums(im.grey)
    params$col.sums = colSums(im.grey)
  }else{
    w = .params$window
    coords = .params$coords
  }
  
  
  # Window expansion factor
  d = round(w * expf)
  
  if(prog) setTxtProgressBar(pb, 75)
  #Pad image
  im.pad =  .padmatrix(im.grey, w, 1)
  
  coords[,c('x','y')] = coords[,c('x','y')]+w
  coords$xl = coords$x - d
  coords$xr = coords$x + d
  coords$yt = coords$y - d
  coords$yb = coords$y + d
  
  loginfo('Fitting bounds...')
  if(prog) setTxtProgressBar(pb, 80)
  coords = .fitRects(coords, im.pad, w)
  
  if(prog) setTxtProgressBar(pb, 84)
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
  
  if(prog) setTxtProgressBar(pb, 90)
  # Save gridded image
  if(!is.null(grid.save) & !.is.ref){
    #imr = drawPeaks(peaks.c=cp.y$peaks, peaks.r=cp.x$peaks, imr)
    imr = .drawRect(coords[,3:6], im.grey)
    save = file.path(grid.save, paste0('gridded_',basename(image.file)))
    loginfo('Saved gridded image to: %s', save)
    
    if(prog) setTxtProgressBar(pb, 93)
    writeJPEG(imr, save)
  }
  
  if(prog) setTxtProgressBar(pb, 98)
  # Save dat file
  if(!is.null(dat.save) & !.is.ref){
    save = file.path(dat.save, paste0(basename(image.file), '.dat'))
    
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
  
  if(prog) setTxtProgressBar(pb, 100)
  if(prog) close(pb)
  return(results)
}

# .getKmeansImage <- function(im, coords){
#   coords = coords[,c('xl', 'xr', 'yt', 'yb')]
#   cent = .getMeans(im)
#   ret = im
#   t = matrix(T, nrow(im), ncol(im))
#   for(i in 1:nrow(coords)){
#     z = as.numeric(coords[i,])
#     spot = im[z[3]:z[4],z[1]:z[2]]
#     im.spot.vec = c(cent, spot)
#     #cent = c(min(im.spot.vec), max(im.spot.vec))
#     cl = kmeans(im.spot.vec, cent)
#     cls = cl$cluster[-(1:2)] - 1
#     spot.bw = matrix(cls, nrow=nrow(spot))
#     ret[z[3]:z[4],z[1]:z[2]] = spot.bw
#     t[z[3]:z[4],z[1]:z[2]] = F
#   }
#   
#   ret[t] = ( ret[t] > 0.7 )+0
#   return(ret)
# }

.threshold <- function(im.grey, nrow, ncol, fast=T, f=1000, pb){
  prog = !is.null(pb)
  ptm <- proc.time()
  if(prog) setTxtProgressBar(pb, 22)
  if(fast){
    loginfo('Running fast background correction')
    im = resize(im.grey, h=f)
    si = round((nrow(im) / nrow) * 1.0)
    loginfo('Opening image with kernel size %s', si)
    kern = makeBrush(.roundOdd(si), 'box')
    if(prog) setTxtProgressBar(pb, 25)
    op = openingGreyScale(im, kern)
    if(prog) setTxtProgressBar(pb, 32)
    op = resize(op, w=nrow(im.grey), h=ncol(im.grey))
    if(prog) setTxtProgressBar(pb, 39)
    im.grey = im.grey - op
    im.grey[im.grey<0]=0
    # Find threshold
    if(prog) setTxtProgressBar(pb, 42)
    thresh = .findOptimalThreshold(resize(im.grey, h=f))
  }else{
    loginfo('Running background correction')
    si = round((nrow(im.grey) / nrow) * 1.0)
    loginfo('Opening image with kernel size %s', si)
    kern = makeBrush(.roundOdd(si), 'box')
    if(prog) setTxtProgressBar(pb, 32)
    im.grey = whiteTopHatGreyScale(im.grey, kern)
    # Find threshold
    if(prog) setTxtProgressBar(pb, 39)
    thresh = .findOptimalThreshold(im.grey)
  }
  
  if(prog) setTxtProgressBar(pb, 48)
  im.grey = (im.grey >= thresh)+0
  
  e = round( (proc.time() - ptm)[[3]], 2)
  loginfo('Thresholding took %s seconds', e)
  return(im.grey)
}

.rmRle <- function(im, p=0.2, margin=1){
  c = p * nrow(im)
  if(margin == 2){
    c = p * ncol(im)
  }
  z = apply(im, MARGIN=margin, function(x){
    r = rle(x)
    any(r$lengths[r$values == 1] > c)
  })
  
  if(margin==1) x = rowSums(im)
  if(margin==2) x = colSums(im)
  
  zh = .splitHalf(z)
  if(sum(zh$left) > 0) x[1:max( which(zh$left) )] = 0
  if(sum(zh$right) > 0) x[min( which(zh$right) + length(zh$left) ) : length(x)] = 0
  
  #x[z] = median(x)
  return(x)
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



gitter.plot <- function(dat, title='', type='heatmap', low='turquoise', mid='black', high='yellow', 
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

gitter.summary <- function(d){
  pf = attr(d, 'format')
  call = attr(d, 'call')
  writeLines(sprintf('###########################\n# gitter V%s data file #\n###########################', .GITTER_VERSION))
  writeLines(sprintf('Function call: %s', deparse(call)))
  writeLines(sprintf('Elapsed time: %s secs', attr(d, 'elapsed')))
  writeLines(sprintf('Plate format: %s x %s (%s)', pf[1], pf[2], prod(pf)))
  writeLines('Colony size statistics:')
  print(summary(d[[3]]))
  writeLines('Dat file (showing 6 rows):')
  print(head(d))
  
}
