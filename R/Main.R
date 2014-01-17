# z = setwd('~/Development/gitter/R')
# #Import packages
# library('EBImage')
# library('jpeg')
# library('tiff')
# library('logging')
# library('parallel')
# library('PET')
# library('ggplot2')
# 
# source('Peaks.R')
# source('Help.R')
# setwd(z)


.GITTER_VERSION = '1.0.3'

.onAttach <- function(lib, pkg, ...) {
#   writeLines(sprintf('gitter version %s - quantification of pinned microbial cultures\n', .GITTER_VERSION))
#   writeLines('Copyright (C) 2014 Omar Wagih\n')
#   writeLines("Type 'gitter.demo()' for a demo, '?gitter' for help ")
#   writeLines("or see http://gitter.ccbr.utoronto.ca for more details")
  packageStartupMessage(sprintf("gitter version %s - quantification of pinned microbial cultures\n
Copyright (C) 2014 Omar Wagih\n
Type 'gitter.demo()' for a demo, '?gitter' for help 
or see http://gitter.ccbr.utoronto.ca for more details", .GITTER_VERSION))
}

.pf = list('1536'=c(32,48),'768'=c(32,48),'384'=c(16,24),'96'=c(8,12))

.readImage <- function(file){
  n = basename(file)
  if(grepl('*.jpg$|*.jpeg$', file, ignore.case=T)){
    return( readJPEG(file) )
  }else if(grepl('*.tiff$', file, ignore.case=T)){
    return( readTIFF(file) )
  }else{
    im = imageData( readImage(file) )
    d = dim(im)
    M = array(NA, dim=d[c(2,1,3)])
    M[,,1] = t(im[,,1])
    if(d[3] > 1){
      M[,,2] = t(im[,,2])
      M[,,3] = t(im[,,3])
    }
    return(M)
  }
  
}



gitter.demo <- function(eg=1){
  if(!eg %in% c(1,2)) stop('Invalid example, please use 1 for a single image or 2 to process an image using a reference image')
  
  if(! interactive() ) stop('Unable to run demo through non-interactive interface!')
  if(eg == 1){ # Process single image
    f = system.file("extdata", "sample.jpg", package="gitter")
    dat = gitter(f, verbose='p')
    p <- plot.gitter(dat, title=sprintf('gitter v%s single image example', .GITTER_VERSION))
    print(p)
    browseURL(f)
    browseURL(file.path(getwd(), paste0('gridded_', basename(f))))
    summary.gitter(dat)
  }
  if(eg == 2){ # Process using reference image
    f = system.file("extdata", "sample_dead.jpg", package="gitter")
    f.ref = system.file("extdata", "sample.jpg", package="gitter")
    gitter.batch(f, f.ref, verbose='p')
    browseURL(f)
    browseURL(file.path(getwd(), paste0('gridded_', basename(f))))
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
    image.files = list.files(image.files, pattern='*.jpg$|*.jpeg$|*.tiff$', full.names=T, ignore.case=T)
    image.files = image.files[!grepl('^gridded_', basename(image.files))] # ignore gridded images
    if(length(image.files) == 0) stop('No images with JPEG or JPG extension found. Images must be JPG format, please convert any non-JPG images to JPG')
  }
  
  z = sapply(image.files, file.exists)
  if(!all(z)) stop(sprintf('Files "%s" do not exist', paste0(image.files[!z], collapse=', ')))
  
  params = NULL
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
                  if(verbose == 'p') cat('\n')
                  e
                })
    
    # If we have an error
    if('error' %in% class(result)){
      failed.plates = c(failed.plates, basename(image.file))
    }
  }
  
  # Save failed plates
  if(length(failed.plates) > 0){
    failed.plates = c(sprintf('# gitter v%s failed images generated on %s', .GITTER_VERSION, format(Sys.time(), "%a %b %d %X %Y")),
                      failed.plates)
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
  
  ptm <- proc.time()
  
  # Set verbose
  logReset()
  if(verbose == 'l') addHandler(writeToConsole, formatter=.defaultFormat)
  
  prog = verbose == 'p' 
  if(!prog) pb <- NULL
  if(prog){
    cat('Processing', basename(image.file), '...\n')
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
  }
  # Read image
  loginfo('Reading image from: %s', image.file)
  if(prog) setTxtProgressBar(pb, 5)
  
  im = .readImage(image.file)
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
  
  # Are we using a reference screen?
  is.ref = is.null(.params)
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
  sum.y = .rmRle(im.grey, p=0.6, 1)
  
  if(prog) setTxtProgressBar(pb, 60)
  loginfo('Computing column sums')
  sum.x = .rmRle(im.grey, p=0.6, 2)
  
  #sum.x = colSums(im.grey)

  if(is.ref){
    # Get peaks of sums
    if(plot) par(mfrow=c(2,1), bty='n', las=1)
    z = nrow*ncol
    loginfo('Getting row peaks...')
    if(prog) setTxtProgressBar(pb, 65)
    cp.y = .colonyPeaks(sum.y, n=nrow, z, plot)
    
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
  d = .roundOdd(d)
  
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
  coords[,1:6][coords[,1:6] < 0] = 1
  coords[,c('xl', 'xr')][coords[,c('xl', 'xr')] > ncol(im.grey)] = ncol(im.grey)
  coords[,c('yt', 'yb')][coords[,c('yt', 'yb')] > nrow(im.grey)] = nrow(im.grey)
  
  
  # Generate rows and columns to be bound to coords below
  rc = expand.grid(1:ncol, 1:nrow)[,c(2,1)]
  names(rc) = c('row', 'col')
  results = cbind(rc, coords)
  
  # Print elapsed time
  elapsed = signif( (proc.time() - ptm)[['elapsed']], 5)
  
  # Rearrange
  results = results[,c('row', 'col', 'size', 'circularity','flags',
                       'x', 'y', 'xl','xr', 'yt', 'yb')]
  
  class(results) = c('gitter', 'data.frame')
  
  if(prog) setTxtProgressBar(pb, 90)
  # Save gridded image
  if(!is.null(grid.save) & !.is.ref){
    #imr = drawPeaks(peaks.c=cp.y$peaks, peaks.r=cp.x$peaks, imr)
    imr = .drawRect(coords[,3:6], im.grey)
    #imr = im.grey
    save = file.path(grid.save, paste0('gridded_',basename(image.file)))
    loginfo('Saved gridded image to: %s', save)
    
    if(prog) setTxtProgressBar(pb, 93)
    #writeJPEG(imr, save,1)
    writeJPEG(imr, save)
  }
  
  if(prog) setTxtProgressBar(pb, 98)
  
  # Save dat file
  if(!is.null(dat.save) & !.is.ref){
    save = file.path(dat.save, paste0(basename(image.file), '.dat'))
    
    results$circularity = round(results$circularity, 4)
    results = results[,1:5]
    
    results = .gitter.write(results, save)
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


.threshold <- function(im.grey, nrow, ncol, fast=T, f=1000, pb){
  prog = !is.null(pb)
  ptm <- proc.time()
  if(prog) setTxtProgressBar(pb, 22)
  if(fast){
    loginfo('Running fast background correction')
    im = resize(im.grey, h=f)
    si = round((nrow(im) / nrow) * 1.5)
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

.xl <- function(z, w){ 
  m = which(z == min(z))
  t = length(z) - m[length(m)]
  if(t < w) t = w
  return(t)
}

.xr <- function(z, w){ 
  t = which.min(z) + 1
  if(t < w) t = w
  return(t)
}
.yt <- .xl
.yb <- .xr


.matBorder <- function(x){
  list(l = x[1:nrow(x),1],
       r = x[1:nrow(x),ncol(x)],
       t = x[1,1:ncol(x)],
       b = x[nrow(x),1:ncol(x)])
}

.spilled <- function(x, f){
  z = .matBorder(x)
  sum( z$l ) > nrow(x)*f | 
    sum( z$r ) > nrow(x)*f | 
    sum( z$t ) > ncol(x)*f | 
    sum( z$b ) > ncol(x)*f
}

.fitRects <- function(coords, im.kmeans, d){
  
  # Minimum border for really small colonies
  # Remove any decimals
  minb = round(d/3)
  
  ret = lapply(1:nrow(coords), function(i){
    # Center of spot
    x = coords$x[i]
    y = coords$y[i]
    
    cent.pixel = im.kmeans[y,x]
    
    # Define expanded rectangle 
    rect = c(coords$xl[i], coords$xr[i], coords$yt[i], coords$yb[i])
    
    spot.bw = im.kmeans[rect[3]:rect[4], rect[1]:rect[2]]
    
    rs = rowSums(spot.bw)
    cs = colSums(spot.bw)
    
    # Crop bw spot and compute new size
    x.rel = x - rect[1]
    y.rel = y - rect[3]
    
    if(cent.pixel == 0){
      z = rep( minb*2, 4)
    }else{
      # Sum rows and cols, split sums in half 
      sp.y = .splitHalf(rs)
      sp.x = .splitHalf(cs)
      # Get minimums relative to the spot and set minimums if too small
      z = c(.xl(sp.x$left, minb), 
            .xr(sp.x$right, minb),
            .yt(sp.y$left, minb), 
            .yb(sp.y$right, minb))
    }
    
#     if(i == 48){
#       print(z)
#       print(c(x.rel, y.rel))
#       
#       par(mfrow=c(1,2))
#       plot(cs, type='l')
#       abline(v=c(x.rel-z[1], x.rel+z[2]), col='red')
#       plot(rs, type='l')
#       abline(v=c(y.rel-z[3], x.rel+z[4]), col='blue')
#     }
    
    rect.rel = c(x.rel - z[1], x.rel + z[2], y.rel - z[3], y.rel + z[4])
    spot.bw = spot.bw[rect.rel[3]:rect.rel[4], rect.rel[1]:rect.rel[2]]
    
    spilled = .spilled(spot.bw, 0.2)
    
    # Compute new spot relative to plate
    rect = c(x - z[1], x + z[2], y - z[3], y + z[4])
    
    c(x, y, rect, sum(spot.bw), .circularity(spot.bw), spilled)
  })
  
  coords = matrix(simplify2array(ret), nrow(coords), byrow=T)
  coords = as.data.frame(coords)
  names(coords) = c('x', 'y', 'xl', 'xr', 'yt', 'yb', 'size', 'circularity', 'spill')
  
  flags = data.frame(spilled = (coords$spill == 1), 
                     lowcirc = (coords$circularity < 0.6 & !is.nan(coords$circularity) & !is.na(coords$circularity)))
  # 1: Colony spill over or edge interference, 2 = Low circularity
  flag.id = c('S', 'C')
  coords$flags = apply(flags, 1, function(r){
    paste0(flag.id[r], collapse=',')
  })
  coords[,names(coords) != 'spill']
  
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


.getFlags <- function(dat){
  s = dat$size / median(dat$size)
  c = dat$circularity
  n = nrow(dat)
  
  fl = c()
  # If 10% of colonies have colonies smaller than 0.1
  if( (sum(s < 0.1) / n) > 0.1){
    fl = append(fl, "1")
  }
  # If empty 
  if( (sum(is.na(c) | c < 0.6) / n) > 0.1 ){
    fl = append(fl, "2")
  }
  return(fl)
}

.flagMap = c("1"="high count of small colony sizes",
             "2"="high count of low colony circularity")
.warningPat = "# Warning possible misgridding: "


.gitter.write <- function(dat, path){
  hd = c(sprintf('# gitter v%s data file generated on %s', .GITTER_VERSION, format(Sys.time(), "%a %b %d %X %Y")))
  id = .getFlags(dat)
  id = id[id %in% names(.flagMap)]
  fl = unname(.flagMap[id])
  if(length(fl) > 0){
    attr(dat, 'warnings') = fl
    fl = sprintf(paste0(.warningPat, "%s"), paste0(fl, collapse=', '))
    hd = append(hd, fl)
  }
  hd = append(hd, '# Flags: S - Colony spill or edge interference, C - Low colony circularity')
  writeLines(hd, path)
  cat('# ', file=path, append=T)
  
  suppressWarnings( write.table(dat, file=path, quote=F, sep='\t', row.names=F, col.names=T, append=T) )
  loginfo('Saved dat file to: %s', path)
  return(dat)
}


gitter.read <- function(path){
  if(is.character(path)){
    dat = read.table(path, stringsAsFactors=F, header=F, sep='\t')
    names(dat) = c('row', 'col', 'size', 'circularity', 'flags')
    
    #Read first 5 lines
    con  <- file(path, open = "r")
    i = 1
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      if(i > 5) break;
      if(grepl(pattern=.warningPat, line)){
        line = gsub(.warningPat, replacement="", line)
        x = strsplit(line, ", ")[[1]]
        x = x[x %in% .flagMap]
        attr(dat, "warnings") = x
      }
    } 
    close(con)
  }else if(!is.data.frame(path)){
    stop('Please enter the filepath to a gitter data file')
  }
  class(dat) = c('gitter', 'data.frame')
  return(dat)
}

plate.warnings <- function(dat){
  if(!is.data.frame(dat) & ! 'gitter' %in% class(dat)) stop('Argument must be a gitter data object')
  return( attr(dat, 'warnings') )
}

plot.gitter <- function(x, title='', type='heatmap', low='turquoise', mid='black', high='yellow', 
                        show.text=F, text.color='white', norm=T, 
                        show.flags=T, flag.color='white', ...){
  dat = x
#   if(is.character(dat)) dat = read.table(dat, stringsAsFactors=F, header=T)
  if(!is.data.frame(dat) & ! 'gitter' %in% class(dat)) stop('Argument must be a gitter data object')
  if(!type %in% c('heatmap', 'bubble')) stop('Invalid plot type. Use "heatmap" or "bubble"')
  
  if(length(dat) > 5 | length(dat) < 3) stop('Invalid number of columns for dat file')
  
  names(dat) = c('r', 'c', 's')
  dat.cs = s = NULL
  
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
  
  if(length(dat) == 5){
    names(dat)[5] = 'flags'
    dat.cs = dat[dat$flags != "",] 
  }
  
  # Round data to 2
  dat$s = round(dat$s, 2)
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
  if(show.flags & !is.null(dat.cs) & nrow(dat.cs) != 0){
    p = p + geom_point(aes(x=c, y=r), data=dat.cs, color=flag.color, size=1)
  }
  
  if(show.text) p = p + geom_text(color = text.color, aes(label=s), size=3)
  return(p)
}

summary.gitter <- function(object, ...){
  d = object
  pf = attr(d, 'format')
  call = attr(d, 'call')
  if(is.null(call)){
    call = "not available"
  }else{
    call = deparse(call)
  }
  
  writeLines(sprintf('# gitter v%s data file #', .GITTER_VERSION))
  writeLines(sprintf('Function call: %s', call))
  writeLines(sprintf('Elapsed time: %s secs', attr(d, 'elapsed')))
  writeLines(sprintf('Plate format: %s x %s (%s)', pf[1], pf[2], prod(pf)))
  writeLines('Colony size statistics:')
  print(summary(d[[3]]))
  writeLines('Dat file (first 6 rows):')
  print(head(d))
  w = attr(object, 'warnings')
  if(! is.null(w)){
    writeLines('Plate warnings:')
    writeLines(paste0(w, collapse=', '))
  }
}
