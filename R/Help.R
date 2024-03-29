# Helper methods for gitter

.setContrast <-function(im, contrast=10){
  c = (100.0 + contrast) / 100.0
  im = ((im-0.5)*c)+0.5
  im[im<0] = 0
  im[im>1] = 1
  return(im)
}

.roundOdd <- function(x){
  x = floor(x)
  if(x %% 2 == 0) x = x+1
  x
}

.roundEven <- function(x){
  2.*round(x/2)
}

#Adds a padding to some matrix mat such that the padding is equal to the value of the nearest cell
#Inputs:
#  mat = matrix to which the padding is added
#  lvl = number of levels (rows/columns) of padding to be added
#	padding = type of padding on the matrix, zero will put zeros as borders, replicate will put the value of the nearest cell
.padmatrix <- function(mat, level, pad.val){
  pc = matrix(pad.val, nrow(mat), level)
  pr = matrix(pad.val, level, ncol(mat) + 2*level)
  t = cbind(pc, mat, pc)
  t = rbind(pr, t, pr)
  return(t)
}

.unpadmatrix <- function(mat, level){
  crem = c(1:level, ncol(mat):(ncol(mat)-level+1))
  rrem = c(1:level, nrow(mat):(nrow(mat)-level+1))
  return(mat[-rrem, -crem])
}

.getMeans <- function(im, low=0.25, high=0.75){
  r = nrow(im)
  c = ncol(im)
  
  d = im[ (low*r):(high*r) , (low*r):(high*r)]
  
  return( c(min(d), max(d)) )
}

# Global thresholding 
.findOptimalThreshold <- function(x, r=2, lim=c(0, 0.4), cap=0.2){
  x[x<lim[1]] = 0
  x[x>lim[2]] = 0
  t = round(mean(x), r)
  
  mu1 = mean(x[x>=t])
  mu2 = mean(x[x<t])
  cmu1 = mu1
  cmu2 = mu2
  i = 1
  while(TRUE){
    if(i > 1){
      cmu1 = mean(x[x>=t])
      cmu2 = mean(x[x<t])
    } 
    if( (i > 1 & cmu1 == mu1 & cmu2 == mu2) | t > 1){
      loginfo('Optimal threshold t = %s', t)
      if(t>cap){
        t = as.numeric(quantile(x, 0.90))
        loginfo('Optimal threshold too high, setting threshold to 90th percentile t = %s', t)
      }
      return(t)
      break;
    }else{
      mu1 = cmu1 
      mu2 = cmu2
      t = round( (mu1 + mu2)/2, r)
      loginfo('Iteration %s, threshold t = %s', i-1, t)
    }
    i = i+1
  }
}

# Rotation methods
# Removed now

.findSlope <- function(bw){
  i <- which(bw==1)
  V <- var(cbind( col(bw)[i], row(bw)[i] ))
  u <- eigen(V)$vectors
  return(u[2,1]/u[1,1])
}

# Circularity functions
# Shift the image in one direction
.s1 <- function(z) cbind(rep(0,nrow(z)), z[,-ncol(z)] )
.s2 <- function(z) cbind(z[,-1], rep(0,nrow(z)) )
.s3 <- function(z) rbind(rep(0,ncol(z)), z[-nrow(z),] )
.s4 <- function(z) rbind(z[-1,], rep(0,ncol(z)) )

.edge <- function(z) z & !(.s1(z)&.s2(z)&.s3(z)&.s4(z))

.perimeter <- function(z) {
  e <- .edge(z)
  ( 
    # horizontal and vertical segments
    sum( e & .s1(e) ) + sum( e & .s2(e) ) + sum( e & .s3(e) ) + sum( e & .s4(e) ) + 
      # diagonal segments
      sqrt(2)*( sum(e & .s1(.s3(e))) + sum(e & .s1(.s4(e))) + sum(e & .s2(.s3(e))) + sum(e & .s2(.s4(e))) )
  ) / 2  # Each segment was counted twice, once for each end
}

.area <- function(z) sum(z)
.circularity <- function(z) 4*pi*.area(z) / .perimeter(z)^2

# Visualization functions

# Drawing on image
.getColorImage <- function(im.grey){
  im.grey = array(im.grey, dim=c(nrow(im.grey), ncol(im.grey), 3))
  return(im.grey)
}

.drawRect <- function(rects, im, color='red', int = 1){
  col.rgb = col2rgb(color)[,1]/255
  if( is(im, 'matrix') )
    im = .getColorImage(im)
  
  for(i in 1:nrow(rects)){
    rect = as.numeric(rects[i,])
    
    im[rect[3]:rect[4],rect[1],1] = col.rgb[1]
    im[rect[3]:rect[4],rect[2],1] = col.rgb[1]
    im[rect[3],rect[1]:rect[2],1] = col.rgb[1]
    im[rect[4],rect[1]:rect[2],1] = col.rgb[1]

    im[rect[3]:rect[4],rect[1],2] = col.rgb[2]
    im[rect[3]:rect[4],rect[2],2] = col.rgb[2]
    im[rect[3],rect[1]:rect[2],2] = col.rgb[2]
    im[rect[4],rect[1]:rect[2],2] = col.rgb[2]
    
    im[rect[3]:rect[4],rect[1],3] = col.rgb[3]
    im[rect[3]:rect[4],rect[2],3] = col.rgb[3]
    im[rect[3],rect[1]:rect[2],3] = col.rgb[3]
    im[rect[4],rect[1]:rect[2],3] = col.rgb[3]
    
  }
  return(im)
}

.drawPeaks <- function(peaks.c, peaks.r, im, channel=1){
  if(length(dim(im)) == 2) im = .getColorImage(im)
  im[,peaks.c,channel] = 1
  im[peaks.r,,channel] = 1
  im[,peaks.c,-channel] = 0
  im[peaks.r,,-channel] = 0
  return(im)
}

.imageSpot <- function(x, title='Spot'){
  x = t(x)
  graphics::image(x, col= gray(0:8 / 8), xaxt = "n", yaxt='n')
}