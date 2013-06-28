
.estimateWindow <- function(y){
  length(y) / (which.max(Re(fft(y-mean(y)))[1:(length(y)/2)]) - 1)
}

.getColonyPeaks2 <- function(x, n, plot=T){
  
  w = .estimateWindow(x) * 1.3
  
  loginfo('Estimated window size: %s', w)
  max.win = (length(x) / n)
  if(w > max.win){
    w = floor(max.win)
    loginfo('Window too large, changing to %s', w)
  }
  
  
  ref = sin( base::seq(0, pi, pi/w) )
  cr = sapply(1:length(x), function(i){
    s = i - (w/2)
    e = i + (w/2)
    if(s < 1) s = 1
    if(e > length(x)) e = length(x)
    
    z = x[s:e]
    cor(ref[1:length(z)], z, method="spearman")
  })
  
  cr[cr < 0.3] = 0
  
  rn = 1:length(cr)
  pe = .getPeaks(cr,floor( (w/1.3) /2))
  rp = rn[pe]
  pv = cr[pe]
  
  if(length(rp) < n)
    stop('Not enough peaks found')
  
  peak.dist = .shift(rp, 1) - rp
  peak.height = .shift(pv, 1) - pv
  
  num.s = (length(rp) - n) + 1
  P = sapply(1: num.s, function(s.i){
    dist = peak.dist[s.i:(s.i + (n-2))]
    w = peak.height[s.i:(s.i + (n-2))]
    sum(log( dnorm(dist, median(dist), sd(dist)) * dnorm(w, median(w), sd(w)) ) )
  })
  
  s = which.max(P)
  peaks = rp[s:(s + (n-1))]
  
  # Compute new delta
  delta = median( .shift(peaks, 1) - peaks, na.rm=T )
  
  if(plot){
    plot(cr, lwd=1, ylab='Sum of pixel intensities', xlab='Index', bty='o', type='l')
    abline(v=rp, col="red")
    abline(v=peaks, col='#2e756d', lwd=2)
    text(peaks, quantile(cr,.1), as.character(1:length(peaks)))
  }
  return(list(peaks=peaks, all.peaks=rp, window=delta/2))
  
}
# 
# .getColonyPeaks <- function(x, n, smoothover=5, plot=T){
#   
#   sm = filter(x,rep(1/smoothover,smoothover),circular=TRUE)
#   
#   smc = sm[ (0.25*length(sm)) : (0.75*length(sm)) ]
#   
#   #tvals = median(smc)
#   tvals = quantile(smc, probs=seq(0.3, 0.7, by=0.05))
#   z = lapply(tvals, function(thresh){
#     cross = which(diff(sign(smc-thresh))!=0)-1
#     c.diff = .shift(cross, 1) - cross
#   })
#   zsd = sapply(z, sd, na.rm=T)
#   z = z[[which.min(zsd)]]
#   thresh = tvals[which.min(zsd)]
#   
#   w = round(median(z, na.rm=T))
#   loginfo('Estimated window size: %s', w)
#   max.win = (length(x) / n)/2
#   if(w > max.win){
#     w = floor(max.win)
#     loginfo('Window too large, changing to %s', w)
#   }
#   
#   rn = 1:length(sm)
#   rp = rn[.getPeaks(sm,w)]
#   
#   if(length(rp) < n)
#     stop('Not enough peaks found')
#   
#   peak.dist = .shift(rp, 1) - rp
#   peak.width = .getPeakRange(rp, sm, w*2, plot=F)
#   peak.width = .shift(peak.width, 1) - peak.width
#   
#   num.s = (length(rp) - n) + 1
#   P = sapply(1: num.s, function(s.i){
#     dist = peak.dist[s.i:(s.i + (n-2))]
#     w = peak.width[s.i:(s.i + (n-2))]
#     sum(log( dnorm(dist, median(dist), sd(dist)) * dnorm(w, median(w), sd(w)) ) )
#   })
#   
#   s = which.max(P)
#   peaks = rp[s:(s + (n-1))]
#   peaks = .fixOuterPeaks(peaks)
#   
#   # Compute new delta
#   delta = median( .shift(peaks, 1) - peaks, na.rm=T )
#   
#   if(plot){
#     plot(sm, lwd=1, ylab='Sum of pixel intensities', xlab='Index', bty='o')
#     abline(h=thresh)
#     #rug(cross)
#     abline(v=rp, col="red")
#     abline(v=peaks, col='#2e756d', lwd=2)
#     text(peaks, quantile(sm,.1), as.character(1:length(peaks)))
#   }
#   return(list(peaks=peaks, all.peaks=rp, window=delta/2))
# }

.getWindow <- function(data, pos, window){
  win.left = window
  win.right = window
  if(pos - window < 0) win.left = pos
  if(pos + window > length(data)) win.right = length(data) - pos
  
  return(data[(pos-win.left):(pos+win.right)])
}

# Shift vector by an offset
.shift <- function(x, offset, na.pad=NA){
  t = (1:length(x)) + offset
  t[t<1] = na.pad
  return(x[ t ])
}


.getPeaks <- function(x, halfWindowSize, type="max") {
  if(type=="min")
    x = 1/x
  
  windowSize <- halfWindowSize * 2 + 1
  windows <- embed(x, windowSize)
  localMaxima <- max.col(windows, "first") == halfWindowSize + 1
  
  return(c(rep(FALSE, halfWindowSize), localMaxima, rep(FALSE, halfWindowSize)))
}

# .getPeakRange <- function(peak.coords, data, win, plot=T){
#   sds = sapply(peak.coords, function(p){
#     d = .getWindow(data, p, win/4)
#     diff(range(d))
#   })
#   return(sds)
# }



.splitHalf <- function(vec){
  t = ceiling(length(vec)/2)
  return(list(left=vec[1:t], right=vec[(t+1):length(vec)]))
}

.fixOuterPeaks <- function(peaks){
  n = length(peaks)
  left = seq(peaks) < n/2
  p1 = peaks[left]
  p2 = peaks[!left]
  
  dist = .shift(peaks, 1) - peaks
  d1 = .shift(p1, 1) - p1
  d2 = p2 - .shift(p2, -1)
  
  s = 4*sd(dist[(0.2*length(dist)):(0.8*length(dist))], na.rm=T)
  m = median(dist, na.rm=T)
  
  f = 1
  for(i in f){
    if(d1[i] < m-s | d1[i] > m+s)
      peaks[1] = round(peaks[2] - m)
    j = length(d2) - (i-1)
    if(d2[j] < m-s | d2[j] > m+s)
      peaks[length(peaks)] = round(peaks[n-i] + m)
  }
  return(peaks)
}
