setwd('~/Development/')
require(devtools)
remove.packages('gitter')

targz = 'gitter_1.1.2.tar.gz'
system(sprintf('rm -rf gitter/build/%s', targz))

# Regenerate Rwd files
document('gitter/')

R_PATH = '/usr/local/bin/R'
# Build 
system(sprintf('%s CMD build gitter', R_PATH))

# Move into build directory
system(sprintf('mv %s gitter/build/%s', targz, targz))

# Install 
system(sprintf('%s CMD INSTALL gitter/build/%s', R_PATH, targz))

#detach("package:gitter", unload=TRUE)
