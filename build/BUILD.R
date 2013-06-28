setwd('~/Development/')

remove.packages('gitter', lib='/Library/Frameworks/R.framework/Versions/Current/Resources/library')

system('rm -rf gitter/build/gitter_1.0.1.tar.gz')
#system('rm -rf gitter/man/*')
system('R CMD build gitter')

system('mv gitter_1.0.1.tar.gz gitter/build/gitter_1.0.1.tar.gz')
#install.packages('gitter_1.0.tar.gz', repos = NULL, type="source")

detach("package:gitter", unload=TRUE)