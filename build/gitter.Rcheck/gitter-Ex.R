pkgname <- "gitter"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "gitter-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('gitter')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("gitter")
### * gitter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gitter
### Title: Process a single plate image
### Aliases: gitter
### Keywords: gitter image process sga single

### ** Examples

# Read sample image
f = system.file("extdata", "sample.jpg", package="gitter")
# Process it
dat = gitter(f)
# View head of the results
head(dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gitter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gitter.batch")
### * gitter.batch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gitter.batch
### Title: Process a batch set of plate images
### Aliases: gitter.batch
### Keywords: batch directory reference

### ** Examples

# Processing image using reference image
# This image would typically fail to process, since its missing several rows
f = system.file("extdata", "sample_dead.jpg", package="gitter")
# We will use this image to successfully process the above image
f.ref = system.file("extdata", "sample.jpg", package="gitter")
# Process
gitter.batch(f, f.ref)

# Remember: output files by default are saved to your working directory



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gitter.batch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gitter.demo")
### * gitter.demo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gitter.demo
### Title: Run a demo of gitter
### Aliases: gitter.demo

### ** Examples

# gitter.demo()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gitter.demo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gitter.read")
### * gitter.read

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gitter.read
### Title: Read in a data file as a 'gitter' data object.
### Aliases: gitter.read
### Keywords: dat file gitter read

### ** Examples

# Get dat file path
f = system.file("extdata", "sample.jpg.dat", package="gitter")
# Read in path as a gitter data object
g = gitter.read(f)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gitter.read", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plate.warnings")
### * plate.warnings

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plate.warnings
### Title: Show any plate-level warnings associated with a 'gitter' data
###   object
### Aliases: plate.warnings
### Keywords: error plate warning

### ** Examples

# dat = gitter("/path/to/image")
# plate.warnings(dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plate.warnings", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.gitter")
### * plot.gitter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.gitter
### Title: Plot a gitter dat file
### Aliases: plot.gitter
### Keywords: bubble display heatmap plot visualize

### ** Examples

f = system.file("extdata", "sample.jpg.dat", package="gitter")
# Read in path as a gitter data object
g = gitter.read(f)
# Plot a heatmap
plot(g, type="heatmap")
# Show a bubble plot
plot(g, type="bubble", low="black", high="red")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.gitter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.gitter")
### * summary.gitter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.gitter
### Title: View the summary of a gitter data file
### Aliases: summary.gitter
### Keywords: error plate warning

### ** Examples

# dat = gitter("/path/to/image")
# summary(dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.gitter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
