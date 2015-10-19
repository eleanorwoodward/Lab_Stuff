## Noah Greenwald

# Notes on how to properly optimize code
# 
# First is to use the syste.time function to figure out how long a given program takes
# User time returns the amount of time the CPU is occupied, System time is the "wall clock" time

system.time(mean(c(1:1000000)))
system.time(sort(rnorm(100000, sd = 10)))

## Great for going line by line, and essentially printing your way to a solution
## However, if your program is too complicated to do line by line analysis of what is hogging all the time,
## the R Profiler can be of help.

output <- Rprof(run.exac(maf))

## The output of rprof is often not readable or useful, so we use the summary function. 
## There are two competing normalizations methods, bytotal or byself. Total is total amount
## of time that a function is present in the stack. Self is this amount, subtracted by the
## amount of time it is open calling other, low level functions. Usually top level functions
## spend their time calling low level functions, and simply act as wrappers. 

sum.stat <- summaryRprof(output)
sum.stat$by.total
sum.stat$by.self
sum.stat$sample.interval
sum.stat$sampling