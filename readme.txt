This R package should be installable using the command 
R INSTALL mgcv in the usual way. 
The R code should mostly be ok with Splus (see contents of R directory), 
although compiled code loading will differ a little (calling should be ok) 
but you'll need to compile the source code in the src directory. e.g. 
gcc -c -o mgcv.o mgcv.c gcv.c matrix.c and then load it with dyn.load("mgcv.o"). 
See mgcv.ps or mgcv.pdf for more detailed documentation.   
