#! /bin/bash
rm .libs/* && rmdir libs
rm *.o *.lo *.la cvRoberts_dns_st
export SUNDIALS_ROOT=/usr/local/src/sundials-3.1.0/cvode_inst
export C_INCLUDE_PATH=${SUNDIALS_ROOT}/include/:$C_INCLUDE_PATH
libtool --mode=compile gcc -g -O -c simu_pipe_functions.c -o simu_pipe_functions.lo -I${SUNDIALS_ROOT}/include/
libtool --mode=compile gcc -g -O -c cvRoberts_dns.c
#export LD_LIBRARY_PATH=$SUNDIALS_ROOT/lib/:$LD_LIBRARY_PATH
#libtool --mode=link gcc -g -O -o libsimu_pipe.la simu_pipe_functions.lo -rpath /usr/local/lib -lm
ar cru libsimu_pipe.a  .libs/simu_pipe_functions.o
ranlib libsimu_pipe.a
gcc -g -O -o cvRoberts_dns_st -static .libs/cvRoberts_dns.o .libs/simu_pipe_functions.o ${SUNDIALS_ROOT}/lib/libsundials_cvode.a ${SUNDIALS_ROOT}/lib/libsundials_nvecserial.a ${SUNDIALS_ROOT}/lib/libsundials_sunmatrixdense.a ${SUNDIALS_ROOT}/lib/libsundials_sunlinsoldense.a -static-libstdc++ -lm

#Results ls -al 
#-rwxr-xr-x 1 fredrik fredrik 1140344 Feb 27 02:51 cvRoberts_dns_st

#At t = 9.9267e-01      y =  1.252777e+05    1.105531e+05
#iout:100000

#Final Statistics:
#nst = 3715765 nfe  = 10324811 nsetups = 4876005 nfeLS = 0      nje = 2712201
#nni = 10324807 ncfn = 1358051 netf = 0      nge = 3815765
 
