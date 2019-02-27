#! /bin/bash
rm .libs/* && rmdir libs
rm *.o *.lo *.la
export SUNDIALS_ROOT=/usr/local/src/sundials-3.1.0/cvode_inst
export C_INCLUDE_PATH=${SUNDIALS_ROOT}/include/:$C_INCLUDE_PATH
libtool --mode=compile gcc -g -O -c simu_pipe_functions.c -o simu_pipe_functions.lo -I${SUNDIALS_ROOT}/include/
libtool --mode=compile gcc -g -O -c cvRoberts_dns.c
libtool --mode=link gcc -g -O -shared -o libsimu_pipe.la simu_pipe_functions.lo -rpath /usr/local/lib -lm
libtool --mode=link gcc -g -O -o cvRoberts_dns_dy cvRoberts_dns.lo .libs/libsimu_pipe.so ${SUNDIALS_ROOT}/lib/libsundials_cvode.so.3 ${SUNDIALS_ROOT}/lib/libsundials_nvecserial.so.3 ${SUNDIALS_ROOT}/lib/libsundials_sunmatrixdense.so.1 ${SUNDIALS_ROOT}/lib/libsundials_sunlinsoldense.so.1
export LD_LIBRARY_PATH=.libs/:${SUNDIALS_ROOT}/lib/:${LD_LIBRARY_PATH}

#Results: ls -al
#-rwxr-xr-x 1 fredrik fredrik   37152 Feb 27 02:53 cvRoberts_dns_dy

#At t = 9.9267e-01      y =  1.252777e+05    1.105531e+05
#iout:100000

#Final Statistics:
#nst = 3715765 nfe  = 10324811 nsetups = 4876005 nfeLS = 0      nje = 2712201
#nni = 10324807 ncfn = 1358051 netf = 0      nge = 3815765
