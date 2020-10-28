#! /bin/bash
rm .libs/* && rmdir libs
rm *.o *.lo *.la cvDirectDemo_ls
export SUNDIALS_ROOT=/home/Utilisateur/sundials-3.1.0/cvode_inst
export C_INCLUDE_PATH=${SUNDIALS_ROOT}/include/:$C_INCLUDE_PATH
#/home/fredrik/Documents/github/modelica-pneumatic/cvDirect_step
#libtool --mode=compile gcc -g -O -c -o sunexportdata.o sunexportdata.cpp 
libtool --mode=compile gcc -g -O -c -o cvDirectDemo_ls.o cvDirectDemo_ls.c -I${SUNDIALS_ROOT}/include/
libtool --mode=link gcc -g -O -o cvDirectDemo_ls.exe .libs/cvDirectDemo_ls.o -I/home/Utilisateur/sundials-3.1.0/cvode_inst/include/ -L${SUNDIALS_ROOT}/lib/ -lsundials_nvecserial -lsundials_cvode -lsundials_sunmatrixdense -lsundials_sunlinsoldense -static-libstdc++ -lm
export LD_LIBRARY_PATH=.libs/:${SUNDIALS_ROOT}/lib/:${LD_LIBRARY_PATH}
#Results ls -al
#-rwxr-xr-x  1 fredrik fredrik 53584 Jan 28 13:20 cvDirectDemo_ls

