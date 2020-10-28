 /******************************************************
 * this must be temporary solution to export scalar data
 * series to disk
 *******************************************************/
#include "sunexportdata.hpp"
#include <stdio.h>
#include <cassert>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <malloc.h>

/********** class IDENT05_IODATA ************/
// constructor 
 SUNEXPORT_DATA::SUNEXPORT_DATA(int size, char* giveFileName) {

    sizeevalsol = size; // for exporting result value
    strncpy(FileName, giveFileName, 14);
 }

//destructor
 SUNEXPORT_DATA::~SUNEXPORT_DATA() {
 }
 
bool SUNEXPORT_DATA::read_extern_output(int dim, li_doubles &li_reals) {
// dim shall be width of spectral fft output: same length as time-val signal

   dpair data;
   double inpdata[2];
   char strInpName[14];

   std::ifstream lh_sundata; 
   strncpy(strInpName, FileName, 8);
   lh_sundata.open(strInpName, std::ifstream::in);  //simplest name from Matlab out.dat (w/o the extension) output
   std::string line;
   int row=1;
   int col;
  // int origin=128;
   data.x = 0.0;
   data.y = 0.0;
 
   while ( (std::getline(lh_sundata, line)) && (row<=256) ){
      std::stringstream ss(line);
      col=0;
      while (ss >> inpdata[col] ) col++; // real part [0] and imag part [1]
      data.x=(double) row;
      data.y=inpdata[0];
      li_reals.push_back(data);
      row++;
   }

   lh_sundata.close();
  return 1;
}


 bool SUNEXPORT_DATA::exportToDisk(li_doubles &li_reals)
{
	// in practice arrayToExport has 6 columns and two of these six shall be selected
  int row;
  dpair data;
  std::ofstream lh_test;
  char strOutName[14];

  strncpy(strOutName, FileName, 8);
  lh_test.open(strOutName, std::ofstream::out); //opening file
  row=0;
  for (auto it=li_reals.begin(); it != li_reals.end(); it++) {
      (data.x)=(*it).x;  // container access, not optimal
      (data.y)=(*it).y;  // container access, not optimal
      lh_test << data.x << "\t";
      lh_test << data.y << "\n";
      row++;
  }
  lh_test.close();

  return 1;  
}

