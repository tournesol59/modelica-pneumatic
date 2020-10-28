#include <string.h>
#include <cassert>
#include <iostream>
#include <list>


typedef struct dpair {
     double x; 
     double y;
  }* p_dpair;

typedef std::list<dpair> li_doubles;

class SUNEXPORT_DATA {
public:
  SUNEXPORT_DATA(int size, char* giveFileName) ;

  ~SUNEXPORT_DATA() ;

  bool read_extern_output(int dim, li_doubles &li_reals) ;

  bool exportToDisk(li_doubles &li_reals) ;

  int sizeevalsol;

  char FileName[14];

};

