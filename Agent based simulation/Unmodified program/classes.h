#include <stdio.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

#include <vector>
using namespace std;

double _drand48(void) ;
void _srand48(int a) ;
void err(char *reason) ;
void err(char *reason, int a) ;
void err(char *reason, double a);
void init();
void end() ;
void reset() ;

extern const int _resol, _bins ;

#ifndef classes_already_defined
#define classes_already_defined

#ifdef __linux
typedef unsigned int DWORD ;
typedef unsigned char BYTE ;
#else 
#include <windows.h>
#endif


struct Cell {
  short unsigned int well ;
  unsigned int gen ; 
  float x,y,z ; // position within the well: (0,0,0) = top centre
};

extern vector <Cell> cells ;

struct Well {
  float food ;
  float a ; 
};

//const unsigned int RESISTANT_PM = 1<<31 ;
const unsigned int DRIVER_PM = 1<<30 ;
const unsigned int L_PM = (1<<30) - 1 ; 

struct Genotype {
  vector <unsigned int> sequence ;
  BYTE no_drivers ;
  int number ; // number of cells of this genotype
  int prev_gen ;
  int index ; // this is set only when saving data
  Genotype(void) ;
  ~Genotype(void) { sequence.clear() ; }
  Genotype(Genotype *mother, int prevg, int no_snp) ;
  float growth(int well) ;
  float death(int well) ;
  float death_during_rep(int w) ;
};

class Hist {
  public:
  int x,n ;
  long long int x2 ;
  Hist() { x=n=0 ; x2=0 ; }
	inline void operator+= ( Hist& H2 ) { x+=H2.x ; x2+=H2.x2 ; n+=H2.n ; }
	inline void operator+= ( int val ) { x+=val ; x2+=val*val ; n++ ; }
  void r() { x=n=0 ; x2=0 ; }
};


#endif

extern vector<Genotype*> genotypes ;
