//
// Program to try and produce the "non-mixed" growth curve using the well-mixed growth data
// Here we model nutrient molecules and their diffusion explicitly
//
//

// v.1.0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "classes.h"
#define __MAIN



//------------parameters-----------------------------------------------------------------------------------
const int nwells=1 ;                // number of compartments (wells)
int wellpops[1];               // no. of cells in each well 


// TESTING DIFF IN SIZE BETWEEN THIS AND diff.cpp
const float w=3.25,  h=11.0,   hole_width=1.5,   hole_depth=3.75;      // size of the well [mm]. Square hole at the top y=(-hole/2,hole/2),z=(-hole,0)


const int L=10000 ;                  // total genome length (it does not have to be realistically large (4.5Mbp for E. coli) 
                                     // but it needs to be larger than the expected number of mutations per cell at the end of the run

const int max_no_drivers=0 ;         // max. no of mutations that cause resistance (driver mutations)
const float gama=0 ;                 // mutation probability per daughter cell. 1e-2 is slightly higher than the base mutation rate (about 2e-3) in the absence of Cipro
                                     // but the actual mutation rate may be even higher in the presence of Cipro
//#### NOT 100% SURE WHAT GAMA IS NOW IVE MODIFIED poisson()

const float diff = 0.8 ;               // bacteria diffusion constant in [mm^2/h] ; typical value for E. coli is 1.5mm^2/h
const float max_growth_rate = 1.0 ;  // growth rate 1/h - 2 is max growth rate for E. coli
const float tmax=20;               // duration of experiment in hours

const float substr_diff = 0.2;     // diffusion constant of nutrient molecules


const float food0 = 1e5 ;          // amount of food per well, bacteria consume food with max rate 1/h. This is unrealistically low at the moment, 
                                     // but I suspect increasing it to a realistic value of 1e7 - 1e8 will change the number of accumulated mutations only as log(food0)
                                     
const float init_bact= 100 ;          // initial number of cells in the first well
const float K = food0+init_bact ;                // carrying capacity. 


const float uni_antibiotic=0, max_antibiotic=0 ;   // Cipro concentrations in ng/ml


// Local volume with which to determine local density 
const float volume_local = 0.1 ;  // small cube to calculate local density mm^3
const float dist_local = pow( volume_local, (1.0/3.0) ) / 2.0 ;   // half the width of the cube
const float volume_well = w * w * h ;


//---------------------------------------------------------------------------------------------------------


//  Well-mixed LB growth data   (length of arrays = 128, hardcoded into growth_function() method)
const float NKs[128] = {0.0, 0.007874015748031496, 0.015748031496062992, 0.023622047244094484, 0.031496062992125984, 0.03937007874015748, 0.047244094488188976, 0.055118110236220465, 0.06299212598425197, 0.07086614173228345, 0.07874015748031495, 0.08661417322834646, 0.09448818897637794, 0.10236220472440945, 0.11023622047244094, 0.11811023622047245, 0.12598425196850394, 0.13385826771653545, 0.14173228346456693, 0.1496062992125984, 0.15748031496062992, 0.16535433070866143, 0.1732283464566929, 0.1811023622047244, 0.1889763779527559, 0.19685039370078738, 0.2047244094488189, 0.21259842519685038, 0.2204724409448819, 0.2283464566929134, 0.23622047244094485, 0.24409448818897636, 0.25196850393700787, 0.25984251968503935, 0.2677165354330709, 0.2755905511811024, 0.28346456692913385, 0.29133858267716534, 0.2992125984251968, 0.30708661417322836, 0.31496062992125984, 0.3228346456692913, 0.33070866141732286, 0.3385826771653543, 0.3464566929133858, 0.3543307086614173, 0.3622047244094488, 0.3700787401574803, 0.3779527559055118, 0.3858267716535433, 0.39370078740157477, 0.40157480314960625, 0.4094488188976378, 0.41732283464566927, 0.42519685039370075, 0.4330708661417323, 0.4409448818897638, 0.4488188976377953, 0.4566929133858268, 0.4645669291338582, 0.4724409448818897, 0.4803149606299212, 0.4881889763779527, 0.4960629921259842, 0.5039370078740157, 0.5118110236220472, 0.5196850393700787, 0.5275590551181102, 0.5354330708661418, 0.5433070866141733, 0.5511811023622047, 0.5590551181102362, 0.5669291338582676, 0.5748031496062992, 0.5826771653543307, 0.5905511811023622, 0.5984251968503936, 0.6062992125984251, 0.6141732283464567, 0.6220472440944882, 0.6299212598425197, 0.6377952755905512, 0.6456692913385826, 0.6535433070866141, 0.6614173228346457, 0.6692913385826771, 0.6771653543307086, 0.68503937007874, 0.6929133858267716, 0.7007874015748031, 0.7086614173228346, 0.7165354330708661, 0.7244094488188976, 0.7322834645669292, 0.7401574803149606, 0.7480314960629921, 0.7559055118110236, 0.763779527559055, 0.7716535433070866, 0.7795275590551181, 0.7874015748031495, 0.795275590551181, 0.8031496062992125, 0.8110236220472441, 0.8188976377952756, 0.8267716535433071, 0.8346456692913385, 0.84251968503937, 0.8503937007874015, 0.8582677165354331, 0.8661417322834646, 0.8740157480314961, 0.8818897637795275, 0.889763779527559, 0.8976377952755906, 0.9055118110236221, 0.9133858267716536, 0.921259842519685, 0.9291338582677164, 0.9370078740157479, 0.9448818897637794, 0.9527559055118109, 0.9606299212598424, 0.968503937007874, 0.9763779527559054, 0.9842519685039369, 0.9921259842519684, 1.0};


const float LB_growth_rates[128] = {2.0, 1.9031812652922926, 1.8265813245206863, 1.7635011585726663, 1.7080276915410892, 1.6575787335152201, 1.6107631120813721, 1.5667325089017647, 1.5249254971440351, 1.4849481927575463, 1.4465118873947458, 1.4093976207682914, 1.3734347263016986, 1.3384871628951769, 1.304444440499734, 1.2712153858195334, 1.238723733808279, 1.206904932381468, 1.1757037765591607, 1.145072623869516, 1.1149700260600632, 1.0853596647963906, 1.0562095132181755, 1.0274911679677254, 0.9991793117650243, 0.9712512773165316, 0.9436866909039032, 0.9164671794225827, 0.8895761285957281, 0.8629984830206581, 0.8367205809152962, 0.8107300181301813, 0.7850155373257476, 0.759566939292014, 0.7343750142888966, 0.7094314920727403, 0.6847290099999233, 0.6602610993075975, 0.6360221904088813, 0.6120076388518786, 0.5882137745314704, 0.5646379778734512, 0.5412787881115484, 0.5181360505507084, 0.49521111198630474, 0.4725070763991393, 0.45002913688961355, 0.4277850048273143, 0.4057854637075078, 0.3840450835933291, 0.36258314261778263, 0.34142481494728866, 0.32060269937530084, 0.30015877734827967, 0.28014689852944674, 0.2606358841920315, 0.24171328936201583, 0.2234897268054033, 0.2061033478628293, 0.1897234746466773, 0.17455135855090048, 0.16081462820839804, 0.14875075819210107, 0.13857550412176903, 0.1304373177875507, 0.12436915329126161, 0.1202590219902531, 0.11785827832197618, 0.11682759738438099, 0.11679950145078165, 0.1174317875622506, 0.11843837538581482, 0.1195985124631672, 0.12075218641205952, 0.12178951986081657, 0.12263914258364501, 0.12325790918328422, 0.12362269555025603, 0.12372421843635174, 0.12356253612704095, 0.12314384303835355, 0.12247821951601835, 0.12157807025434818, 0.12045705251905794, 0.1191293502439622, 0.11760919156398184, 0.11591053760389598, 0.11404689198760082, 0.11203119587626333, 0.10987578419480296, 0.10759238640493368, 0.10519216069489586, 0.1026857544670921, 0.10008338700039136, 0.09739495248785307, 0.09463014354642953, 0.09179859693260198, 0.08891006468533032, 0.08597461530169025, 0.08300287079586806, 0.08000628643797443, 0.0769974802436361, 0.07399061816285943, 0.07100185709096928, 0.06804983908292031, 0.06515621298973427, 0.06234612908644776, 0.05964860188114751, 0.05709656097053268, 0.05472631139323961, 0.052576024270212025, 0.05068283172383112, 0.04907820613542173, 0.04778167233964922, 0.046793548427132875, 0.046088117628753186, 0.04560893369670621, 0.04526743210865565, 0.04494472450957251, 0.0444950261864247, 0.043748273227268306, 0.042509190458121474, 0.040549655441472016, 0.037589142323385866, 0.033250582060310126, 0.02695134980973683, 0.01754722247245933, 0.0};






//-----------------------------------------------------------------------------------------------------------


//what to determine experimentally specifically for MG1655:
//- no. of bacteria per well in 384 well plate with 0.15% agarose
//- birth vs food
//- migration rate
//- mutation rate vs CIPR concentration


char *NUM ; 
int RAND ;
double tt=0, tt_at_start ;
int start_clock ;
int sample=0 ;
vector<Cell> cells ;

// Model food molecules explicitly
vector<Cell> foodmols;

FILE *times ; 
char *timesbuffer ;
Well *well ; 




//----------------------------------- OPERATING SYSTEM STUFF ----------------------------------------------
#if defined __linux
#include <unistd.h>    
typedef unsigned int DWORD ;
int memory_taken() // return memory available in MB
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available memory in MB
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#elif defined __APPLE__
typedef unsigned int DWORD ;
int memory_taken()
{
  return 0 ; // not implemented
}
#else
#include <windows.h>
#include <psapi.h>
int memory_taken() // return memory available in MB
{
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (int) (info.WorkingSetSize/(1<<20));
}
#endif

void err(char *reason)
{
  cout <<reason<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif  
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, char *a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}
//----------------------------------- OPERATING SYSTEM STUFF ---------------------------------------------






static long long unsigned int _x=0x000100010001LL, _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  // works only on compilers with long long int!
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}



void _srand48(int a) { _x=a ; }




// Modified to use half the mutation rate to determine #SNPs in each daughter cell. ############   IS THIS RIGHT??
int poisson(float g)  // generates k from P(k)=exp(-gamma) gamma^k / k!
{
  const double l=exp(-(g/2.0)) ;
  double p=1. ;
  int k=0 ;
  do {
    k++ ;
    p*=_drand48() ;
  } while (p > l) ;
  return k - 1 ;
}





double gauss()
{
 static int iset=0;
 static float gset;
 float gasdev,v1,v2,r,fac;
 if (iset==0)
 {
  do
  {
	v1=2*_drand48()-1.0 ; v2=2*_drand48()-1.0 ;
	r=v1*v1+v2*v2;
  } while (r>1) ;
  fac=sqrt(-2.0*log(r)/r) ;
  gset=v1*fac ; gasdev=v2*fac ; iset=1 ;
 }
 else
 {
  gasdev=gset ; iset=0 ;
 }
 return (gasdev) ;
}




float Genotype::growth(int well) 
{
  return 1.0 ;  
//Barteks exp. data already has max at 2 contained in growth_function()
}



float Genotype::death(int w) 
{
  return 0 ;
}



float Genotype::death_during_rep(int w)
{
  float d=0 ;
  float mic=20 ; if (no_drivers>0) mic=100 ;   // MIC of the wild type is 20ng/ml, MIC of the mutant is assumed to be 100ng/ml. 
                                               // In reality there is a spectrum of mutants with different MICs
//  mic+=20*no_drivers ;
  
  d+=SQR(well[w].a/mic) ;
  if (d>1) d=1 ;
  return d ;  
}




// Modified to incorporate Bartek's experimental LB-growth data mapped to LOCAL SUBSTRATE CONCENTRATION
float growth_function( float converted_loc_dens )
{
  // There is now no food, just a growth rate that depends on N/K
  // where N is #cells in well, K is carrying capacity of well
  float g=-9.0 ;


  float NvK = converted_loc_dens ;


  int entry = 0;
  for( int i=0; i<128; i++ ){
    if( NKs[i] > NvK ){ 
      entry = i ;
      i = 200;  //break loop??
    }
  }

  if( NvK <= 0 ){ g = max_growth_rate ; }
  else if( NvK >= 1 ) { g = 0; }   // >= since we have no strict upper limit due to migration between wells (... implies unlimited "food" now)
  else{

    if( entry==0 ) {
       cout<<"NvK="<<NvK<<"  entry="<<entry<<" "<<endl;
       for( int gbg=0; gbg<nwells; gbg++ ){
       cout<<wellpops[gbg]<< " ";
       }
       cout<<endl;
       cout << "ERROR: ENTRY=0 IN growth_function"; exit (EXIT_FAILURE) ; 
    }

    g = (LB_growth_rates[entry-1]+LB_growth_rates[entry]) / 2.0 ;

    if( g<0 ) { 
        cout<<"ERROR, NEGATIVE GROWTH RATE in growth_function()"<<endl;  
        cout<<"g = "<<g<<"  entry = "<<entry<<"  NvK = "<<NvK<<endl;
        exit (EXIT_FAILURE) ;
    }

  }

   return g ;
}







Genotype::Genotype(void) 
{ 
  number=1 ; no_drivers=0 ; sequence.clear() ; prev_gen=-1 ;
}






Genotype::Genotype(Genotype *mother, int prevg, int no_snp) {
//
// This method still never checks whether the generated genotype already exists though, we just create a new entry every time.

  prev_gen=prevg ;  // i.e. prev_gen of this NEW genotype is the current genotype, cell[i].gen, fed into the method  
  sequence=mother->sequence ; 
  no_drivers=mother->no_drivers; 
  number=1;                     // but genotype might already exist (not important for data we are outputting)
  int snps_implemented = 0 ;
  while ( snps_implemented < no_snp ){
    
    int k ;
    int j= int(_drand48()*L) ; 

    for (k=0; k<sequence.size(); k++) if (sequence[k]==j) break ;      // code ignores back-mutations

    if (k==sequence.size()) {
      sequence.push_back(j) ;
      if (j<max_no_drivers) no_drivers++ ;
      snps_implemented++ ; 
    }

  }// -------- CAREFUL - DONT WANT INFINITE LOOP WHEN GENOME IS FULL. (only problem if we don't set L large enough)


}















vector<Genotype*> genotypes ;





// --------------------------------  SAVE AND SORT METHODS --------------------------------------------------



void quicksort2(float *n, int *nums, int lower, int upper)
{
	int i, m, temp ;
  float pivot, tempd;
	
	if (lower < upper)
	{
		SWAPD(n[lower], n[(upper + lower) / 2]); SWAP(nums[lower], nums[(upper + lower) / 2]);
		pivot = n[lower];
		m = lower;
		for (i = lower + 1; i <= upper; i++)
			if (n[i] > pivot)
			{
				m++;
				SWAPD(n[m], n[i]); SWAP(nums[m], nums[i]);
			}
		SWAPD(n[lower], n[m]); SWAP(nums[lower], nums[m]);
		quicksort2(n, nums, lower, m - 1);
		quicksort2(n, nums, m + 1, upper);
	}
}








//######################

void save_well_pops()
{
  int *nperw=new int[nwells] ;
  int i,j;
  for (i=0;i<nwells;i++) nperw[i]=0 ;
  for (i=0;i<cells.size();i++) nperw[cells[i].well]++ ;
  char name[256] ;
  sprintf(name,"%s/nwell.dat",NUM) ;
  FILE *f=fopen(name,"a") ;
  fprintf(f, "%lf ", tt);
  for (i=0;i<nwells;i++) fprintf(f,"%d ",nperw[i]) ; fprintf(f,"\n");
  fclose(f) ;
  delete [] nperw ;

}


//####################



//######################

void save_ODs()
{
   //  print out ODs as measured in diff.cpp / experiments
   // do 2 cross-sections (both cubiods) 1mm^2 and 0.25mm^2
  //int *nperw=new int[nwells] ;
  //int *ODs_perw_1 = new int[nwells] ;
  int *ODs_perw_2 = new int[nwells] ;

  for (int i=0;i<nwells;i++){
    //ODs_perw_1[i] = 0 ; 
    ODs_perw_2[i] = 0;
  }
  //for (i=0;i<nwells;i++) nperw[i]=0;
  //for (i=0;i<cells.size();i++) nperw[cells[i].well]++ ;

  for (int j=0; j<cells.size(); j++){
    //if( fabs(cells[j].x)<0.125 && fabs(cells[j].y)<0.125 ) {
    //    ODs_perw_1[ cells[j].well ]++ ;
    //}
    if( fabs(cells[j].x)<0.5 && fabs(cells[j].y)<0.5 ) {
        ODs_perw_2[ cells[j].well ]++ ;
    }
  }

  //char name[256] ;
  //sprintf(name,"%s/ODs_wells_0.25.dat",NUM) ;
  //FILE *f=fopen(name,"a") ;
  //fprintf(f, "%lf ", tt);
  //for (int z=0;z<nwells;z++) fprintf(f,"%d ", ODs_perw_1[z] ) ; fprintf(f,"\n");
  //fclose(f) ;

  char name2[256] ;
  sprintf(name2,"%s/ODs_wells_1.0.dat",NUM) ;
  FILE *f2=fopen(name2,"a") ;
  fprintf(f2, "%lf ", tt);
  for (int u=0;u<nwells;u++) fprintf(f2,"%d ", ODs_perw_2[u] ) ; fprintf(f2,"\n");
  fclose(f2) ;

  //delete [] nperw ;
  //delete [] ODs_perw_1 ;
  delete [] ODs_perw_2 ; 
}

//####################








void save_data()
{
  int *nperw=new int[nwells] ;
  int i,j;
  for (i=0;i<nwells;i++) nperw[i]=0 ;
  for (i=0;i<cells.size();i++) nperw[cells[i].well]++ ;
  char name[256] ;
  sprintf(name,"%s/state.dat",NUM) ;
  FILE *f=fopen(name,"a") ;
  for (i=0;i<nwells;i++) fprintf(f,"%d ",nperw[i]) ; fprintf(f,"\n");
  fclose(f) ;
  delete [] nperw ;

  float cutoff=0.1 ;
  int ntot=cells.size() ;
  int cells_drv=0 ;
  double drv_per_cell=0, pms_per_cell=0 ;

  int *snp_no=new int[L] ; // array of SNPs abundances
  for (i=0;i<L;i++) { snp_no[i]=0 ; }
  for (i=0;i<genotypes.size();i++) {
    if (genotypes[i]!=NULL && genotypes[i]->number>0) {
      for (int j=0;j<genotypes[i]->sequence.size();j++) snp_no[(genotypes[i]->sequence[j])]+=genotypes[i]->number ;      
    }
  }  
  int snps_det=0 ;
  for (i=0;i<L;i++) if (snp_no[i]>cutoff*ntot) snps_det++ ;
  delete [] snp_no ;

  for (i=0;i<ntot;i++) {
    Genotype *g=genotypes[cells[i].gen] ; if (g==NULL) err("g=NULL)") ;
    pms_per_cell+=g->sequence.size() ;
    if (g->no_drivers>0) {
      cells_drv++ ; drv_per_cell+=g->no_drivers ; 
    }
  }
  drv_per_cell/=ntot ; pms_per_cell/=ntot ;

  fprintf(times,"%lf %d %d  ",tt,ntot,genotypes.size()) ; //sqrt(raver2-raver*raver)) ;
  fprintf(times,"%d  %lf ",cells_drv,drv_per_cell) ;
  fprintf(times,"%d %f\n",memory_taken(),float(1.*(clock()-start_clock)/CLOCKS_PER_SEC)) ;
  fflush(times) ; 
}











void save_positions(char *name)
{
  FILE *data=fopen(name,"w") ;
  int di=cells.size()/100000 ; if (di<1) di=1 ; // do not save all cells if too many
  for (int i=0;i<cells.size();i+=di) {
    fprintf(data,"%d %f %f %f\n",i,cells[i].x+(w+1)*cells[i].well,cells[i].y,cells[i].z) ; 
  }
  fclose(data) ;  
}
// -----------------------------------------  SAVE AND SORT METHODS ----------------------------------------------------







int where_is_bug(float x, float y, float z)
{
  if (x>=-w/2 && x<w/2 && y>=-w/2 && y<w/2 && z<=0 && z>-h) return 0 ; // inside the well
  if (x>=w/2 && fabs(y)<hole_width/2 && z<=0 && z>-hole_depth) return 1 ; // going to the right well
  if (x<-w/2 && fabs(y)<hole_width/2 && z<=0 && z>-hole_depth) return 2 ; // going to the left well  
  return -1 ; // outside the well
}






int main_proc()
{
  int i,j,k,n,mol,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;
  double tt_old_2=tt ; 


  for(int it=0;;it++) {      // main loop 
#ifdef PAUSE_WHEN_MEMORY_LOW
    timeout++ ; if (timeout>1000000) {
      timeout=0 ; //printf("xxx") ;
      while (freemem()<PAUSE_WHEN_MEMORY_LOW) { sleep(1) ; } //printf("%d\n",freemem()) ; }
    }    
#endif


    double max_death_rate=max_growth_rate ;

    double tot_rate=(max_growth_rate+max_death_rate) ;

    double dt=-log(1-_drand48())/tot_rate ;

    tt +=  dt/( cells.size()+foodmols.size()  ) ; 


    float prob = _drand48();


    if( prob < ( 1.0*foodmols.size()/ (1.0*foodmols.size()+cells.size())  ) ){
      // Choose and move food molecule

      mol=_drand48()*foodmols.size();

        // DIFFUSION OF FOOD MOLECULE 
        float dx,dy,dz ;
        int ok ;
        do { 
          float sigma=sqrt(  substr_diff*dt ) ;
          dx=gauss()*sigma ;
          dy=gauss()*sigma ;
          dz=gauss()*sigma ;
          ok=where_is_bug( foodmols[mol].x+dx, foodmols[mol].y+dy, foodmols[mol].z+dz) ;
        } while (ok<0) ;
        if (ok==0) { // cell within the same well
        foodmols[mol].x+=dx ; foodmols[mol].y+=dy ; foodmols[mol].z+=dz ;
        } else if (ok==1 && foodmols[mol].well<nwells-1) { // cell in another well to the right
          foodmols[mol].x+=dx-w ; foodmols[mol].y+=dy ; foodmols[mol].z+=dz ;
          foodmols[mol].well++ ;     
        } else if (ok==2 && foodmols[mol].well>0) { // cell in another well to the left
          foodmols[mol].x+=dx+w ; foodmols[mol].y+=dy ; foodmols[mol].z+=dz ;      
          foodmols[mol].well-- ;  
        }
    }
    else{
      // Choose, move, reproduce cell

      n=_drand48()*cells.size() ;

      int foodmoltoremove;

      // migration
        float dx,dy,dz ;
        int ok ;
        do { 
          float sigma=sqrt(diff*dt) ;
          dx=gauss()*sigma ;
          dy=gauss()*sigma ;
          dz=gauss()*sigma ;
          ok=where_is_bug(cells[n].x+dx,cells[n].y+dy,cells[n].z+dz) ;
        } while (ok<0) ;
        if (ok==0) { // cell within the same well
        cells[n].x+=dx ; cells[n].y+=dy ; cells[n].z+=dz ;
        } else if (ok==1 && cells[n].well<nwells-1) { // cell in another well to the right
          cells[n].x+=dx-w ; cells[n].y+=dy ; cells[n].z+=dz ;
          cells[n].well++ ;     
          // track wellpops:
          wellpops[ cells[n].well ]++;  wellpops[ cells[n].well-1 ]--;
        } else if (ok==2 && cells[n].well>0) { // cell in another well to the left
          cells[n].x+=dx+w ; cells[n].y+=dy ; cells[n].z+=dz ;      
          cells[n].well-- ;  
          // track movement
          wellpops[ cells[n].well ]++;  wellpops[ cells[n].well+1 ]--;
        }


        double q=_drand48()*(max_growth_rate+max_death_rate), br,dr ;
         int mode=0 ; 


    //  We use the "local density" to determine growth rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        int cur_well = cells[n].well ; 
        float cur_x = cells[n].x ; float cur_y = cells[n].y ; float cur_z = cells[n].z ;

       // Get local density of food molecules
       int count_loc = 0;
       bool pickedmol = false;
       for( int fm=0; fm<foodmols.size(); fm++ ){

            if( foodmols[fm].well==cur_well ){
                 if( foodmols[fm].x < cur_x+dist_local  &&  foodmols[fm].x > cur_x-dist_local ){
                    if( foodmols[fm].y < cur_y+dist_local  &&  foodmols[fm].y > cur_y-dist_local ){
                        if( foodmols[fm].z < cur_z+dist_local  &&  foodmols[fm].z > cur_z-dist_local ){
                            count_loc++ ; 
                            if(!pickedmol){
                                // just remove first food molecule we find in neighbourhood; still random
                                foodmoltoremove = fm;
                                pickedmol = true;
                            }
                        }
                    }
                }
            }
        }


        // find correct local volume (i.e. ONLY INCLUDING WELL INTERIOR)
     float zout, yout, xout;
      if( cur_z > (-1.0*dist_local) ){
         zout = cur_z + dist_local ;
     }
     else{  zout = 0 ; }
     if( fabs(cur_y) + dist_local > w/2.0 ){
            yout = fabs(cur_y) + dist_local - w/2.0 ;
      }
        else{ yout = 0 ; }
     if( fabs(cur_x) + dist_local > w/2.0 ){
            xout = fabs(cur_x) + dist_local - w/2.0 ;
      }
     else{ xout = 0 ; }

        float correct_local_volume = (2.0*dist_local-zout)*(2.0*dist_local-yout)*(2.0*dist_local-xout) ;

        //float density_loc = (1.0*count_loc)/volume_local ; 
        float density_loc = (1.0*count_loc)/correct_local_volume ; 

        float converted_loc_dens = density_loc * ( volume_well / (1.0*K) ) ; // converts Barteks N/K to density -- ## IS THIS CORRECT FOR SUBSTRATE MOLECULES?!?!


     if( correct_local_volume < 0 ){ // THIS ARISES FROM BUGS BEING OUTSIDE THE WELLS... THIS NEEDS FIXED PROPERLY. DONT NOW HOW IT OCCURS

            //cout << "ZZZ  correct_local_volume = " << correct_local_volume << "\n";
            //cout << "xout, yout, zout = " << xout << " " << yout << " " << zout << "\n" ;
            //cout << "cur_(x,y,z) = " << cur_x << ", " << cur_y << ", " << cur_z << "\n" ;
            //cout << "current well = " << cur_well << "\n" ;
         //exit (6) ;

         br = 0 ;
     }
        else{
            br=genotypes[cells[n].gen]->growth(cells[n].well) * growth_function( 1.0-converted_loc_dens ) ;      // Use (1-local density of nutrient) here!!!! 
     }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (br<0) br=0 ;
        dr=genotypes[cells[n].gen]->death(cells[n].well) ;  // death in volume
        if (q<br) mode=0 ;
        else if (q<br+dr) mode=1 ;
        else mode=2 ; 


        if (mode==0) { // reproduction (make 2 new cells, kill mother cell)
          if (_drand48()>genotypes[cells[n].gen]->death_during_rep(cells[n].well)) { // replication

            int no_SNPs=poisson(gama) ; // newly produced cell mutants
            Cell c ; c.well=cells[n].well ; c.x=cells[n].x ; c.y=cells[n].y ; c.z=cells[n].z ;

            int no_SNPs_2 = poisson( gama ) ; // newly produced cell mutants
            Cell c_2 ; c_2.well=cells[n].well ; c_2.x=cells[n].x ; c_2.y=cells[n].y ; c_2.z=cells[n].z ;


            //add 1 to well population (mother is deleted later so add 2)
            wellpops[ cells[n].well ]+=2;


            if (no_SNPs>0) { 
              c.gen=genotypes.size() ; 
              genotypes.push_back(new Genotype(genotypes[cells[n].gen], cells[n].gen, no_SNPs)) ; // mutate 
            } else { 
              c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
            }
            cells.push_back(c) ;

            if (no_SNPs_2 > 0) { 
              c_2.gen = genotypes.size() ;   
              genotypes.push_back(new Genotype(genotypes[cells[n].gen], cells[n].gen, no_SNPs_2 ) ) ; // mutate 
            } else { 
              c_2.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
            }
            cells.push_back(c_2) ;
    
            // ####   We have produced 2 new cells above, hence now must kill mother cell ##############
            mode=1 ;
            // #########################################################################################


            // If bacterium replicates, must remove 1 food molecule from local neighbourhood
            if (foodmoltoremove != foodmols.size()-1) { 
            foodmols[ foodmoltoremove ] = foodmols[foodmols.size()-1] ;
            }
            foodmols.pop_back() ; 
//            bool removedfood = false;
//            for( int fm=0; fm<foodmols.size(); fm++ ){
//                if( foodmols[fm].well==cur_well && !removedfood ){
//                     if( foodmols[fm].x < cur_x+dist_local  &&  foodmols[fm].x > cur_x-dist_local ){
 //                       if( foodmols[fm].y < cur_y+dist_local  &&  foodmols[fm].y > cur_y-dist_local ){
  //                          if( foodmols[fm].z < cur_z+dist_local  &&  foodmols[fm].z > cur_z-dist_local ){
 //                               // delete food molecule
  //                              if (fm!=foodmols.size()-1) { 
  //                              foodmols[fm]=foodmols[foodmols.size()-1] ;
 //                               }
  //                              foodmols.pop_back() ; 
  //                              removedfood=true;
 //                               fm=foodmols.size()+10;
 //                           }
  //                      }
 //                   }
 //               }
  //          }
            



    
          } else { // death during replication
            mode=1 ;
          }
        }
      
        // now we implement death
        if (mode==1) {
          
          wellpops[ cells[n].well ]--;

          genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
            delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
          }
          if (n!=cells.size()-1) { 
            cells[n]=cells[cells.size()-1] ;
          }
          cells.pop_back() ; 
        }
    
   
    }//end else{ cell }






    ntot=cells.size() ; 
    int ntot_2=0;
    for( int dsa=0; dsa<nwells; dsa++){ ntot_2 += wellpops[ dsa ] ; }
//    if (it%1000000==0) cout <<tt<<": "<<ntot<<endl ;
    if (it%1000000==0) cout <<tt<<": "<<ntot<<"   "<<ntot_2<<endl ;

    double wait_time=1 ;                                       
    if (tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
//    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }


// --- added to produce Fisher plots -----------------------------------------
    double wait_time_2=0.005 ;
    if (tt > tt_old_2+wait_time_2 ) { tt_old_2=tt ; save_well_pops();  save_ODs(); }
// ---------------------------------------------------------------------------


    if (cells.size()==0) return 1 ; 
    if (tt>tmax) return 2 ; // return after prescribed time  

  }

}















void reset() 
{
//  cout <<"reset\n"  ;
  tt=0 ; //max_growth_rate=1 ;
  for (int i=0;i<genotypes.size();i++) if (genotypes[i]!=NULL) delete genotypes[i] ;
  genotypes.clear() ; genotypes.push_back(new Genotype) ;  
  for (int i=0; i<nwells; i++ ){ wellpops[i]=0; }
  cells.clear() ; 


  // Place cells randomly throughout the first well
  for (int i=0;i<init_bact;i++) { 
    Cell c ; 
    c.well=0 ; c.gen=0 ;   
    c.x = w * _drand48() * (0.99) - (w*0.99/2.0) ;
    c.y = w * _drand48() * (0.99) - (w*0.99/2.0) ;
    /////c.x = w * _drand48() * (0.4286) - (w*0.4286/2.0) ; //leave 1mm gaps in x and y direction
    /////c.y = w * _drand48() * (0.4286) - (w*0.4286/2.0) ;
    c.z = h * _drand48() * (-0.15) ;
    //
    //c.x = w * _drand48() * (0.3) - (w*0.3/2.0) ;
    //c.y = w * _drand48() * (0.3) - (w*0.3/2.0) ;
    //c.z = h * _drand48() * (-0.3) - (h*0.3) ;
    cells.push_back(c) ; 
// debugged
    genotypes[0]->number++ ;
    wellpops[0]++;
  }


  for (int i=0; i<food0; i++){
    Cell c ; 
    c.well=0 ; c.gen=0 ;   
    c.x = w * _drand48() * (0.99) - (w*0.99/2.0) ;
    c.y = w * _drand48() * (0.99) - (w*0.99/2.0) ;
    c.z = h * _drand48() * (-0.99) ;

    foodmols.push_back(c) ; 
  }




  for (int i=0;i<nwells;i++) {
    well[i].a=uni_antibiotic+pow(max_antibiotic,1.*(i+1)/nwells) ;  // CHANGED SO THAT WELL 0 HAS ZERO ANTIBIOTICS
  }



  // erase output buffer for "times"
#if defined __linux
  times->_IO_write_ptr = times->_IO_write_base ;
#elif defined __APPLE__
  // not defined yet
#else
  times->_ptr = times->_base ; // this operates on elements of _iobuf and is specific to Windows GNU C++
#endif
}





void init()
{
  well=new Well[nwells] ;
  char txt[256] ;
  sprintf(txt,"mkdir %s",NUM) ; system(txt) ;
  sprintf(txt,"%s/%s_%d.dat",NUM,NUM,RAND) ; times=fopen(txt,"w") ;
  sprintf(txt,"%s/state.dat",NUM) ; FILE *f=fopen(txt,"w") ; fclose(f) ;

  timesbuffer=new char[(1<<16)] ;
  setvbuf (times , timesbuffer , _IOFBF , (1<<16));   // this is to prevent saving data if no fflush is attempted 
                                                      // (this e.g. allows one to discard N<256)
  start_clock=clock() ;
}





void end() {
  fclose(times) ; 
}








int main(int argc, char *argv[])
{
  int nsam ;
  if (argc!=4) { err(" Error:: arguments needed: name, no_samples, RAND. Program terminated. \n"); } 
  else { 
    NUM=argv[1] ;
    nsam=atoi(argv[2]) ;
    RAND=atoi(argv[3]) ;
  }
  cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<endl ;
  _srand48(RAND) ;
  init();

  // PRINT PARAMETERS TO FILE----------------------------------
  char params[256] ;
  sprintf( params, "%s/parameters.txt", NUM);
  FILE *pff = fopen( params,"w") ;
  fprintf(pff, "nwells = %d\n", nwells);
  fprintf(pff, "L = %d\n", L);
  fprintf(pff, "max_no_drivers = %d\n", max_no_drivers);
  fprintf(pff, "gama = %f\n", gama);
  fprintf(pff, "max_growth_rate = %f\n", max_growth_rate);
  fprintf(pff, "tmax = %f\n", tmax);
  fprintf(pff, "uni_antibiotic = %f\n", uni_antibiotic);
  fprintf(pff, "max_antibiotic = %f\n", max_antibiotic);
  fprintf(pff, "hole_width = %f, hole_depth = %f\n", hole_width, hole_depth ) ;
  fprintf(pff, "well_width = %f, well_depth = %f\n", w, h ) ;
  fprintf(pff, "volume of \"local cube\" = %f", volume_local ) ;
  fprintf(pff, "K = %f", K ) ;
  fprintf(pff, "init_bact = %f", init_bact ) ;
  fprintf(pff, "food0 = %f", food0 ) ;
  fprintf(pff, "bacteria diff = %f\n", diff);
  fprintf(pff, "substr diff = %f\n", substr_diff);

  fclose( pff ); 
  //-------------------------------------------



  char name[256] ;
  sprintf(name,"%s/each_run.dat",NUM) ;
  FILE *er=fopen(name,"w") ; fclose(er) ;
  for (sample=0;sample<nsam;sample++) { 

    reset() ;


    if( tt==0 ) { save_ODs(); save_well_pops(); }


    int s=0 ; while (main_proc()==1) { s++ ; reset() ; } // initial growth until max size is reached


    if (s>0) printf("resetted %d times\n",s) ;
    fflush(stdout) ;
    save_data() ; 
    
    sprintf(name,"%s/each_run.dat",NUM) ;
    FILE *er=fopen(name,"a") ;
    fprintf(er,"%d\t%d %lf\n",sample,s,tt) ;
    fclose(er) ;

    sprintf(name,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(name) ;


   
  } 
  end() ;
	return 0 ;
}
