// this program simulates growth and evolution of a bacterial population invading a series of compartments
// the compartments are connected in a chain and have dimensions exactly as in my experiments

// v.1.0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "classes.h"
#define __MAIN

//------------parameters-----------
const int nwells=24 ; // number of compartments (wells)
const float w=3.5,h=10,hole=1 ; // size of the well [mm]. Square hole at the top y=(-hole/2,hole/2),z=(-hole,0)

const int L=10000 ; // total genome length (it does not have to be realistically large (4.5Mbp for E. coli) 
                    //but it needs to be larger than the expected number of mutations per cell at the end of the run
const int max_no_drivers=10 ; // max. no of mutations that cause resistance (driver mutations)
const float gama=1e-2 ; // mutation probability per daughter cell. 1e-2 is slightly higher than the base mutation rate (about 2e-3) in the absence of Cipro
                        // but the actual mutation rate may be even higher in the presence of Cipro

const float diff=1.5 ; // bacteria diffusion constant in [mm^2/h] ; typical value for E. coli is 1.5mm^2/h
const float max_growth_rate=2 ; // growth rate 1/h - 2 is max growth rate for E. coli
const float tmax=3*24 ; // duration of experiment in hours

const float init_bact=100 ; // initial number of cells in the first well
const float food0 = 1e5 ; // amount of food per well, bacteria consume food with max rate 1/h. This is unrealistically low at the moment, 
                          //but I suspect increasing it to a realistic value of 1e7 - 1e8 will change the number of accumulated mutations only as log(food0)
                          // so it is not very important to have a realistic value here
const float uni_antibiotic=0, max_antibiotic=50 ; // Cipro concentrations in ng/ml


//--------------------


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
FILE *times ; 
char *timesbuffer ;
Well *well ; 

#if defined __linux
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


static long long unsigned int _x=0x000100010001LL, _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  // works only on compilers with long long int!
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}

void _srand48(int a) { _x=a ; }


int poisson(void)  // generates k from P(k)=exp(-gamma) gamma^k / k!
{
  const double l=exp(-gama) ;
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
//  float gr=1 ; 
//  if (gr<0) gr=0 ;
  return max_growth_rate ;
}

float Genotype::death(int w) 
{
/*  float d=0.001 ;
//  d+=0.01*(1-well[w].food/food0) ;
  float mic=20 ; if (no_drivers>0) mic=1000 ;
  d+=SQR(well[w].a/mic) ;
  if (d>1) d=1 ;
  return d ; */
  return 0 ;
}

float Genotype::death_during_rep(int w)
{
  float d=0 ;
  float mic=20 ; if (no_drivers>0) mic=100 ; // MIC of the wild type is 20ng/ml, MIC of the mutant is assumed to be 100ng/ml. 
                                              // In reality there is a spectrum of mutants with different MICs
//  mic+=20*no_drivers ;
  
  d+=SQR(well[w].a/mic) ;
  if (d>1) d=1 ;
  return d ;  
}

float growth_function(float q)
{
  return q ; //BW q*q*q ; 
}


//double max_growth_rate ;

Genotype::Genotype(void) 
{ 
  number=1 ; no_drivers=0 ; sequence.clear() ; prev_gen=-1 ;
}


Genotype::Genotype(Genotype *mother, int prevg, int no_snp) { 
  prev_gen=prevg ;
  sequence=mother->sequence ; no_drivers=mother->no_drivers; 
  for (int i=0;i<no_snp;i++) {
    int k,j=int(_drand48()*L) ; 
    for (k=0;k<sequence.size();k++) if (sequence[k]==j) break ;
    if (k==sequence.size()) {
      sequence.push_back(j) ;
      if (j<max_no_drivers) no_drivers++ ; 
    }
  }
  number=1 ;
}

vector<Genotype*> genotypes ;


void save_genotypes(char *name)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d  %d  %d  %d\t",i, g->prev_gen,g->no_drivers, g->number) ;
      for (int j=0;j<g->sequence.size();j++) fprintf(data," %u",g->sequence[j]) ; 
      fprintf(data,"\n") ;
    } 
  }

  fclose(data) ;  
}

void save_most_abund_gens(char *name, int *most_abund)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *gg=genotypes[i] ;
    if (gg!=NULL && gg->number>0) {
      int r=0,g=0,b=0 ;
      for (int j=0;j<gg->sequence.size();j++) {
        if ((gg->sequence[j]&L_PM)==most_abund[0]) r=1 ;
        if ((gg->sequence[j]&L_PM)==most_abund[1]) g=1 ; 
        if ((gg->sequence[j]&L_PM)==most_abund[2]) b=1 ;
      }
      if (r || g || b) fprintf(data,"%d %d %d\t%d\n",r,g,b,gg->index) ;
    }
  }
  fclose(data) ;  

}

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

void save_snps(char *name,int *n, int total, int mode, int *most_abund) 
{
  const float cutoff=0.01 ;
  FILE *f=fopen(name,"w") ;
  if (f==NULL) err(name) ;
  int i, j, nsnps=0, nsnpsc=0;
  for (i=0;i<L;i++) { 
    if (n[i]>(1e-4)*total) nsnps++ ;
    if (1.*n[i]/total>cutoff) nsnpsc++ ;
//    cout <<n[i]<<" ";
  }
  float *abund=new float[nsnps], tempd ;
  int *num=new int[nsnps], temp ;
  nsnps=0 ;
  for (i=0;i<L;i++) if (n[i]>(1e-4)*total) { num[nsnps]=i ; abund[nsnps]=float(1.*n[i]/total)*(1+0.000001*i/L) ; nsnps++ ; }
  quicksort2(abund,num,0,nsnps-1) ;

  if (mode) {
    for (i=0;i<nsnps;i++) fprintf(f,"%d %d %f\n",i,num[i],abund[i]) ;
  } else {
    for (i=0;i<nsnps;i++) if (abund[i]>cutoff || i<100) fprintf(f,"%d %d %f\n",i,num[i],abund[i]) ;
  }
  if (most_abund!=NULL) { // store first 100 most abundant PMs
    for (i=0;i<MIN(100,nsnps);i++) most_abund[i]=num[i] ; 
  }
  delete [] abund ; delete [] num ;
  fclose(f) ;
}


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
      for (int j=0;j<genotypes[i]->sequence.size();j++) snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
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


int where_is_bug(float x, float y, float z)
{
  if (x>=-w/2 && x<w/2 && y>=-w/2 && y<w/2 && z<=0 && z>-h) return 0 ; // inside the well
  if (x>=w/2 && fabs(y)<hole/2 && z<=0 && z>-hole) return 1 ; // going to the right well
  if (x<-w/2 && fabs(y)<hole/2 && z<=0 && z>-hole) return 2 ; // going to the left well  
  return -1 ; // outside the well
}

int main_proc()
{
  int i,j,k,n,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;

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
    tt+=dt/cells.size() ; 
    n=_drand48()*cells.size() ;

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
      cells[n].well++ ;     // cout <<cells[n].well ; 
    } else if (ok==2 && cells[n].well>0) { // cell in another well to the left
      cells[n].x+=dx+w ; cells[n].y+=dy ; cells[n].z+=dz ;      
      cells[n].well-- ;  //cout <<cells[n].well ; 
    }


    double q=_drand48()*(max_growth_rate+max_death_rate), br,dr ; //,ml,mr  ;
    int mode=0 ; 
    br=genotypes[cells[n].gen]->growth(cells[n].well) * growth_function(well[cells[n].well].food/food0) ; 
    if (br<0) br=0 ;
    dr=genotypes[cells[n].gen]->death(cells[n].well) ;  // death in volume
//    if (cells[n].well>0) ml=migr ; else ml=0 ;
//    if (cells[n].well<nwells-1) mr=migr ; else mr=0 ;
    if (q<br) mode=0 ;
    else if (q<br+dr) mode=1 ;
    else mode=2 ; //if (q<br+dr+ml) mode=2 ;
//    else if (q<br+dr+ml+mr) mode=3 ; 
//    else mode=4 ;
//    cout <<mode ; 


    if (mode==0) { // reproduction
      if (_drand48()>genotypes[cells[n].gen]->death_during_rep(cells[n].well)) { // replication
        int no_SNPs=poisson() ; // newly produced cell mutants
        Cell c ; c.well=cells[n].well ; c.x=cells[n].x ; c.y=cells[n].y ; c.z=cells[n].z ;
        well[c.well].food-=1.0 ;
  //      if (c.well<nwells-1 && _drand48()<dt*migr) c.well++ ;
  //      if (c.well>0 && _drand48()<dt*migr) c.well-- ;
        if (no_SNPs>0) { 
          c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // mutate 
        } else { 
          c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
        }
        cells.push_back(c) ; 
      } else { // death during replication
        mode=1 ;
      }
    }
  
// now we implement death
    if (mode==1) {
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      if (n!=cells.size()-1) { 
        cells[n]=cells[cells.size()-1] ;
      }
      cells.pop_back() ; 
    }


//    cout <<dt<<endl ;
    ntot=cells.size() ; if (it%1000000==0) cout <<tt<<": "<<ntot<<endl ;
//    cout <<ntot<<": "<<dx<<endl ;

/*    if (it%100==0) cout <<ntot<<endl ;    
    if (it%1000==0) {
      char name[256] ;
      sprintf(name,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(name) ;
    }*/


    double wait_time=1 ; 
    if (tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
//    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }

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
  cells.clear() ; 
  Cell c ; 
  for (int i=0;i<init_bact;i++) { c.well=0 ; c.gen=0 ; c.x=c.y=c.z=0 ; cells.push_back(c) ; }
  for (int i=0;i<nwells;i++) {
    well[i].food=food0 ;
    well[i].a=uni_antibiotic+pow(max_antibiotic,1.*i/nwells) ;
//    cout <<well[i].a<<" " ;
  }
//  drivers.clear() ;
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
  setvbuf (times , timesbuffer , _IOFBF , (1<<16));  // this is to prevent saving data if no fflush is attempted 
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
  char name[256] ;
  sprintf(name,"%s/each_run.dat",NUM) ;
  FILE *er=fopen(name,"w") ; fclose(er) ;
  for (sample=0;sample<nsam;sample++) { 
    reset() ;
    int s=0 ; while (main_proc()==1) { s++ ; reset() ; } // initial growth until max size is reached
    if (s>0) printf("resetted %d times\n",s) ;
    fflush(stdout) ;
    save_data() ; 
    
    sprintf(name,"%s/each_run.dat",NUM) ;
    FILE *er=fopen(name,"a") ;
    fprintf(er,"%d\t%d %lf\n",sample,s,tt) ;
    fclose(er) ;

    sprintf(name,"%s/genotypes_%d_%d.dat",NUM,RAND,sample) ; save_genotypes(name) ;

    // save some more data

    sprintf(name,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(name) ;

    
    int *snp_no=new int[L] ; // array of SNPs abundances
    for (int i=0;i<L;i++) { snp_no[i]=0 ; }
    for (int i=0;i<genotypes.size();i++) {
      if (genotypes[i]!=NULL && genotypes[i]->number>0) 
        for (int j=0;j<genotypes[i]->sequence.size();j++) {
          snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
//          if (((genotypes[i]->sequence[j])&DRIVER_PM)) snp_drivers[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;
        }
    }

    printf("saving PMs...\n") ;
    int most_abund[100] ;
    sprintf(name,"%s/all_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_no,cells.size(),1,most_abund) ;
//    if (driver_adv>0 || driver_migr_adv>0) { printf("saving driver PMs...\n") ; sprintf(name,"%s/drv_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_drivers,max_size,0,NULL) ; }
    
// save SNPs in different wells
    int wells_to_save[7]={0,3,7,11,15,19,23} ;
    for (int w=0;w<7;w++) {  
      int ninwell=0 ;
      for (int i=0;i<L;i++) { snp_no[i]=0 ; }
      for (int i=0;i<cells.size();i++) {
        if (cells[i].well==wells_to_save[w]) {
          Genotype *g=genotypes[cells[i].gen] ; ninwell++ ;
        //if (genotypes[i]!=NULL && genotypes[i]->number>0) 
          for (int j=0;j<g->sequence.size();j++) snp_no[((g->sequence[j])&L_PM)]++ ;      
        }
      }
      sprintf(name,"%s/well_%d_PMs_%d_%d.dat",NUM,wells_to_save[w],RAND,sample) ; save_snps(name,snp_no,ninwell,1,most_abund) ;
    }
    
    delete [] snp_no ;
  } 
  end() ;
	return 0 ;
}
