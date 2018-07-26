// #######################################################################################################################
//
//  CIPROFLOXACIN-DEPENDENT MOTILITY:  (i) hard cut-off => no motility > f*MIC where f is some arbitrary cutoff fraction.
//
//  Non-motile cells cannot mutate (i.e. filaments that grow, consume nutrients and contribute to carrying capacity)
//
// #######################################################################################################################



// v.1.0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "classes_strain.h"
#include <algorithm>


#include <sstream>
#include <map>
#include <fstream>

#define __MAIN




//------------parameters-----------------------------------------------------------------------------------

const float motilityCutoff = 0.15 ;



const int nwells=24 ;                  // number of compartments (wells)
int wellpops[24];                    // no. of cells in each well 

const float w=3.25,  h=11.0,   hole_width=1.5,   hole_depth=3.75 ;      // size of the well [mm]. Square hole at the top y=(-hole/2,hole/2),z=(-hole,0)


//const int L=10000 ;                   // total genome length (it does not have to be realistically large (4.5Mbp for E. coli) 
                                        // but it needs to be larger than the expected number of mutations per cell at the end of the run


const int L = 5;
const int max_no_drivers = 5 ;         // max. no of mutations that cause resistance (driver mutations)
const float gama=1e-3 ;                // mutation probability per daughter cell. 1e-2 is slightly higher than the base mutation rate (about 2e-3) in the absence of Cipro
                                       // but the actual mutation rate may be even higher in the presence of Cipro
//#### NOT 100% SURE WHAT GAMA IS NOW IVE MODIFIED poisson()

const float diff = 1.4 ;                 // bacteria diffusion constant in [mm^2/h] ; typical value for E. coli is 1.5mm^2/h
const float max_growth_rate = 1.0 ;      // growth rate 1/h - 2 is max growth rate for E. coli
const float tmax=3*24;                   // duration of experiment in hours

const float init_bact = 1000 ;          // initial number of cells in the first well
const float K = 1e5 ;                   // carrying capacity. Have removed "food".

//const float food0 = 1e6 ;            // NO LONGER USED, remove this from code

const float uni_antibiotic=0, max_antibiotic=300;   // Cipro concentrations in ng/ml

const float MIC_WT = 32 ;


// Define probabilities for each mutation in Marcusson2009:
//   (note: this is the probability that GIVEN A MUTATION IS GOING TO HAPPEN, it is in one of the driver positions. HOW TO ESTIMATE THESE?!?)
const float p_g0 = 1.0/1e3 ;     //gyrA1
const float p_g1 = p_g0 ;        //gyrA2
const float p_g2 = p_g0 ;        //parC
const float p_g3 = 10*p_g0 ;    //marR
const float p_g4 = 10*p_g0;     //acrR
 

// Mutation rate 'mu' depends on Cipro concentration, 'a', with form:  mu = gama * [ 1 + ( X.a/MIC_WT )^Y  ], with some upper limit.
const float X = 15.0;   const float Y = 1.5;    //parameters governing form of mut-rate dependence on antibiotic.
const int max_increase = 200;                   //i.e. max_increase*gama is max possible mutation rate per replication



const int num_strains = 32 ;  // WT + 28 Marcusson2009 strains + 3 not included in Marcusson
//---------------------------------------------------------------------------------------------------------


//  LB growth data   (length of arrays = 132, hardcoded into growth_function() method)

const float NKs[132] = {0.0, 0.0058823529411764705, 0.011764705882352941, 0.01764705882352941, 0.023529411764705882, 0.029411764705882356, 0.03529411764705883, 0.041176470588235294, 0.047058823529411764, 0.052941176470588235, 0.058823529411764705, 0.06470588235294118, 0.07058823529411765, 0.07647058823529412, 0.0823529411764706, 0.08823529411764708, 0.09411764705882353, 0.1, 0.10588235294117648, 0.11176470588235295, 0.11764705882352942, 0.1235294117647059, 0.12941176470588237, 0.13529411764705884, 0.14117647058823532, 0.14705882352941177, 0.15294117647058825, 0.15882352941176472, 0.1647058823529412, 0.17058823529411768, 0.17647058823529413, 0.1823529411764706, 0.18823529411764706, 0.19411764705882353, 0.2, 0.2058823529411765, 0.21176470588235297, 0.21764705882352942, 0.2235294117647059, 0.22941176470588237, 0.23529411764705885, 0.24117647058823533, 0.2470588235294118, 0.2529411764705882, 0.25882352941176473, 0.2647058823529412, 0.2705882352941177, 0.27647058823529413, 0.28235294117647064, 0.2882352941176471, 0.29411764705882354, 0.3, 0.3058823529411765, 0.31176470588235294, 0.31764705882352945, 0.3235294117647059, 0.3294117647058824, 0.33529411764705885, 0.34117647058823536, 0.34705882352941175, 0.35294117647058826, 0.3588235294117647, 0.3647058823529412, 0.37058823529411766, 0.3764705882352941, 0.3823529411764706, 0.38823529411764707, 0.3941176470588236, 0.4, 0.4058823529411765, 0.411764705882353, 0.4176470588235295, 0.4235294117647059, 0.4294117647058824, 0.43529411764705883, 0.4411764705882353, 0.4470588235294118, 0.45294117647058824, 0.45882352941176474, 0.4647058823529412, 0.4705882352941177, 0.47647058823529415, 0.48235294117647065, 0.4882352941176471, 0.4941176470588236, 0.5, 0.5058823529411764, 0.511764705882353, 0.5176470588235295, 0.5235294117647059, 0.5294117647058824, 0.5352941176470588, 0.5411764705882354, 0.5470588235294118, 0.5529411764705883, 0.5588235294117647, 0.5647058823529413, 0.5705882352941176, 0.5764705882352942, 0.5823529411764706, 0.5882352941176471, 0.5941176470588235, 0.6, 0.6058823529411765, 0.611764705882353, 0.6176470588235294, 0.6235294117647059, 0.6294117647058824, 0.6352941176470589, 0.6411764705882353, 0.6470588235294118, 0.6529411764705884, 0.6588235294117648, 0.6647058823529413, 0.6705882352941177, 0.6764705882352943, 0.6823529411764707, 0.6882352941176471, 0.6941176470588235, 0.7, 0.7058823529411765, 0.711764705882353, 0.7176470588235294, 0.7235294117647059, 0.7294117647058824, 0.7352941176470589, 0.7411764705882353, 0.7470588235294118, 0.7529411764705882, 0.7588235294117648, 0.7647058823529412, 1.0};



const float LB_growth_rates[132] = {2.0, 1.4515383528820198, 1.3268456392990229, 1.2460908592325395, 1.1846451596215506, 1.134283973750877, 1.0911905551573338, 1.0532647233941137, 1.0192207364737311, 0.9882112281288045, 0.9596466297437259, 0.9330993012989593, 0.9082486451646022, 0.8848477738334221, 0.8627022830879314, 0.8416561953338396, 0.8215823409458175, 0.8023755902452877, 0.7839479753985176, 0.7662251000781128, 0.7491434478301104, 0.732648330996518, 0.7166923048293633, 0.7012339251501226, 0.6862367635732859, 0.6716686184846697, 0.6575008766551941, 0.6437079920911613, 0.6302670570802816, 0.6171574464378226, 0.6043605203866059, 0.5918593747889519, 0.5796386299112848, 0.5676842507675355, 0.5559833935142047, 0.5445242734708579, 0.5332960511964979, 0.5222887337240966, 0.5114930885863537, 0.5009005686879714, 0.49050324641787374, 0.48029375566719057, 0.4702652406396036, 0.46041131052055717, 0.4507259992192211, 0.44120372951843595, 0.4318392810682613, 0.42262776174216915, 0.41356458194455936, 0.4046454315166416, 0.3958662589368562, 0.3872232525534866, 0.3787128236223371, 0.37033159095227813, 0.3620763669870753, 0.3539441451738354, 0.3459320884872817, 0.3380375189953642, 0.33025790836581714, 0.3225908692255455, 0.3150341472954447, 0.3075856142326207, 0.3002432611202417, 0.2930051925525486, 0.2858696212690258, 0.2788348632975086, 0.27189933357120677, 0.2650615419892983, 0.25832008989503097, 0.2516736669491774, 0.245121048380296, 0.23866109259666474, 0.2322927391478957, 0.22601500702730692, 0.21982699330902147, 0.2137278721165939, 0.20771689392275633, 0.20179338518257317, 0.1959567483050643, 0.190206461971, 0.18454208180734244, 0.17896324143143022, 0.17346965388069063, 0.16806111344620409, 0.16273749793085315, 0.15749877135495985, 0.1523449871340351, 0.14727629175442222, 0.14229292897283247, 0.1373952445647053, 0.1325836916434338, 0.127858836566989, 0.12322136543960888, 0.11867209120243125, 0.11421196128692858, 0.10984206577647826, 0.10556364598170644, 0.10137810328118886, 0.09728700800635194, 0.09329210805319035, 0.08939533677799223, 0.08559881957299086, 0.08190487831467703, 0.07831603262609003, 0.07483499659181396, 0.07146466921154675, 0.06820811648543784, 0.0650685426155766, 0.062049247428196244, 0.059153566843654026, 0.056384793155643355, 0.05374607217598019, 0.051240275138533097, 0.0488698448281974, 0.046636617864559866, 0.044541628470238084, 0.042584903227772364, 0.04076526081233089, 0.03908013467394313, 0.03752543905436482, 0.03609549843832685, 0.03478305676227284, 0.03357937537289606, 0.032474418764680796, 0.03145711635618678, 0.0305156792519155, 0.02963794507643564, 0.028811722645465793, 0.028025111338640003, 0.027266776402217897, 0.02652616933704952, 0.0};
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



//int poisson(void)  // generates k from P(k)=exp(-gamma) gamma^k / k!
//{
//  const double l=exp(-gama) ;
//  double p=1. ;
//  int k=0 ;
//  do {
//    k++ ;
//    p*=_drand48() ;
//  } while (p > l) ;
//  return k - 1 ;
//}




// Modified to use half the mutation rate to determine #SNPs in each daughter cell. 
int poisson(float g)  // generates k from P(k)=exp(-g) g^k / k!
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







//###############################################################################################################################################

float Genotype::growth(int well) 
{
  //float growth = max_growth_rate ;
  float growth = 1.0 ;                   // This is now 1.0 since Bartek's experimental LB data begins at 2.0

  bool g0=false, g1=false, g2=false, g3=false, g4=false;
  if( find( sequence.begin(), sequence.end(), 0) != sequence.end()  )  {  g0=true;  }  // gyrA1
  if( find( sequence.begin(), sequence.end(), 1) != sequence.end()  )  {  g1=true;  }  // gyrA2
  if( find( sequence.begin(), sequence.end(), 2) != sequence.end()  )  {  g2=true;  }  // parC
  if( find( sequence.begin(), sequence.end(), 3) != sequence.end()  )  {  g3=true;  }  // marR
  if( find( sequence.begin(), sequence.end(), 4) != sequence.end()  )  {  g4=true;  }  // acrR

  // Check for all 28 types listed in Marcusson2009:
  //  CAREFUL - NEED TO SET EXACTLY CORRECT NUMBER OF DRIVER MUTATIONS NOW
  if ( no_drivers == 0  ){
    growth = max_growth_rate ;
  }
  else if ( no_drivers == 1 ){
    if      ( g4==true )                                           { growth = 0.91  ; }  // LM351
    else if ( g3==true )                                           { growth = 0.83  ; }  // LM202
    else if ( g2==true )                                           { growth = 0.99  ; }  // LM792
    else if ( g1==true )                                           { growth = 0.99  ; }  // LM534
    else if ( g0==true )                                           { growth = 1.01  ; }  // LM378
  }
  else if ( no_drivers == 2 ){
    if      ( g3==true && g4==true )                               { growth = 0.82  ; }  // LM367
    else if ( g1==true && g4==true )                               { growth = 0.92  ; }  // LM592
    else if ( g1==true && g3==true )                               { growth = 0.83  ; }  // LM538
    else if ( g1==true && g2==true )                               { growth = 1.02  ; }  // LM1124
    else if ( g0==true && g4==true )                               { growth = 0.95  ; }  // LM647
    else if ( g0==true && g3==true )                               { growth = 0.86  ; }  // LM421
    else if ( g0==true && g2==true )                               { growth = 0.98  ; }  // LM862
    else if ( g0==true && g1==true )                               { growth = 0.97  ; }  // LM625
    // extra
    else if ( g2==true && g3==true )                               { growth = 0.83  ; }  // parC+LM202
    else if ( g2==true && g4==true )                               { growth = 0.91  ; }  // parC+LM351
  }
  else if ( no_drivers == 3 ){
    if      ( g1==true && g3==true && g4==true )                   { growth = 0.60  ; }  // LM595
    else if ( g0==true && g3==true && g4==true )                   { growth = 0.78  ; }  // LM709
    else if ( g1==true && g2==true && g4==true )                   { growth = 0.95  ; }  // LM1125
    else if ( g1==true && g2==true && g3==true )                   { growth = 0.84  ; }  // LM882
    else if ( g0==true && g2==true && g4==true )                   { growth = 0.92  ; }  // LM873
    else if ( g0==true && g2==true && g3==true )                   { growth = 0.86  ; }  // LM871
    else if ( g0==true && g1==true && g4==true )                   { growth = 0.93  ; }  // LM691
    else if ( g0==true && g1==true && g3==true )                   { growth = 0.79  ; }  // LM695
    else if ( g0==true && g1==true && g2==true )                   { growth = 1.01  ; }  // LM693
    // extra
    else if ( g2==true && g3==true && g4==true )                   { growth = 0.82  ; }  // parC+LM367
  }
  else if ( no_drivers == 4 ){
    if      ( g1==true && g2==true && g3==true && g4==true )       { growth = 0.72  ; }  // LM878
    else if ( g0==true && g2==true && g3==true && g4==true )       { growth = 0.71  ; }  // LM875
    else if ( g0==true && g1==true && g2==true && g4==true )       { growth = 0.94  ; }  // LM703
    else if ( g0==true && g1==true && g2==true && g3==true )       { growth = 0.89  ; }  // LM707
    else if ( g0==true && g1==true && g3==true && g4==true )       { growth = 0.66  ; }  // LM701
  }
  else if ( no_drivers == 5 ){
    if( g0==true && g1==true && g2==true && g3==true && g4==true ) { growth = 0.68 ; }  // LM705
  }

  return growth ;
}

//###############################################################################################################################################



float Genotype::death_during_rep(int w)      // MAYBE WE SHOULD MULTIPLE ALL MICs BY 2 TO BE IN AGREEMENT WITH BARTEK'S EXPERIMENTS?
{

  float d=0 ;
  float mic = 0.5 * MIC_WT ; //(this will be x2 later in method)

  bool g0=false, g1=false, g2=false, g3=false, g4=false;
  if( find( sequence.begin(), sequence.end(), 0) != sequence.end()  )  {  g0=true;  }  // gyrA1
  if( find( sequence.begin(), sequence.end(), 1) != sequence.end()  )  {  g1=true;  }  // gyrA2
  if( find( sequence.begin(), sequence.end(), 2) != sequence.end()  )  {  g2=true;  }  // parC
  if( find( sequence.begin(), sequence.end(), 3) != sequence.end()  )  {  g3=true;  }  // marR
  if( find( sequence.begin(), sequence.end(), 4) != sequence.end()  )  {  g4=true;  }  // acrR

  // Check for all 28 types listed in Marcusson2009:
  //  CAREFUL - NEED TO SET EXACTLY CORRECT NUMBER OF DRIVER MUTATIONS NOW
  if ( no_drivers == 0  ){
    // do nothing
  }
  else if ( no_drivers == 1 ){
    if      ( g4==true )                                           { mic = 47 ; }  // LM351
    else if ( g3==true )                                           { mic = 32 ; }  // LM202
    else if ( g2==true )                                           { mic = 16 ; }  // LM792
    else if ( g1==true )                                           { mic = 250 ; }  // LM534
    else if ( g0==true )                                           { mic = 380 ; }  // LM378
  }
  else if ( no_drivers == 2 ){
    if      ( g3==true && g4==true )                               { mic = 125 ; }  // LM367
    else if ( g1==true && g4==true )                               { mic = 380 ; }  // LM592
    else if ( g1==true && g3==true )                               { mic = 1000 ; }  // LM538
    else if ( g1==true && g2==true )                               { mic = 380 ; }  // LM1124
    else if ( g0==true && g4==true )                               { mic = 500 ; }  // LM647
    else if ( g0==true && g3==true )                               { mic = 1000 ; }  // LM421
    else if ( g0==true && g2==true )                               { mic = 1000 ; }  // LM862
    else if ( g0==true && g1==true )                               { mic = 380 ; }  // LM625
    // extra
    else if ( g2==true && g3==true )                               { mic = 32 ; }  // parC+LM202
    else if ( g2==true && g4==true )                               { mic = 47 ; }  // parC+LM351
  }
  else if ( no_drivers == 3 ){
    if      ( g1==true && g3==true && g4==true )                   { mic = 1500 ; }  // LM595
    else if ( g0==true && g3==true && g4==true )                   { mic = 1500 ; }  // LM709
    else if ( g1==true && g2==true && g4==true )                   { mic = 500 ; }  // LM1125
    else if ( g1==true && g2==true && g3==true )                   { mic = 750 ; }  // LM882
    else if ( g0==true && g2==true && g4==true )                   { mic = 3000 ; }  // LM873
    else if ( g0==true && g2==true && g3==true )                   { mic = 6000 ; }  // LM871
    else if ( g0==true && g1==true && g4==true )                   { mic = 750 ; }  // LM691
    else if ( g0==true && g1==true && g3==true )                   { mic = 750 ; }  // LM695
    else if ( g0==true && g1==true && g2==true )                   { mic = 32000 ; }  // LM693
    // extra
    else if ( g2==true && g3==true && g4==true )                   { mic = 125 ; }  // parC+LM367
  }
  else if ( no_drivers == 4 ){
    if      ( g1==true && g2==true && g3==true && g4==true )       { mic = 2000 ; }  // LM878
    else if ( g0==true && g2==true && g3==true && g4==true )       { mic = 8000 ; }  // LM875
    else if ( g0==true && g1==true && g2==true && g4==true )       { mic = 32000 ; }  // LM703
    else if ( g0==true && g1==true && g2==true && g3==true )       { mic = 32000 ; }  // LM707
    else if ( g0==true && g1==true && g3==true && g4==true )       { mic = 2000 ; }  // LM701
  }
  else if ( no_drivers == 5 ){
    if( g0==true && g1==true && g2==true && g3==true && g4==true ) { mic = 32000 ; }  // LM705
  }


  // ------ MUTLIPLY ALL MARCUSSON2009 VALUES BY 2 to better agrre with Bartek's experiments------
  mic = mic * 2;
  // ---------------------------------------------------------------------------------------------


  d = SQR(well[w].a/mic) ;

  //if (d>1) d=0.98 ;
  if (d>1) d=1.0 ;
   

  return d ;  
}


//################################################################################################################################################





float Genotype::getMIC()      // MAYBE WE SHOULD MULTIPLE ALL MICs BY 2 TO BE IN AGREEMENT WITH BARTEK'S EXPERIMENTS?
{

  float mic = 0.5 * MIC_WT ; //(this will be x2 later in method)

  bool g0=false, g1=false, g2=false, g3=false, g4=false;
  if( find( sequence.begin(), sequence.end(), 0) != sequence.end()  )  {  g0=true;  }  // gyrA1
  if( find( sequence.begin(), sequence.end(), 1) != sequence.end()  )  {  g1=true;  }  // gyrA2
  if( find( sequence.begin(), sequence.end(), 2) != sequence.end()  )  {  g2=true;  }  // parC
  if( find( sequence.begin(), sequence.end(), 3) != sequence.end()  )  {  g3=true;  }  // marR
  if( find( sequence.begin(), sequence.end(), 4) != sequence.end()  )  {  g4=true;  }  // acrR

  // Check for all 28 types listed in Marcusson2009:
  //  CAREFUL - NEED TO SET EXACTLY CORRECT NUMBER OF DRIVER MUTATIONS NOW
  if ( no_drivers == 0  ){
    // do nothing
  }
  else if ( no_drivers == 1 ){
    if      ( g4 )                                           { mic = 47 ; }  // LM351
    else if ( g3 )                                           { mic = 32 ; }  // LM202
    else if ( g2 )                                           { mic = 16 ; }  // LM792
    else if ( g1 )                                           { mic = 250 ; }  // LM534
    else if ( g0 )                                           { mic = 380 ; }  // LM378
  }
  else if ( no_drivers == 2 ){
    if      ( g3 && g4 )                               { mic = 125 ; }  // LM367
    else if ( g1 && g4 )                               { mic = 380 ; }  // LM592
    else if ( g1 && g3 )                               { mic = 1000 ; }  // LM538
    else if ( g1 && g2 )                               { mic = 380 ; }  // LM1124
    else if ( g0 && g4 )                               { mic = 500 ; }  // LM647
    else if ( g0 && g3 )                               { mic = 1000 ; }  // LM421
    else if ( g0 && g2 )                               { mic = 1000 ; }  // LM862
    else if ( g0 && g1 )                               { mic = 380 ; }  // LM625
    // extra
    else if ( g2 && g3 )                               { mic = 32 ; }  // parC+LM202
    else if ( g2 && g4 )                               { mic = 47 ; }  // parC+LM351
  }
  else if ( no_drivers == 3 ){
    if      ( g1==true && g3==true && g4==true )                   { mic = 1500 ; }  // LM595
    else if ( g0==true && g3==true && g4==true )                   { mic = 1500 ; }  // LM709
    else if ( g1==true && g2==true && g4==true )                   { mic = 500 ; }  // LM1125
    else if ( g1==true && g2==true && g3==true )                   { mic = 750 ; }  // LM882
    else if ( g0==true && g2==true && g4==true )                   { mic = 3000 ; }  // LM873
    else if ( g0==true && g2==true && g3==true )                   { mic = 6000 ; }  // LM871
    else if ( g0==true && g1==true && g4==true )                   { mic = 750 ; }  // LM691
    else if ( g0==true && g1==true && g3==true )                   { mic = 750 ; }  // LM695
    else if ( g0==true && g1==true && g2==true )                   { mic = 32000 ; }  // LM693
    // extra
    else if ( g2==true && g3==true && g4==true )                   { mic = 125 ; }  // parC+LM367
  }
  else if ( no_drivers == 4 ){
    if      ( g1==true && g2==true && g3==true && g4==true )       { mic = 2000 ; }  // LM878
    else if ( g0==true && g2==true && g3==true && g4==true )       { mic = 8000 ; }  // LM875
    else if ( g0==true && g1==true && g2==true && g4==true )       { mic = 32000 ; }  // LM703
    else if ( g0==true && g1==true && g2==true && g3==true )       { mic = 32000 ; }  // LM707
    else if ( g0==true && g1==true && g3==true && g4==true )       { mic = 2000 ; }  // LM701
  }
  else if ( no_drivers == 5 ){
    if( g0==true && g1==true && g2==true && g3==true && g4==true ) { mic = 32000 ; }  // LM705
  }


  // ------ MUTLIPLY ALL MARCUSSON2009 VALUES BY 2 to better agrre with Bartek's experiments------
  mic = mic * 2;
  // ---------------------------------------------------------------------------------------------

  return mic ;  
}


//################################################################################################################################################













// METHOD TO FIND STRAIN OF BACTERIA (LABELLED AS INTEGER)
int Genotype::getStrain( )     
{

  int str = -9 ;

  bool g0=false, g1=false, g2=false, g3=false, g4=false;
  if( find( sequence.begin(), sequence.end(), 0) != sequence.end()  )  {  g0=true;  }  // gyrA1
  if( find( sequence.begin(), sequence.end(), 1) != sequence.end()  )  {  g1=true;  }  // gyrA2
  if( find( sequence.begin(), sequence.end(), 2) != sequence.end()  )  {  g2=true;  }  // parC
  if( find( sequence.begin(), sequence.end(), 3) != sequence.end()  )  {  g3=true;  }  // marR
  if( find( sequence.begin(), sequence.end(), 4) != sequence.end()  )  {  g4=true;  }  // acrR


  if ( no_drivers == 0  ){
    str = 0; 
  }
  else if ( no_drivers == 1 ){
    if      ( g4==true )                                           { str = 5 ; }  // LM351
    else if ( g3==true )                                           { str = 4 ; }  // LM202
    else if ( g2==true )                                           { str = 3 ; }  // LM792
    else if ( g1==true )                                           { str = 2 ; }  // LM534
    else if ( g0==true )                                           { str = 1 ; }  // LM378
  }
  else if ( no_drivers == 2 ){
    if      ( g3==true && g4==true )                               { str = 13 ; }  // LM367
    else if ( g1==true && g4==true )                               { str = 12 ; }  // LM592
    else if ( g1==true && g3==true )                               { str = 11 ; }  // LM538
    else if ( g1==true && g2==true )                               { str = 10 ; }  // LM1124
    else if ( g0==true && g4==true )                               { str = 9 ; }  // LM647
    else if ( g0==true && g3==true )                               { str = 8 ; }  // LM421
    else if ( g0==true && g2==true )                               { str = 7 ; }  // LM862
    else if ( g0==true && g1==true )                               { str = 6 ; }  // LM625
    // extra
    else if ( g2==true && g3==true )                               { str = 29 ; }  // parC+LM202
    else if ( g2==true && g4==true )                               { str = 30 ; }  // parC+LM351
  }
  else if ( no_drivers == 3 ){
    if      ( g1==true && g3==true && g4==true )                   { str = 22 ; }  // LM595
    else if ( g0==true && g3==true && g4==true )                   { str = 21 ; }  // LM709
    else if ( g1==true && g2==true && g4==true )                   { str = 20 ; }  // LM1125
    else if ( g1==true && g2==true && g3==true )                   { str = 19 ; }  // LM882
    else if ( g0==true && g2==true && g4==true )                   { str = 18 ; }  // LM873
    else if ( g0==true && g2==true && g3==true )                   { str = 17 ; }  // LM871
    else if ( g0==true && g1==true && g4==true )                   { str = 16 ; }  // LM691
    else if ( g0==true && g1==true && g3==true )                   { str = 15 ; }  // LM695
    else if ( g0==true && g1==true && g2==true )                   { str = 14 ; }  // LM693
    // extra
    else if ( g2==true && g3==true && g4==true )                   { str = 31 ; }  // parC+LM367
  }
  else if ( no_drivers == 4 ){
    if      ( g1==true && g2==true && g3==true && g4==true )       { str = 27 ; }  // LM878
    else if ( g0==true && g2==true && g3==true && g4==true )       { str = 26 ; }  // LM875
    else if ( g0==true && g1==true && g2==true && g4==true )       { str = 25 ; }  // LM703
    else if ( g0==true && g1==true && g2==true && g3==true )       { str = 24 ; }  // LM707
    else if ( g0==true && g1==true && g3==true && g4==true )       { str = 23 ; }  // LM701
  }
  else if ( no_drivers == 5 ){
    if( g0==true && g1==true && g2==true && g3==true && g4==true ) { str = 28; }  // LM705
  }


  if( str == -9 ){ 
    cout<<"Strain in getStrain() is -9 : ERROR "<<endl;  
    cout<<g0<<" "<<g1<<" "<<g2<<" "<<g3<<" "<<g4<<endl ; 
    cout<<"num_drivers = "<<no_drivers<<endl ;
    exit (EXIT_FAILURE) ;
  }

  return str ;  
}


//################################################################################################################################################






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





//float growth_function(float q)
//{
//  return q ; //BW q*q*q ; 
//}



// Modified to incorporate Bartek's experimental LB-growth data
float growth_function( int wl )
{
  // There is now no food, just a growth rate that depends on N/K
  // where N is #cells in well, K is carrying capacity of well
  float g=-9.0 ;

  float NvK = (float)wellpops[ wl ] / (float)K ;
  int entry = 0;
  for( int i=0; i<132; i++ ){
    if( NKs[i] > NvK ){ 
      entry = i ;
      i = 200;  //break loop??
    }
  }

  if( NvK == 0 ){ g = max_growth_rate ; }
  else if( NvK >= 1 ) { g = 0; }   // >= since we have no strict upper limit due to migration between wells (... implies unlimited "food" now)
  else{

    if( entry==0 ) {
       cout<<"NvK="<<NvK<<"  entry="<<entry<<" "<<endl;
       for( int gbg=0; gbg<nwells; gbg++){
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
  number=1 ; no_drivers=0 ; sequence.clear() ; prev_gen=-1 ; first_well = 0 ;
}




//------------------------------------------------------------------------------------------------------------
// MODIFIED THIS METHOD SO THAT DRIVER MUTATIONS OCCUR WITH DIFFERENT PROBABILITIES
//
Genotype::Genotype(Genotype *mother, int prevg, int no_snp, int well_num) {

  prev_gen=prevg ;  // i.e. prev_gen of this NEW genotype is the current genotype, cell[i].gen, fed into the method  
  sequence=mother->sequence ; 
  no_drivers=mother->no_drivers; 
  number=1;                     
  first_well = well_num ;  


  int snps_implemented = 0 ;  
  while ( snps_implemented < no_snp ){

    //try driver mutations
    float prob_driver = _drand48();
    int site_to_mutate = -9;

    if     ( prob_driver <  p_g0 ) { site_to_mutate = 0 ; }
    else if( prob_driver < (p_g0 + p_g1)  ) { site_to_mutate = 1 ; }
    else if( prob_driver < (p_g0 + p_g1 + p_g2)  ) { site_to_mutate = 2 ; }
    else if( prob_driver < (p_g0 + p_g1 + p_g2 + p_g3)  ) { site_to_mutate = 3 ; }
    else if( prob_driver < (p_g0 + p_g1 + p_g2 + p_g3 + p_g4)  ) { site_to_mutate = 4 ; }
 
    //if not one of these, just choose another *gene*(?) at random
    if( site_to_mutate >= 0 ){
      int k ;
      for (k=0; k<sequence.size(); k++) if (sequence[k]==site_to_mutate) break ;
      if  ( k==sequence.size() ) { 
        sequence.push_back( site_to_mutate ) ;  
        snps_implemented++ ;  
        no_drivers++ ; 
      }
    }
    else if( site_to_mutate == -9 ){     

      snps_implemented++ ; 

// DO NOTHING SINCE WE DON'T TRACK PMs ANYMORE /////////////////////
//      // choose other random gene
//      int k ;  int j= int(_drand48()*L) ; 
//      for (k=0; k<sequence.size(); k++) if (sequence[k]==j) break ;
//      if (k==sequence.size()  &&  j>=max_no_drivers) {
//        sequence.push_back(j) ;
//        snps_implemented++ ; 
//      }  /////////////////////////////////////////////////////////

    }
    else{
      // start process again
      // snps_implemented++ ;  // uncomment this is I am biasing the driver mutations?!??
    }

  }


}
//------------------------------------------------------------------------------------------------------------------














vector<Genotype*> genotypes ;





// --------------------------------  SAVE AND SORT METHODS --------------------------------------------------
void save_genotypes(char *name, int cutoff)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if ( g!=NULL &&  g->number > cutoff ) {
      fprintf(data,"%d  %d  %d  %d  %d\t", i, g->prev_gen, g->first_well, g->no_drivers, g->number) ;
      for (int j=0;j<g->sequence.size();j++) fprintf(data," %u",g->sequence[j]) ; 
      fprintf(data,"\n") ;
    } 
  }
  fclose(data) ;  
}




// Print all genotypes and their abundance in a specific well --------  THIS IS RIDICULOUSLY SLOW. WHAT AM I DOING WRONG?!
void save_genotypes_per_well(char *name, int cutoff, int well_num)
{
  FILE *data=fopen(name,"w") ;

  vector < unsigned int > gen_list ;
  for (int i=0; i<cells.size(); i++) {
    if( cells[i].well == well_num ){
      int geno = cells[i].gen ;
      //if( geno >= 0 ) {    gen_num_in_well[ geno ]++;  }
      //gen_num_in_well[ geno ]++;
      //if( find( gen_list.begin(), gen_list.end(), geno ) == gen_list.end()   &&  geno >= 0 )  {  
      if( find( gen_list.begin(), gen_list.end(), geno ) == gen_list.end()  )  {  
        gen_list.push_back( geno ) ;  
      }
    }
  }
 
  int *geno_num_in_well = new int [ gen_list.size() ] ;                  // PROBLEM WITH MEMORY!??!
  for( int xx=0; xx<gen_list.size(); xx++ ) { geno_num_in_well[xx] = 0 ; }

  // count cells with each genotype in this well
  for( int gt=0; gt<gen_list.size(); gt++ ){
    for (int i=0; i<cells.size(); i++ ){
      if (cells[i].well == well_num  &&  cells[i].gen==gen_list[gt] ) {
        geno_num_in_well[ gt ]++;
      }
    }
  }

  // print genotypes present in well above the given cutoff  
  for( int a=0; a<gen_list.size(); a++){
    Genotype *gg = genotypes[ gen_list[a] ] ;
    if( gg!=NULL && geno_num_in_well[a] > cutoff ){
      //fprintf(data,"%d  %d  %d  %d  %d\t", gen_list[a], g->prev_gen, g->first_well, g->no_drivers, g->number) ;
      fprintf(data,"%d  %d  %d  %d  %d\t", gen_list[a], gg->prev_gen, gg->first_well, gg->no_drivers, geno_num_in_well[a] ) ;
      for (int j=0;j<gg->sequence.size();j++)   fprintf(data," %u",gg->sequence[j]) ;
      fprintf(data,"\n") ;  
    }
  }
  fclose(data) ;
  delete [] geno_num_in_well ; 
 
}



/**
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
**/




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
    if (n[i]>(5e-3)*total) nsnps++ ;
    if (1.*n[i]/total>cutoff) nsnpsc++ ;
//    cout <<n[i]<<" ";
  }
  float *abund=new float[nsnps], tempd ;
  int *num=new int[nsnps], temp ;
  nsnps=0 ;
  for (i=0;i<L;i++) if (n[i]>(5e-3)*total) { num[nsnps]=i ; abund[nsnps]=float(1.*n[i]/total)*(1+0.000001*i/L) ; nsnps++ ; }
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




// Save number of each strain in given well
void save_strains(char *name, int wl) 
{

  int *nstrain=new int[num_strains] ;
  for (int j=0; j<num_strains ;  j++)   nstrain[j]=0 ;

  FILE *f=fopen(name,"w") ;
  if (f==NULL) err(name) ;

  // Check each cell in well "wl" and save it's strain
  for(int i=0; i<cells.size(); i++){
    
    if( cells[i].well == wl ){     
      // check what strain it is
      int str = genotypes[cells[i].gen]->getStrain( ) ;    
      nstrain[ str ] += 1 ;
    }
  }


  for (int s=0; s<num_strains; s++)  fprintf(f,"%d %f\n", s, float(nstrain[s])/K   ) ;


  fclose(f) ;
  delete [] nstrain ;
}









//####################

void save_ODs()
{
   //  print out ODs as measured in diff.cpp / experiments
   //  (both cubiods of area 1mm^2) 

  int *ODs_perw_2 = new int[nwells] ;

  for (int i=0;i<nwells;i++){
    ODs_perw_2[i] = 0;
  }

  for (int j=0; j<cells.size(); j++){
    if( fabs(cells[j].x)<0.5 && fabs(cells[j].y)<0.5 ) {
        ODs_perw_2[ cells[j].well ]++ ;
    }
  }

  char name2[256] ;
  sprintf(name2,"%s/ODs_wells_1.0.dat",NUM) ;
  FILE *f2=fopen(name2,"a") ;
  fprintf(f2, "%lf ", tt);
  for (int u=0;u<nwells;u++) fprintf(f2,"%d ", ODs_perw_2[u] ) ; fprintf(f2,"\n");
  fclose(f2) ;

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
      for (int j=0;j<genotypes[i]->sequence.size();j++)     snp_no[  (genotypes[i]->sequence[j]) & L_PM  ]+=genotypes[i]->number ;      
    }
  }  
  int snps_det=0 ;
  for (i=0;i<L;i++) if (snp_no[i]>cutoff*ntot) snps_det++ ;
  delete [] snp_no ;

  for (i=0;i<ntot;i++) {
    Genotype *g=genotypes[cells[i].gen] ; if (g==NULL) err("g=NULL") ;
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








// Print out genomic sequences in the final well THAT WAS REACHED BY BACTERIA
void save_sequences(char *name, int well_num)
{
  // Map to hold each sequence
  map <string, int> trajectory_counter;

  //FILE *data=fopen(name,"w") ;

  ofstream myfile (name, ofstream::out);
  //myfile.open();


  for (int i=0; i<cells.size(); i++) {
    if( cells[i].well == well_num ){
      //int geno = cells[i].gen ;
      Genotype *gg = genotypes[ cells[i].gen ] ;  //if (g==NULL) err("g=NULL") ;
      //CreateMap
      string gen_lab ;
      for (int j=0;j<gg->sequence.size();j++){
          int bp = gg->sequence[j] ;
          string bp_s;
          ostringstream convert;
          convert << bp ;
          bp_s = convert.str() ;
          gen_lab = gen_lab + " " + bp_s ;
      }
      
      trajectory_counter[ gen_lab ]++ ;
    }
  }

  // Print each key and value pair
  map<string,int>::iterator iter;   
  for( iter = trajectory_counter.begin(); iter != trajectory_counter.end(); iter++ ) {
    //cout << "gen: " << iter->first << ", count: " << iter->second << endl;
    //fprintf(data, "%d "+iter->first+"\n",  iter->second );
    if (myfile.is_open()){
      myfile << iter->second; 
      myfile << iter->first;   myfile << "\n";
    }
    else cout << "Unable to open file";
  }

  myfile.close();

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
  int i,j,k,n,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;
//######
  double tt_old_2=tt ; 
//######

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


// ONLY ALLOW MOTILITY IF [CIPRO] < 0.5*MIC -------------------------------------
//
// STEP FUNCTION:
    float well_cipro = well[ cells[n].well  ].a ;
    float thisMIC = genotypes[cells[n].gen]->getMIC() ;
    ///cout << thisMIC ;
    ///exit (EXIT_FAILURE) ; 
    int mpf = 1;                                             
    if( well_cipro > motilityCutoff * thisMIC ){ mpf = 0; }                // DEFINE MIC-MOTILITY CUTOFF HERE
//
//
// migration
    float dx,dy,dz ;
    int ok ;
    do { 
      float sigma=sqrt(diff*dt) ;
      dx=gauss()*sigma * mpf ;
      dy=gauss()*sigma * mpf ;
      dz=gauss()*sigma * mpf ;
      ok=where_is_bug(cells[n].x+dx,cells[n].y+dy,cells[n].z+dz) ;
    } while (ok<0) ;
    if (ok==0) { // cell within the same well
      cells[n].x+=dx ; cells[n].y+=dy ; cells[n].z+=dz ;
    } else if (ok==1 && cells[n].well<nwells-1) { // cell in another well to the right
      cells[n].x+=dx-w ; cells[n].y+=dy ; cells[n].z+=dz ;
      cells[n].well++ ;     // cout <<cells[n].well ; 
      // track wellpops:
      wellpops[ cells[n].well ]++;  wellpops[ cells[n].well-1 ]--;
    } else if (ok==2 && cells[n].well>0) { // cell in another well to the left
      cells[n].x+=dx+w ; cells[n].y+=dy ; cells[n].z+=dz ;      
      cells[n].well-- ;  //cout <<cells[n].well ; 
      // track movement
      wellpops[ cells[n].well ]++;  wellpops[ cells[n].well+1 ]--;
    }
//---------------------------------------------------------------------------------


    double q=_drand48()*(max_growth_rate+max_death_rate), br,dr ; //,ml,mr  ;
    int mode=0 ; 
//
    //br=genotypes[cells[n].gen]->growth(cells[n].well) * growth_function(well[cells[n].well].food/food0) ; 
    br=genotypes[cells[n].gen]->growth(cells[n].well) * growth_function( cells[n].well ) ;                 //EDITED FOR LB-GROWTH DATA
//  
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


    if (mode==0) { // reproduction (make 2 new cells, kill mother cell)
      if (_drand48()>genotypes[cells[n].gen]->death_during_rep(cells[n].well)) { // replication

      
        // STEVE- MODIFIED TO INCREASE MUTATION PROB WITH [ANTIBIOTIC] --------------------    
        float antib = well[  cells[n].well  ].a ;

        float mu = gama * ( 1 + pow( X*antib/MIC_WT, Y)  ) ;
        if( mu > gama*max_increase ) { mu = gama*max_increase ; }

        if( mpf==0 ){  mu = 0 ; }           //  <<--------------------------- IF NONMOTILE (mpf==0) DO NOT ALLOW MUTATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!


        int no_SNPs=poisson( mu ) ; // newly produced cell mutants

        Cell c ; c.well=cells[n].well ; c.x=cells[n].x ; c.y=cells[n].y ; c.z=cells[n].z ;

        int no_SNPs_2 = poisson( mu ) ; // newly produced cell mutants
        Cell c_2 ; c_2.well=cells[n].well ; c_2.x=cells[n].x ; c_2.y=cells[n].y ; c_2.z=cells[n].z ;
        // --------------------------------------------------------------------------------

        //add 1 to well population (mother is deleted later so add 2)
        wellpops[ cells[n].well ]+=2;
        //
        //well[c.well].food -= 1.0 ; // Use 1 food per bacteria  NO LONGER USE FOOD BUT N/K DEPENDENCE OF BARTEKS LB DATA

        if (no_SNPs>0) { 
          c.gen=genotypes.size() ; 
          genotypes.push_back(new Genotype( genotypes[cells[n].gen], cells[n].gen, no_SNPs, cells[n].well)  ) ; // mutate 
        } else { 
          c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
        }
        cells.push_back(c) ;

        if (no_SNPs_2 > 0) { 
          c_2.gen = genotypes.size() ;   
          genotypes.push_back(new Genotype( genotypes[cells[n].gen], cells[n].gen, no_SNPs_2, cells[n].well)  ) ; // mutate 
        } else { 
          c_2.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
        }
        cells.push_back(c_2) ;

        // ####   We have produced 2 new cells above, hence now must kill mother cell ##############
        mode=1 ;
        // #########################################################################################

      } else { // death during replication
        mode=1 ;
      }
    }
  
// now we implement death
    if (mode==1) {
      //
      wellpops[ cells[n].well ]--;
      //
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      if (n!=cells.size()-1) { 
        cells[n]=cells[cells.size()-1] ;
      }
      cells.pop_back() ; 
    }



    ntot=cells.size() ;//if (it%1000000==0) cout <<tt<<": "<<ntot<<endl ;
    int ntot_2=0;
    for( int dsa=0; dsa<nwells; dsa++){ ntot_2 += wellpops[ dsa ] ; }
    if (it%1000000==0) cout <<tt<<": "<<ntot<<"   "<<ntot_2<<endl ;

// --- added to produce Fisher plots -----------------------------------------
    double wait_time_2=0.1 ;
    if (tt > tt_old_2+wait_time_2 ) { tt_old_2=tt ; save_well_pops();  save_ODs();}
// ---------------------------------------------------------------------------


    int total_food_remaining = -10;
    double wait_time=1 ;                                       
    if (tt>tt_old+wait_time) { 
       tt_old=tt ; 
       save_data(); 
    
       //total_food_remaining = 0;
      // for( int i=0; i < nwells; i++ ){
      //   total_food_remaining += well[i].food ;
      // }
    }
    // If no food left, end experiment
    //if (total_food_remaining == 0) return 2;

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
    c.z = h * _drand48() * (-0.99) ;
    cells.push_back(c) ; 

    // THE SEGMENTATION FAULT !!
    genotypes[0]->number++ ;
    wellpops[0]++;
  }

  for (int i=0;i<nwells;i++) {
    well[i].food=0 ;
    //well[i].a=uni_antibiotic+pow(max_antibiotic,1.*i/nwells) ;
    well[i].a=uni_antibiotic+pow(max_antibiotic,1.*(i+1)/nwells) ;  // CHANGED SO THAT WELL 0 HAS ZERO ANTIBIOTICS

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
  sprintf(txt,"mkdir %s",NUM) ; system(txt) ;    // what is this system(txt) doing?
  //
  sprintf(txt,"mkdir %s/allPMs",NUM) ; system(txt) ;
  sprintf(txt,"mkdir %s/wellFiles",NUM) ; system(txt) ;
  sprintf(txt,"mkdir %s/strainFiles",NUM) ; system(txt) ;
  //
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

  if( p_g0 + p_g1 + p_g2 + p_g3 + p_g4 > 1.0 ){ cout<<"(p_g0+...+p_g4)>1. Exiting.\n"; exit (EXIT_FAILURE); }


  // PRINT PARAMETERS TO FILE----------------------------------
  char params[256] ;
  sprintf( params, "%s/parameters.txt", NUM);
  FILE *pff = fopen( params,"w") ;
  fprintf(pff, "motilityCutoff = %f\n", motilityCutoff );
  fprintf(pff, "nwells = %d\n", nwells);
  fprintf(pff, "L = %d\n", L);
  fprintf(pff, "max_no_drivers = %d\n", max_no_drivers);
  fprintf(pff, "gama = %f\n", gama);
  fprintf(pff, "diff = %f\n", diff);
  fprintf(pff, "max_growth_rate = %f\n", max_growth_rate);
  fprintf(pff, "tmax = %f\n", tmax);
  fprintf(pff, "init_bact = %f\n", init_bact);
  //fprintf(pff, "food0 = %f\n", food0);
  fprintf(pff, "uni_antibiotic = %f\n", uni_antibiotic);
  fprintf(pff, "max_antibiotic = %f\n", max_antibiotic);
  fprintf(pff, "MIC_WT = %f\n", MIC_WT) ;
  fprintf(pff, "hole_width = %f, hole_depth = %f\n", hole_width, hole_depth ) ;
  fprintf(pff, "well_width = %f, well_depth = %f\n", w, h ) ;
  fprintf(pff, "Parameters for cipro-dep mutation rate:\n") ;
  fprintf(pff, "X = %f, Y = %f, max_increase = %d\n", X, Y, max_increase ) ;
  fprintf(pff, "P_gyrA1 = %f, P_marR = %f", p_g0, p_g3 ) ;
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




// ===================== DISABLED SINCE FOR K=1e8 FILE SIZE MIGHT BE HUGE ========================
//    sprintf(name,"%s/genotypes_ALL_%d_%d.dat",NUM,RAND,sample) ; save_genotypes(name, 0) ;
//
//    sprintf(name,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(name) ;
// ===============================================================================================


    
    int *snp_no=new int[L] ; // array of SNPs abundances
    for (int i=0;i<L;i++) { snp_no[i]=0 ; }
    for (int i=0;i<genotypes.size();i++) {
      if (genotypes[i]!=NULL && genotypes[i]->number>0) 
        for (int j=0;j<genotypes[i]->sequence.size();j++) {
          snp_no[ (genotypes[i]->sequence[j]) & L_PM   ]+=genotypes[i]->number ;      

        }
    }

    printf("saving PMs...\n") ;
    int most_abund[100] ;
    sprintf(name,"%s/allPMs/all_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_no,cells.size(),1,most_abund) ;

    
// save SNPs in different wells
    //int wells_to_save[10]={0,3,5,7,11,13,15,19,21,23} ;
    int wells_to_save[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23} ;

    for (int w=0;  w<24;   w++) {  
      int ninwell=0 ;
      for (int i=0;i<L;i++) { snp_no[i]=0 ; }
      for (int i=0;i<cells.size();i++) {
        if (cells[i].well==wells_to_save[w]) {
          Genotype *g=genotypes[cells[i].gen] ; ninwell++ ;
          //for (int j=0;j<g->sequence.size();j++) snp_no[(g->sequence[j])]++ ;
          for (int j=0;j<g->sequence.size();j++)     snp_no[ (g->sequence[j]) & L_PM   ]++ ;            
        }
      }
      sprintf(name,"%s/wellFiles/well_%d_PMs_%d_%d.dat",NUM,wells_to_save[w],RAND,sample) ; save_snps(name,snp_no,ninwell,1,most_abund) ;
    }
    delete [] snp_no ;


// save strains present in each well
    for (int w=0;  w<nwells;   w++) {  
      sprintf(name,"%s/strainFiles/strains_w%d___%d_%d.dat",NUM,wells_to_save[w],RAND,sample) ;     
      save_strains(name, wells_to_save[w]) ;
    }

// save all trajectories for cells present in each well if wellpop > 10% of carrying capaity
    for (int w=0;  w<nwells;   w++) {  
      if( wellpops[w] > K/10.0 ){
        sprintf(name,"%s/sequences_(well%d)_%d_%d.dat",NUM,w,RAND,sample) ;  
        save_sequences(name,w) ;   
      }
    }






  } 
 




  // Print antibiotic concentration profile 
  char aname[256] ;
  sprintf( aname, "%s/antibiotic.dat", NUM);
  FILE *ff=fopen(aname,"w") ;
  for (int i=0; i<nwells; i++){
      fprintf(ff, "%d %f\n", i, well[i].a);
  }
  fclose( ff );



  end() ;
	return 0 ;
}
