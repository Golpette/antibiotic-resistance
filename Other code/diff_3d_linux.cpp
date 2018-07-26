#include <stdio.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#define NUM "test"
const int nn=15 ; // no of wells (5 is enough)
const float wellw=3.5,wellh=11.5,                      wellsep=1.0,  holew=1.5,   holeh=4.0 ; // sizes of a single well in [mm]
const float lx=nn*wellw+wellsep*(nn-1), ly=wellw, lz=wellh ; // sizes in [mm]
const float dd=2e-4 ; // diffusion constant in mm^2/s
const float dt=20; // temporal resolution [s]
const float dx=0.25 ; // spatial resolution [mm]

int wx=lx/dx, wy=ly/dx, wz=lz/dx ; 

void err(char *txt)
{
  printf("%s\n",txt) ;
  //exit(0) ;
}
  

double drand48(void)  // works only on compilers with long long int!
{
  static long long int x=0x000100010001LL, mn=0x0005deece66dLL, dod=0xbLL ;
  x=mn*x+dod ; x&=0xffffffffffffLL ;
 return (x/281474976710656.0) ;
}


void init();
void main_proc();
void end() ;

double tt ;

int main()
{
  init() ;
  do { main_proc() ; } while (tt<3*24*3600) ; // how long to simulate
  end() ; 
}



double ***c, ***dc ;
FILE *out ;




void init()
{
  int i,j,k,n;
  
  c=new double**[wx] ; if (c==NULL) err("cannot allocate p[][][]") ;
  dc=new double**[wx] ; if (dc==NULL) err("cannot allocate p[][][]") ;
  for (i=0;i<wx;i++) { 
    c[i]=new double*[wy] ; if (c[i]==NULL) err("out of memory") ;
    dc[i]=new double*[wy] ; if (dc[i]==NULL) err("out of memory") ;
    for (j=0;j<wy;j++) { 
      c[i][j]=new double[wz] ; if (c[i][j]==NULL) err("out of memory") ; 
      dc[i][j]=new double[wz] ; if (dc[i][j]==NULL) err("out of memory") ; 
    }
  }    


  for (i=0;i<wx;i++) 
    for (j=0;j<wy;j++) 
      for (k=0;k<wz;k++) c[i][j][k]=1e-12 ; 

// initial condition:
  for (i=1/dx;i<(wellw-1)/dx;i++)                                       // IC: so ignore i={0,1,2,3} and last 4 dx's (same for dy)
    for (j=1/dx;j<(wellw-1)/dx;j++)                                     //   i.e. IC: no bacteria 1mm from each wall in x-y
      for (k=0;k<wz;k++) c[i][j][k]=1 ;

  for (i=0;i<wx;i++) for (j=0;j<wy;j++) c[i][j][0]=c[i][j][wz-1]=-1 ;   // STEVE: Value of (-1) signifies nothing can be here? Thickness of actual wall?? 
  for (i=0;i<wx;i++) for (j=0;j<wz;j++) c[i][0][j]=c[i][wy-1][j]=-1 ;   // this means well is dx*2 (0.5mm) smaller in all directions?? This is why density at tt=0 is
  for (i=0;i<wy;i++) for (j=0;j<wz;j++) c[0][i][j]=c[wx-1][i][j]=-1 ;   // not 1 (but actually 44/46) since first and last dx pieces in z-axis are empty.


// make separating walls
  for (n=0;n<nn-1;n++) 
    for (  i=n*(wellw+wellsep)/dx + wellw/dx;   i<n*(wellw+wellsep)/dx+(wellw+wellsep)/dx;    i++) 
      for (k=0;k<wz;k++) 
        for (j=0;j<wy;j++) 
          if (j<wy/2-holew/(2*dx) || j>wy/2+holew/(2*dx) || k<wz-holeh/dx  )   c[i][j][k]=-1 ;



  char name[256] ;
  //string hole_width = to_string( holew ) ;  
  //string hole_height = to_string( holeh ) ; 
  //sprintf(name, "ODs_%f_%f.dat", hole_width, hole_height) ;
  sprintf(name, "ODs_%f_%f_wellsep_%f.dat", holew, holeh, wellsep) ;
  out=fopen(name,"w") ;
  //out=fopen("ODs.dat","w") ; 

    
  tt=0 ;    

}




void end() {
  fclose(out) ; 
}

//-----------------------------------------------------

void main_proc()
{
  int i,j,k,n,l;
  const int kx[6]={-1,1,0,0,0,0},
            ky[6]={0,0,-1,1,0,0},
            kz[6]={0,0,0,0,-1,1};
  


  if( tt==0 ){
    // save ODs  AT START ~~~~~
    fprintf(out,"%f\t",tt/3600.) ;
    for (i=0;i<wx;i++) {
      float cc=0 ;
      n=0 ;
      for (k=0;k<wz;k++) if (c[i][wy/2][k]>=0) { cc+=c[i][wy/2][k] ; n++ ; }                      
      if ((i%int((wellw+wellsep)/dx))==int(0.5*wellw/dx)) fprintf(out,"%lf ",cc/(wz-2)) ;      
    }
    fprintf(out,"\n") ;
  }//SC ~~~~~~~~~~~~~~~~~~~~~



  //int maxst=18 ; 
  int maxst=100 ; 
  for (int st=0;st<maxst;st++) {      
    for (i=1;i<wx-1;i++) 
      for (j=1;j<wy-1;j++)
        for (k=1;k<wz-1;k++) {
          int nn=0 ;
          dc[i][j][k]=0 ;
          for (n=0;n<6;n++) 
            if (c[i+kx[n]][j+ky[n]][k+kz[n]]>=0) {
              nn++ ; 
              dc[i][j][k]+=c[i+kx[n]][j+ky[n]][k+kz[n]] ;
            } 
          dc[i][j][k]-=nn*c[i][j][k] ;  
          dc[i][j][k]*=dd/(dx*dx) ;
        }
    for (i=1;i<wx-1;i++) 
      for (j=1;j<wy-1;j++)
        for (k=1;k<wz-1;k++) 
          if (c[i][j][k]>=0) c[i][j][k]+=dt*dc[i][j][k] ;
    tt+=dt ;
  }




  // save ODs        
  fprintf(out,"%f\t",tt/3600.) ;
  for (i=0;i<wx;i++) {
    float cc=0 ;
    n=0 ;
    for (k=0;k<wz;k++) if (c[i][wy/2][k]>=0) { cc+=c[i][wy/2][k] ; n++ ; }       
    if ((i%int((wellw+wellsep)/dx))==int(0.5*wellw/dx)) fprintf(out,"%lf ",cc/(wz-2)) ;   //ok, so initial OD is (46-2)*1 - i.e. (46-2) is our normalisation factor
  }
  fprintf(out,"\n") ;



}
