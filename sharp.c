//#include "grid/quadtree.h"
//#include "axi.h"

scalar cs[], css[];
face vector fs[], fss[];
int pic=0;

#define SigR 1
#define kp 0.5
//#define CDT 0.0002

#include "slipboundary/2nd_drop_m.h"

#include "slipboundary/embed3D_concen.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"


//#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "navier-stokes/concen.h"
//#include "henry.h"
//#include "view.h"

double bubble_radius = 2.;
double box_size = 20.;
//double conc_liq1 = 0, conc_gas1 = 1.;
double end_time[4] = {7.5,3,3,2};	
	
int MAXLEVEL, dcase = 1;

int main (int argc, char **argv)
{
  size (box_size);
  MAXLEVEL = 9;

  N = 512;

//  CFL=0.3;
  rho1 = 1.;
  rho2 = 0.01;
  TOLERANCE = 1e-4;
 
//   phi.y.inverse = true;
    
 

  for (dcase = 1; dcase <= 4; dcase++) {
    switch (dcase) {
    case 1:
      f.sigma = 10.0;	
  // f.sigma = 100.0;
  //   G.x = - 2.5;
    G.x = 0.;	
      break;
    case 2:
      f.sigma = 10.0;	
      G.x = - 7.8125;	
      break;
    case 3:
      f.sigma = 1.0;	
      G.x = - 10.0;
      break;
    case 4:
      f.sigma = 10.0;	
      G.x = - 7.8125;
      break;
    default:
      fprintf (stderr, "Error: must specify case\n");
      exit (1);
    }

  mu1 = 0.08944;
  mu2 = 0.08944/20.;
    
    fprintf (stderr, "\n\n# case %d\n", dcase);
    run();
  }
  
  }


event init (t = 0)
{
//  refine (sq(2.*bubble_radius) - sq(x - box_size*0.2) - sq(y-10) > 0 &&
//	  level < MAXLEVEL);
  
  
    vertex scalar phi2[], phi_inverse[];
  foreach_vertex(){
     phi2[] = sq(x-box_size*0.2) + sq(y-10) - sq(bubble_radius);
    phi_inverse[] = -(sq(x-box_size*0.2) + sq(y-10) - sq(bubble_radius));
    }
  boundary ({phi2,phi_inverse});
  fractions (phi2, cs, fs);

  foreach(){
  f[]=cs[];
  css[]=1.-cs[];
  }
  foreach_face()
  fss.x[]=1.-fs.x[];
  boundary((scalar*){f,css,fss});
  
   
  foreach(){
    phi.x[] = 0.;
     phi.y[] = 0.;
   if (cs[]<1)
   phi.y[] = 1.;
   if (cs[]>0)
    phi.x[] = 0.;
    c[]=1.-cs[];
    }
    
    boundary((scalar*){phi,c});
    
    foreach_face(){
    sigmae.x[]=0.08944*fs.x[];
    sigmae2.x[]=0.08944*(1.-fs.x[]);
    }
    boundary((scalar*){sigmae,sigmae2});
  
  
    char outfile11[100];
      sprintf(outfile11, "phiy-%d.ppm",pic);
  FILE * fp_interface2 = fopen (outfile11, "w");
  output_ppm (fs.x, fp_interface2, linear=true, min = 0., max = 1., n=1024);
  
    fclose(fp_interface2);
  pic++;

}

//#if TREE
//event adapt (i++) {
//  adapt_wavelet ({f, u}, (double[]){0.01,0.01,0.01},
//		 maxlevel = MAXLEVEL, minlevel=6); 
//}
//#endif
	
	

event pictures (t+=0.1)
//event pictures (i++)
{

    char outfile11[100];
      sprintf(outfile11, "phiy-%d.ppm",pic);
  FILE * fp_interface2 = fopen (outfile11, "w");
  output_ppm (phi.y, fp_interface2, linear=true, min = 0, max = 1., n=512);
 
    fclose(fp_interface2);
    
    
    
        char outfile12[100];
      sprintf(outfile12, "phix-%d.ppm",pic);
  FILE * fp_interface3 = fopen (outfile12, "w");
  output_ppm (phi.x, fp_interface3, linear=true, min = 0, max = 1., n=512);
 
    fclose(fp_interface3);
    
    
    
           char outfile13[100];
      sprintf(outfile13, "f-%d.ppm",pic);
  FILE * fp_interface4 = fopen (outfile13, "w");
  output_ppm (fs.y, fp_interface4, linear=true, min = 0, max = 1., n=512);
 
    fclose(fp_interface4);
   
   
 //  dump("1");

pic++;
}


event extract (t +=0.2)
{	
char dumpname[100];
      sprintf(dumpname, "dump-%g",t);
dump(file=dumpname);

}

event pictures2 (t = end)
{

}
