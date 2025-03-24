#include "poisson.h"
#include "fracface_avg.h"
#include "navier-stokes/concen_solver_huang.h"
//----无量纲单位B0=sqrt(NN*rho*u/segmae/L)=sqrt(NN*M*U/(segmae*L^4))~0.2T

mgstats mgphi; 
vector phi[];
face vector sigmae[],sigmae2[];
scalar phib[];
scalar c[];
scalar * stracers = {c};
//scalar * stracers2 = {phi.y};

attribute {
  double D1, D2;  //D1 for f = 1, D2 for f = 0
  scalar phi1, phi2, phi3, phi4; // private
}


/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

//#if TREE
//event defaults (i = 0)
//{
//  for (scalar s in stracers) {
//#if EMBED
//    s.refine = s.prolongation = refine_embed_linear;
//#else
//    s.refine  = refine_linear;
//#endif
 ///   s.restriction = restriction_volume_average;
 // }
//}
//#endif // TREE

/**
## Advection

To avoid numerical diffusion through the interface we use the [VOF
tracer transport scheme](/src/vof.h) for the temporary fields
$\phi_1$ and $\phi_2$, see section 3.2 of [Farsoiya et al.,
2021](#farsoiya2021). */

static scalar * phi_tracers = NULL;

event vof (i++)
{

  phi_tracers = f.tracers;
  for (scalar c in stracers) {
    scalar phi1 = new scalar, phi2 = new scalar;
    c.phi1 = phi1, c.phi2 = phi2;
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    phi2.inverse = true;
    
    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c \frac{\alpha f}{\alpha f + (1 - f)}
    $$
    $$
    \phi_2 = c \frac{1 - f}{\alpha f + (1 - f)}
    $$
    */
		  
    foreach() {
      phi1[] = f[]*phi.x[];
      phi2[] = (1-f[])*phi.x[];
  //     phi1[] = phi.x[];
 //     phi2[] = phi.y[];
    }
    boundary ({phi1, phi2});
  }
  
  
  
    for (scalar c in stracers) {
    scalar phi3 = new scalar, phi4 = new scalar;
    c.phi3 = phi3, c.phi4 = phi4;
    scalar_clone (phi3, c);
    scalar_clone (phi4, c);
    phi4.inverse = true;
    
    f.tracers = list_append (f.tracers, phi3);
    f.tracers = list_append (f.tracers, phi4);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c \frac{\alpha f}{\alpha f + (1 - f)}
    $$
    $$
    \phi_2 = c \frac{1 - f}{\alpha f + (1 - f)}
    $$
    */ 
		  
    foreach() {
      phi3[] = f[]*phi.y[];
      phi4[] = (1-f[])*phi.y[];
  //     phi1[] = phi.x[];
 //     phi2[] = phi.y[];
    }
    boundary ({phi3, phi4});
  }
  
  
  
  //vof_advection({f},i);
  
  
  



}




event tracer_diffusion (i++){ 

  free (f.tracers);
  f.tracers = phi_tracers;
  for (scalar c in stracers) {

    /**
    The advected concentration is computed from $\phi_1$ and $\phi_2$ as
    $$
    c = \phi_1 + \phi_2
    $$
    and these fields are then discarded. */
    
    scalar phi1 = c.phi1, phi2 = c.phi2;
    scalar phi3 = c.phi3, phi4 = c.phi4;
    foreach() {
      phi.x[] = phi1[]+phi2[];
      phi.y[] = phi3[]+phi4[];
      if (cs[]<1 && cs[]>0){
      
    //  phi.x[] = phi1[]/f[];
    //  phi.y[] = phi2[]/(1-f[]);      
      
      }
      
    }
    delete ({phi1, phi2, phi3, phi4});
     boundary ({phi.x, phi.y});
    }


foreach(){
cs[]=f[];
css[]=1.-cs[];
}

face_fraction_C (f,fs);
foreach_face()
fss.x[]=1.-fs.x[];

boundary((scalar*){cs,css,fs,fss});

    foreach_face(){
//    sigmae.x[]=0.08944*fs.x[];
//    sigmae2.x[]=0.08944*(1.-fs.x[]);
   sigmae.x[]=0.08944*fs.x[];
    sigmae2.x[]=0.08944*(1.-fs.x[]);
    }
    boundary((scalar*){sigmae,sigmae2});


mgphi = Phi_jump(phi, sigmae, sigmae2, phi.x, phi.y, cs, dt);


}


