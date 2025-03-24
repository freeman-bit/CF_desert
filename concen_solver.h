#include "poisson.h"


mgstats mg_solve2 (struct MGSolve p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
//    embedded_velocity ();
//    embedded_velocity2 ();
//	embedded_phi();
 // fprintf (stderr, "pass\n");
  
    mg_cycle (p.a, res, da, p.relax, p.data,
	      s.nrelax,
	      p.minlevel,
	      grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = p.minlevel;
  
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr, 
	     "WARNING: convergence for %s not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", v.name,
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}






struct Phi_jump {
  vector u;
  face vector mu;
  face vector mu2; 
  scalar rhoj1;
  scalar rhoj2; 
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
  double (* embed_flux_phi) (Point, scalar, scalar,scalar,scalar, vector, vector, vector, double *, int);
};


# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif

// Note how the relaxation function uses "naive" gradients i.e. not
// the face_gradient_* macros.

static void Phi_relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
  struct Phi_jump * p = (struct Phi_jump *) data;
   (const) scalar rho = p->rho;
   double dt = p->dt;
  (const) face vector mu = p->mu;
  (const) face vector mu2 = p->mu2;
  vector u = vector(a[0]), r = vector(b[0]);

   double (* embed_flux_phi) (Point, scalar, scalar,scalar,scalar, vector, vector, vector, double *, int) = p->embed_flux_phi;
  
  // embedded_phi();

   
  foreach_level_or_leaf (l) {
    double avgmu = 0.;
    foreach_dimension()
      avgmu += mu.x[] + mu.x[1];
    avgmu = dt*avgmu + SEPS;
    
    double avgmu2 = 0.;
    foreach_dimension()
      avgmu2 += mu2.x[] + mu2.x[1];
    avgmu2 = dt*avgmu2 + SEPS;    
    
  
      double c = 0.;
      scalar s1 = u.x;
      scalar s2 = u.y;  
      int flag=1;      
      double d = embed_flux_phi (point, s1, s2, cs, css, fs, fss, mu, &c, flag);

      double a = 0.;
      foreach_dimension()
	a += mu.x[1]*s1[1] + mu.x[]*s1[-1];
      u.x[] = (dt*a + (r.x[] - dt*c)*sq(Delta))/
	(sq(Delta)*(rho[]*lambda.x + dt*d) + avgmu);
    
    
    
      double c2 = 0.;
      flag=0;
      double d2 =  embed_flux_phi (point, s2, s1, css, cs, fss, fs, mu2, &c2, flag);
      double a2 = 0.;
      foreach_dimension()
	a2 += mu2.x[1]*s2[1] + mu2.x[]*s2[-1];
      u.y[] = (dt*a2 + (r.y[] - dt*c2)*sq(Delta))/
	(sq(Delta)*((1.-rho[])*lambda.y + dt*d2) + avgmu2);   
  }
  
#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

static double Phi_residual_diffusion (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Phi_jump * p = (struct Phi_jump *) data;
  (const) face vector mu = p->mu;
  (const) face vector mu2 = p->mu2;
  (const) scalar rho = p->rho;
  double dt = p->dt;
double (* embed_flux_phi) (Point, scalar, scalar,scalar,scalar, vector, vector, vector, double *, int)= p->embed_flux_phi;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;

  /* conservative coarse/fine discretisation (2nd order) */
#if TREE
    scalar s1 = u.x;
    scalar s2 = u.y;   
    face vector g[];
    foreach_face()
      g.x[] = mu.x[]*face_gradient_x (s1, 0);
    boundary_flux ({g});
    foreach (reduction(max:maxres)) {
      double a = 0.;
      foreach_dimension()
	a += g.x[] - g.x[1];
      res.x[] = r.x[] - rho[]*lambda.x*u.x[] - dt*a/Delta;
      int flag = 1;
	double c, d =  embed_flux_phi (point, u.x, u.y, cs, css, fs, fss, mu, &c, flag);
	res.x[] -= dt*(c + d*u.x[]);
      
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  ///////
  
  

    face vector g2[];
    foreach_face()
      g2.x[] = mu2.x[]*face_gradient_x (s2, 0);
    boundary_flux ({g2});
    foreach (reduction(max:maxres)) {
      double a2 = 0.;
      foreach_dimension()
	a2 += g2.x[] - g2.x[1];
      res.y[] = r.y[] - (1.-rho[])*lambda.y*u.y[] - dt*a2/Delta;
      int flag = 0;
	double c2, d2 = embed_flux_phi (point, u.y, u.x, css, cs, fss, fs, mu2, &c2, flag);
	res.y[] -= dt*(c2 + d2*u.y[]);     
      if (fabs (res.y[]) > maxres)
	maxres = fabs (res.y[]);
    }
  
  
  //////
  
  
  boundary (resl);
  
#else


    /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres))
  
      scalar s = u.x;
      scalar s2 = u.y;
      double a = 0.;
      foreach_dimension()
	a += mu.x[0]*face_gradient_x (s1, 0) - mu.x[1]*face_gradient_x (s1, 1);
      res.x[] = r.x[] - rho[]*lambda.x*u.x[] - dt*a/Delta;
	int flag = 1;
	double c, d =  embed_flux_phi (point, u.x, u.y, cs, css, fs, fss, mu, &c, flag);
	res.x[] -= dt*(c + d*u.x[]);
      
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
	
	
 
      double a2 = 0.;
      foreach_dimension()
	a2 += mu2.x[0]*face_gradient_x (s2, 0) - mu2.x[1]*face_gradient_x (s2, 1);
      res.y[] = r.y[] - (1.-rho[])*lambda.x*u.y[] - dt*a2/Delta;
  	flag=0;
	double c2, d2 =  embed_flux_phi (point, u.y, u.x, css, cs, fss, fs, mu2, &c2, flag);
	res.y[] -= dt*(c2 + d2*u.y[]);
      
      if (fabs (res.y[]) > maxres)
	maxres = fabs (res.y[]);	

  
#endif
  return maxres;
}

#undef lambda

double TOLERANCE_PHI = 0; // default to TOLERANCE

trace
mgstats Phi_jump (struct Phi_jump p)
{
  vector u = p.u, r[];
  scalar rhoj1 = p.rhoj1;
  scalar rhoj2 = p.rhoj2; 
  scalar rho = p.rho; 
  foreach(){
      r.x[] = rho[]*rhoj1[];
      r.y[] = (1.-rho[])*rhoj2[];
      
  //    r.x[] = rhoj1[];
  //    r.y[] = rhoj2[];
  //   r.x[] = 0.;
  //   r.y[] = 0.;
      }
     boundary((scalar*){r});

  face vector mu = p.mu;
  face vector mu2 = p.mu2;
  restriction ((scalar*){mu,mu2,rhoj1,rhoj2});

  p.embed_flux_phi = embed_flux_phi;
  return mg_solve2 ((scalar *){u}, (scalar *){r},
		   Phi_residual_diffusion, Phi_relax_diffusion, &p, p.nrelax, p.res,
		   minlevel = 1, // fixme: because of root level
                                  // BGHOSTS = 2 bug on trees
		   tolerance = TOLERANCE_PHI);
}
