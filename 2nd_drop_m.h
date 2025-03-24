
#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
///////////////////////////////////////////////////////////////////////////////////////////计算v1
foreach_dimension()
 double embedded_velocity1_normal_x (Point point, scalar s,        
                                              scalar cs, face vector fs, coord n, coord p) 
{
 foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
return v[0];
}

double embedded_velocity1_normal (Point point, scalar s,        
                                              scalar cs, face vector fs, coord n, coord p)  
{  
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return embedded_velocity1_normal_x (point, s, cs, fs, n, p) ;  
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return embedded_velocity1_normal_x (point, s, cs, fs, n, p) ;         
                                          
  }
  else if (fabs(n.y) >= fabs(n.z))
    return embedded_velocity1_normal_y (point, s, cs, fs, n, p) ;       
                                          
  return embedded_velocity1_normal_z (point, s, cs, fs, n, p) ;       
                                           
#endif // dimension == 3
  return nodata;
}


////////////////////////////////////////////////////////////////////////////////////////////////////计算v2

foreach_dimension()
 double embedded_velocity2_normal_x (Point point, scalar s,        
                                              scalar cs, face vector fs, coord n, coord p)        
{
 foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
return v[1];
}

double embedded_velocity2_normal (Point point, scalar s,        
                                              scalar cs, face vector fs, coord n, coord p)  
{  
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return embedded_velocity2_normal_x (point, s, cs, fs, n, p) ;  
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return embedded_velocity2_normal_x (point, s, cs, fs, n, p) ;         
                                          
  }
  else if (fabs(n.y) >= fabs(n.z))
    return embedded_velocity2_normal_y (point, s, cs, fs, n, p) ;       
                                          
  return embedded_velocity2_normal_z (point, s, cs, fs, n, p) ;       
                                           
#endif // dimension == 3
  return nodata;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////计算d1，d2（已*Delta）
foreach_dimension()
double distance1_x (Point point, scalar cs, face vector fs, coord n, coord p) 
{
     foreach_dimension()
    n.x = - n.x;
    double d[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      }
  return d[0]*Delta;

}
foreach_dimension()
double distance2_x (Point point, scalar cs, face vector fs, coord n, coord p) 
{
     foreach_dimension()
    n.x = - n.x;
  double d[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      }
  return d[1]*Delta;
  
}

double distance1 (Point point,scalar cs, face vector fs, coord n,coord p) 
{ 
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return distance1_x (point,cs,fs,n,p);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return distance1_x (point,cs,fs,n,p);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return distance1_y (point,cs,fs,n,p);
  return distance1_z (point,cs,fs,n,p);
#endif // dimension == 3
  return nodata;
}

double distance2 (Point point,scalar cs, face vector fs,coord n,coord p) 
{ 
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return distance2_x (point,cs,fs,n,p);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return distance2_x (point,cs,fs,n,p);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return distance2_y (point,cs,fs,n,p);
  return distance2_z (point,cs,fs,n,p);
#endif // dimension == 3
  return nodata;
}






