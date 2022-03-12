
/**
# A tsunami in mediterranean sea

A clone of Popinet's  tsunami  code (see links)
 
 

## Solver setup

The following headers specify that we use spherical coordinates anthe [Saint-Venant solver](/src/saint-venant.h) together with
[(dynamic) terrain reconstruction](/src/terrain.h). We will use [inputs](/src/input.h) only
when restarting a simulation from a snapshot. */

#include "spherical.h"
#include "saint-venant.h"
#include "terrain.h"
#include "input.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 7 //8   
#define MINLEVEL 4 //5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define HMAXE    5e-2 // error on maximum free surface elevation (5 cm)


int main()
{
  /**
  Here we setup the domain geometry. We choose to use metre as length
  unit, so we set the radius of the Earth (required for the [spherical
  coordinates](/src/spherical.h)) in metres. The *x* and *y*
  coordinates are longitude and latitude in degrees, so we set the
  size of the box *L0* and the coordinates of the lower-left corner
  *(X0,Y0)* in degrees.
   the domain is 20 degrees squared,centered on 10,36 longitude,latitude
   */

  Radius = 6371220.;
  size (20);
  origin (0, 30);

  /**
  *G* is the acceleration of gravity required by the Saint-Venant
  solver. This is the only dimensional parameter. We rescale it so that
  time is in minutes. */

  // acceleration of gravity in m/min^2
  G = 9.81*sq(60.);

  /**
  When using a quadtree (i.e. adaptive) discretisation, we want to start
  with the coarsest grid, otherwise we directly refine to the maximum
  level. Note that *1 << n* is C for $2^n$. */


  N = 1 << MAXLEVEL;


  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */
  run();
}

/**
We declare and allocate another scalar field which will be used to
store the maximum wave elevation reached over time. */

scalar hmax[];

/**
## Boundary conditions

We set the normal velocity component on the left, right and bottom
boundaries to a "radiation condition" with a reference sealevel of
zero. The top boundary is always "dry" in this example so can be left
alone. Note that the sign is important and needs to reflect the
orientation of the boundary. Any way, in the mediterranean sea, it is  not a problem
 */

u.n[left]   = - radiation(0);
u.n[right]  = + radiation(0);
u.n[bottom] = - radiation(0);


/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two *#if...#else* branches selecting
whether the simulation is being run on an (adaptive) quadtree or a
(static) Cartesian grid.

We want to adapt according to two criteria: an estimate of the error
on the free surface position -- to track the wave in time -- and an
estimate of the error on the maximum wave height *hmax* -- to make
sure that the final maximum wave height field is properly resolved.

We first define a temporary field (in the
[automatic variable](http://en.wikipedia.org/wiki/Automatic_variable)
*η*) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int adapt() {
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  /**
  We can now use wavelet adaptation on the list of scalars *{η,hmax}*
  with thresholds *{ETAE,HMAXE}*. The compiler is not clever enough yet
  and needs to be told explicitly that this is a list of *double*s,
  hence the *(double[])*
  [type casting](http://en.wikipedia.org/wiki/Type_conversion). 
  
  The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({eta, hmax}, (double[]){ETAE,HMAXE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
}

/**
## Initial conditions

We first specify the terrain database to use to reconstruct the
topography $z_b$. This KDT database needs to be built beforehand. See the
[*xyz2kdt* manual](http://gfs.sourceforge.net/wiki/index.php/Xyz2kdt)
for explanations on how to do this.

 
 [https://www.ngdc.noaa.gov/mgg/global/etopo2.html]()
 
 The horizontal grid spacing is 2-minutes of latitude and longitude (1 minute of latitude = 1.853 km at the Equator). The vertical precision is 1 meter.
 


The next line tells the Saint-Venant solver to conserve water surface
elevation rather than volume when adapting the mesh. This is important
for tsunamis since most of the domain will be close to "lake-at-rest"
balance. */


event init (i = 0)
{
  

  terrain (zb, "./topoT", NULL);
  conserve_elevation();

  /**
   The initial still water surface is at $z=0$ (The horizontal datum is WGS-84, the vertical datum is Mean Sea Level) so that the water depth
  $h$ is... */

  foreach()
    h[] = max(0., - zb[]);
  boundary ({h});

  /**
  For real applications, the initial deformation is given by the Okada model.
   
   [https://www.clawpack.org/geoclaw/Okada.html]()
   
   
   we replace it by a instantaneous variation of water of 5 m in a circle of radius .5 degree (1/2 NMile)
    */
    
    foreach()
    h[] = h[] +  7*(sqrt(sq(x-5.5)+sq(y-40.26)) < .6) ;
    boundary ({h});
    
   
}

/**
## Outputs

### At each timestep

We output simple summary statistics for *h* and *u.x* on standard
error. */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min,  s.max, s.sum, n.rms, n.max, dt);

  /**
  We also use a simple implicit scheme to implement quadratic bottom
  friction i.e.
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=10^{-4}$. */
  
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
      boundary ({u.x,u.y});
  }
}

/**
 
![Maximum wave elevation (metres) reached over 10 hours.](tsunami/plot.png)

### Movies

This is done every minute (*t++*). The static variable *fp* is *NULL*
when the simulation starts and is kept between calls (that is what
*static* means). The first time the event is called we set *fp* to a
*ppm2mpeg* pipe. This will convert the stream of PPM images into an
mpeg video. 

We use the *mask* option of *output_ppm()* to mask out the dry
topography. Any part of the image for which *m[]* is negative
(i.e. for which *etam[] < zb[]*) will be masked out. */


event photo ( i = 0 ){
  double zbmin = 0,zbmax = 0;
  foreach(){
    if ( zb[] < zbmin )
      zbmin = zb[];
    if ( zb[] > zbmax )
      zbmax = zb[];
  }
  fprintf(stderr,"  z max  = %lf z min %lf \n", zbmax,zbmin);
  static FILE * fzb = fopen ("topo.ppm", "w");
  output_ppm (zb, fzb, min = zbmin, max = 0, n = 1 << MAXLEVEL, linear = true);
  static FILE * fzt = fopen ("terre.ppm", "w");
  output_ppm (zb, fzt, min = 0, max = 500, n = 1 << MAXLEVEL, linear = true);
}

event movies (t++;t<300) {
  //static FILE * fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam,  file = "eta.mp4", mask = m, min = -3, max = 3, n = 512, linear = true);

  /**
  After completion this will give the following animation
  
  ![[Animation](tsunami/eta.mpg) of the wave elevation. Dark blue is -2 metres
  and less. Dark red is +2 metres and more.](tsunami/eta.png)
  
  We also use the *box* option to only output a subset of the domain
  (defined by the lower-left, upper-right coordinates). */
  
  //static FILE * fp2 = popen ("ppm2mpeg > eta-zoom.mpg", "w");
  output_ppm (etam, file =  "eta-zoom.mp4",   min = -2, max = 2, n = 512, linear = true,
	      box = {{3,41},{8,45}});

  /**
  ![[Animation](tsunami/eta-zoom.mpg) of the wave elevation. Dark blue is 
  -2 metres and less. Dark red is +2 metres and more.](tsunami/eta-zoom.png)
  
  And repeat the operation for the level of refinement...*/

 // static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l,  file = "level.mp4" , min = MINLEVEL, max = MAXLEVEL, n = 512);

  /**
  ![[Animation](tsunami/level.mpg) of the level of refinement. Dark blue is 5
  and dark red is 10.](tsunami/level.png)
    */
  
}

/**

### Tide gauges

We define a list of file names, locations and descriptions and use the
*output_gauges()* function to output timeseries (for each timestep) of
$\eta$ for each location. */

Gauge gauges[] = {
  // file   lon      lat 
  {"Nice.ga", 7.21,  42.05},//7.21,  43.65
  {"Montpellier.ga", 3.93,  42.0},//3.93,  43.53
  {"Perpignan.ga", 3.04,  42.70},//3.04,  42.70
  {"Cannes.ga", 7.03,    43.54},//7.03,    43.54
  {"Marseille.ga", 5.36,    43.30},//5.36,    43.30
  {"Calvi.ga", 8.77,    42.57},//8.77,    42.57
  {"Bastia.ga", 9.5,    42.70},//9.5,    42.70
  {"Livorno.ga", 10.29,    43.54},//10.29,    43.54
  {"Napples.ga", 14.27,    40.83},//14.27,    40.83
  {"Tunis.ga", 11.9,    37.26},//10.34,    36.84
  {"Epicentre.ga", 5.5, 40.26},//11.9, 37.26
  {NULL}
};
event gauges1 (i++) output_gauges (gauges,{eta});
/**
 
 
### Snapshots

Every 60 minutes, the $h$, $z_b$ and *hmax* fields are interpolated
bilinearly onto a *n x n* regular grid and written on standard
output. */


event gauges1 (i++) output_gauges (gauges, {eta});
/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();
    
/**
 
## RUN
  
 
 
 IMPORTANT : il faut rajouter le chemin de basilisk/src/kdt/kdt.o quand on compile.
 
~~~bash
 
 qcc -O2 -options... tsutunis.c  ../../src/kdt/kdt.o -lm
~~~
 
example with *fopenmp*
 
 
~~~bash

 qcc -fopenmp  -g -O2 -DTRASH=1 -Wall  tsutunis.c -o tsutunis /Users/pyl/basilisk/src/kdt/kdt.o -lm
~~~
 
with *qcc*

~~~bash
 
 qcc tsutunis.c  ~/basilisk/src/kdt/kdt.o -lm -o tsutunis
~~~

 
You may have to process the topo with
~~~bash
 
 kdt2kdt -v ./topoT
 
~~~
 
  on raspberry, you have to update the topo with *kdt*, then compile

~~~bash
 
 ./gerris/gerris-snapshot-131206/modules/kdt/kdt2kdt /home/pi/Desktop/tsunami_tunis/topoT

 qcc tsutunis.c /home/pi/basilisk/src/kdt/kdt.o -lm -o tsutunis 

 
 

 
~~~ 

# links
 
 * basilsk tutorial
 
 * tsunami lisbonne 1755
 
 * tsunami 2004 Indian Ocean earthquake and tsunami
 
 * tsunami 2011 Tōhoku earthquake and tsunami
 
 * tsunami [https://en.wikipedia.org/wiki/List_of_tsunamis]()
 
 * marees bretagne
 
 * [http://gfs.sourceforge.net/wiki/index.php/Installing_from_source]()
 
 
# bib

 
 
 */

