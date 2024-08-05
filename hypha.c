/* Title: Flow of a drop through a single hypha branch
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# Version 0.1
# Updated: Aug 5, 2024
*/

// 1 is drop, 2 is film and 3 is air

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "src-local/three-phase-nonCoalescing-elastic.h"
#include "src-local/log-conform-elastic.h"
#include "tension.h"

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define AErr (1e-3) // error tolerance in Conformation tensor
#define MINlevel 4 // minimum level
#define tsnap (0.01)


// boundary conditions
u.n[right] = neumann(0.0);
p[right] = dirichlet(0.0);
u.n[left] = neumann(0.0);
p[left] = dirichlet(0.0);

int MAXlevel;
double tmax;

// moving Drop is assumed Newtonian for now
double Ohd; 
// hypha is modelled as Kelvin-Voigt solid
double Ohf, hf, Ec; // De is \infty
// cytoplasm is assumed Newtonian as well
double Ohc = 1e-3;
#define gap(x,y,hf,width,x0,c0) (y - ((1e0+c0) + 0.5 * (hf - (1e0+c0)) * (1 + tanh((x-x0) / width))))

double Bond; // driving force
double Ldomain;


int main(int argc, char const *argv[]) {

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  MAXlevel = 10; //atoi(argv[1]);
  tmax = 2; //atof(argv[2]);

  // Drop
  Ohd = 1e0; // <0.000816/sqrt(816*0.017*0.00075) = 0.008>

  // Hypha
  Ohf = 1e0;
  hf = 0.5; //atof(argv[3]); // ratio of the gap thickness to the drop radius, far awaay from the drop.
  Ec = 1e-1; //atof(argv[4]); // Elasto-capillary number: 1e-4 (very soft) to 1e3 (very stiff)
  
  Bond = 1e1; // Bond number: we keep the driving fixed

  Ldomain = 16.0; // Dimension of the domain: should be large enough to get a steady solution to drop velocity.

  fprintf(ferr, "Level %d tmax %g. Ohd %g, Ohf %3.2e, hf %3.2f, Ec %3.2f, Bo %3.2f, De infty \n", MAXlevel, tmax, Ohd, Ohf, hf, Ec, Bond);

  L0=Ldomain;
  X0=-2e0; Y0=0.0;
  init_grid (1 << (MINlevel));

  // drop
  rho1 = 1e0; mu1 = Ohd; G1 = 0.;
  // Hypha
  rho2 = 1e0; mu2 = Ohf; G2 = Ec;
  // cytoplasm
  rho3 = 1e-1; mu3 = Ohc; G3 = 0.;

  f1.sigma = 1.0; // drop-cytoplasm interfacial tension
  f2.sigma = 1e-1; // hypha-cytoplasm interfacial tension

  run();

}

event acceleration(i++){
  face vector av = a;
  foreach_face(x){
    av.x[] += Bond*f1[];
  }
}

event init(t = 0){
  if (!restore (file = "restart")) {
    double width = 2.5e-1; // width of the tanh function
    double clearance = 0.01; // clearance from the drop inside the cytoplasm at t = 0
    double x0tanh = 0.70; // midpoint of the tanh function

    refine (y < 1.2 && level < MAXlevel);
    
    fraction(f1, sq(1e0) - sq(y) - sq(x));
    fraction(f2, gap(x,y,hf,width,x0tanh,clearance));
  }
}

event adapt(i++){
  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);

  adapt_wavelet ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, conform_p.x.x, conform_p.x.y, conform_p.y.y},
  (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, AErr, AErr, AErr}, 
  MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  dump (file = "restart");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke), reduction(+:vcm), reduction(+:wt)){
    ke += (0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    vcm += (f1[]*u.y[])*sq(Delta);
    wt += f1[]*sq(Delta);
  }
  if (wt > 0.0) vcm /= wt;
  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. Ohd %g, Ohf %3.2e, hf %3.2f, Ec %3.2f, Bo %3.2f, De infty \n", MAXlevel, tmax, Ohd, Ohf, hf, Ec, Bond);
      fprintf (fp, "i dt t ke vcm\n");
    } else {
      fp = fopen ("log", "a");
    }
    fprintf (fp, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
    fclose(fp);
    fprintf (ferr, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
  }
  assert(ke > -1e-10);
  // assert(ke < 1e2);
  // dump(file = "dumpTest");
}
