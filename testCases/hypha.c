/* Title: Flow of a drop through a single hypha branch
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# Version 1.0
# Updated: Aug 11, 2024

# Change log: 
* Version 1.0, Aug 11, 2024: allowing the stresses to relax. This allows for the cytoplasm as well as the drop to be viscoelastic.
*/

// 1 is drop, 2 is film and 3 is air

// #include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "src-local/three-phase-nonCoalescing-viscoelastic.h"
#include "src-local/log-conform-viscoelastic.h"
#include "tension.h"
#include "src-local/reduced-three-phase-nonCoalescing.h"

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define AErr (1e-3) // error tolerance in Conformation tensor
#define MINlevel 4 // minimum level
#define tsnap (1e-1)
#define tsnap2 (1e-3)


int MAXlevel;
double tmax;

// moving Drop is assumed Newtonian for now
double Ohd, RhoR_dc, Ec_d, De_d; 

// hypha is modelled as Kelvin-Voigt solid
double RhoR_hc, Ohf, hf, Ec_h, De_h; 

// cytoplasm is assumed Newtonian as well
double Ohc, Ec_c, De_c;

#define gap(x,y,hf,width,x0,c0) (y - ((c0+1e0) + 0.5 * (hf - (c0+1e0)) * (1 + tanh(sq(x-x0) / width))))

double Bond; // driving force
double Ldomain;


int main(int argc, char const *argv[]) {

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  MAXlevel = 9; //atoi(argv[1]);
  tmax = 1e2; //atof(argv[2]);

  // Drop
  Ohd = 1e0; // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  RhoR_dc = 1.2e0;
  Ec_d = 0.0; //atof(argv[3]);
  De_d = 0.0; //atof(argv[4]);

  // Hypha
  Ohf = 1e0;
  hf = 0.90; //atof(argv[5]); // ratio of the gap thickness to the drop radius, far awaay from the drop.
  Ec_h = 1e0; //atof(argv[6]); // Elasto-capillary number: 1e-4 (very soft) to 1e3 (very stiff)
  De_h = 1e30; //atof(argv[7]);
  RhoR_hc = 1e0; // density ratio of hypha to cytoplasm

  // cytoplasm
  Ohc = 1e-2;
  Ec_c = 0.0; //atof(argv[8]);
  De_c = 0.0; //atof(argv[9]);
  
  Bond = 1e0; // Bond number: we keep the driving fixed

  Ldomain = 16.0; // Dimension of the domain: should be large enough to get a steady solution to drop velocity.

  fprintf(ferr, "Level %d tmax %g. Ohd %3.2f, Ec_d %3.2f, De_d %3.2e, Ohc %3.2f, Ec_c %3.2f, De_c %3.2e, Ohf %3.2f, Ec_h %3.2f, De_h %4.3e, hf %3.2f, Bo %3.2f\n", MAXlevel, tmax, Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf, Bond);

  L0=Ldomain;
  X0=-4.0; Y0=0.0;
  init_grid (1 << (MINlevel));
  periodic(right);

  // drop
  rho1 = RhoR_dc; mu1 = Ohd; G1 = Ec_d; lambda1 = De_d;
  
  // Hypha
  rho2 = RhoR_hc; mu2 = Ohf; G2 = Ec_h; lambda2 = De_h;

  // cytoplasm
  rho3 = 1e0; mu3 = Ohc; G3 = Ec_c; lambda3 = De_c;

  Bf1.x = Bond;
  Bf2.x = Bond;

  f1.sigma = 1e0; // drop-cytoplasm interfacial tension
  f2.sigma = 1e0; // hypha-cytoplasm interfacial tension

  run();

}

event init(t = 0){
  if (!restore (file = "restart")) {
    double width = 2e0; // width of the tanh function
    double clearance = 0.20; // clearance from the drop inside the cytoplasm at t = 0
    double x0tanh = 0.0; // midpoint of the tanh function

    refine (y < 1.2 && level < MAXlevel);
    
    fraction(f1, sq(1e0) - sq(y) - sq(x/1.5));
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

  unrefine(x > L0-1e0);
}

// Outputs
event writingFiles (t = 0, t += tsnap; t <= tmax+tsnap) {
  dump (file = "restart");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0, t += tsnap2; t <= tmax+tsnap) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke), reduction(+:vcm), reduction(+:wt)){
    ke += (0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    vcm += (f1[]*u.x[])*sq(Delta);
    wt += f1[]*sq(Delta);
  }
  if (wt > 0.0) vcm /= wt;
  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. Ohd %3.2f, Ec_d %3.2f, De_d %3.2e, Ohc %3.2f, Ec_c %3.2f, De_c %3.2e, Ohf %3.2f, Ec_h %3.2f, De_h %4.3e, hf %3.2f, Bo %3.2f\n", MAXlevel, tmax, Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf, Bond);
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
