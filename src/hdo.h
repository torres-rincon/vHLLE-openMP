#include <deque>
#include "cll.h"

class Cell;
class Fluid;
class EoS;
class TransportCoeff;
class Particle;

struct NSThermo {
 double p;
 double T;
 double mub;
 double s;
 double etaS;
};

struct NSCenterState {
 double e;
 double p;
 double nb;
 double nq;
 double ns;
 double vx;
 double vy;
 double vz;
};

// this class implements the hydrodynamic evolution and
// contains the hydrodynamic algorithm
class Hydro {
private:
 Fluid *f;
 EoS *eos;
 TransportCoeff *trcoeff;
 double dt, tau;  // dt: timestep, tau: current value of the proper time
 double t;  // time in Cartesian frame
 double tau_z;    // effective value of the proper time used in 1/tau factors in
                  // the fluxes. Used to increase the accuracy
 bool vorticityOn = false;
 bool cartesian;

public:
 Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0, double _dt, bool _cartesian);
 ~Hydro();
 void enableVorticity(); // enable vorticity
 void setDtau(double deltaTau);  // change the timestep
 double getDtau() { return dt; }  // return current value of timestep
 void setFluid(Fluid *_f) { f = _f; }
 Fluid* getFluid()  const { return f; }

 // HLLE (ideal)flux between two neighbouring cells in a given direction
 // mode: PREDICT = used in predictor step; calculates fluxes for dt/2
 // CORRECT = used in corrector step, calculates fluxes based on predicted
 // half-step quantities
 bool hlle_flux_value(Cell *left, Cell *right, int direction, int mode, double flux[7]);
 void hlle_flux(Cell *left, Cell *right, int direction, int mode);
 void accumulate_hlle_fluxes(int ix, int iy, int iz, int mode);
 void accumulate_hlle_fluxes_buffered(int mode);
 // viscous flux \delta F
 bool visc_flux_value(Cell *left, Cell *right, int direction, double flux[4]);
 void visc_flux(Cell *left, Cell *right, int direction);
 void accumulate_visc_fluxes(int ix, int iy, int iz);
 void accumulate_visc_fluxes_buffered();
 // viscous source step for a given cell (ix,iy,iz)
 void visc_source_step(int ix, int iy, int iz);
 void source(double tau, double x, double y, double z, double Q[7],
             double S[7]);
 // ideal source step for a given cell (ix,iy,iz)
 void source_step(int ix, int iy, int iz, int mode);
 // shear stress tensor and bulk pressure in Navier-Stokes (NS) limit
 // plus \partial_\mu u^\nu matrix (dmu) and
 // expansion scalar \partial_mu u^\mu (du)
 // for a given cell (ix,iy,iz)
 bool NSquant(int ix, int iy, int iz, double pi[4][4], double &Pi,
              double dmu[4][4], std::unique_ptr<Matrix2D> &dbeta, double &du,
              const NSCenterState &centerState, NSThermo &thermo);
 // sets the values of shear stress/bulk pressure in NS limit in all hydro grid
 void setNSvalues();
 // advances numerical solution for shear/bulk in a whole grid over one
 // timestep
 void ISformal();
 // advances numerical solution for Q (including ideal and viscous fluxes and
 // source terms) over one timestep
 void performStep(void);
 // gets the current proper time
 inline double getTau() const { return tau; }
 // gets the current time in cartesian coordinates
 inline double getTime() { return t; }
 // adds sources from incoming particles into the hydro
 void addParticles(std::deque<Particle>* particles);
};
