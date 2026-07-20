/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <unistd.h>
#include <memory>
#include "hdo.h"
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "eos.h"
#include "trancoeff.h"
#include "particle.h"

using namespace std;

namespace {
double hydroElapsedSeconds(std::chrono::steady_clock::time_point start) {
 return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
}

// extern bool debugRiemann ;

namespace {
constexpr int PI_INDEX_44[4][4] = {{0, 1, 3, 6},
                                   {1, 2, 4, 7},
                                   {3, 4, 5, 8},
                                   {6, 7, 8, 9}};
}

double sign(double x) {
 if (x > 0)
  return 1.;
 else if (x < 0.)
  return -1.;
 else
  return 0.;
}

// this version contains NO PRE-ADVECTION for the IS solution

// enable this to use formal solution for the relaxation part of
// Israel-Stewart equations (not recommended)
//#define FORMAL_SOLUTION
// else: use 1st order finite difference update

Hydro::Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0,
             double _dt, bool _cartesian) {
 eos = _eos;
 trcoeff = _trcoeff;
 f = _f;
 dt = _dt;
 cartesian = _cartesian;
 if (cartesian) {
  t = _t0;
  tau = 1.0;
 }
 else {
  tau = _t0;
 }
}

Hydro::~Hydro() {}

void Hydro::enableVorticity() {
  //enable vorticity in hydro and fluid
  vorticityOn = true;
  f->enableVorticity();
}

void Hydro::setDtau(double deltaTau) {
 dt = deltaTau;
 if (dt > f->getDx() / 2. ||
     dt > f->getDy() / 2. /*|| dt>tau*f->getDz()/2. */) {
  cout << "too big delta_tau " << dt << "  " << f->getDx() << "  " << f->getDy()
       << "  " << tau * f->getDz() << endl;
  exit(1);
 }
}

bool Hydro::hlle_flux_value(Cell *left, Cell *right, int direction, int mode,
                            double flux[7]) {
 // for all variables, suffix "l" = left state, "r" = right state
 // with respect to the cell boundary
 double el, er, pl, pr, nbl, nql, nsl, nbr, nqr, nsr, vxl, vxr, vyl, vyr, vzl,
     vzr, bl = 0., br = 0., csb, vb, El, Er, dx = 0.;
 double Ftl = 0., Fxl = 0., Fyl = 0., Fzl = 0., Fbl = 0., Fql = 0., Fsl = 0.,
        Ftr = 0., Fxr = 0., Fyr = 0., Fzr = 0., Fbr = 0., Fqr = 0., Fsr = 0.;
 double U1l, U2l, U3l, U4l, Ubl, Uql, Usl, U1r, U2r, U3r, U4r, Ubr, Uqr, Usr;
 const double dta = mode == 0 ? dt / 2. : dt;
 double tauFactor = 1.0;  // fluxes are also multiplied by tau, 1 for Cartesian
 if (left->getM(direction) < 1. && right->getM(direction) < 1.) return false;
 if (mode == PREDICT) {
  // get primitive quantities from Q_{i+} at previous timestep
  left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                        direction);
  // ... and Q_{(i+1)-}
  right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                        direction);
  El = (el + pl) / (1 - vxl * vxl - vyl * vyl - vzl * vzl);
  Er = (er + pr) / (1 - vxr * vxr - vyr * vyr - vzr * vzr);
  if (!cartesian) {
   tauFactor = tau + 0.25 * dt;
  }
 } else {
  // use half-step updated Q's for corrector step
  left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                         direction);
  right->getPrimVarHLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                         direction);
  El = (el + pl) / (1 - vxl * vxl - vyl * vyl - vzl * vzl);
  Er = (er + pr) / (1 - vxr * vxr - vyr * vyr - vzr * vzr);
  if (!cartesian) {
   tauFactor = tau + 0.5 * dt;
  }
 }

 if (el < 0.) {
  el = 0.;
  pl = 0.;
 }
 if (er < 0.) {
  er = 0.;
  pr = 0.;
 }

 if (el > 1e10) {
  cout << "e>1e10; debug info below:\n";
  left->Dump(tau);
  // debugRiemann = true ;
  if (mode == PREDICT)
   left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                         direction);
  else
   left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                          direction);
  // debugRiemann = false ;
  exit(0);
 }

 // skip the procedure for two empty cells
 if (el == 0. && er == 0.) return false;
 if (pr < 0.) {
  cout << "Negative pressure" << endl;
  left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                        direction);
  right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                        direction);
 }

 double gammal = 1. / sqrt(1 - vxl * vxl - vyl * vyl - vzl * vzl);
 double gammar = 1. / sqrt(1 - vxr * vxr - vyr * vyr - vzr * vzr);
 U1l = gammal * gammal * (el + pl) * vxl;
 U2l = gammal * gammal * (el + pl) * vyl;
 U3l = gammal * gammal * (el + pl) * vzl;
 U4l = gammal * gammal * (el + pl) - pl;
 Ubl = gammal * nbl;
 Uql = gammal * nql;
 Usl = gammal * nsl;

 U1r = gammar * gammar * (er + pr) * vxr;
 U2r = gammar * gammar * (er + pr) * vyr;
 U3r = gammar * gammar * (er + pr) * vzr;
 U4r = gammar * gammar * (er + pr) - pr;
 Ubr = gammar * nbr;
 Uqr = gammar * nqr;
 Usr = gammar * nsr;

 if (direction == X_) {
  Ftl = U4l * vxl + pl * vxl;
  Fxl = U1l * vxl + pl;
  Fyl = U2l * vxl;
  Fzl = U3l * vxl;
  Fbl = Ubl * vxl;
  Fql = Uql * vxl;
  Fsl = Usl * vxl;

  Ftr = U4r * vxr + pr * vxr;
  Fxr = U1r * vxr + pr;
  Fyr = U2r * vxr;
  Fzr = U3r * vxr;
  Fbr = Ubr * vxr;
  Fqr = Uqr * vxr;
  Fsr = Usr * vxr;

  // for the case of constant c_s only
  csb = sqrt(eos->cs2() +
             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vxl - vxr, 2));
  vb = (sqrt(El) * vxl + sqrt(Er) * vxr) / (sqrt(El) + sqrt(Er));
  bl = min(0., min((vb - csb) / (1 - vb * csb),
                   (vxl - eos->cs()) / (1 - vxl * eos->cs())));
  br = max(0., max((vb + csb) / (1 + vb * csb),
                   (vxr + eos->cs()) / (1 + vxr * eos->cs())));

  dx = f->getDx();

  // bl or br in the case of boundary with vacuum
  if (el == 0.) bl = -1.;
  if (er == 0.) br = 1.;
 }
 if (direction == Y_) {
  Ftl = U4l * vyl + pl * vyl;
  Fxl = U1l * vyl;
  Fyl = U2l * vyl + pl;
  Fzl = U3l * vyl;
  Fbl = Ubl * vyl;
  Fql = Uql * vyl;
  Fsl = Usl * vyl;

  Ftr = U4r * vyr + pr * vyr;
  Fxr = U1r * vyr;
  Fyr = U2r * vyr + pr;
  Fzr = U3r * vyr;
  Fbr = Ubr * vyr;
  Fqr = Uqr * vyr;
  Fsr = Usr * vyr;

  // for the case of constant c_s only
  csb = sqrt(eos->cs2() +
             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vyl - vyr, 2));
  vb = (sqrt(El) * vyl + sqrt(Er) * vyr) / (sqrt(El) + sqrt(Er));
  bl = min(0., min((vb - csb) / (1 - vb * csb),
                   (vyl - eos->cs()) / (1 - vyl * eos->cs())));
  br = max(0., max((vb + csb) / (1 + vb * csb),
                   (vyr + eos->cs()) / (1 + vyr * eos->cs())));

  dx = f->getDy();

  // bl or br in the case of boundary with vacuum
  if (el == 0.) bl = -1.;
  if (er == 0.) br = 1.;
 }
 if (direction == Z_) {
  double tau1 = tauFactor;
  Ftl = U4l * vzl / tau1 + pl * vzl / tau1;
  Fxl = U1l * vzl / tau1;
  Fyl = U2l * vzl / tau1;
  Fzl = U3l * vzl / tau1 + pl / tau1;
  Fbl = Ubl * vzl / tau1;
  Fql = Uql * vzl / tau1;
  Fsl = Usl * vzl / tau1;

  Ftr = U4r * vzr / tau1 + pr * vzr / tau1;
  Fxr = U1r * vzr / tau1;
  Fyr = U2r * vzr / tau1;
  Fzr = U3r * vzr / tau1 + pr / tau1;
  Fbr = Ubr * vzr / tau1;
  Fqr = Uqr * vzr / tau1;
  Fsr = Usr * vzr / tau1;

  // for the case of constant c_s only
  // factor 1/tau accounts for eta-coordinate

  // different estimate
  csb = sqrt(eos->cs2() +
             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vzl - vzr, 2));
  vb = (sqrt(El) * vzl + sqrt(Er) * vzr) / (sqrt(El) + sqrt(Er));
  bl = 1. / tau * min(0., min((vb - csb) / (1 - vb * csb),
                              (vzl - eos->cs()) / (1 - vzl * eos->cs())));
  br = 1. / tau * max(0., max((vb + csb) / (1 + vb * csb),
                              (vzr + eos->cs()) / (1 + vzr * eos->cs())));

  dx = f->getDz();

  // bl or br in the case of boundary with vacuum
  if (el == 0.) bl = -1. / tau;
  if (er == 0.) br = 1. / tau;
 }

 if (bl == 0. && br == 0.) return false;

 // finally, HLLE formula for the fluxes
 flux[T_] = tauFactor * dta / dx *
            (-bl * br * (U4l - U4r) + br * Ftl - bl * Ftr) / (-bl + br);
 flux[X_] = tauFactor * dta / dx *
            (-bl * br * (U1l - U1r) + br * Fxl - bl * Fxr) / (-bl + br);
 flux[Y_] = tauFactor * dta / dx *
            (-bl * br * (U2l - U2r) + br * Fyl - bl * Fyr) / (-bl + br);
 flux[Z_] = tauFactor * dta / dx *
            (-bl * br * (U3l - U3r) + br * Fzl - bl * Fzr) / (-bl + br);
 flux[NB_] = tauFactor * dta / dx *
             (-bl * br * (Ubl - Ubr) + br * Fbl - bl * Fbr) / (-bl + br);
 flux[NQ_] = tauFactor * dta / dx *
             (-bl * br * (Uql - Uqr) + br * Fql - bl * Fqr) / (-bl + br);
 flux[NS_] = tauFactor * dta / dx *
             (-bl * br * (Usl - Usr) + br * Fsl - bl * Fsr) / (-bl + br);

 if (flux[NB_] != flux[NB_]) {  // if things failed
  cout << "---- error in hlle_flux: f_nb undefined!\n";
  cout << setw(12) << U4l << setw(12) << U1l << setw(12) << U2l << setw(12)
       << U3l << endl;
  cout << setw(12) << U4r << setw(12) << U1r << setw(12) << U2r << setw(12)
       << U3r << endl;
  cout << setw(12) << Ubl << setw(12) << Uql << setw(12) << Usl << endl;
  cout << setw(12) << Ubr << setw(12) << Uqr << setw(12) << Usr << endl;
  cout << setw(12) << Ftl << setw(12) << Fxl << setw(12) << Fyl << setw(12)
       << Fzl << endl;
  cout << setw(12) << Ftr << setw(12) << Fxr << setw(12) << Fyr << setw(12)
       << Fzr << endl;
  exit(1);
 }

 return true;
}

void Hydro::hlle_flux(Cell *left, Cell *right, int direction, int mode) {
 double flux[7];
 if (!hlle_flux_value(left, right, direction, mode, flux)) return;

 // update the cumulative fluxes in both neighbouring cells
 left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
               -flux[NQ_], -flux[NS_]);
 right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
                flux[NS_]);
}

void Hydro::accumulate_hlle_fluxes(int ix, int iy, int iz, int mode) {
 Cell *cell = f->getCell(ix, iy, iz);
 double flux[7];

 if (ix > 0 && hlle_flux_value(f->getCell(ix - 1, iy, iz), cell, X_, mode, flux)) {
  cell->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
                flux[NS_]);
 }
 if (ix < f->getNX() - 1 &&
     hlle_flux_value(cell, f->getCell(ix + 1, iy, iz), X_, mode, flux)) {
  cell->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
                -flux[NQ_], -flux[NS_]);
 }
 if (iy > 0 && hlle_flux_value(f->getCell(ix, iy - 1, iz), cell, Y_, mode, flux)) {
  cell->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
                flux[NS_]);
 }
 if (iy < f->getNY() - 1 &&
     hlle_flux_value(cell, f->getCell(ix, iy + 1, iz), Y_, mode, flux)) {
  cell->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
                -flux[NQ_], -flux[NS_]);
 }
 if (iz > 0 && hlle_flux_value(f->getCell(ix, iy, iz - 1), cell, Z_, mode, flux)) {
  cell->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
                flux[NS_]);
 }
 if (iz < f->getNZ() - 1 &&
     hlle_flux_value(cell, f->getCell(ix, iy, iz + 1), Z_, mode, flux)) {
  cell->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
                -flux[NQ_], -flux[NS_]);
 }
}

void Hydro::accumulate_hlle_fluxes_buffered(int mode) {
 const int nx = f->getNX();
 const int ny = f->getNY();
 const int nz = f->getNZ();

#ifdef _OPENMP
#pragma omp parallel
#endif
 {
  {
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
   for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
     for (int ix = 0; ix < nx - 1; ix++) {
      double flux[7];
      Cell *left = f->getCell(ix, iy, iz);
      Cell *right = f->getCell(ix + 1, iy, iz);
      if (hlle_flux_value(left, right, X_, mode, flux)) {
       left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
                     -flux[NQ_], -flux[NS_]);
       right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_],
                      flux[NQ_], flux[NS_]);
      }
     }
    }
  }
  {
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
   for (int iz = 0; iz < nz; iz++)
    for (int ix = 0; ix < nx; ix++) {
     for (int iy = 0; iy < ny - 1; iy++) {
      double flux[7];
      Cell *left = f->getCell(ix, iy, iz);
      Cell *right = f->getCell(ix, iy + 1, iz);
      if (hlle_flux_value(left, right, Y_, mode, flux)) {
       left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
                     -flux[NQ_], -flux[NS_]);
       right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_],
                      flux[NQ_], flux[NS_]);
      }
     }
    }
  }
  {
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
   for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
     for (int iz = 0; iz < nz - 1; iz++) {
      double flux[7];
      Cell *left = f->getCell(ix, iy, iz);
      Cell *right = f->getCell(ix, iy, iz + 1);
      if (hlle_flux_value(left, right, Z_, mode, flux)) {
       left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
                     -flux[NQ_], -flux[NS_]);
       right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_],
                      flux[NQ_], flux[NS_]);
      }
     }
    }
  }
 }
}

void Hydro::source(double tau1, double x, double y, double z, double Q[7],
                   double S[7]) {
 if (cartesian) {
  // geometrical source term is zero in Cartesian frame
  for (int i = 0; i < 7; i++) S[i] = 0.0;
 }
 else {
  double _Q[7], e, p, nb, nq, ns, vx, vy, vz;
  for (int i = 0; i < 7; i++) _Q[i] = Q[i] / tau1;  // no tau factor in  _Q
  transformPV(eos, _Q, e, p, nb, nq, ns, vx, vy, vz);
  S[T_] = -_Q[T_] * vz * vz - p * (1. + vz * vz);
  S[X_] = 0.;
  S[Y_] = 0.;
  S[Z_] = -_Q[Z_];
  S[NB_] = 0.;
  S[NQ_] = 0.;
  S[NS_] = 0.;
 }
}

void Hydro::source_step(int ix, int iy, int iz, int mode) {
 if (cartesian) return;

 double _dt;
 if (mode == PREDICT)
  _dt = dt / 2.;
 else
  _dt = dt;

 double tau1;
 double Q[7];
 double k[7];

 Cell *c = f->getCell(ix, iy, iz);

 if (mode == PREDICT) {
  c->getQ(Q);
  tau1 = (cartesian) ? t : tau;
 } else {
  c->getQh(Q);
  tau1 = (cartesian) ? t : tau;
  tau1 = tau1 + 0.5 * dt;
 }
 if (Q[T_] == 0.0 && Q[Z_] == 0.0) return;

 source(tau1, 0.0, 0.0, 0.0, Q, k);
 for (int i = 0; i < 7; i++) k[i] *= _dt;

 if (k[NB_] != k[NB_]) {  // something failed
  cout << "---- error in source_step: k_nb undefined!\n";
  cout << setw(12) << k[0] << setw(12) << k[1] << setw(12) << k[2] << setw(12)
       << k[3] << endl;
  cout << setw(12) << k[4] << setw(12) << k[5] << setw(12) << k[6] << endl;
  exit(1);
 }
 c->addFlux(k[T_], k[X_], k[Y_], k[Z_], k[NB_], k[NQ_], k[NS_]);
}

void Hydro::visc_source_step(int ix, int iy, int iz) {
 if (cartesian) return;

 double e, p, nb, nq, ns, vx, vy, vz;
 double uuu[4];
 double k[7];

 Cell *c = f->getCell(ix, iy, iz);
 double Qh[7];
 c->getQh(Qh);
 if (Qh[T_] <= 0.0) return;

 c->getPrimVarHCenter(eos, tau - dt / 2., e, p, nb, nq, ns, vx, vy, vz);  // TODO Cartesian
 if (e <= 0.) return;
 uuu[0] = 1. / sqrt(1. - vx * vx - vy * vy - vz * vz);
 uuu[1] = uuu[0] * vx;
 uuu[2] = uuu[0] * vy;
 uuu[3] = uuu[0] * vz;

 k[T_] = -c->getpiH(3, 3) + c->getPiH() * (-1.0 - uuu[3] * uuu[3]);
 k[X_] = 0.;
 k[Y_] = 0.;
 k[Z_] = -(c->getpiH(0, 3) + c->getPiH() * uuu[0] * uuu[3]);
 for (int i = 0; i < 4; i++) k[i] *= dt;
 c->addFlux(k[T_], k[X_], k[Y_], k[Z_], 0., 0., 0.);
}

// for the procedure below, the following approximations are used:
// dv/d(tau) = v^{t+dt}_ideal - v^{t}
// dv/dx_i ~ v^{x+dx}-v{x-dx},
// which makes sense after non-viscous step
bool Hydro::NSquant(int ix, int iy, int iz, double pi[4][4], double &Pi,
                    double dmu[4][4], unique_ptr<Matrix2D> &dbeta, double &du,
                    const NSCenterState &centerState, NSThermo &thermo) {
 const double VMIN = 1e-2;
 const double UDIFF = 3.0;
 double e0, e1, p, nb, nq, ns, vx1, vy1, vz1, vx0, vy0, vz0, vxH, vyH, vzH;
 double T, T0, T1, mub, muq, mus;
 double ut0, ux0, uy0, uz0, ut1, ux1, uy1, uz1;
 //	double dmu [4][4] ; // \partial_\mu u^\nu matrix
 // coordinates: 0=tau, 1=x, 2=y, 3=eta
 double uuu[4];         // the 4-velocity
 double gmunu[4][4] = {{1, 0, 0, 0},
                       {0, -1, 0, 0},
                       {0, 0, -1, 0},
                       {0, 0, 0, -1}};  // omit 1/tau^2 in g^{eta,eta}
 Cell *c = f->getCell(ix, iy, iz);
 double dx = f->getDx(), dy = f->getDy(), dz = f->getDz();
 // check if the cell is next to vacuum from +-x, +-y side:
 if (f->getCell(ix + 1, iy, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy + 1, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix - 1, iy, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy - 1, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy, iz + 1)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy, iz - 1)->getMaxM() <= 0.9 ||

     f->getCell(ix + 1, iy + 1, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix + 1, iy - 1, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix - 1, iy + 1, iz)->getMaxM() <= 0.9 ||
     f->getCell(ix - 1, iy - 1, iz)->getMaxM() <= 0.9 ||
     
     f->getCell(ix + 1, iy, iz + 1)->getMaxM() <= 0.9 ||
     f->getCell(ix + 1, iy, iz - 1)->getMaxM() <= 0.9 ||
     f->getCell(ix - 1, iy, iz + 1)->getMaxM() <= 0.9 ||
     f->getCell(ix - 1, iy, iz - 1)->getMaxM() <= 0.9 ||
     
     f->getCell(ix, iy + 1, iz + 1)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy + 1, iz - 1)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy - 1, iz + 1)->getMaxM() <= 0.9 ||
     f->getCell(ix, iy - 1, iz - 1)->getMaxM() <= 0.9) {
  
  for (int i = 0; i < 4; i++)
   for (int j = 0; j < 4; j++) {
    pi[i][j] = 0.;
    dmu[i][j] = 0.;
  }
  Pi = du = 0.;
  return false;
 }
 // calculation of \partial_\mu u^\nu matrix
 // mu=first index, nu=second index
 // centered differences with respect to the values at (it+1/2, ix, iy, iz)
 // d_tau u^\mu
 if (cartesian) {
  c->getPrimVarPrev(eos, 1.0, e0, p, nb, nq, ns, vx0, vy0, vz0);
  c->getPrimVar(eos, 1.0, e1, p, nb, nq, ns, vx1, vy1, vz1);
 }
 else {
  c->getPrimVarPrev(eos, tau - dt, e0, p, nb, nq, ns, vx0, vy0, vz0);
  c->getPrimVar(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
 }
 e1 = centerState.e;
 p = centerState.p;
 nb = centerState.nb;
 nq = centerState.nq;
 ns = centerState.ns;
 vxH = centerState.vx;
 vyH = centerState.vy;
 vzH = centerState.vz;
 double tauPlusHalf = (cartesian) ? 1.0 : (tau + 0.5 * dt);
 double tauMinusHalf = (cartesian) ? 1.0 : (tau - 0.5 * dt);
 if (vorticityOn) {
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
 }
 //############## get transport coefficients
 double etaS, zetaS;
 double s = eos->s(e1, nb, nq, ns);  // entropy density in the current cell
 eos->eos(e1, nb, nq, ns, T, mub, muq, mus, p);
 trcoeff->getEta(e1, nb, T, mub,  s,  p, etaS, zetaS);
 thermo.p = p;
 thermo.T = T;
 thermo.mub = mub;
 thermo.s = s;
 thermo.etaS = etaS;
 //##############
 // if(e1<0.00004) s=0. ; // negative pressure due to pi^zz for small e
 ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
 ux0 = ut0 * vx0;
 uy0 = ut0 * vy0;
 uz0 = ut0 * vz0;
 ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
 ux1 = ut1 * vx1;
 uy1 = ut1 * vy1;
 uz1 = ut1 * vz1;
 uuu[0] = 1.0 / sqrt(1.0 - vxH * vxH - vyH * vyH - vzH * vzH);
 uuu[1] = uuu[0] * vxH;
 uuu[2] = uuu[0] * vyH;
 uuu[3] = uuu[0] * vzH;

 dmu[0][0] = (ut1 * ut1 - ut0 * ut0) / 2. / uuu[0] / dt;
 dmu[0][1] = (ux1 * ux1 - ux0 * ux0) / 2. / uuu[1] / dt;
 dmu[0][2] = (uy1 * uy1 - uy0 * uy0) / 2. / uuu[2] / dt;
 dmu[0][3] = (uz1 * uz1 - uz0 * uz0) / 2. / uuu[3] / dt;
 if (vorticityOn) {
  (*dbeta)[0][0] = (ut1 / T1 - ut0 / T0) / dt;
  (*dbeta)[0][1] = (ux1 / T1 - ux0 / T0) / dt;
  (*dbeta)[0][2] = (uy1 / T1 - uy0 / T0) / dt;
  (*dbeta)[0][3] = (uz1 / T1 - uz0 / T0) / dt;
  if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.) {
    (*dbeta)[0][0] = (*dbeta)[0][1] = (*dbeta)[0][2] = (*dbeta)[0][3] = 0.;
  }
 }
 if (fabs(0.5 * (ut1 + ut0) / ut1) > UDIFF) dmu[0][0] = (ut1 - ut0) / dt;
 if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / ux1) > UDIFF)
  dmu[0][1] = (ux1 - ux0) / dt;
 if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uy1) > UDIFF)
  dmu[0][2] = (uy1 - uy0) / dt;
 if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uz1) > UDIFF)
  dmu[0][3] = (uz1 - uz0) / dt;
 if (e1 <= 0. || e0 <= 0.) {  // matter-vacuum
  dmu[0][0] = dmu[0][1] = dmu[0][2] = dmu[0][3] = 0.;
 }
 // d_x u^\mu
 f->getCell(ix + 1, iy, iz)
     ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
 if (vorticityOn) {
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
 }
 f->getCell(ix - 1, iy, iz)
     ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
 if (vorticityOn) {
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
 }
 if (e1 > 0. && e0 > 0.) {
  ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
  ux0 = ut0 * vx0;
  uy0 = ut0 * vy0;
  uz0 = ut0 * vz0;
  ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
  ux1 = ut1 * vx1;
  uy1 = ut1 * vy1;
  uz1 = ut1 * vz1;
  dmu[1][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dx;
  dmu[1][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dx;
  dmu[1][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dx;
  dmu[1][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dx;
  if(vorticityOn) {
    (*dbeta)[1][0] = 0.5 * (ut1 / T1 - ut0 / T0) / dx;
    (*dbeta)[1][1] = 0.5 * (ux1 / T1 - ux0 / T0) / dx;
    (*dbeta)[1][2] = 0.5 * (uy1 / T1 - uy0 / T0) / dx;
    (*dbeta)[1][3] = 0.5 * (uz1 / T1 - uz0 / T0) / dx;
    if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.) {
      (*dbeta)[1][0] = (*dbeta)[1][1] = (*dbeta)[1][2] = (*dbeta)[1][3] = 0.;
    }
  }
  if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
   dmu[1][0] = 0.5 * (ut1 - ut0) / dx;
  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
   dmu[1][1] = 0.5 * (ux1 - ux0) / dx;
  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
   dmu[1][2] = 0.5 * (uy1 - uy0) / dx;
  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
   dmu[1][3] = 0.5 * (uz1 - uz0) / dx;
 } else {  // matter-vacuum
  dmu[1][0] = dmu[1][1] = dmu[1][2] = dmu[1][3] = 0.;
 }
 if (fabs(dmu[1][3]) > 1e+10)
  cout << "dmu[1][3]:  " << uz1 << "  " << uz0 << "  " << uuu[3] << endl;
 // d_y u^\mu
 f->getCell(ix, iy + 1, iz)
     ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
 if (vorticityOn) {
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
 }
 f->getCell(ix, iy - 1, iz)
     ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
 if (vorticityOn) {
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
 }
 if (e1 > 0. && e0 > 0.) {
  ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
  ux0 = ut0 * vx0;
  uy0 = ut0 * vy0;
  uz0 = ut0 * vz0;
  ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
  ux1 = ut1 * vx1;
  uy1 = ut1 * vy1;
  uz1 = ut1 * vz1;
  dmu[2][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dy;
  dmu[2][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dy;
  dmu[2][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dy;
  dmu[2][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dy;
  if (vorticityOn) {
    (*dbeta)[2][0] = 0.5 * (ut1 / T1 - ut0 / T0) / dy;
    (*dbeta)[2][1] = 0.5 * (ux1 / T1 - ux0 / T0) / dy;
    (*dbeta)[2][2] = 0.5 * (uy1 / T1 - uy0 / T0) / dy;
    (*dbeta)[2][3] = 0.5 * (uz1 / T1 - uz0 / T0) / dy;
    if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.) {
      (*dbeta)[2][0] = (*dbeta)[2][1] = (*dbeta)[2][2] = (*dbeta)[2][3] = 0.;
    }
  }
  if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
   dmu[2][0] = 0.5 * (ut1 - ut0) / dy;
  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
   dmu[2][1] = 0.5 * (ux1 - ux0) / dy;
  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
   dmu[2][2] = 0.5 * (uy1 - uy0) / dy;
  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
   dmu[2][3] = 0.5 * (uz1 - uz0) / dy;
 } else {  // matter-vacuum
  dmu[2][0] = dmu[2][1] = dmu[2][2] = dmu[2][3] = 0.;
 }
 // d_z u^\mu
 f->getCell(ix, iy, iz + 1)
     ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
 if (vorticityOn) {
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
 }
 f->getCell(ix, iy, iz - 1)
     ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
 if (vorticityOn) {
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
 }
 if (e1 > 0. && e0 > 0.) {
  ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
  ux0 = ut0 * vx0;
  uy0 = ut0 * vy0;
  uz0 = ut0 * vz0;
  ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
  ux1 = ut1 * vx1;
  uy1 = ut1 * vy1;
  uz1 = ut1 * vz1;
  dmu[3][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dz / tauPlusHalf;
  dmu[3][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dz / tauPlusHalf;
  dmu[3][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dz / tauPlusHalf;
  dmu[3][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dz / tauPlusHalf;
  if (vorticityOn) {
    (*dbeta)[3][0] = 0.5 * (ut1 / T1 - ut0 / T0) / (dz * tauPlusHalf);
    (*dbeta)[3][1] = 0.5 * (ux1 / T1 - ux0 / T0) / (dz * tauPlusHalf);
    (*dbeta)[3][2] = 0.5 * (uy1 / T1 - uy0 / T0) / (dz * tauPlusHalf);
    (*dbeta)[3][3] = 0.5 * (uz1 / T1 - uz0 / T0) / (dz * tauPlusHalf);
    if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.){
      (*dbeta)[3][0] = (*dbeta)[3][1] = (*dbeta)[3][2] = (*dbeta)[3][3] = 0.;
 }
  }
  if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
   dmu[3][0] = 0.5 * (ut1 - ut0) / dz / tauPlusHalf;
  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
   dmu[3][1] = 0.5 * (ux1 - ux0) / dz / tauPlusHalf;
  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
   dmu[3][2] = 0.5 * (uy1 - uy0) / dz / tauPlusHalf;
  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
   dmu[3][3] = 0.5 * (uz1 - uz0) / dz / tauPlusHalf;
 } else {  // matter-vacuum
  dmu[3][0] = dmu[3][1] = dmu[3][2] = dmu[3][3] = 0.;
 }
 // additional terms from Christoffel symbols :)
 if (!cartesian) {
  dmu[3][0] += uuu[3] / tauMinusHalf;
  dmu[3][3] += uuu[0] / tauMinusHalf;
  if(vorticityOn && T>0.){
   (*dbeta)[3][0] += uuu[3] / (T * tauMinusHalf);
   (*dbeta)[3][3] += uuu[0] / (T * tauMinusHalf);
  }
 }
 // calculating sigma[mu][nu]
 for (int i = 0; i < 4; i++)
  for (int j = 0; j < 4; j++) {
   pi[i][j] = 0.0;
   for (int k = 0; k < 4; k++)
    for (int l = 0; l < 4; l++) {
     double Z = 0.0;
     if (j == l)
      Z += 0.5 * (gmunu[i][k] - uuu[i] * uuu[k]);
     if (i == l)
      Z += 0.5 * (gmunu[j][k] - uuu[j] * uuu[k]);
     if (k == l)
      Z -= (gmunu[i][j] - uuu[i] * uuu[j]) / 3.0;
     pi[i][j] += Z * dmu[k][l] * 2.0 * etaS * s / 5.068;
    }
  }
 Pi = -zetaS * s * (dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3]) /
      5.068;  // fm^{-4} --> GeV/fm^3
 du = dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3];
 //--------- debug part: NaN/inf check, trace check, diag check, transversality
 // check
 for (int i = 0; i < 4; i++)
  for (int j = 0; j < 4; j++) {
   if (pi[i][j] != 0. && fabs(1.0 - pi[j][i] / pi[i][j]) > 1e-10)
    cout << "non-diag: " << pi[i][j] << "  " << pi[j][i] << endl;
   if (std::isinf(pi[i][j]) || std::isnan(pi[i][j])) {
    cout << "hydro:NSquant: inf/nan i " << i << " j " << j << endl;
    exit(1);
   }
  }
 return true;
}

void Hydro::setNSvalues() {
 double e, p, nb, nq, ns, vx, vy, vz, piNS[4][4], PiNS;
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
    if (e <= 0.) continue;
    // NSquant(ix, iy, iz, piNS, PiNS, dmu, du) ;
    //############## set NS values assuming initial zero flow + Bjorken z
    // flow
    double T, mub, muq, mus;
    double etaS, zetaS;
    double s = eos->s(e, nb, nq, ns);  // entropy density in the current cell
    eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
    trcoeff->getEta(e, nb, T, mub,  s,  p, etaS, zetaS);
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) piNS[i][j] = 0.0;  // reset piNS
    piNS[1][1] = piNS[2][2] = 2.0 / 3.0 * etaS * s / tau / 5.068;
    piNS[3][3] = -2.0 * piNS[1][1];
    PiNS = 0.0;
    for (int i = 0; i < 4; i++)
     for (int j = 0; j <= i; j++) c->setpi(i, j, piNS[i][j]);
    c->setPi(PiNS);
   }
 cout << "setNS done\n";
}

void Hydro::ISformal() {
 const bool timingEnabled = std::getenv("VHLLE_TIMING") != nullptr;
 auto timingStart = std::chrono::steady_clock::now();
 double timingRelaxSource = 0.0;
 double timingAdvection = 0.0;
 double timingNSquantCpuSum = 0.0;
 long long timingNSquantCalls = 0;
 long long timingNSquantThermoHits = 0;
 const double gmumu[4] = {1., -1., -1., -1.};
 double tauMinusHalf = (cartesian) ? 1.0 : (tau - 0.5 * dt);
 double tauMinusDt = (cartesian) ? 1.0 : (tau - dt);
 
#ifdef _OPENMP
#pragma omp parallel
#endif
 {
 double timingNSquantCpuSumLocal = 0.0;
 long long timingNSquantCallsLocal = 0;
 long long timingNSquantThermoHitsLocal = 0;
 if (timingEnabled) {
#ifdef _OPENMP
#pragma omp single
#endif
  {
   timingStart = std::chrono::steady_clock::now();
  }
 }
 // loop #1 (relaxation+source terms)
#ifdef _OPENMP
#pragma omp for collapse(3) schedule(static)
#endif
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    double e, p, nb, nq, ns, vx, vy, vz, T, mub, muq, mus;
    double piNS[4][4], sigNS[4][4], PiNS, dmu[4][4], du;
    std::unique_ptr<Matrix2D> dbeta = vorticityOn
       ? std::make_unique<Matrix2D>(Matrix2D(4, std::vector<double>(4, 0.0)))
       : nullptr;
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarPrev(eos, tauMinusDt, e, p, nb, nq, ns, vx, vy,
                         vz);  // instead of getPrimVar()
    if (e <= 0.) {             // empty cell?
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
       c->setpiH0(i, j, 0.0);
       c->setpi0(i, j, 0.0);
      }
     c->setPiH0(0.0);
     c->setPi0(0.0);
     if (vorticityOn) {
      c->resetDbeta();
     }
    } else {  // non-empty cell
     // 1) relaxation(pi)+source(pi) terms for half-step
     double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
     double u[4];
     u[0] = gamma;
     u[1] = u[0] * vx;
     u[2] = u[0] * vy;
     u[3] = u[0] * vz;
     // source term  + tau*delta_Q_i/delta_tau
     double flux[4];
     for (int i = 0; i < 4; i++)
      flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * u[0] * u[i]);
     flux[0] += -tauMinusDt * c->getPi();
     c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
     c->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vx, vy, vz);
     gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
     u[0] = gamma;
     u[1] = u[0] * vx;
     u[2] = u[0] * vy;
     u[3] = u[0] * vz;
	     // now calculating viscous terms in NS limit
	     NSThermo nsThermo;
	     const NSCenterState centerState = {e, p, nb, nq, ns, vx, vy, vz};
	     bool hasNSThermo;
	     if (timingEnabled) {
	      const auto timingNSquantStart = std::chrono::steady_clock::now();
	      hasNSThermo = NSquant(ix, iy, iz, piNS, PiNS, dmu, dbeta, du,
	                            centerState, nsThermo);
	      timingNSquantCpuSumLocal += hydroElapsedSeconds(timingNSquantStart);
	      timingNSquantCallsLocal++;
	      if (hasNSThermo) timingNSquantThermoHitsLocal++;
	     } else {
	      hasNSThermo = NSquant(ix, iy, iz, piNS, PiNS, dmu, dbeta, du,
	                            centerState, nsThermo);
	     }
	     if (vorticityOn) {
	      c->setDbeta(dbeta);
	     }
     double etaS, s;
     if (hasNSThermo) {
      p = nsThermo.p;
      T = nsThermo.T;
      mub = nsThermo.mub;
      s = nsThermo.s;
      etaS = nsThermo.etaS;
     } else {
      eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
      double zetaS;
      s = eos->s(e, nb, nq, ns);
      trcoeff->getEta(e, nb, T, mub,  s,  p, etaS, zetaS);
     }
     const double eta = etaS * s;
     // auxiliary variable sigmaNS = piNS / (2*eta), 
     // mainly to protect against division by zero in the eta=0 case.
     for(int i=0; i<4; i++)
     for(int j=0; j<4; j++) {
      sigNS[i][j] = 0.5 * piNS[i][j] / eta * 5.068;
      if(eta<=0.0) sigNS[i][j] = 0.0;
     }
     //############# get relaxation times
     double taupi, tauPi;
     trcoeff->getTau(e, nb, T, mub, s, p, taupi, tauPi);
     double deltapipi, taupipi, lambdapiPi, phi7, delPiPi, lamPipi;
     trcoeff->getOther(e, nb, nq, ns, deltapipi, taupipi, lambdapiPi, phi7);
     phi7 = phi7/taupi;  // dividing by tau_pi here, to avoid NaNs when tau_pi==0
     if(taupi<0.5*dt)
      deltapipi = taupipi = lambdapiPi = phi7 = 0.0;
     trcoeff->getOtherBulk(e, nb, nq, ns, delPiPi, lamPipi);
     if(tauPi<0.5*dt)
      delPiPi = lamPipi = 0.0;
     //#############
	     double piCell[4][4];
	     for (int i = 0; i < 4; i++)
	      for (int j = 0; j < 4; j++)
	       piCell[i][j] = c->getpi(i, j);
	     const double PiCell = c->getPi();
	     double Delta[10];
	     // relaxation term, piH,PiH-->half-step
	     for (int i = 0; i < 4; i++)
	      for (int j = 0; j <= i; j++) {
	       const int ij = PI_INDEX_44[i][j];
	       Delta[ij] = -u[i] * u[j];
	       if (i == j) Delta[ij] += gmumu[i];
#ifdef FORMAL_SOLUTION
	       c->setpiH0(i, j, (piCell[i][j] - piNS[i][j]) *
	                                exp(-dt / 2.0 / gamma / taupi) +
	                            piNS[i][j]);
#else
	      if(taupi>0.5*dt)
	       c->setpiH0(i, j, piCell[i][j] -
	          (piCell[i][j] - piNS[i][j]) * dt / 2.0 / gamma / taupi);
	      else
	       c->setpiH0(i, j, piNS[i][j]);
#endif
	      }
#ifdef FORMAL_SOLUTION
	     c->setPiH0((PiCell - PiNS) * exp(-dt / 2.0 / gamma / tauPi) + PiNS);
#else
	    if(tauPi>0.5*dt)
	     c->setPiH0(PiCell - (PiCell - PiNS) * dt / 2.0 / gamma / tauPi);
	    else
	     c->setPiH0(PiNS);
#endif
	     // sources from Christoffel symbols from \dot pi_munu - only in tau-eta coordinate frame
	     if (!cartesian) {
	      double tau1 = tau - dt * 0.75;
	      c->addpiH0(0, 0,
	                -2. * vz * piCell[0][3] / tau1 * dt / 2.);  // *gamma/gamma
	      c->addpiH0(3, 3, -(2. * vz / tau1 * piCell[0][3]) * dt / 2.);
	      c->addpiH0(
	         3, 0,
	         -(vz / tau1 * piCell[0][0] + vz / tau1 * piCell[3][3]) * dt / 2.);
	      c->addpiH0(1, 0, -vz / tau1 * piCell[1][3] * dt / 2.);
	      c->addpiH0(2, 0, -vz / tau1 * piCell[2][3] * dt / 2.);
	      c->addpiH0(3, 1, -(vz / tau1 * piCell[0][1]) * dt / 2.);
	      c->addpiH0(3, 2, -(vz / tau1 * piCell[0][2]) * dt / 2.);
	     }
     
     // source from full IS equations (see  draft for the description)
	     for (int i = 0; i < 4; i++)
	      for (int j = 0; j <= i; j++) {
	       const int ij = PI_INDEX_44[i][j];
	       // now transversality and cross terms
	       c->addpiH0(i, j, (- deltapipi * piCell[i][j] * du +
	         lambdapiPi * PiCell * sigNS[i][j]) / gamma * 0.5 * dt);
	       for (int k = 0; k < 4; k++) {
	        // parts of terms with one internal summation index
	        c->addpiH0(i, j, (phi7 * piCell[i][k] * piCell[j][k] * gmumu[k] - taupipi * 0.5 * (piCell[i][k] * sigNS[j][k] * gmumu[k] + piCell[j][k] * sigNS[i][k] * gmumu[k])) / gamma * 0.5 * dt);
	        // parts of terms with two internal summation indexes
	        for (int l = 0; l < 4; l++){
	         c->addpiH0(i, j, (-piCell[i][k] * u[j] - piCell[j][k] * u[i]) * u[l] * dmu[l][k] * gmumu[k] / gamma * 0.5 * dt
	          - 1. / 3. * Delta[ij] * piCell[k][l] * ( phi7 * piCell[k][l] - taupipi * sigNS[k][l]) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt);
	        }
	       }
	      }
	     for(int k = 0; k < 4; k++)
	     for(int l = 0; l < 4; l++) {
		      c->addPiH0(lamPipi * piCell[k][l] * sigNS[k][l] * gmumu[k] * gmumu[l] / gamma * 0.5 * dt);
	     }
		     c->addPiH0(-delPiPi * PiCell * du / gamma * 0.5 * dt);
	     double piH0Cell[4][4];
	     for (int i = 0; i < 4; i++)
	      for (int j = 0; j < 4; j++)
	       piH0Cell[i][j] = c->getpiH0(i, j);
	     const double PiH0Cell = c->getPiH0();
	     // 1) relaxation(piH)+source(piH) terms for full-step
	     for (int i = 0; i < 4; i++)
	      for (int j = 0; j <= i; j++) {
#ifdef FORMAL_SOLUTION
	       c->setpi0(i, j,
	                 (piCell[i][j] - piNS[i][j]) * exp(-dt / gamma / taupi) +
	                     piNS[i][j]);
#else
	      if(taupi>0.5*dt)
	       c->setpi0(i, j, piCell[i][j] -
	          (piH0Cell[i][j] - piNS[i][j]) * dt / gamma / taupi);
	      else
	       c->setpi0(i, j, piNS[i][j]);
#endif
	      }
#ifdef FORMAL_SOLUTION
	     c->setPi0((PiCell - PiNS) * exp(-dt / gamma / tauPi) + PiNS);
#else
	    if(tauPi>0.5*dt)
	     c->setPi0(PiCell - (PiH0Cell - PiNS) * dt / gamma / tauPi);
	    else
	     c->setPi0(PiNS);
#endif
	     if (!cartesian) {
	      double tau1 = tau - dt * 0.5;
	      c->addpi0(0, 0, -2. * vz / tau1 * piH0Cell[0][3] * dt);  // *gamma/gamma
	      c->addpi0(3, 3, -(2. * vz / tau1 * piH0Cell[0][3]) * dt);
	      c->addpi0(
	         3, 0,
	         -(vz / tau1 * piH0Cell[0][0] + vz / tau1 * piH0Cell[3][3]) * dt);
	      c->addpi0(1, 0, -vz / tau1 * piH0Cell[1][3] * dt);
	      c->addpi0(2, 0, -vz / tau1 * piH0Cell[2][3] * dt);
	      c->addpi0(3, 1, -(vz / tau1 * piH0Cell[0][1]) * dt);
	      c->addpi0(3, 2, -(vz / tau1 * piH0Cell[0][2]) * dt);
	     }
	     
	     // source from full IS equations (see draft for the description)
	     for (int i = 0; i < 4; i++)
	      for (int j = 0; j <= i; j++) {
	       const int ij = PI_INDEX_44[i][j];
	       // now transversality and cross terms
	       c->addpi0(i, j, (- deltapipi * piH0Cell[i][j] * du +
	         lambdapiPi * PiH0Cell * sigNS[i][j]) / gamma * dt);
	       for (int k = 0; k < 4; k++) {
	        // parts of terms with one internal summation index
	        c->addpi0(i, j, (phi7 * piCell[i][k] * piH0Cell[j][k] * gmumu[k] - taupipi * 0.5 * (piH0Cell[i][k] * sigNS[j][k] * gmumu[k] + piH0Cell[j][k] * sigNS[i][k] * gmumu[k])) / gamma * dt);
	        for (int l = 0; l < 4; l++){
	         c->addpi0(i, j, ((-piH0Cell[i][k] * u[j] - piH0Cell[j][k] * u[i]) * u[l] * dmu[l][k] * gmumu[k]
	          - 1. / 3. * Delta[ij] * piH0Cell[k][l] * ( phi7 * piH0Cell[k][l] - taupipi * sigNS[k][l]) * gmumu[k] * gmumu[l]) / gamma * dt);
	        }
	       }
	      }
	     for (int k = 0; k < 4; k++)
	     for (int l = 0; l < 4; l++) {
		      c->addPi0(lamPipi * piH0Cell[k][l] * sigNS[k][l] * gmumu[k] * gmumu[l] / gamma * dt);
	     }
		     c->addPi0(-delPiPi * PiH0Cell * du / gamma * dt);
	    }  // end non-empty cell
	   }   // end loop #1
 if (timingEnabled) {
#ifdef _OPENMP
#pragma omp critical
#endif
  {
   timingNSquantCpuSum += timingNSquantCpuSumLocal;
   timingNSquantCalls += timingNSquantCallsLocal;
   timingNSquantThermoHits += timingNSquantThermoHitsLocal;
  }
 }
 if (timingEnabled) {
#ifdef _OPENMP
#pragma omp single
#endif
  {
   timingRelaxSource = hydroElapsedSeconds(timingStart);
   timingStart = std::chrono::steady_clock::now();
  }
 }
	
	 // 3) -- advection ---
#ifdef _OPENMP
#pragma omp for collapse(3) schedule(static)
#endif
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    double e, p, nb, nq, ns, vx, vy, vz;
    double pi[4][4], piH[4][4], Pi, PiH;
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vx, vy,
                         vz);  // getPrimVar() before
    if (e <= 0.) continue;
    const double dx = f->getDx();
    const double dy = f->getDy();
    const double dz = f->getDz();
    double xm = -vx * dt / dx;
    double ym = -vy * dt / dy;
    double zm = -vz * dt / dz / tauMinusHalf;
    double xmH = -vx * dt / dx / 2.0;
    double ymH = -vy * dt / dy / 2.0;
    double zmH = -vz * dt / dz / tauMinusHalf / 2.0;
    const int sx = (xm > 0.) - (xm < 0.);
    const int sy = (ym > 0.) - (ym < 0.);
    const int sz = (zm > 0.) - (zm < 0.);
    double wx[2] = {(1. - fabs(xm)), fabs(xm)};
    double wy[2] = {(1. - fabs(ym)), fabs(ym)};
    double wz[2] = {(1. - fabs(zm)), fabs(zm)};
    double wxH[2] = {(1. - fabs(xmH)), fabs(xmH)};
    double wyH[2] = {(1. - fabs(ymH)), fabs(ymH)};
    double wzH[2] = {(1. - fabs(zmH)), fabs(zmH)};
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) {
      pi[i][j] = piH[i][j] = 0.0;
     }
    Pi = PiH = 0.0;
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       // pi,Pi-->full step, piH,PiH-->half-step
       Cell *c1 = f->getCell(ix + jx * sx, iy + jy * sy, iz + jz * sz);
       const double w = wx[jx] * wy[jy] * wz[jz];
       const double wH = wxH[jx] * wyH[jy] * wzH[jz];
       const double *pi0 = c1->getpi0Array();
       const double *piH0 = c1->getpiH0Array();
       for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
         const int ij = PI_INDEX_44[i][j];
         pi[i][j] += w * pi0[ij];
         piH[i][j] += wH * piH0[ij];
        }
       Pi += w * c1->getPi0();
       PiH += wH * c1->getPiH0();
      }
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < i; j++) {
      pi[j][i] = pi[i][j];
      piH[j][i] = piH[i][j];
     }
    //--------- debug part: trace check, diag check, transversality check
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) {
      if (pi[i][j] != 0. && fabs(1.0 - pi[j][i] / pi[i][j]) > 1e-10)
       cout << "non-diag: " << pi[i][j] << "  " << pi[j][i] << endl;
     }
    //------ end debug
    //======= hydro applicability check (viscous corrections limiter):
    // double maxT0 = max(fabs((e+p)*vx*vx/(1.-vx*vx-vy*vy-vz*vz)+p),
    //   fabs((e+p)*vy*vy/(1.-vx*vx-vy*vy-vz*vz)+p)) ;
    double maxT0 = max((e + p) / (1. - vx * vx - vy * vy - vz * vz) - p,
                       (e + p) * (vx * vx + vy * vy + vz * vz) /
                               (1. - vx * vx - vy * vy - vz * vz) +
                           p);
    // double maxpi = max(fabs(pi[1][1]),fabs(pi[2][2])) ;
    double maxpi = 0.;
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++)
      if (fabs(pi[i][j]) > maxpi) maxpi = fabs(pi[i][j]);
    bool rescaled = false;
    if (maxT0 / maxpi < 1.0) {
     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
       pi[i][j] = 0.1 * pi[i][j] * maxT0 / maxpi;
       piH[i][j] = 0.1 * piH[i][j] * maxT0 / maxpi;
      }
     rescaled = true;
    }
    if (fabs(Pi) > p) {
     if (Pi != 0.) Pi = 1 * Pi / fabs(Pi) * p;     //modified .1 -> 1
     if (PiH != 0.) PiH = 1 * PiH / fabs(PiH) * p; //modified .1 -> 1
     rescaled = true;
    }
    if (rescaled)
     c->setViscCorrCutFlag(maxT0 / maxpi);
    else
     c->setViscCorrCutFlag(1.);
    // updating to the new values
    for (int i = 0; i < 4; i++)
     for (int j = 0; j <= i; j++) {
      c->setpi(i, j, pi[i][j]);
      c->setpiH(i, j, piH[i][j]);
     }
    c->setPi(Pi);
    c->setPiH(PiH);
    // source term  - (tau+dt)*delta_Q_(i+1)/delta_tau
    double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    double u[4];
    u[0] = gamma;
    u[1] = u[0] * vx;
    u[2] = u[0] * vy;
    u[3] = u[0] * vz;
    double flux[4];
    for (int i = 0; i < 4; i++)
     flux[i] = -tau * (c->getpi(0, i) + c->getPi() * u[0] * u[i]);
    flux[0] += tau * c->getPi();
	    c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
	   }  // advection loop (all cells)
 if (timingEnabled) {
#ifdef _OPENMP
#pragma omp single
#endif
  {
   timingAdvection = hydroElapsedSeconds(timingStart);
  }
 }
	 }
 if (timingEnabled) {
  std::cout << "Timing ISformal detail [sec]:"
            << " relaxSource=" << timingRelaxSource
            << " advection=" << timingAdvection
            << " total=" << timingRelaxSource + timingAdvection
            << std::endl;
  std::cout << "Timing ISformal relaxSource detail [sec]:"
            << " nsquantCpuSum=" << timingNSquantCpuSum
            << " nsquantCalls=" << timingNSquantCalls
            << " nsquantThermoHits=" << timingNSquantThermoHits
            << std::endl;
 }
}

// this procedure explicitly uses T_==0, X_==1, Y_==2, Z_==3
bool Hydro::visc_flux_value(Cell *left, Cell *right, int direction,
                            double flux[4]) {
 int ind2 = 0;
 double dxa = 0.;
 double tauMinusHalf = (cartesian) ? 1.0 : (tau - 0.5 * dt);
 double tauPlusHalf = (cartesian) ? 1.0 : (tau + 0.5 * dt);
 // exit if noth cells are not full with matter
 if (left->getM(direction) < 1. && right->getM(direction) < 1.) return false;

 if (direction == X_)
  dxa = f->getDx();
 else if (direction == Y_)
  dxa = f->getDy();
 else if (direction == Z_)
  dxa = f->getDz() * tauPlusHalf;
 double e, p, nb, nq, ns, vxl, vyl, vzl, vxr, vyr, vzr;
 // we need to know the velocities at both cell centers at (n+1/2) in order to
 // interpolate to
 // get the value at the interface
 left->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vxl, vyl, vzl);
 right->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vxr, vyr, vzr);
 vxl = 0.5 * (vxl + vxr);
 vyl = 0.5 * (vyl + vyr);
 vzl = 0.5 * (vzl + vzr);
 double v = sqrt(vxl * vxl + vyl * vyl + vzl * vzl);
 if (v > 1.) {
  vxl = 0.99 * vxl / v;
  vyl = 0.99 * vyl / v;
  vzl = 0.99 * vzl / v;
 }
 double gamma = 1. / sqrt(1. - v * v);
 double uuu[4] = {gamma, gamma * vxl, gamma * vyl, gamma * vzl};
 double gmumu[4] = {1., -1., -1., -1.};
 if (direction == X_)
  ind2 = 1;
 else if (direction == Y_)
  ind2 = 2;
 else if (direction == Z_)
  ind2 = 3;
 for (int ind1 = 0; ind1 < 4; ind1++) {
  flux[ind1] = 0.5 * (left->getpiH(ind1, ind2) + right->getpiH(ind1, ind2));
  if (ind1 == ind2)
   flux[ind1] += -0.5 * (left->getPiH() + right->getPiH()) *
                 gmumu[ind1];  // gmunu is diagonal
  flux[ind1] +=
      0.5 * (left->getPiH() + right->getPiH()) * uuu[ind1] * uuu[ind2];
 }
 for (int i = 0; i < 4; i++) flux[i] = flux[i] * tauMinusHalf * dt / dxa;
 return true;
}

void Hydro::visc_flux(Cell *left, Cell *right, int direction) {
 double flux[4];
 if (!visc_flux_value(left, right, direction, flux)) return;

 left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
 right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
}

void Hydro::accumulate_visc_fluxes(int ix, int iy, int iz) {
 Cell *cell = f->getCell(ix, iy, iz);
 double flux[4];

 if (ix > 0 && visc_flux_value(f->getCell(ix - 1, iy, iz), cell, X_, flux)) {
  cell->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
 }
 if (ix < f->getNX() - 1 &&
     visc_flux_value(cell, f->getCell(ix + 1, iy, iz), X_, flux)) {
  cell->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
 }
 if (iy > 0 && visc_flux_value(f->getCell(ix, iy - 1, iz), cell, Y_, flux)) {
  cell->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
 }
 if (iy < f->getNY() - 1 &&
     visc_flux_value(cell, f->getCell(ix, iy + 1, iz), Y_, flux)) {
  cell->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
 }
 if (iz > 0 && visc_flux_value(f->getCell(ix, iy, iz - 1), cell, Z_, flux)) {
  cell->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
 }
 if (iz < f->getNZ() - 1 &&
     visc_flux_value(cell, f->getCell(ix, iy, iz + 1), Z_, flux)) {
  cell->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
 }
}

void Hydro::accumulate_visc_fluxes_buffered() {
 const int nx = f->getNX();
 const int ny = f->getNY();
 const int nz = f->getNZ();

#ifdef _OPENMP
#pragma omp parallel
#endif
 {
  {
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
   for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
     for (int ix = 0; ix < nx - 1; ix++) {
      double flux[4];
      Cell *left = f->getCell(ix, iy, iz);
      Cell *right = f->getCell(ix + 1, iy, iz);
      if (visc_flux_value(left, right, X_, flux)) {
       left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
       right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
      }
     }
    }
  }
  {
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
   for (int iz = 0; iz < nz; iz++)
    for (int ix = 0; ix < nx; ix++) {
     for (int iy = 0; iy < ny - 1; iy++) {
      double flux[4];
      Cell *left = f->getCell(ix, iy, iz);
      Cell *right = f->getCell(ix, iy + 1, iz);
      if (visc_flux_value(left, right, Y_, flux)) {
       left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
       right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
      }
     }
    }
  }
  {
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
   for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
     for (int iz = 0; iz < nz - 1; iz++) {
      double flux[4];
      Cell *left = f->getCell(ix, iy, iz);
      Cell *right = f->getCell(ix, iy, iz + 1);
      if (visc_flux_value(left, right, Z_, flux)) {
       left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
       right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
      }
     }
    }
  }
 }
}

void Hydro::performStep(void) {
 // debugRiemann = false ; // turn off debug output
 const bool timingEnabled = std::getenv("VHLLE_TIMING") != nullptr;
 const auto timingTotalStart = std::chrono::steady_clock::now();
 auto timingPhaseStart = timingTotalStart;
 double timingUpdateM = 0.0, timingResetPredict = 0.0, timingHllePredict = 0.0;
 double timingSourcePredict = 0.0, timingHlleCorrect = 0.0, timingSourceCorrect = 0.0;
 double timingCorrectImag = 0.0, timingISformal = 0.0, timingViscFlux = 0.0;
 double timingViscSource = 0.0, timingCorrectImagFull = 0.0;
 const auto finishTimingPhase = [&](double &accumulator) {
  if (timingEnabled) {
   accumulator += hydroElapsedSeconds(timingPhaseStart);
   timingPhaseStart = std::chrono::steady_clock::now();
  }
 };

 f->updateM(tau, dt);
 finishTimingPhase(timingUpdateM);

 tau_z = dt / 2. / log(1 + dt / 2. / tau);

 //-----PREDICTOR-ideal
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->saveQprev();
    c->clearFlux();
   }
 finishTimingPhase(timingResetPredict);
 accumulate_hlle_fluxes_buffered(PREDICT);
 finishTimingPhase(timingHllePredict);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    Cell *c = f->getCell(ix, iy, iz);
    source_step(ix, iy, iz, PREDICT);
    c->updateQtoQhByFlux();
    c->clearFlux();
   }
 finishTimingPhase(timingSourcePredict);

 //----CORRECTOR-ideal

 tau_z = dt / log(1 + dt / tau);
 accumulate_hlle_fluxes_buffered(CORRECT);
 finishTimingPhase(timingHlleCorrect);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    Cell *c = f->getCell(ix, iy, iz);
    source_step(ix, iy, iz, CORRECT);
    c->updateByFlux();
    c->clearFlux();
   }
 finishTimingPhase(timingSourceCorrect);
 if (cartesian) {
  t += dt;
 }
 else {
  tau += dt;
 }
 
 f->correctImagCells();
 finishTimingPhase(timingCorrectImag);

 //===== viscous part ======
 if (trcoeff->isViscous()) {
  ISformal();  // evolution of viscous quantities according to IS equations
  finishTimingPhase(timingISformal);

  accumulate_visc_fluxes_buffered();
  finishTimingPhase(timingViscFlux);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++)
    for (int ix = 0; ix < f->getNX(); ix++) {
     visc_source_step(ix, iy, iz);
     f->getCell(ix, iy, iz)->updateByViscFlux();
     f->getCell(ix, iy, iz)->clearFlux();
    }
  finishTimingPhase(timingViscSource);
 } else {  // end viscous part
 }
 //==== finishing work ====
 f->correctImagCellsFull();
 finishTimingPhase(timingCorrectImagFull);
 if (timingEnabled) {
  std::cout << "Timing performStep detail [sec]:"
            << " updateM=" << timingUpdateM
            << " resetPredict=" << timingResetPredict
            << " hllePredict=" << timingHllePredict
            << " sourcePredict=" << timingSourcePredict
            << " hlleCorrect=" << timingHlleCorrect
            << " sourceCorrect=" << timingSourceCorrect
            << " correctImag=" << timingCorrectImag
            << " ISformal=" << timingISformal
            << " viscFlux=" << timingViscFlux
            << " viscSource=" << timingViscSource
            << " correctImagFull=" << timingCorrectImagFull
            << " total=" << hydroElapsedSeconds(timingTotalStart)
            << std::endl;
 }
}

void Hydro::addParticles(deque<Particle>* particles) {
 //==== particles coming in ====
 double particle_t = particles->front().getT();
 // a quick fix for the code to compile in tau-eta but this fn is not intended to run in tau-eta
 double current_t = (cartesian) ? t : tau;
 while (particle_t < current_t) {
   if (particles->size() > 0) {
    Particle particleToInject = particles->front();
    //cout << "particle at t, e:\n";
    //cout << particleToInject.getT() << " " << particleToInject.getE() << endl;
    f->addParticle(particleToInject);
    particles->pop_front();
    //cout << "particles in queue: " << particles->size() << endl;
    if (particles->size() > 0) particle_t = particles->front().getT();
    else particle_t = 1000.;
   } 
 }
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    f->getCell(ix, iy, iz)->updateByParticleSource();
    f->getCell(ix, iy, iz)->clearParticleSource();
   } 
}
