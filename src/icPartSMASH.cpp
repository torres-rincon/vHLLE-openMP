#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cfloat>
#include <vector>

#include "eos.h"
#include "eoChiral.h"
#include "rmn.h"
#include "fld.h"
#include "icPartSMASH.h"
#include "colour.h"

using namespace std;

namespace {
double elapsedSeconds(std::chrono::steady_clock::time_point start) {
 return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
}

IcPartSMASH::IcPartSMASH(Fluid* f, const char* filename, double _Rgt, double _Rgz,
                         int _smoothingType) {

 std::string tau0_warning = "Warning! tau0 will not be taken from the configuration file.\n"
                            "         It will be automatically set from SMASH initial state input file.\n";
 std::cout << yellow << tau0_warning << reset;
 
 nx = f->getNX();
 ny = f->getNY();
 nz = f->getNZ();
 dx = f->getDx();
 dy = f->getDy();
 dz = f->getDz();
 xmin = f->getX(0);
 xmax = f->getX(nx - 1);
 ymin = f->getY(0);
 ymax = f->getY(ny - 1);
 zmin = f->getZ(0);
 zmax = f->getZ(nz - 1);

 Rgx = _Rgt;
 Rgy = _Rgt;
 Rgz = _Rgz;
 nsmoothx = (int)(4.0 * Rgx / dx);  // smoothly distribute to +- this many cells
 nsmoothy = (int)(4.0 * Rgy / dy);
 nsmoothz = (int)(1.5 * Rgz / dz);

 isKernelInvariant = _smoothingType;

 T00 = new double**[nx];
 T0x = new double**[nx];
 T0y = new double**[nx];
 T0z = new double**[nx];
 QB = new double**[nx];
 QE = new double**[nx];
 QS = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  T00[ix] = new double*[ny];
  T0x[ix] = new double*[ny];
  T0y[ix] = new double*[ny];
  T0z[ix] = new double*[ny];
  QB[ix] = new double*[ny];
  QE[ix] = new double*[ny];
  QS[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   T00[ix][iy] = new double[nz];
   T0x[ix][iy] = new double[nz];
   T0y[ix][iy] = new double[nz];
   T0z[ix][iy] = new double[nz];
   QB[ix][iy] = new double[nz];
   QE[ix][iy] = new double[nz];
   QS[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    T00[ix][iy][iz] = 0.0;
    T0x[ix][iy][iz] = 0.0;
    T0y[ix][iy][iz] = 0.0;
    T0z[ix][iy][iz] = 0.0;
    QB[ix][iy][iz] = 0.0;
    QE[ix][iy][iz] = 0.0;
    QS[ix][iy][iz] = 0.0;
   }
  }
 }
#ifdef TSHIFT
 tau0 += tshift;
#endif
 // ---- read the events
 nevents = 0;
 ifstream fin(filename);
 if (!fin.good()) {
  std::cout << "I/O error with " << filename << endl;
  exit(1);
 }
 int np = 0;  // particle counter
 std::string line;
 std::istringstream instream;
 bool isBaryonIncluded {false};
 std::string quantity[14];
 
 double E_smash {0.0}; 
 int B_smash {0};
 int Q_smash {0};
 int S_smash {0};

 // Read first two lines of the header
 // to decide what type of input file 
 // SMASH 3.0: includes baryon number and strangeness
 getline(fin, line);
 getline(fin, line);
 instream.str(line);
 instream.seekg(0);
 instream.clear();
 instream >> quantity[0] >> quantity[1] >> quantity[2] >> quantity[3] >>
   quantity[4] >> quantity[5] >> quantity[6] >> quantity[7] >> quantity[8] >>
   quantity[9] >> quantity[10] >> quantity[11] >> quantity[12] >> quantity[13];
 if (quantity[11] != "") 
  isBaryonIncluded = true;
  
 // Read subsequent lines
 while (!fin.eof()) {
  getline(fin, line);
  instream.str(line);
  instream.seekg(0);
  instream.clear();
  // Read line
  // if n_B is included: read in all quantities, incl. nB and nS
  if (isBaryonIncluded)
   instream >> Tau_val >> X_val >> Y_val >> Eta_val >> Mt_val >> Px_val >>
               Py_val >> Rap_val >> Id_val >> Charge_val >> Baryon_val >> 
               Strangeness_val;
  // if not included: nB from PDG code, nS = 0
  else {
   instream >> Tau_val >> X_val >> Y_val >> Eta_val >> Mt_val >> Px_val >>
               Py_val >> Rap_val >> Id_val >> Charge_val;
   Baryon_val = 0;
   for (auto &code : PDG_Codes_Baryons) {
        if (Id_val == code) {
          // PDG code > 0: baryon, else anti baryon
          (code > 0) ? Baryon_val = 1 : Baryon_val = -1;
        }
   }
   Strangeness_val = 0;
  }
  // Fill arrays
  if (!instream.fail()) {
   tau0 = Tau_val;
   Tau.push_back(Tau_val);
   X.push_back(X_val);
   Y.push_back(Y_val);
   Eta.push_back(Eta_val);
   Mt.push_back(Mt_val);
   Px.push_back(Px_val);
   Py.push_back(Py_val);
   Rap.push_back(Rap_val);
   Id.push_back(Id_val);
   Charge.push_back(Charge_val);
   Baryon.push_back(Baryon_val);
   Strangeness.push_back(Strangeness_val);

#ifdef TSHIFT
   Eta[np] = TMath::ATanH(Tau[np] * sinh(Eta[np]) /
                          (Tau[np] * cosh(Eta[np]) + tshift));
   Tau[np] += tshift;
#endif
   E_smash += Mt_val * cosh(Rap_val);
   B_smash += Baryon_val;
   Q_smash += Charge_val;
   S_smash += Strangeness_val;

   np++;
  }
  else if (np > 0) {
   // cout<<"readF14:instream: failure reading data\n" ;
   // cout<<"stream = "<<instream.str()<<endl ;
   if (nevents % 100 == 0) {
    std::cout << "event = " << nevents << "  np = " << np << "\r";
    std::cout << flush;
   }
   makeSmoothTable(np);
   np = 0;

   // Clear arrays for next event
   Tau.clear();
   X.clear();
   Y.clear();
   Eta.clear();
   Mt.clear();
   Px.clear();
   Py.clear();
   Rap.clear();
   Id.clear();
   Charge.clear();
   Baryon.clear();
   Strangeness.clear();

   nevents++;
   // if(nevents>10000) return ;
  }
 }
 if (nevents > 1) {
  std::cout << "++ Warning: loaded " << nevents << "  initial SMASH events\n";
  std::cout << "Running vHLLE on averaged initial state from SMASH\n";
 }
 std::cout << "particle E = " << E_smash/ nevents << "  Nbar = " << B_smash / nevents 
      << "  Ncharge = " << Q_smash / nevents << " Ns = " << S_smash / nevents <<
      "Rgt = " << _Rgt << "Rgz = " << _Rgz << std::endl;
}

IcPartSMASH::~IcPartSMASH() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] T00[ix][iy];
   delete[] T0x[ix][iy];
   delete[] T0y[ix][iy];
   delete[] T0z[ix][iy];
   delete[] QB[ix][iy];
   delete[] QE[ix][iy];
   delete[] QS[ix][iy];
  }
  delete[] T00[ix];
  delete[] T0x[ix];
  delete[] T0y[ix];
  delete[] T0z[ix];
  delete[] QB[ix];
  delete[] QE[ix];
  delete[] QS[ix];
 }
 delete[] T00;
 delete[] T0x;
 delete[] T0y;
 delete[] T0z;
 delete[] QB;
 delete[] QE;
 delete[] QS;
}

// parameter smoothingType: 0 (default) for kernel contracted in eta, 1 for invariant kernel
// switches on isKernelInvariant
// if true only Rgz parameter is needed to smooth out particles
void IcPartSMASH::makeSmoothTable(int npart) {
 const bool timingEnabled = std::getenv("VHLLE_TIMING") != nullptr;
 double timingKernel = 0.0;
 double timingDeposit = 0.0;
 long long timingCells = 0;
 struct KernelCell {
  int ix;
  int iy;
  int iz;
  double weight;
  double coshWeight;
  double sinhWeight;
 };
 struct ParticleKernelData {
  double normGauss;
  std::vector<KernelCell> cells;
 };
 const std::size_t maxKernelCells = static_cast<std::size_t>(2 * nsmoothx + 1)
                                    * static_cast<std::size_t>(2 * nsmoothy + 1)
                                    * static_cast<std::size_t>(2 * nsmoothz + 1);
 const int kernelBlockSize = 32;
 std::vector<ParticleKernelData> blockData(kernelBlockSize);
 for (ParticleKernelData& particleData : blockData) {
  particleData.normGauss = 0.0;
  particleData.cells.reserve(maxKernelCells);
 }

 for (int blockStart = 0; blockStart < npart; blockStart += kernelBlockSize) {
  const int blockEnd = std::min(blockStart + kernelBlockSize, npart);
  const int blockCount = blockEnd - blockStart;
  auto timingStart = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int blockIndex = 0; blockIndex < blockCount; blockIndex++) {
   const int ip = blockStart + blockIndex;
   ParticleKernelData& particleData = blockData[blockIndex];
   std::vector<KernelCell>& kernelCells = particleData.cells;
   kernelCells.clear();
   particleData.normGauss = 0.0;

   int ixc = (int)round((X[ip] - xmin) / dx);
   int iyc = (int)round((Y[ip] - ymin) / dy);
   int izc = (int)round((Eta[ip] - zmin) / dz);
   // finding the norm
   double norm_gauss {0.0};
   const double gammaz = cosh(Rap[ip] - Eta[ip]);
   for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++) {
    for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++) {
     for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++) {
      // check if within the grid minus ghost cells
      if (ix >= 2 && ix < (nx - 2) && iy >= 2 && iy < (ny - 2)
         && iz >= 2 && iz < (nz - 2)) {
       const double xdiff = X[ip] - (xmin + ix * dx);
       const double ydiff = Y[ip] - (ymin + iy * dy);
       const double zdiff = Eta[ip] - (zmin + iz * dz);
       spatialVector rdiff {xdiff, ydiff, zdiff};
       double kernel;
       
       // calculate norm of the kernel
       if (isKernelInvariant) {
         velocityVector velocity = velocityHyperbolic(Mt[ip], Px[ip], Py[ip], Rap[ip], Eta[ip], zdiff, tau0);
         kernel = smoothingKernelInvariant(rdiff, velocity, Rgz, tau0);
       }
       else {
         kernel = smoothingKernel(rdiff, gammaz, tau0, Rgx, Rgy, Rgz);
       }
       norm_gauss += kernel;
       const double etaMomentumDiff = Rap[ip] - Eta[ip] + zdiff;
       kernelCells.push_back({ix, iy, iz, kernel, cosh(etaMomentumDiff), sinh(etaMomentumDiff)});
      }
     }
    }
   }
   particleData.normGauss = norm_gauss;
   for (KernelCell& cell : kernelCells) {
    double weight = cell.weight / norm_gauss;
    if (weight != weight || std::abs(weight) > DBL_MAX) {
     weight = 0.0;
    }
    cell.weight = weight;
   }
  }
  if (timingEnabled) {
   timingKernel += elapsedSeconds(timingStart);
  }

  timingStart = std::chrono::steady_clock::now();
  for (int blockIndex = 0; blockIndex < blockCount; blockIndex++) {
   const int ip = blockStart + blockIndex;
   const ParticleKernelData& particleData = blockData[blockIndex];
   if (timingEnabled) {
    timingCells += static_cast<long long>(particleData.cells.size());
   }

   for (const KernelCell& cell : particleData.cells) {
      const int ix = cell.ix;
      const int iy = cell.iy;
      const int iz = cell.iz;
      const double weight = cell.weight;
      
      // Add components of energy-momentum tensor
      T00[ix][iy][iz] += weight * Mt[ip] * cell.coshWeight;
      T0x[ix][iy][iz] += weight * Px[ip];
      T0y[ix][iy][iz] += weight * Py[ip];
      T0z[ix][iy][iz] += weight * Mt[ip] * cell.sinhWeight;

      // Add baryon number
      QB[ix][iy][iz] += weight * Baryon[ip];
      
      // Add electrice charge
      QE[ix][iy][iz] += weight * Charge[ip];

      // Add strangeness
      QS[ix][iy][iz] += weight * Strangeness[ip];
   }
  }
  if (timingEnabled) {
   timingDeposit += elapsedSeconds(timingStart);
  }
 }  // end particle block loop
 if (timingEnabled) {
  std::cout << "Timing SMASH fluidization [sec]: particles=" << npart
            << " cells=" << timingCells
            << " kernel=" << timingKernel
            << " deposit=" << timingDeposit
            << " total=" << timingKernel + timingDeposit << std::endl;
 }
}

void IcPartSMASH::setIC(Fluid* f, EoS* eos) {
 const bool timingEnabled = std::getenv("VHLLE_TIMING") != nullptr;
 auto timingStart = std::chrono::steady_clock::now();
 double E {0.0}, Px {0.0}, Py {0.0}, Pz {0.0};
 double Nb {0.0}, Nq {0.0}, Ns {0.0}, S {0.0};
 double centerX {0.0}, centerZ {0.0}, centerQT {0.0}, centerQZ {0.0};
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static) reduction(+:E,Px,Py,Pz,Nb,Nq,Ns,S)
#endif
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    double Q[7], e, p, nb, nq, ns, vx, vy, vz;
    Q[T_] = T00[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
    Q[X_] = T0x[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[Y_] = T0y[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[Z_] = T0z[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[NB_] = QB[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[NQ_] = QE[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[NS_] = QS[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    if (ix == nx / 2 && iy == ny / 2 && iz == nz / 2) {
     centerX = xmin + ix * dx;
     centerZ = zmin + iz * dz;
     centerQT = Q[T_];
     centerQZ = Q[Z_];
    }
    transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz);
    if (e < 1e-7 || fabs(f->getX(ix)) > 10. || fabs(f->getY(iy)) > 10. ||
        fabs(f->getZ(iz)) > 5.) {
     e = nb = nq = 0.0;
     vx = vy = vz = 0.0;
    }
    Cell* c = f->getCell(ix, iy, iz);
    c->setPrimVar(eos, tau0, e, nb, nq, ns, vx, vy, vz);
    if (e > 0.) c->setAllM(1.);
    const double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    double u[4] = {gamma, gamma * vx, gamma * vy, gamma * vz};
    double eta = zmin + iz * dz;
    E += tau0 * ((e + p) * u[0] * (u[0] * cosh(eta) + u[3] * sinh(eta)) -
                 p * cosh(eta)) *
         dx * dy * dz;
    Pz += tau0 * ((e + p) * u[0] * (u[0] * sinh(eta) + u[3] * cosh(eta)) -
                  p * sinh(eta)) *
          dx * dy * dz;
    Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
    Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
    Nb += tau0 * nb * u[0] * dx * dy * dz;
    Nq += tau0 * nq * u[0] * dx * dy * dz;
    Ns += tau0 * ns * u[0] * dx * dy * dz;
    S += tau0 * eos->s(e, nb, nq, ns) * u[0] * dx * dy * dz;
   }
 if (timingEnabled) {
  std::cout << "Timing SMASH setIC [sec]: cells="
            << static_cast<long long>(nx) * ny * nz
            << " total=" << elapsedSeconds(timingStart) << std::endl;
 }
 std::cout << "IC SMASH, center: " << centerX << "  " << centerZ
      << "  " << centerQT << "  " << centerQZ << endl;
 std::cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb << "  Ncharge = " << Nq
      << " Ns = " << Ns << std::endl
      << "  Px = " << Px << "  Py = " << Py << std::endl;
 std::cout << "initial_entropy S_ini = " << S << std::endl;
}

velocityVector velocityHyperbolic(double _mt, double _px, double _py, double _y, double _eta, double _etaDiff, double _tau){

    double _p0 = _mt * cosh(_y);
    double _pz = _mt * sinh(_y);
    double mass = sqrt(_p0 * _p0 - _px * _px - _py *_py - _pz * _pz);
    double pEta = _mt * sinh(_y - _eta + _etaDiff);
    
    double vx = _px / mass;
    double vy = _py / mass;
    double vz = pEta / mass;
    velocityVector v {vx, vy, vz};
    
    return v;
}

// kernel contracted only in eta
double smoothingKernel(spatialVector _r, double _gammaz, double _tau, 
  double _Rgx, double _Rgy, double _Rgz){
    
    double rxSquared =  _r.xdiff * _r.xdiff / _Rgx / _Rgx;
    double rySquared =  _r.ydiff * _r.ydiff / _Rgy / _Rgy;
    // z = eta
    // contracted in eta, tau factor from Milne metric
    double rzSquared = _r.zdiff * _r.zdiff / _Rgz / _Rgz;
    rzSquared *= _gammaz * _gammaz * _tau * _tau; 
    
    if (std::abs(_Rgz) < 1e-5){
      return exp(-rxSquared - rySquared);
    } // jun19 
    else {
      return exp(-rxSquared - rySquared - rzSquared);
    }
}

// invariant kernel as in SMASH
double smoothingKernelInvariant(spatialVector _r, velocityVector _v, double _R, 
  double _tau){

    // z = eta, tau factors from Milne metric
    double rSquared = _r.xdiff * _r.xdiff + _r.ydiff * _r.ydiff + 
      _r.zdiff * _r.zdiff * _tau * _tau;
    double ruScalar = _r.xdiff * _v.vx + _r.ydiff * _v.vy + 
      _r.zdiff *  _v.vz * _tau * _tau;
    
    return exp( (-rSquared - ruScalar * ruScalar) / _R / _R);
}
