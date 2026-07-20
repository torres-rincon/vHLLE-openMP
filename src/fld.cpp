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
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <array>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "cll.h"
#include "eos.h"
#include "trancoeff.h"
#include "cornelius.h"
#include "colour.h"
#include "particle.h"

#define OUTPI

// change to hadron EoS (e.g. Laine) to calculate v,T,mu at the surface
#define SWAP_EOS

using namespace std;

namespace output{  // a namespace containing all the output streams
  ofstream fkw, fkw_dim, fxvisc, fyvisc, fdiagvisc, fx,
     fy, fdiag, fz, faniz, f2d, ffreeze, fbeta;

  constexpr std::size_t fileBufferSize = 1 << 20;
  using FileBuffer = std::array<char, fileBufferSize>;
  FileBuffer fxviscBuffer, fyviscBuffer, fdiagviscBuffer, fxBuffer,
     fyBuffer, fdiagBuffer, fzBuffer, fanizBuffer, f2dBuffer,
     ffreezeBuffer, fbetaBuffer;
  constexpr std::size_t maxFreezeoutQueuedBuffers = 4;
  std::mutex freezeoutWriterMutex;
  std::condition_variable freezeoutWriterReady;
  std::condition_variable freezeoutWriterSpace;
  std::deque<std::string> freezeoutWriteQueue;
  std::thread freezeoutWriterThread;
  bool freezeoutWriterRunning = false;
  bool freezeoutWriterStopping = false;

  void setBuffer(ofstream &stream, FileBuffer &buffer) {
    stream.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
  }

  void appendFreezeoutField(std::string &buffer, double value) {
    char field[64];
    const int written = std::snprintf(field, sizeof(field), "%24.15g", value);
    if (written > 0) buffer.append(field, static_cast<std::size_t>(written));
  }

  void freezeoutWriterLoop() {
    std::unique_lock<std::mutex> lock(freezeoutWriterMutex);
    for (;;) {
      freezeoutWriterReady.wait(lock, [] {
        return freezeoutWriterStopping || !freezeoutWriteQueue.empty();
      });
      if (freezeoutWriterStopping && freezeoutWriteQueue.empty()) break;

      std::string buffer = std::move(freezeoutWriteQueue.front());
      freezeoutWriteQueue.pop_front();
      freezeoutWriterSpace.notify_one();

      lock.unlock();
      ffreeze.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      lock.lock();
    }
  }

  void startFreezeoutWriter() {
    std::lock_guard<std::mutex> lock(freezeoutWriterMutex);
    if (freezeoutWriterRunning) return;
    freezeoutWriteQueue.clear();
    freezeoutWriterStopping = false;
    freezeoutWriterRunning = true;
    freezeoutWriterThread = std::thread(freezeoutWriterLoop);
  }

  void writeFreezeout(std::string &&buffer) {
    if (buffer.empty()) return;

    std::unique_lock<std::mutex> lock(freezeoutWriterMutex);
    if (!freezeoutWriterRunning || freezeoutWriterStopping) {
      lock.unlock();
      ffreeze.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      return;
    }

    freezeoutWriterSpace.wait(lock, [] {
      return freezeoutWriteQueue.size() < maxFreezeoutQueuedBuffers ||
             freezeoutWriterStopping;
    });
    if (freezeoutWriterStopping) {
      lock.unlock();
      ffreeze.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      return;
    }
    freezeoutWriteQueue.emplace_back(std::move(buffer));
    lock.unlock();
    freezeoutWriterReady.notify_one();
  }

  void finishFreezeoutWriter() {
    {
      std::lock_guard<std::mutex> lock(freezeoutWriterMutex);
      if (!freezeoutWriterRunning) {
        ffreeze.flush();
        return;
      }
      freezeoutWriterStopping = true;
    }
    freezeoutWriterReady.notify_one();
    freezeoutWriterSpace.notify_all();
    if (freezeoutWriterThread.joinable()) freezeoutWriterThread.join();
    {
      std::lock_guard<std::mutex> lock(freezeoutWriterMutex);
      freezeoutWriterRunning = false;
      freezeoutWriterStopping = false;
    }
    ffreeze.flush();
  }
}

// returns the velocities in cartesian coordinates, fireball rest frame.
// Y=longitudinal rapidity of fluid
// Cartesian mode: Y returns vz, longitudinal velocity
void Fluid::getCMFvariables(Cell *c, double tau, double &e, double &nb,
                            double &nq, double &ns, double &vx, double &vy,
                            double &Y) {
 double p, vz;
 tau = (cartesian) ? 1.0 : tau;
 c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
 if (cartesian) {
  Y = vz;
 }
 else {
  double eta = getZ(c->getZ());
  //	Y = eta + TMath::ATanH(vz) ;
  Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
  vx = vx * cosh(Y - eta) / cosh(Y);
  vy = vy * cosh(Y - eta) / cosh(Y);
 }
}

Fluid::Fluid(EoS *_eos, EoS *_eosH, TransportCoeff *_trcoeff, int _nx, int _ny,
             int _nz, double _minx, double _maxx, double _miny, double _maxy,
             double _minz, double _maxz, double _dt, double eCrit, bool _cartesian) {
 eos = _eos;
 eosH = _eosH;
 trcoeff = _trcoeff;
 nx = _nx;
 ny = _ny;
 nz = _nz;
 minx = _minx;
 maxx = _maxx;
 miny = _miny;
 maxy = _maxy;
 minz = _minz;
 maxz = _maxz;
 dx = (maxx - minx) / (nx - 1);
 dy = (maxy - miny) / (ny - 1);
 dz = (maxz - minz) / (nz - 1);
 dt = _dt;
 cartesian = _cartesian;

 cell = new Cell[nx * ny * nz];

 cell0 = new Cell;
 cell0->setPrimVar(eos, 1.0, 0., 0., 0., 0., 0., 0.,
                   0.);  // tau is not important here, since *0
 cell0->setAllM(0.);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    getCell(ix, iy, iz)->setPrev(X_, getCell(ix - 1, iy, iz));
    getCell(ix, iy, iz)->setNext(X_, getCell(ix + 1, iy, iz));
    getCell(ix, iy, iz)->setPrev(Y_, getCell(ix, iy - 1, iz));
    getCell(ix, iy, iz)->setNext(Y_, getCell(ix, iy + 1, iz));
    getCell(ix, iy, iz)->setPrev(Z_, getCell(ix, iy, iz - 1));
    getCell(ix, iy, iz)->setNext(Z_, getCell(ix, iy, iz + 1));
    getCell(ix, iy, iz)->setPos(ix, iy, iz);
   }

 output_nt = 0;
 output_nx = 0;
 output_ny = 0;

 //---- Cornelius init
 double arrayDx[4] = {dt, dx, dy, dz};
 cornelius = new Cornelius;
 cornelius->init(4, eCrit, arrayDx);
 ecrit = eCrit;
 vEff = 0.;
 EtotSurf = 0.0;
}

Fluid::~Fluid() {
 output::finishFreezeoutWriter();
 delete[] cell;
 delete cell0;
}

void Fluid::initOutput(const char *dir, double tau0, bool hsOnly) {
 // hsOnly (default false):
 // if true only the hypersurface output is initialized
 freezeoutOnlyOutput = hsOnly;
 std::string outfreeze = dir;
 bool return_mkdir = std::filesystem::create_directory(outfreeze);
 cout << "mkdir returns: " << return_mkdir << endl;
 outfreeze.append("/freezeout.dat");
 checkOutputDirectory(outfreeze);
 outfreeze.append(".unfinished");
 checkOutputDirectory(outfreeze);
 output::setBuffer(output::ffreeze, output::ffreezeBuffer);
 output::ffreeze.open(outfreeze.c_str());
 output::startFreezeoutWriter();

 // initialize vorticity output if enabled
 if (vorticityOn) {
  string outbeta = dir;
  outbeta.append("/beta.dat");
  output::setBuffer(output::fbeta, output::fbetaBuffer);
  output::fbeta.open(outbeta.c_str());
 }

 if (!hsOnly) {
  string outx = dir;
  outx.append("/outx.dat");
  string outxvisc = dir;
  outxvisc.append("/outx.visc.dat");
  string outyvisc = dir;
  outyvisc.append("/outy.visc.dat");
  string outdiagvisc = dir;
  outdiagvisc.append("/diag.visc.dat");
  string outy = dir;
  outy.append("/outy.dat");
  string outdiag = dir;
  outdiag.append("/outdiag.dat");
  string outz = dir;
  outz.append("/outz.dat");
  string outaniz = dir;
  outaniz.append("/out.aniz.dat");
  string out2d = dir;
  out2d.append("/out2D.dat");
  output::setBuffer(output::fx, output::fxBuffer);
  output::fx.open(outx.c_str());
  output::setBuffer(output::fy, output::fyBuffer);
  output::fy.open(outy.c_str());
  output::setBuffer(output::fz, output::fzBuffer);
  output::fz.open(outz.c_str());
  output::setBuffer(output::fdiag, output::fdiagBuffer);
  output::fdiag.open(outdiag.c_str());
  output::setBuffer(output::f2d, output::f2dBuffer);
  output::f2d.open(out2d.c_str());
  output::setBuffer(output::fxvisc, output::fxviscBuffer);
  output::fxvisc.open(outxvisc.c_str());
  output::setBuffer(output::fyvisc, output::fyviscBuffer);
  output::fyvisc.open(outyvisc.c_str());
  output::setBuffer(output::fdiagvisc, output::fdiagviscBuffer);
  output::fdiagvisc.open(outdiagvisc.c_str());
  output::setBuffer(output::faniz, output::fanizBuffer);
  output::faniz.open(outaniz.c_str());
  //################################################################
  // important remark. for correct diagonal output, nx=ny must hold.
  //################################################################
  outputGnuplot(tau0);
  output::faniz << "#  tau  <<v_T>>  e_p  e'_p  (to compare with SongHeinz)\n";
 }
}

void Fluid::printDbetaHeader() {
  // print header of beta.dat if vorticity is enabled
  if (num_corona_cells == -1){
    //throw error that header cannot be printed due to num_corona_cells not set
    throw std::runtime_error("Error: num_corona_cells not set, cannot print header of beta.dat");
  } else {
    output::fbeta << "# The derivatives of β are given in Cartesian coordinates and include a factor of 1/2, such that ∂ₘβₙ=1/2 * ∂ₘ(uₙ/T)" << '\n';
    output::fbeta << "# Number of corona cells: " << num_corona_cells << '\n';
    output::fbeta << "#  τ  x  y  η  dΣ[0]  dΣ[1]  dΣ[2]  dΣ[3]  "
                  << "u[0]  u[1]  u[2]  u[3]  T  μB  μQ  μS  "
                  << "∂₀β₀  ∂₀β₁  ∂₀β₂  ∂₀β₃  ∂₁β₀  ∂₁β₁  ∂₁β₂  ∂₁β₃  "
                  << "∂₂β₀  ∂₂β₁  ∂₂β₂  ∂₂β₃  ∂₃β₀  ∂₃β₁  ∂₃β₂  ∂₃β₃  ϵ" << '\n';
  }
}

void Fluid::renameOutput(const char *dir) {
  // renames hypersurface output in case of clean exit
  output::finishFreezeoutWriter();
  std::string directory = dir;
  std::string oldFile = directory + "/freezeout.dat.unfinished";
  std::string newFile = directory + "/freezeout.dat";
  if (std::rename(oldFile.c_str(), newFile.c_str()) != 0)
		perror("Error renaming freezeout.dat.unfinished!");
}

void Fluid::checkOutputDirectory(std::string freezeoutFile) {
  // remove old freezeout.dat(.unfinished) file
  bool isFilePresent = std::filesystem::exists(freezeoutFile);
  std::string filename = std::filesystem::path(freezeoutFile).filename();
  if (isFilePresent) {
    std::string file_warning = "Warning! A '" + filename +
                               "' is present in your output directory.\n" +
                               "         It will be deleted automatically.\n";
    std::cout << yellow << file_warning << reset;
    bool isDeleted = std::filesystem::remove(freezeoutFile);
  }
}

void Fluid::correctImagCells(void) {
#ifdef _OPENMP
#pragma omp parallel
#endif
 {
 // Z
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   double Q[7], Qh[7];
   // left boundary
   getCell(ix, iy, 2)->getQ(Q);
   getCell(ix, iy, 1)->setQ(Q);
   getCell(ix, iy, 0)->setQ(Q);

   getCell(ix, iy, 2)->getQh(Qh);
   getCell(ix, iy, 1)->setQh(Qh);
   getCell(ix, iy, 0)->setQh(Qh);
   // right boundary
   getCell(ix, iy, nz - 3)->getQ(Q);
   getCell(ix, iy, nz - 2)->setQ(Q);
   getCell(ix, iy, nz - 1)->setQ(Q);

   getCell(ix, iy, nz - 3)->getQh(Qh);
   getCell(ix, iy, nz - 2)->setQh(Qh);
   getCell(ix, iy, nz - 1)->setQh(Qh);
  }
 // Y
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
 for (int ix = 0; ix < nx; ix++)
  for (int iz = 0; iz < nz; iz++) {
   double Q[7], Qh[7];
   // left boundary
   getCell(ix, 2, iz)->getQ(Q);
   getCell(ix, 1, iz)->setQ(Q);
   getCell(ix, 0, iz)->setQ(Q);

   getCell(ix, 2, iz)->getQh(Qh);
   getCell(ix, 1, iz)->setQh(Qh);
   getCell(ix, 0, iz)->setQh(Qh);
   // right boundary
   getCell(ix, ny - 3, iz)->getQ(Q);
   getCell(ix, ny - 2, iz)->setQ(Q);
   getCell(ix, ny - 1, iz)->setQ(Q);

   getCell(ix, ny - 3, iz)->getQh(Qh);
   getCell(ix, ny - 2, iz)->setQh(Qh);
   getCell(ix, ny - 1, iz)->setQh(Qh);
  }
 // X
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
 for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
   double Q[7], Qh[7];
   // left boundary
   getCell(2, iy, iz)->getQ(Q);
   getCell(1, iy, iz)->setQ(Q);
   getCell(0, iy, iz)->setQ(Q);

   getCell(2, iy, iz)->getQh(Qh);
   getCell(1, iy, iz)->setQh(Qh);
   getCell(0, iy, iz)->setQh(Qh);
   // right boundary
   getCell(nx - 3, iy, iz)->getQ(Q);
   getCell(nx - 2, iy, iz)->setQ(Q);
   getCell(nx - 1, iy, iz)->setQ(Q);

   getCell(nx - 3, iy, iz)->getQh(Qh);
   getCell(nx - 2, iy, iz)->setQh(Qh);
   getCell(nx - 1, iy, iz)->setQh(Qh);
  }
 }
}

void Fluid::correctImagCellsFull(void) {
#ifdef _OPENMP
#pragma omp parallel
#endif
 {
 // Z
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   double Q[7], _pi[4][4], _Pi;
   // left boundary
   getCell(ix, iy, 2)->getQ(Q);
   getCell(ix, iy, 1)->setQ(Q);
   getCell(ix, iy, 0)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, iy, 2)->getpi(i, j);
   _Pi = getCell(ix, iy, 2)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, iy, 0)->setpi(i, j, _pi[i][j]);
     getCell(ix, iy, 1)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, iy, 0)->setPi(_Pi);
   getCell(ix, iy, 1)->setPi(_Pi);
   // right boundary
   getCell(ix, iy, nz - 3)->getQ(Q);
   getCell(ix, iy, nz - 2)->setQ(Q);
   getCell(ix, iy, nz - 1)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = getCell(ix, iy, nz - 3)->getpi(i, j);
   _Pi = getCell(ix, iy, nz - 3)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, iy, nz - 2)->setpi(i, j, _pi[i][j]);
     getCell(ix, iy, nz - 1)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, iy, nz - 2)->setPi(_Pi);
   getCell(ix, iy, nz - 1)->setPi(_Pi);
  }
 // Y
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
 for (int ix = 0; ix < nx; ix++)
  for (int iz = 0; iz < nz; iz++) {
   double Q[7], _pi[4][4], _Pi;
   // left boundary
   getCell(ix, 2, iz)->getQ(Q);
   getCell(ix, 1, iz)->setQ(Q);
   getCell(ix, 0, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, 2, iz)->getpi(i, j);
   _Pi = getCell(ix, 2, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, 0, iz)->setpi(i, j, _pi[i][j]);
     getCell(ix, 1, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, 0, iz)->setPi(_Pi);
   getCell(ix, 1, iz)->setPi(_Pi);
   // right boundary
   getCell(ix, ny - 3, iz)->getQ(Q);
   getCell(ix, ny - 2, iz)->setQ(Q);
   getCell(ix, ny - 1, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = getCell(ix, ny - 3, iz)->getpi(i, j);
   _Pi = getCell(ix, ny - 3, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, ny - 2, iz)->setpi(i, j, _pi[i][j]);
     getCell(ix, ny - 1, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, ny - 2, iz)->setPi(_Pi);
   getCell(ix, ny - 1, iz)->setPi(_Pi);
  }
 // X
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
 for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
   double Q[7], _pi[4][4], _Pi;
   // left boundary
   getCell(2, iy, iz)->getQ(Q);
   getCell(1, iy, iz)->setQ(Q);
   getCell(0, iy, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = getCell(2, iy, iz)->getpi(i, j);
   _Pi = getCell(2, iy, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(0, iy, iz)->setpi(i, j, _pi[i][j]);
     getCell(1, iy, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(0, iy, iz)->setPi(_Pi);
   getCell(1, iy, iz)->setPi(_Pi);
   // right boundary
   getCell(nx - 3, iy, iz)->getQ(Q);
   getCell(nx - 2, iy, iz)->setQ(Q);
   getCell(nx - 1, iy, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = getCell(nx - 3, iy, iz)->getpi(i, j);
   _Pi = getCell(nx - 3, iy, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(nx - 2, iy, iz)->setpi(i, j, _pi[i][j]);
     getCell(nx - 1, iy, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(nx - 2, iy, iz)->setPi(_Pi);
   getCell(nx - 1, iy, iz)->setPi(_Pi);
  }
 }
}

void Fluid::updateM(double tau, double dt) {
#ifdef _OPENMP
#pragma omp parallel
#endif
 {
#ifdef _OPENMP
#pragma omp for collapse(3) schedule(static)
#endif
 for (int ix = 0; ix < getNX(); ix++)
  for (int iy = 0; iy < getNY(); iy++)
   for (int iz = 0; iz < getNZ(); iz++) {
    Cell *c = getCell(ix, iy, iz);
    c->setDM(X_, 0.);
    c->setDM(Y_, 0.);
    c->setDM(Z_, 0.);
    if (getCell(ix, iy, iz)->getMaxM() < 1.) {
     if (getCell(ix + 1, iy, iz)->getM(X_) >= 1. ||
         getCell(ix - 1, iy, iz)->getM(X_) >= 1.)
      c->setDM(X_, dt / dx);
     if (getCell(ix, iy + 1, iz)->getM(Y_) >= 1. ||
         getCell(ix, iy - 1, iz)->getM(Y_) >= 1.)
      c->setDM(Y_, dt / dy);
     if (getCell(ix, iy, iz + 1)->getM(Z_) >= 1. ||
         getCell(ix, iy, iz - 1)->getM(Z_) >= 1.)
      c->setDM(Z_, dt / dz / tau);

     if (c->getDM(X_) == 0. && c->getDM(Y_) == 0. && c->getDM(Z_) == 0.) {
      int diagCase = 0;
      
      if (getCell(ix + 1, iy + 1, iz)->getMaxM() >= 1. ||
          getCell(ix + 1, iy - 1, iz)->getMaxM() >= 1. ||
          getCell(ix - 1, iy + 1, iz)->getMaxM() >= 1. ||
          getCell(ix - 1, iy - 1, iz)->getMaxM() >= 1.)
        diagCase = 1;
          
      if (getCell(ix, iy + 1, iz + 1)->getMaxM() >= 1. ||
          getCell(ix, iy + 1, iz - 1)->getMaxM() >= 1. ||
          getCell(ix, iy - 1, iz + 1)->getMaxM() >= 1. ||
          getCell(ix, iy - 1, iz - 1)->getMaxM() >= 1.)
        diagCase = 2;
          
      if (getCell(ix + 1, iy, iz + 1)->getMaxM() >= 1. ||
          getCell(ix + 1, iy, iz - 1)->getMaxM() >= 1. ||
          getCell(ix - 1, iy, iz + 1)->getMaxM() >= 1. ||
          getCell(ix - 1, iy, iz - 1)->getMaxM() >= 1.)
        diagCase = 3;

      switch(diagCase)
      {
       case 1:
        c->setDM(X_, 0.707 * dt / dx);
        c->setDM(Y_, 0.707 * dt / dy);
        break;
       case 2:
        c->setDM(Y_, 0.707 * dt / dy);
        c->setDM(Z_, 0.707 * dt / dz / tau);
        break;
       case 3:
        c->setDM(X_, 0.707 * dt / dx);
        c->setDM(Z_, 0.707 * dt / dz / tau);
        break;
        default:
         break;
      }
     }
	   }  // if
	   }

#ifdef _OPENMP
#pragma omp for collapse(3) schedule(static)
#endif
 for (int ix = 0; ix < getNX(); ix++)
  for (int iy = 0; iy < getNY(); iy++)
   for (int iz = 0; iz < getNZ(); iz++) {
    Cell *c = getCell(ix, iy, iz);
    c->addM(X_, c->getDM(X_));
    c->addM(Y_, c->getDM(Y_));
    c->addM(Z_, c->getDM(Z_));
   }
 }
}

void Fluid::outputGnuplot(double tau) {
 // in Cartesian frame:
 // time t is passed as tau parameter. 
 // in such case getCMFvariables() re-sets tau=1 internally
 double e, p, nb, nq, ns, T, mub, muq, mus, vx, vy, vz;

 // X direction
 for (int ix = 0; ix < nx; ix++) {
  double x = getX(ix);
  Cell *c = getCell(ix, ny / 2, nz / 2);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
  output::fx << setw(14) << tau << setw(14) << x << setw(14) << vx << setw(14) << vy
        << setw(14) << e << setw(14) << nb << setw(14) << T << setw(14) << mub;
  output::fx << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1) << setw(14)
        << c->getpi(0, 2);
  output::fx << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1) << setw(14)
        << c->getpi(1, 2);
  output::fx << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2) << setw(14)
        << c->getpi(2, 3);
  output::fx << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
        << c->getViscCorrCutFlag() << endl;
 }
 output::fx << endl;

 // Y direction
 for (int iy = 0; iy < ny; iy++) {
  double y = getY(iy);
  Cell *c = getCell(nx / 2, iy, nz / 2);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
  output::fy << setw(14) << tau << setw(14) << y << setw(14) << vy << setw(14) << vx
        << setw(14) << e << setw(14) << nb << setw(14) << T << setw(14) << mub;
  output::fy << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1) << setw(14)
        << c->getpi(0, 2);
  output::fy << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1) << setw(14)
        << c->getpi(1, 2);
  output::fy << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2) << setw(14)
        << c->getpi(2, 3);
  output::fy << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
        << c->getViscCorrCutFlag() << endl;
 }
 output::fy << endl;

 // diagonal
 for (int ix = 0; ix < nx; ix++) {
  double x = getY(ix);
  Cell *c = getCell(ix, ix, nz / 2);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
  output::fdiag << setw(14) << tau << setw(14) << sqrt(2.) * x << setw(14) << vx
           << setw(14) << vy << setw(14) << e << setw(14) << nb << setw(14) << T
           << setw(14) << mub << endl;
  output::fdiag << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
           << setw(14) << c->getpi(0, 2);
  output::fdiag << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
           << setw(14) << c->getpi(1, 2);
  output::fdiag << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
           << setw(14) << c->getpi(2, 3);
  output::fdiag << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
           << c->getViscCorrCutFlag() << endl;
 }
 output::fdiag << endl;

 // Z direction
 for (int iz = 0; iz < nz; iz++) {
  double z = getZ(iz);
  Cell *c = getCell(nx / 2, ny / 2, iz);
  getCMFvariables(getCell(nx / 2, ny / 2, iz), tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
  output::fz << setw(14) << tau << setw(14) << z << setw(14) << vz << setw(14) << vx
        << setw(14) << e << setw(14) << nb << setw(14) << T << setw(14) << mub;
  output::fz << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1) << setw(14)
        << c->getpi(0, 2);
  output::fz << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1) << setw(14)
        << c->getpi(1, 2);
  output::fz << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2) << setw(14)
        << c->getpi(2, 3);
  output::fz << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
        << c->getViscCorrCutFlag() << endl;
 }
 output::fz << endl;
}

// input: geom. rapidity + velocities in Bjorken frame, --> output: velocities
// in lab.frame
void transformToLab(double eta, double &vx, double &vy, double &vz) {
 const double Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
 vx = vx * cosh(Y - eta) / cosh(Y);
 vy = vy * cosh(Y - eta) / cosh(Y);
 vz = tanh(Y);
}

// return value: the number of surface elements reconstructed at the current timestep
int Fluid::outputSurface(double tau, bool extendFO) {
 static double EtotSurfPos = 0.0, EtotSurfPosVisc = 0.0, EtotSurfNeg = 0.0;
 static double nbTotSurf = 0.0, nbTotSurfPos = 0.0, nbTotSurfNeg = 0.0;
 double e, p, vx, vy, vz;
 double E = 0., Ecore = 0., Efull = 0., S = 0., Px = 0., vt_num = 0., vt_den = 0.,
        vxvy_num = 0., vxvy_den = 0., pi0x_num = 0., pi0x_den = 0.,
        txxyy_num = 0., txxyy_den = 0., Nb1 = 0., Nb2 = 0., Nbcore = 0.,
        eps_p = 0.;
 const double gmumu[4] = {1., -1., -1., -1.};
 int nelements = 0, nsusp = 0;  // all surface elements and suspicious ones
 int nCoreCells = 0,
     nCoreCutCells = 0;  // cells with e>eCrit and cells with cut visc.corr.
#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 const bool writeDiagnostics = !freezeoutOnlyOutput;
 std::string freezeoutOutput;
 freezeoutOutput.reserve(1024 * 1024);
 struct SurfaceThermoCell {
  double e, nb, nq, ns, vx, vy, vz, p, T, s;
 };
 struct SurfaceDiagnosticCell {
  double E, Ecore, Efull, S, Px, vt_num, vt_den, vxvy_num, vxvy_den,
         txxyy_num, txxyy_den, pi0x_num, pi0x_den, Nb1, Nb2, Nbcore;
  int nCoreCells, nCoreCutCells;
 };
 const auto surfaceEnergyIndex = [this](int ix, int iy, int iz) {
  return (ix * ny + iy) * nz + iz;
 };
 const int interiorNx = max(0, nx - 4);
 const int interiorNy = max(0, ny - 4);
 const int interiorNz = max(0, nz - 4);
 const auto surfaceThermoIndex = [interiorNy, interiorNz](int ix, int iy, int iz) {
  return ((ix - 2) * interiorNy + (iy - 2)) * interiorNz + (iz - 2);
 };
 std::vector<double> surfaceEnergyNow(nx * ny * nz);
 std::vector<double> surfaceEnergyPrev(nx * ny * nz);
 std::vector<SurfaceThermoCell> surfaceThermo;
 std::vector<SurfaceDiagnosticCell> surfaceDiagnostic;
 std::vector<double> surfaceCoshInt;
 std::vector<double> surfaceSinhInt;
 if (writeDiagnostics) {
  surfaceThermo.resize(interiorNx * interiorNy * interiorNz);
  surfaceDiagnostic.resize(interiorNx * interiorNy * interiorNz);
 }
 if (writeDiagnostics && !cartesian) {
  surfaceCoshInt.resize(nz);
  surfaceSinhInt.resize(nz);
  for (int iz = 0; iz < nz; iz++) {
   const double etaValue = getZ(iz);
   surfaceCoshInt[iz] = (sinh(etaValue + 0.5 * dz) - sinh(etaValue - 0.5 * dz)) / dz;
   surfaceSinhInt[iz] = (cosh(etaValue + 0.5 * dz) - cosh(etaValue - 0.5 * dz)) / dz;
  }
 }
 if (writeDiagnostics) {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int ix = 0; ix < nx; ix++)
   for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
     double e_now, e_prev;
     double p_now, nb_now, nq_now, ns_now, vx_now, vy_now, vz_now;
     double p_prev, nb_prev, nq_prev, ns_prev, vx_prev, vy_prev, vz_prev;
     Cell *cc = getCell(ix, iy, iz);
     if (cartesian) {
      cc->getPrimVar(eos, 1.0, e_now, p_now, nb_now, nq_now, ns_now, vx_now, vy_now, vz_now);
      cc->getPrimVarPrev(eos, 1.0, e_prev, p_prev, nb_prev, nq_prev, ns_prev, vx_prev, vy_prev, vz_prev);
     }
     else {
      cc->getPrimVar(eos, tau, e_now, p_now, nb_now, nq_now, ns_now, vx_now, vy_now, vz_now);
      cc->getPrimVarPrev(eos, tau - dt, e_prev, p_prev, nb_prev, nq_prev, ns_prev, vx_prev, vy_prev, vz_prev);
     }
     const int index = surfaceEnergyIndex(ix, iy, iz);
     surfaceEnergyNow[index] = e_now;
     surfaceEnergyPrev[index] = e_prev;
     if (ix >= 2 && ix < nx - 2 && iy >= 2 && iy < ny - 2 && iz >= 2 && iz < nz - 2) {
      double vx_cache = vx_now;
      double vy_cache = vy_now;
      double vz_cache = vz_now;
      if (!cartesian) {
       const double etaValue = getZ(iz);
       const double rapidity = etaValue + 0.5 * log((1. + vz_cache) / (1. - vz_cache));
       vx_cache = vx_cache * cosh(rapidity - etaValue) / cosh(rapidity);
       vy_cache = vy_cache * cosh(rapidity - etaValue) / cosh(rapidity);
       vz_cache = rapidity;
      }
      double T_cache, mub_cache, muq_cache, mus_cache, p_cache;
      eos->eos(e_now, nb_now, nq_now, ns_now, T_cache, mub_cache, muq_cache, mus_cache, p_cache);
      const double s_cache = eos->s(e_now, nb_now, nq_now, ns_now);
      const int thermoIndex = surfaceThermoIndex(ix, iy, iz);
      surfaceThermo[thermoIndex].e = e_now;
      surfaceThermo[thermoIndex].nb = nb_now;
      surfaceThermo[thermoIndex].nq = nq_now;
      surfaceThermo[thermoIndex].ns = ns_now;
      surfaceThermo[thermoIndex].vx = vx_cache;
      surfaceThermo[thermoIndex].vy = vy_cache;
      surfaceThermo[thermoIndex].vz = vz_cache;
      surfaceThermo[thermoIndex].p = p_cache;
      surfaceThermo[thermoIndex].T = T_cache;
      surfaceThermo[thermoIndex].s = s_cache;

      double Q_diag[7];
      cc->getQ(Q_diag);
      SurfaceDiagnosticCell diagnostic{};
      if (cartesian) {
       const double denominator = 1. - vx_cache * vx_cache - vy_cache * vy_cache - vz_cache * vz_cache;
       diagnostic.E = (e_now + p_cache) / denominator - p_cache;
       if (e_now > ecrit) {
        diagnostic.Ecore = diagnostic.E;
        diagnostic.Nbcore = Q_diag[NB_];
        diagnostic.nCoreCells = 1;
        if (cc->getViscCorrCutFlag() < 0.9) diagnostic.nCoreCutCells = 1;
       }
       diagnostic.Nb1 = Q_diag[NB_];
       diagnostic.Nb2 = nb_now / sqrt(denominator);
       diagnostic.Efull = diagnostic.E;
       if (trcoeff->isViscous()) diagnostic.Efull += cc->getpi(0, 0);
       double deltas = 0.;
       if (trcoeff->isViscous())
        for (int i = 0; i < 4; i++)
         for (int j = 0; j < 4; j++)
          deltas += pow(cc->getpi(i, j), 2) * gmumu[i] * gmumu[j];
       if (T_cache > 0.02) {
        double s_value = s_cache;
        s_value += 1.5 * deltas / ((e_now + p_cache) * T_cache);
        diagnostic.S = s_value / sqrt(denominator);
       }
      }
      else {
       const double cosh_int = surfaceCoshInt[iz];
       const double sinh_int = surfaceSinhInt[iz];
       const double tanh_vz = tanh(vz_cache);
       const double denominator = 1. - vx_cache * vx_cache - vy_cache * vy_cache - tanh_vz * tanh_vz;
       const double flow_factor = cosh_int - tanh_vz * sinh_int;
       diagnostic.E = tau * (e_now + p_cache) / denominator * flow_factor -
                      tau * p_cache * cosh_int;
       diagnostic.Nb1 = Q_diag[NB_];
       diagnostic.Nb2 = tau * nb_now * flow_factor / sqrt(denominator);
       diagnostic.Efull = tau * (e_now + p_cache + cc->getPi()) / denominator * flow_factor -
                          tau * (p_cache + cc->getPi()) * cosh_int;
       if (trcoeff->isViscous())
        diagnostic.Efull +=
           tau * cc->getpi(0, 0) * cosh_int + tau * cc->getpi(0, 3) * sinh_int;
       if (e_now > ecrit) {
        diagnostic.nCoreCells = 1;
        if (cc->getViscCorrCutFlag() < 0.9) diagnostic.nCoreCutCells = 1;
       }
       double deltas = 0.;
       if (trcoeff->isViscous())
        for (int i = 0; i < 4; i++)
         for (int j = 0; j < 4; j++)
          deltas += pow(cc->getpi(i, j), 2) * gmumu[i] * gmumu[j];
       if (T_cache > 0.02) {
        double s_value = s_cache;
        s_value += 1.5 * deltas / ((e_now + p_cache) * T_cache);
        diagnostic.S = tau * s_value * flow_factor / sqrt(denominator);
       }
       diagnostic.Px = tau * (e_now + p_cache) * vx_cache / denominator;
       if (iz > nz/2-3 and iz < nz/2+3) {
        diagnostic.vxvy_num = e_now * (fabs(vx_cache) - fabs(vy_cache));
        diagnostic.vxvy_den = e_now;
        diagnostic.vt_den = e_now / sqrt(1. - vx_cache * vx_cache - vy_cache * vy_cache);
        diagnostic.vt_num = e_now / sqrt(1. - vx_cache * vx_cache - vy_cache * vy_cache) *
                            sqrt(vx_cache * vx_cache + vy_cache * vy_cache);
        diagnostic.txxyy_num = (e_now + p_cache) / denominator *
                               (vx_cache * vx_cache - vy_cache * vy_cache);
        diagnostic.txxyy_den = (e_now + p_cache) / denominator *
                                   (vx_cache * vx_cache + vy_cache * vy_cache) +
                               2. * p_cache;
       }
       diagnostic.pi0x_num = e_now / denominator * fabs(cc->getpi(0, 1));
       diagnostic.pi0x_den = e_now / denominator;
      }
      surfaceDiagnostic[thermoIndex] = diagnostic;
     }
    }
 }
 else {
  const double tauNow = cartesian ? 1.0 : tau;
  const double tauPrev = cartesian ? 1.0 : tau - dt;
  std::vector<unsigned char> upperEnergyAbove(nx * ny * nz);
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int ix = 0; ix < nx; ix++)
   for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
     double QNow[7], QPrev[7];
     Cell *cc = getCell(ix, iy, iz);
     cc->getQ(QNow);
     cc->getQprev(QPrev);
     const int index = surfaceEnergyIndex(ix, iy, iz);
     upperEnergyAbove[index] =
        (QNow[T_] / tauNow >= ecrit || QPrev[T_] / tauPrev >= ecrit) ? 1 : 0;
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int ix = 2; ix < nx - 1; ix++)
   for (int iy = 2; iy < ny - 1; iy++)
    for (int iz = 2; iz < nz - 1; iz++) {
     bool needsExactEnergy = false;
     for (int jx = max(2, ix - 1); jx <= min(nx - 2, ix + 1) && !needsExactEnergy; jx++)
      for (int jy = max(2, iy - 1); jy <= min(ny - 2, iy + 1) && !needsExactEnergy; jy++)
       for (int jz = max(2, iz - 1); jz <= min(nz - 2, iz + 1); jz++) {
        if (upperEnergyAbove[surfaceEnergyIndex(jx, jy, jz)]) {
         needsExactEnergy = true;
         break;
        }
       }
     if (!needsExactEnergy) continue;

     double e_now, e_prev;
     double p_now, nb_now, nq_now, ns_now, vx_now, vy_now, vz_now;
     double p_prev, nb_prev, nq_prev, ns_prev, vx_prev, vy_prev, vz_prev;
     Cell *cc = getCell(ix, iy, iz);
     if (cartesian) {
      cc->getPrimVar(eos, 1.0, e_now, p_now, nb_now, nq_now, ns_now, vx_now, vy_now, vz_now);
      cc->getPrimVarPrev(eos, 1.0, e_prev, p_prev, nb_prev, nq_prev, ns_prev, vx_prev, vy_prev, vz_prev);
     }
     else {
      cc->getPrimVar(eos, tau, e_now, p_now, nb_now, nq_now, ns_now, vx_now, vy_now, vz_now);
      cc->getPrimVarPrev(eos, tau - dt, e_prev, p_prev, nb_prev, nq_prev, ns_prev, vx_prev, vy_prev, vz_prev);
     }
     const int index = surfaceEnergyIndex(ix, iy, iz);
     surfaceEnergyNow[index] = e_now;
     surfaceEnergyPrev[index] = e_prev;
    }
 }
 if (writeDiagnostics) {
  output::f2d << '\n';
  for (int ix = 2; ix < nx - 2; ix++)
   for (int iy = 2; iy < ny - 2; iy++)
    for (int iz = 2; iz < nz - 2; iz++) {
     const SurfaceThermoCell &thermo = surfaceThermo[surfaceThermoIndex(ix, iy, iz)];
     e = thermo.e;
     vx = thermo.vx;
     vy = thermo.vy;
     vz = thermo.vz;
     p = thermo.p;
     const SurfaceDiagnosticCell &diagnostic =
        surfaceDiagnostic[surfaceThermoIndex(ix, iy, iz)];
     E += diagnostic.E;
     Ecore += diagnostic.Ecore;
     Efull += diagnostic.Efull;
     S += diagnostic.S;
     Px += diagnostic.Px;
     vt_num += diagnostic.vt_num;
     vt_den += diagnostic.vt_den;
     vxvy_num += diagnostic.vxvy_num;
     vxvy_den += diagnostic.vxvy_den;
     txxyy_num += diagnostic.txxyy_num;
     txxyy_den += diagnostic.txxyy_den;
     pi0x_num += diagnostic.pi0x_num;
     pi0x_den += diagnostic.pi0x_den;
     Nb1 += diagnostic.Nb1;
     Nb2 += diagnostic.Nb2;
     Nbcore += diagnostic.Nbcore;
     nCoreCells += diagnostic.nCoreCells;
     nCoreCutCells += diagnostic.nCoreCutCells;
     //---- inf check
     if (std::isinf(E)) {
      cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
          << setw(14) << vy << setw(14) << vz << endl;
      exit(1);
     }
    }
 }
 struct SurfaceElementContribution {
  double dVEff = 0.0;
  double dEtotSurf = 0.0;
  double nbDVEff = 0.0;
  double dEtotSurfVisc = 0.0;
  bool flowingOut = false;
  bool suspicious = false;
  std::string warningOutput;
 };
 struct SurfaceLineResult {
  std::string freezeoutOutput;
  std::string betaOutput;
  std::string warningOutput;
  std::vector<SurfaceElementContribution> contributions;
  int nelements = 0;
 };
 const auto surfaceLineIndex = [interiorNy](int ix, int iy) {
  return (ix - 2) * interiorNy + (iy - 2);
 };
 std::vector<SurfaceLineResult> surfaceLines(interiorNx * interiorNy);
#ifdef _OPENMP
#pragma omp parallel
#endif
 {
  double corneliusDx[4] = {dt, dx, dy, dz};
  Cornelius localCornelius;
  localCornelius.init(4, ecrit, corneliusDx);
  double ccubeValues[2][2][2][2];
  double *ccubeLevel3[2][2][2];
  double **ccubeLevel2[2][2];
  double ***ccubeLevel1[2];
  for (int i1 = 0; i1 < 2; i1++) {
   ccubeLevel1[i1] = ccubeLevel2[i1];
   for (int i2 = 0; i2 < 2; i2++) {
    ccubeLevel2[i1][i2] = ccubeLevel3[i1][i2];
    for (int i3 = 0; i3 < 2; i3++) {
     ccubeLevel3[i1][i2][i3] = ccubeValues[i1][i2][i3];
    }
   }
  }
  double ****ccube = ccubeLevel1;
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
  for (int ix = 2; ix < nx - 2; ix++)
   for (int iy = 2; iy < ny - 2; iy++) {
    SurfaceLineResult lineResult;
    lineResult.freezeoutOutput.reserve(4096);
    for (int iz = 2; iz < nz - 2; iz++) {
     int surfaceAbove = 0;
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++) {
        const int index = surfaceEnergyIndex(ix + jx, iy + jy, iz + jz);
        const double energyNow = surfaceEnergyNow[index];
        const double energyPrev = surfaceEnergyPrev[index];
        ccube[1][jx][jy][jz] = energyNow;
        ccube[0][jx][jy][jz] = energyPrev;
        if (energyNow >= ecrit) surfaceAbove++;
        if (energyPrev >= ecrit) surfaceAbove++;
       }
     int Nsegm = 0;
     if (surfaceAbove != 0 && surfaceAbove != 16) {
      localCornelius.find_surface_4d(ccube);
      Nsegm = localCornelius.get_Nelements();
     }
     if (Nsegm > 0) {
      double QCube[2][2][2][2][7];
      double piSquare[2][2][2][10], PiSquare[2][2][2];
      std::unique_ptr<Block3D> dbetaBlock =
              vorticityOn ? std::make_unique<Block3D>(
                                2, std::vector<std::vector<Matrix2D>>(
                                       2, std::vector<Matrix2D>(
                                              2, Matrix2D{{0.0, 0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0, 0.0}})))
                          : nullptr;
      for (int jx = 0; jx < 2; jx++)
       for (int jy = 0; jy < 2; jy++)
        for (int jz = 0; jz < 2; jz++) {
         Cell *cc = getCell(ix + jx, iy + jy, iz + jz);
         cc->getQ(QCube[1][jx][jy][jz]);
         cc->getQprev(QCube[0][jx][jy][jz]);
         for (int ii = 0; ii < 4; ii++) {
          for (int jj = 0; jj <= ii; jj++) {
           piSquare[jx][jy][jz][index44(ii, jj)] = cc->getpi(ii, jj);
          }
         }
         PiSquare[jx][jy][jz] = cc->getPi();
         if(vorticityOn) {
          if(!dbetaBlock) {
            std::runtime_error("dbetaBlock is a nullptr");
          }
          for(int column = 0 ; column < 4 ; column++) {
           for(int row = 0 ; row < 4 ; row++) {
            (*dbetaBlock)[jx][jy][jz][column][row] = cc -> getDbeta(column, row);
           }
          }
         }
        }
      for (int isegm = 0; isegm < Nsegm; isegm++) {
       SurfaceElementContribution contribution;
       if(vorticityOn) {
        output::appendFreezeoutField(lineResult.betaOutput, tau + localCornelius.get_centroid_elem(isegm, 0));
        output::appendFreezeoutField(lineResult.betaOutput, getX(ix) + localCornelius.get_centroid_elem(isegm, 1));
        output::appendFreezeoutField(lineResult.betaOutput, getY(iy) + localCornelius.get_centroid_elem(isegm, 2));
        output::appendFreezeoutField(lineResult.betaOutput, getZ(iz) + localCornelius.get_centroid_elem(isegm, 3));
       }
       double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
              musC = 0., piC[10] = {0.}, PiC = 0., nbC = 0.,
              nqC = 0.;
       double QC[7] = {0.};
       double eC = 0., pC = 0.;
       double wCenT[2] = {1. - localCornelius.get_centroid_elem(isegm, 0) / dt,
                          localCornelius.get_centroid_elem(isegm, 0) / dt};
       double wCenX[2] = {1. - localCornelius.get_centroid_elem(isegm, 1) / dx,
                          localCornelius.get_centroid_elem(isegm, 1) / dx};
       double wCenY[2] = {1. - localCornelius.get_centroid_elem(isegm, 2) / dy,
                          localCornelius.get_centroid_elem(isegm, 2) / dy};
       double wCenZ[2] = {1. - localCornelius.get_centroid_elem(isegm, 3) / dz,
                          localCornelius.get_centroid_elem(isegm, 3) / dz};
       for (int jt = 0; jt < 2; jt++)
        for (int jx = 0; jx < 2; jx++)
         for (int jy = 0; jy < 2; jy++)
          for (int jz = 0; jz < 2; jz++)
           for (int i = 0; i < 7; i++) {
            QC[i] += QCube[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                     wCenY[jy] * wCenZ[jz];
           }
       if (!cartesian) {
        for (int i = 0; i < 7; i++)
         QC[i] = QC[i] / (tau + localCornelius.get_centroid_elem(isegm, 0));
       }
       double _ns = 0.0;
       transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
       eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
       if (TC > 0.4 || fabs(mubC) > 0.85) {
        std::ostringstream warning;
        warning << "#### Error (surface): high T/mu_b (T=" << TC << "/mu_b=" << mubC << ") #### FREEZEOUT \n";
        if (writeDiagnostics)
         contribution.warningOutput = warning.str();
        else
         lineResult.warningOutput.append(warning.str());
       }
       contribution.suspicious = eC > ecrit * 2.0 || eC < ecrit * 0.5;
       std::unique_ptr<Matrix2D> dbetaInterpolated = vorticityOn
          ? std::make_unique<Matrix2D>(Matrix2D(4, std::vector<double>(4, 0.0)))
          : nullptr;
       for (int jx = 0; jx < 2; jx++)
        for (int jy = 0; jy < 2; jy++)
         for (int jz = 0; jz < 2; jz++) {
          for (int ii = 0; ii < 10; ii++) {
           piC[ii] += piSquare[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
          }
          PiC += PiSquare[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
          if (vorticityOn) {
           if(!dbetaInterpolated) {
             std::runtime_error("dbetaInterpolated is a nullptr");
           }
           for(int column = 0 ; column < 4 ; column++) {
            for(int row = 0 ; row < 4 ; row++) {
             (*dbetaInterpolated)[column][row] +=
             (*dbetaBlock)[jx][jy][jz][column][row] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
            }
           }
          }
         }
       double etaC = 0.0;
       double ch = 0.0, sh = 0.0;
       if (!cartesian){
        etaC = getZ(iz) + localCornelius.get_centroid_elem(isegm, 3);
        transformToLab(etaC, vxC, vyC, vzC);
       }
       double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
       if (v2C > 1.) {
        vxC *= sqrt(0.99 / v2C);
        vyC *= sqrt(0.99 / v2C);
        vzC *= sqrt(0.99 / v2C);
        v2C = 0.99;
       }
       double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);
       double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
       double dsigma[4];
       if (cartesian) {
        for(int ii=0; ii<4; ii++)
         dsigma[ii] = localCornelius.get_normal_elem(isegm, ii);
       }
       else {
        const double tauC = tau + localCornelius.get_centroid_elem(isegm, 0);
        ch = cosh(etaC);
        sh = sinh(etaC);
        dsigma[0] = tauC * (ch * localCornelius.get_normal_elem(isegm, 0) -
                           sh / tauC * localCornelius.get_normal_elem(isegm, 3));
        dsigma[3] = tauC * (-sh * localCornelius.get_normal_elem(isegm, 0) +
                           ch / tauC * localCornelius.get_normal_elem(isegm, 3));
        dsigma[1] = tauC * localCornelius.get_normal_elem(isegm, 1);
        dsigma[2] = tauC * localCornelius.get_normal_elem(isegm, 2);
       }
       double dVEff = 0.0, dsigma2 = 0.0;
       for (int ii = 0; ii < 4; ii++) {
        dVEff += dsigma[ii] * uC[ii];
        dsigma2 += gmumu[ii] * dsigma[ii] *  dsigma[ii];
       }
       const double dEtotSurf = (eC + pC) * uC[0] * dVEff - pC * dsigma[0];
       contribution.dVEff = dVEff;
       contribution.dEtotSurf = dEtotSurf;
       contribution.nbDVEff = nbC * dVEff;
       contribution.flowingOut = (dsigma2>0. and dsigma[0]>0.) or (dsigma2<0. and dEtotSurf>0.);
       if(contribution.flowingOut) {
        output::appendFreezeoutField(lineResult.freezeoutOutput, tau + localCornelius.get_centroid_elem(isegm, 0));
        output::appendFreezeoutField(lineResult.freezeoutOutput, getX(ix) + localCornelius.get_centroid_elem(isegm, 1));
        output::appendFreezeoutField(lineResult.freezeoutOutput, getY(iy) + localCornelius.get_centroid_elem(isegm, 2));
        output::appendFreezeoutField(lineResult.freezeoutOutput, getZ(iz) + localCornelius.get_centroid_elem(isegm, 3));
        for (int ii = 0; ii < 4; ii++) output::appendFreezeoutField(lineResult.freezeoutOutput, dsigma[ii]);
        for (int ii = 0; ii < 4; ii++) output::appendFreezeoutField(lineResult.freezeoutOutput, uC[ii]);
        output::appendFreezeoutField(lineResult.freezeoutOutput, TC);
        output::appendFreezeoutField(lineResult.freezeoutOutput, mubC);
        output::appendFreezeoutField(lineResult.freezeoutOutput, muqC);
        output::appendFreezeoutField(lineResult.freezeoutOutput, musC);
        if (vorticityOn) {
         for (int ii = 0; ii < 4; ii++) output::appendFreezeoutField(lineResult.betaOutput, dsigma[ii]);
         for (int ii = 0; ii < 4; ii++) output::appendFreezeoutField(lineResult.betaOutput, uC[ii]);
         output::appendFreezeoutField(lineResult.betaOutput, TC);
         output::appendFreezeoutField(lineResult.betaOutput, mubC);
         output::appendFreezeoutField(lineResult.betaOutput, muqC);
         output::appendFreezeoutField(lineResult.betaOutput, musC);
        }
#ifdef OUTPI
        double picart[10] = {0.};
        if (trcoeff->isViscous()) {
         if (cartesian) {
          for(int ii=0; ii<4; ii++)
           for(int jj=0; jj<=ii; jj++)
             picart[index44(ii, jj)] = piC[index44(ii, jj)];
         }
         else {
          picart[index44(0, 0)] = ch * ch * piC[index44(0, 0)] +
                                           2. * ch * sh * piC[index44(0, 3)] +
                                           sh * sh * piC[index44(3, 3)];
          picart[index44(0, 1)] =
             ch * piC[index44(0, 1)] + sh * piC[index44(3, 1)];
          picart[index44(0, 2)] =
             ch * piC[index44(0, 2)] + sh * piC[index44(3, 2)];
          picart[index44(0, 3)] =
             ch * sh * (piC[index44(0, 0)] + piC[index44(3, 3)]) +
             (ch * ch + sh * sh) * piC[index44(0, 3)];
          picart[index44(1, 1)] = piC[index44(1, 1)];
          picart[index44(1, 2)] = piC[index44(1, 2)];
          picart[index44(1, 3)] =
             sh * piC[index44(0, 1)] + ch * piC[index44(3, 1)];
          picart[index44(2, 2)] = piC[index44(2, 2)];
          picart[index44(2, 3)] =
             sh * piC[index44(0, 2)] + ch * piC[index44(3, 2)];
          picart[index44(3, 3)] = sh * sh * piC[index44(0, 0)] +
                                           ch * ch * piC[index44(3, 3)] +
                                           2. * sh * ch * piC[index44(0, 3)];
         }
        }
        for (int ii = 0; ii < 10; ii++) output::appendFreezeoutField(lineResult.freezeoutOutput, picart[ii]);
        output::appendFreezeoutField(lineResult.freezeoutOutput, PiC);
        if (extendFO) {
         output::appendFreezeoutField(lineResult.freezeoutOutput, eC);
         output::appendFreezeoutField(lineResult.freezeoutOutput, nbC);
         lineResult.freezeoutOutput.push_back('\n');
        }
        else {
         lineResult.freezeoutOutput.push_back('\n');
        }
#else
        output::appendFreezeoutField(lineResult.freezeoutOutput, dVEff);
        lineResult.freezeoutOutput.push_back('\n');
#endif
        std::unique_ptr<Matrix2D> jacobian = vorticityOn
         ? std::make_unique<Matrix2D>(Matrix2D{
            {1., 0., 0., 0.},
            {0., 1., 0., 0.},
            {0., 0., 1., 0.},
            {0., 0., 0., 1.}})
         : nullptr;
        if (!cartesian) {
         jacobian = vorticityOn
         ? std::make_unique<Matrix2D>(Matrix2D{
            {ch, 0., 0., -sh},
            {0., 1., 0., 0.},
            {0., 0., 1., 0.},
            {-sh, 0., 0., ch}})
         : nullptr;
        }
        std::unique_ptr<Matrix2D> dbetaCartesian = vorticityOn
          ? std::make_unique<Matrix2D>(Matrix2D(4, std::vector<double>(4, 0.0)))
          : nullptr;
        if(vorticityOn) {
         if (!dbetaCartesian || !jacobian) {
           std::runtime_error("dbetaCartesian and/or jacobian is a nullptr");
         }
         for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
           for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
             (*dbetaCartesian)[i][j] += (*jacobian)[i][k] * (*jacobian)[j][l]
                                      * (*dbetaInterpolated)[k][l] * gmumu[l];
            }
           }
          }
         }
         for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
           output::appendFreezeoutField(lineResult.betaOutput, (*dbetaCartesian)[i][j]);
          }
         }
         output::appendFreezeoutField(lineResult.betaOutput, eC);
         lineResult.betaOutput.push_back('\n');
        }
#ifdef OUTPI
        double dEsurfVisc = 0.;
        for (int i = 0; i < 4; i++)
         dEsurfVisc += picart[index44(0, i)] * dsigma[i];
        if (writeDiagnostics) contribution.dEtotSurfVisc = dEtotSurf + dEsurfVisc;
#endif
       }
       lineResult.nelements++;
       if (writeDiagnostics) lineResult.contributions.push_back(std::move(contribution));
      }
     }
    }
    surfaceLines[surfaceLineIndex(ix, iy)] = std::move(lineResult);
   }
 }
 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++) {
   SurfaceLineResult &lineResult = surfaceLines[surfaceLineIndex(ix, iy)];
   if (writeDiagnostics) {
    for (const SurfaceElementContribution &contribution : lineResult.contributions) {
     nelements++;
     if (!contribution.warningOutput.empty()) cout << contribution.warningOutput;
     if (contribution.suspicious) nsusp++;
     vEff += contribution.dVEff;
     EtotSurf += contribution.dEtotSurf;
     nbTotSurf += contribution.nbDVEff;
     if(contribution.flowingOut) {
      EtotSurfPos += contribution.dEtotSurf;
      nbTotSurfPos += contribution.nbDVEff;
      EtotSurfPosVisc += contribution.dEtotSurfVisc;
     } else {
      EtotSurfNeg += contribution.dEtotSurf;
      nbTotSurfNeg += contribution.nbDVEff;
     }
    }
   } else {
   nelements += lineResult.nelements;
   if (!lineResult.warningOutput.empty()) cout << lineResult.warningOutput;
  }
  freezeoutOutput.append(lineResult.freezeoutOutput);
  if (!lineResult.betaOutput.empty()) {
   output::fbeta.write(lineResult.betaOutput.data(), lineResult.betaOutput.size());
  }
 }
 output::writeFreezeout(std::move(freezeoutOutput));
 if (writeDiagnostics) {
 E = E * dx * dy * dz;
 Ecore = Ecore * dx * dy * dz;
 Efull = Efull * dx * dy * dz;
 S = S * dx * dy * dz;
 Nb1 *= dx * dy * dz;
 Nb2 *= dx * dy * dz;
 Nbcore *= dx * dy * dz;
 eps_p = txxyy_num / txxyy_den ;
 output::faniz << setw(12) << tau << setw(14) << vt_num / vt_den << setw(14)
           << vxvy_num / vxvy_den << setw(14) << pi0x_num / pi0x_den << '\n';
// cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13)
//      << nbTotSurf << setw(13) << S << setw(13) << EtotSurf
//      << setw(10) << nelements << setw(10) << nsusp
//      << setw(13) << (float)(nCoreCutCells) / (float)(nCoreCells) << endl;
 // alternative print-out:
 cout << setw(10) << tau << setw(13) << Efull << setw(13) << Nb1 << setw(13)
      << EtotSurf<< setw(13) << EtotSurfPos
      << setw(13) << EtotSurfNeg << setw(13) << nbTotSurf
      << setw(13) << nbTotSurfPos << endl; // trim output line << setw(13) << nbTotSurfNeg << endl;
 cout << setw(10) << "core " << setw(13) << tau << setw(13) << Ecore
      << setw(13) << Nbcore << setw(13) << Nb2 << endl;
 }
#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 return nelements;
}

void Fluid::outputCorona(double tau, bool extendFO) {
 static double nbSurf = 0.0;
 double e, p, nb, nq, ns, T, mub, muq, mus, vx, vy, vz, Q[7];
 double E = 0., Efull = 0., S = 0., Px = 0., vt_num = 0., vt_den = 0.,
        vxvy_num = 0., vxvy_den = 0., pi0x_num = 0., pi0x_den = 0.,
        txxyy_num = 0., txxyy_den = 0., Nb1 = 0., Nb2 = 0.;
 double eta = 0;
 double ch {0.}, sh {0.};
 double etaC {0.};
 int nelements = 0;
 const double gmumu[4] = {1., -1., -1., -1.};
 const bool writeDiagnostics = !freezeoutOnlyOutput;
 std::string coronaFreezeoutOutput;
 coronaFreezeoutOutput.reserve(1024 * 1024);

#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 if (writeDiagnostics) output::f2d << '\n';
 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = 2; iz < nz - 2; iz++) {
    Cell *c = getCell(ix, iy, iz);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    c->getQ(Q);
    eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
    double s = eos->s(e, nb, nq, ns);
    eta = getZ(iz);
    if (cartesian) {
     E += (e + p) / (1. - vx*vx - vy*vy - vz*vz) - p;
     Nb1 += Q[NB_];
     Nb2 += nb / sqrt(1. - vx*vx - vy*vy - vz*vz);
     //---- inf check
     if (std::isinf(E)) {
      cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
          << setw(14) << vy << setw(14) << vz << endl;
      exit(1);
     }
     //--------------
     Efull += (e + p) / (1. - vx*vx - vy*vy - vz*vz) - p;
     if (trcoeff->isViscous())
      Efull += c->getpi(0, 0);
     // -- noneq. corrections to entropy flux
     double deltas = 0.;
     if (trcoeff->isViscous())
      for (int i = 0; i < 4; i++)
       for (int j = 0; j < 4; j++)
        deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
     if (T > 0.02) {
      s += 1.5 * deltas / ((e + p) * T);
      S += s / sqrt(1. - vx*vx - vy*vy - vz*vz);
     }
    }
    else {
     const double cosh_int = (sinh(eta + 0.5 * dz) - sinh(eta - 0.5 * dz)) / dz;
     const double sinh_int = (cosh(eta + 0.5 * dz) - cosh(eta - 0.5 * dz)) / dz;
     E += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
             (cosh_int - tanh(vz) * sinh_int) -
         tau * p * cosh_int;
     Nb1 += Q[NB_];
     Nb2 += tau * nb * (cosh_int - tanh(vz) * sinh_int) /
           sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
     //---- inf check
     if (std::isinf(E)) {
      cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
          << setw(14) << vy << setw(14) << vz << endl;
      exit(1);
     }
     //--------------
     Efull += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (cosh(eta) - tanh(vz) * sinh(eta)) -
             tau * p * cosh(eta);
     if (trcoeff->isViscous())
      Efull +=
         tau * c->getpi(0, 0) * cosh(eta) + tau * c->getpi(0, 3) * sinh(eta);
     // -- noneq. corrections to entropy flux
     double deltas = 0.;
     if (trcoeff->isViscous())
      for (int i = 0; i < 4; i++)
       for (int j = 0; j < 4; j++)
        deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
     if (T > 0.02) {
      s += 1.5 * deltas / ((e + p) * T);
      S += tau * s * (cosh_int - tanh(vz) * sinh_int) /
          sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
     }
     Px += tau * (e + p) * vx / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
     vt_num += e / sqrt(1. - vx * vx - vy * vy) * sqrt(vx * vx + vy * vy);
     vt_den += e / sqrt(1. - vx * vx - vy * vy);
     vxvy_num += e * (fabs(vx) - fabs(vy));
     vxvy_den += e;
     txxyy_num += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (vx * vx - vy * vy);
     txxyy_den += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (vx * vx + vy * vy) +
                 2. * p;
     pi0x_num += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                fabs(c->getpi(0, 1));
     pi0x_den += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)); 
    }
    
    //----- Cornelius stuff
    bool isCorona = true, isTail = true;
    double QCube[2][2][2][7];
    double piSquare[2][2][2][10], PiSquare[2][2][2];
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       double _p, _nb, _nq, _ns, _vx, _vy, _vz;
       Cell *cc = getCell(ix + jx, iy + jy, iz + jz);
       if (cartesian) {
        cc->getPrimVar(eos, 1.0, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
       }
       else {
        cc->getPrimVar(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
       }
       cc->getQ(QCube[jx][jy][jz]);
       if (e > ecrit) isCorona = false;
       // ---- get viscous tensor
       for (int ii = 0; ii < 4; ii++)
        for (int jj = 0; jj <= ii; jj++)
         piSquare[jx][jy][jz][index44(ii, jj)] = cc->getpi(ii, jj);
       PiSquare[jx][jy][jz] = cc->getPi();
      }

    // ---- interpolation procedure
    double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
           musC = 0., piC[10]={0.}, PiC = 0., nbC = 0.,
           nqC = 0.;  // values at the centre, to be interpolated
    double QC[7] = {0.};
    double eC = 0., pC = 0.;
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++)
       for (int i = 0; i < 7; i++) {
        QC[i] += QCube[jx][jy][jz][i] * 0.125;
       }
    if (!cartesian) 
     for (int i = 0; i < 7; i++) QC[i] = QC[i] / tau;
    double _ns = 0.0;
    transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
    if (eC >= 0.01) isTail = false;

    if (isCorona && !isTail) {
      nelements++;
      output::appendFreezeoutField(coronaFreezeoutOutput, tau);
      output::appendFreezeoutField(coronaFreezeoutOutput, getX(ix) + 0.5 * dx);
      output::appendFreezeoutField(coronaFreezeoutOutput, getY(iy) + 0.5 * dy);
      output::appendFreezeoutField(coronaFreezeoutOutput, getZ(iz) + 0.5 * dz);
     eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
     if (TC > 0.4 || fabs(mubC) > 0.99) {
      cout << "#### Error (surface): high T/mu_b ####\n";
     }
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++) {
        for (int ii = 0; ii < 10; ii++)
         piC[ii] += piSquare[jx][jy][jz][ii] * 0.125;
        PiC += PiSquare[jx][jy][jz] * 0.125;
       }
     if (!cartesian) {
      etaC = getZ(iz) + 0.5 * dz;
      transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
     }

     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     
     double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);
     double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
     double dsigma[4];
     
     if (cartesian) {
      dsigma[0] = dx * dy * dz;
      dsigma[1] = dsigma[2] = dsigma[3] = 0.0; 
     }
     else {
      const double tauC = tau;
      
      // ---- transform dsigma to lab.frame :
      ch = cosh(etaC);
      sh = sinh(etaC);
      dsigma[0] = tauC * (ch * dx * dy * dz);
      dsigma[3] = tauC * (-sh * dx * dy * dz);
      dsigma[1] = 0.0;
      dsigma[2] = 0.0;
     }
     
     double dVEff = 0.0;
     for (int ii = 0; ii < 4; ii++)
      dVEff += dsigma[ii] * uC[ii];  // normalize for Delta eta=1
     vEff += dVEff;
     for (int ii = 0; ii < 4; ii++) output::appendFreezeoutField(coronaFreezeoutOutput, dsigma[ii]);
     for (int ii = 0; ii < 4; ii++) output::appendFreezeoutField(coronaFreezeoutOutput, uC[ii]);
     output::appendFreezeoutField(coronaFreezeoutOutput, TC);
     output::appendFreezeoutField(coronaFreezeoutOutput, mubC);
     output::appendFreezeoutField(coronaFreezeoutOutput, muqC);
     output::appendFreezeoutField(coronaFreezeoutOutput, musC);
     #ifdef OUTPI
     double picart[10] = {0.};
     if (trcoeff->isViscous()) {
      if (cartesian) {
       for(int ii=0; ii<4; ii++)
        for(int jj=0; jj<=ii; jj++)
         picart[index44(ii, jj)] = piC[index44(ii, jj)];
      }
      else {
       /*pi00*/ picart[index44(0, 0)] = ch * ch * piC[index44(0, 0)] +
                                        2. * ch * sh * piC[index44(0, 3)] +
                                        sh * sh * piC[index44(3, 3)];
       /*pi01*/ picart[index44(0, 1)] =
           ch * piC[index44(0, 1)] + sh * piC[index44(3, 1)];
       /*pi02*/ picart[index44(0, 2)] =
           ch * piC[index44(0, 2)] + sh * piC[index44(3, 2)];
       /*pi03*/ picart[index44(0, 3)] =
           ch * sh * (piC[index44(0, 0)] + piC[index44(3, 3)]) +
           (ch * ch + sh * sh) * piC[index44(0, 3)];
       /*pi11*/ picart[index44(1, 1)] = piC[index44(1, 1)];
       /*pi12*/ picart[index44(1, 2)] = piC[index44(1, 2)];
       /*pi13*/ picart[index44(1, 3)] =
           sh * piC[index44(0, 1)] + ch * piC[index44(3, 1)];
       /*pi22*/ picart[index44(2, 2)] = piC[index44(2, 2)];
       /*pi23*/ picart[index44(2, 3)] =
           sh * piC[index44(0, 2)] + ch * piC[index44(3, 2)];
       /*pi33*/ picart[index44(3, 3)] = sh * sh * piC[index44(0, 0)] +
                                        ch * ch * piC[index44(3, 3)] +
                                        2. * sh * ch * piC[index44(0, 3)];
     }
    }
     for (int ii = 0; ii < 10; ii++) output::appendFreezeoutField(coronaFreezeoutOutput, picart[ii]);
     output::appendFreezeoutField(coronaFreezeoutOutput, PiC);
     if (extendFO) {
	      output::appendFreezeoutField(coronaFreezeoutOutput, eC);
	      output::appendFreezeoutField(coronaFreezeoutOutput, nbC);
	      coronaFreezeoutOutput.push_back('\n');
     }
     else {
	      coronaFreezeoutOutput.push_back('\n');
     }
#else
	     output::appendFreezeoutField(coronaFreezeoutOutput, dVEff);
	     coronaFreezeoutOutput.push_back('\n');
#endif
     double dEsurfVisc = 0.;
     for (int i = 0; i < 4; i++)
      dEsurfVisc += picart[index44(0, i)] * dsigma[i];
     EtotSurf += (eC + pC) * uC[0] * dVEff - pC * dsigma[0] + dEsurfVisc;
     nbSurf += nbC * dVEff;
    }
    //----- end Cornelius
   }
 output::writeFreezeout(std::move(coronaFreezeoutOutput));
 // Set number of corona cells in Fluid. This is needed for the
 // header of beta.dat in case vorticity is used
 num_corona_cells = nelements;

 E = E * dx * dy * dz;
 Efull = Efull * dx * dy * dz;
 S = S * dx * dy * dz;
 Nb1 *= dx * dy * dz;
 Nb2 *= dx * dy * dz;
	 if (writeDiagnostics) {
	  output::faniz << setw(12) << tau << setw(14) << vt_num / vt_den << setw(14)
	            << vxvy_num / vxvy_den << setw(14) << pi0x_num / pi0x_den << '\n';
	 }
 cout << setw(10) << "tau" << setw(13) << "E" << setw(13) << "Efull" << setw(13)
      << "Nb" << setw(13) << "Sfull" << setw(13) << "EtotSurf" << setw(13) << "elements" << setw(10)
      << "susp." << setw(13) << "\%cut" << endl;
 cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13)
      << nbSurf << setw(13) << S << setw(13) << EtotSurf << endl;
#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 cout << "corona elements : " << nelements << endl;
}


void Fluid::InitialAnisotropies(double tau0) {
 double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;

 double xcm = 0.0, xcm_nom = 0.0, xcm_denom = 0.0 ;
 double ycm = 0.0, ycm_nom = 0.0, ycm_denom = 0.0 ;

 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   for (int iz = 0; iz < nz; iz++) {
    if(fabs(getZ(iz)) < 0.5) {
     Cell* c = getCell(ix, iy, iz);
     double x = getX(ix) ;
     double y = getY(iy) ;
     getCMFvariables(c, tau0, e, nb, nq, ns, vx, vy, vz);
     xcm_nom += x * e ;
     xcm_denom += e ;
     ycm_nom += y * e ;
     ycm_denom += e ;
    }
   }
  }
 }

 xcm = xcm_nom / xcm_denom ;
 ycm = ycm_nom / ycm_denom ;

 double eps2 = 0.0, eps2_nom_real = 0.0, eps2_nom_imag = 0.0, eps2_denom = 0.0 ;
 double eps3 = 0.0, eps3_nom_real = 0.0, eps3_nom_imag = 0.0, eps3_denom = 0.0 ;

 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   for (int iz = 0; iz < nz; iz++) {
    if(fabs(getZ(iz)) < 0.5) {
     Cell* c = getCell(ix, iy, iz);
     double x = getX(ix) ;
     double y = getY(iy) ;
     x = x - xcm ;
     y = y - ycm ;
     double r = sqrt(x*x + y*y) ;
     double phi = atan2(y, x) ;
     getCMFvariables(c, tau0, e, nb, nq, ns, vx, vy, vz);
     eps2_denom += pow(r, 2) * e ;
     eps2_nom_real += pow(r, 2) * cos(2 * phi) * e ;
     eps2_nom_imag += pow(r, 2) * sin(2 * phi) * e ;
     eps3_denom += pow(r, 3) * e ;
     eps3_nom_real += pow(r, 3) * cos(3 * phi) * e ;
     eps3_nom_imag += pow(r, 3) * sin(3 * phi) * e ;
    }
   }
  }
 }

 eps2 = sqrt(pow(eps2_nom_real, 2) + pow(eps2_nom_imag, 2)) / eps2_denom ;
 eps3 = sqrt(pow(eps3_nom_real, 2) + pow(eps3_nom_imag, 2)) / eps3_denom ;

 cout << "epsilon2 = " << eps2 << endl;
 cout << "epsilon3 = " << eps3 << endl;

 double xcm_eta[nz], xcm_nom_eta[nz], xcm_denom_eta[nz];
 double ycm_eta[nz], ycm_nom_eta[nz], ycm_denom_eta[nz];

 for (int iz = 0; iz < nz; iz++) {
  xcm_eta[iz] = 0;
  xcm_nom_eta[iz] = 0;
  xcm_denom_eta[iz] = 0;
  ycm_eta[iz] = 0;
  ycm_nom_eta[iz] = 0;
  ycm_denom_eta[iz] = 0;
  for (int ix = 0; ix < nx; ix++) {
   for (int iy = 0; iy < ny; iy++) {
    Cell* c = getCell(ix, iy, iz);
    double x = getX(ix) ;
    double y = getY(iy) ;
    getCMFvariables(c, tau0, e, nb, nq, ns, vx, vy, vz);
    xcm_nom_eta[iz] += x * e ;
    xcm_denom_eta[iz] += e ;
    ycm_nom_eta[iz] += y * e ;
    ycm_denom_eta[iz] += e ;
   }
  }
  xcm_eta[iz] = xcm_nom_eta[iz] / xcm_denom_eta[iz];
  ycm_eta[iz] = ycm_nom_eta[iz] / ycm_denom_eta[iz];
 }

 for (int iz = 0; iz < nz; iz++) {
  eps2 = 0.0;
  eps2_nom_real = 0.0;
  eps2_nom_imag = 0.0;
  eps2_denom = 0.0;
  for (int ix = 0; ix < nx; ix++) {
   for (int iy = 0; iy < ny; iy++) {
     Cell* c = getCell(ix, iy, iz);
     double x = getX(ix) ;
     double y = getY(iy) ;
     x = x - xcm_eta[iz] ;
     y = y - ycm_eta[iz] ;
     double r = sqrt(x*x + y*y) ;
     double phi = atan2(y, x) ;
     getCMFvariables(c, tau0, e, nb, nq, ns, vx, vy, vz);
     eps2_denom += pow(r, 2) * e ;
     eps2_nom_real += pow(r, 2) * cos(2 * phi) * e ;
     eps2_nom_imag += pow(r, 2) * sin(2 * phi) * e ;
   }
  }
  eps2 = sqrt(pow(eps2_nom_real, 2) + pow(eps2_nom_imag, 2)) / eps2_denom ;
  double eta = getZ(iz) ;
  cout << eta << " " << eps2 << endl;
 }

 exit(1) ;
}

void Fluid::addParticle(Particle _particle) {
 double source[7] = {0.};
 double dv = dx * dy * dz;
 // where to smooth the particle out
 int ixc = _particle.getIxc();
 int smoothx = _particle.getNsmoothX();
 int iyc = _particle.getIyc();
 int smoothy = _particle.getNsmoothY();
 int izc = _particle.getIzc();
 int smoothz = _particle.getNsmoothZ();
 const double scale = _particle.getScale();
 
 for (int ix = ixc - smoothx; ix < ixc + smoothx + 1; ix++) 
  for (int iy = iyc - smoothy; iy < iyc + smoothy + 1; iy++) 
   for (int iz = izc - smoothz; iz < izc + smoothz + 1; iz++) 
     if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
      
      const double xdiff = _particle.getX() - (minx + ix * dx);
      const double ydiff = _particle.getY() - (miny + iy * dy);
      const double zdiff = _particle.getZ() - (minz + iz * dz);

      double weight = _particle.getWeight(xdiff, ydiff, zdiff);
      source[0] = _particle.getE() * weight * scale / dv;
      source[1] = _particle.getPx() * weight * scale / dv;
      source[2] = _particle.getPy() * weight * scale / dv;
      source[3] = _particle.getPz() * weight * scale / dv;
      source[4] = _particle.getB() * weight * scale / dv;
      source[5] = _particle.getQ() * weight * scale / dv;
      source[6] = _particle.getS() * weight * scale / dv;
      if(isnan(source[0]) or isinf(source[0])) {
        cout << "Fluid::addParticle: NaN/inf\n" ;
      }
      getCell(ix,iy,iz)->addParticleSource(source);
 }
}
