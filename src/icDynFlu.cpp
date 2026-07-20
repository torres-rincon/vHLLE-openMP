#include <cstdlib>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cfloat>
#include <vector>
#include <deque>
#include <algorithm>  // For std::sort

#include "eos.h"
#include "eoChiral.h"
#include "rmn.h"
#include "fld.h"
#include "icPartSMASH.h"
#include "icDynFlu.h"
#include "colour.h"
#include "particle.h"

using namespace std;

IcDynFlu::IcDynFlu(Fluid* f, const char* filename, double gaussian_sigma,
                         deque<Particle>* particles) {
 // works only in CARTESIAN coordinates
 nx = f->getNX();
 ny = f->getNY();
 nz = f->getNZ();
 
 dx = f->getDx();
 dy = f->getDy();
 dz = f->getDz();

 R = sqrt(2)*gaussian_sigma;

 // helper variables for read-out
 double x1, x2, x3, x4, x5, x6, x7;
 double x8, x9;
 std::vector<Particle> all_particles;
 all_particles.clear();
 all_particles.reserve(10000);

 // initialize the grid of conserved quantities
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

 //cout << "I am in the class \n";                          
 // ---- read the events
 nevents = 0;
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 int np = 0;  // particle counter
 string line;
 istringstream instream;
 getline(fin, line);
 getline(fin, line);
 instream.str(line);
 instream.seekg(0);
 instream.clear();
 while (true) {
  getline(fin, line);
  if( fin.eof() ) break;
  instream.str(line);
  instream.seekg(0);
  instream.clear();
  if (line[0] == '#') {
    if (line.size() >= 5 && line.substr(line.size() - 5) == "start") {
      nevents++;
      np = 0;
    }
    continue;
  }
  // Read line
  if (!instream.fail()) {
  instream >> T_val >> X_val >> Y_val >> Z_val >> M_val >> E_val >> Px_val >>
              Py_val >> Pz_val >> Id_val >> x1 >> Charge_val >> x2 >> x3 >>
              x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> Baryon_val >> Strangeness_val;
  
  if(abs(Id_val)>100) { // exclude photons, W and Z bosons and Higgs
    Particle particleIn(f, R, Baryon_val, Charge_val, Strangeness_val,
    T_val, X_val, Y_val, Z_val, E_val, Px_val, Py_val, Pz_val, Id_val, nevents-1);
    all_particles.push_back(particleIn);
    np++;
    //cout << np << " " << particleIn.getT() << " " << particleIn.getE() << endl;
  }
  }
  else if (np > 0) {
   if (nevents % 100 == 0) {
    cout << "event = " << nevents << "  np = " << np << "\r";
    cout << flush;
   }
  }
 }
 // sort particles by time
 std::sort(all_particles.begin(), all_particles.end(),
   [](const Particle &a, const Particle &b) -> bool { return a.getT() < b.getT(); });
 cout << "IcPartSMASH: earliest particle at t = " << all_particles[0].getT() << endl;
 cout << "IcPartSMASH:   latest particle at t = " << all_particles[all_particles.size()-1].getT() << endl;
 cout << "total number of particles read = " << all_particles.size() << endl;
 // copy the contents of the vector to the queue
 for (auto &particle : all_particles) {
   particle.setScale(1./nevents);
   particles->push_back(particle);
 }
 if (nevents > 1)
  cout << "++ Warning: loaded " << nevents << "  initial SMASH events\n";
 }

IcDynFlu::~IcDynFlu() {
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

void IcDynFlu::setIC(Fluid* f, EoS* eos, deque<Particle>* particles, double &timeInit,
     int min_particles_FO, double &timeInitFO) {
  // works only in CARTESIAN coordinates
  double Q[7], e, p, nb, nq, ns, vx, vy, vz;
  double weight;

  timeInit = particles->front().getT();
  timeInitFO = particles->at(min_particles_FO*nevents).getT();
  
  cout << "Hydro starting at " << timeInit << endl;
  cout << "Freezeout starting at " << timeInitFO << endl;

  double t = timeInit;
  // pick particles that left SMASH the earliest
  while (t <= timeInit && particles->size()>0) {
    Particle particleToSmooth = particles->front();
    int ixc = particleToSmooth.getIxc();
    int smoothx = particleToSmooth.getNsmoothX();
    int iyc = particleToSmooth.getIyc();
    int smoothy = particleToSmooth.getNsmoothY();
    int izc = particleToSmooth.getIzc();
    int smoothz = particleToSmooth.getNsmoothZ();
    const double scale = particleToSmooth.getScale();
      
    for (int ix = ixc - smoothx; ix < ixc + smoothx + 1; ix++) 
     for (int iy = iyc - smoothy; iy < iyc + smoothy + 1; iy++) 
      for (int iz = izc - smoothz; iz < izc + smoothz + 1; iz++) 
      if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
        
        weight = particleToSmooth.getWeight(ix, iy, iz, f, R);
        
        T00[ix][iy][iz] += particleToSmooth.getE() * weight * scale;
        T0x[ix][iy][iz] += particleToSmooth.getPx() * weight * scale;
        T0y[ix][iy][iz] += particleToSmooth.getPy() * weight * scale;
        T0z[ix][iy][iz] += particleToSmooth.getPz() * weight * scale;
        QB[ix][iy][iz] += particleToSmooth.getB() * weight * scale;
        QE[ix][iy][iz] += particleToSmooth.getQ() * weight * scale;
        QS[ix][iy][iz] += particleToSmooth.getS() * weight * scale;
  }
  
  particles->pop_front();
  t = particles->front().getT(); 
  } // all particles for IC should be in now

  // initialize cell values
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++) {
        Q[T_] = T00[ix][iy][iz] / dx / dy / dz;
        Q[X_] = T0x[ix][iy][iz] / dx / dy / dz;
        Q[Y_] = T0y[ix][iy][iz] / dx / dy / dz;
        Q[Z_] = T0z[ix][iy][iz] / dx / dy / dz;
        Q[NB_] = QB[ix][iy][iz] / dx / dy / dz;
        Q[NQ_] = QE[ix][iy][iz] / dx / dy / dz;
        Q[NS_] = QS[ix][iy][iz] / dx / dy / dz;
        transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz);

        Cell* c = f->getCell(ix, iy, iz);
        c->setPrimVar(eos, 1.0, e, nb, nq, ns, vx, vy, vz);
        c->getQ(Q);
        if (e > 0.) 
          c->setAllM(1.);
      }
 // fluid initialized + particles queue awaits
}

// dynamical IC: set up IC with 1st particles, others stay in queue
void outputCoronaParticles(std::deque<Particle>* particles, std::string outputDir) { 
 // works only in CARTESIAN coordinates

 // sort particles by event number
  std::vector<Particle> tempVec;
    while (!particles->empty()) {
        tempVec.push_back(particles->front());
        particles->pop_front();
    }
  std::sort(tempVec.begin(), tempVec.end(),
   [](const Particle &a, const Particle &b) -> bool { return a.getEventNo() < b.getEventNo(); }); 
  for (const auto& item : tempVec) {
        particles->push_back(item);
    }

 
  int particle_number = particles->size();
  int n_part = 0;
  int n_event = 0;
    
    if (particle_number > 0) {
      string filename = outputDir.c_str();
      filename.append("/particle_lists.oscar");
      constexpr std::size_t particleOutputBufferSize = 1 << 20;
      std::array<char, particleOutputBufferSize> outputBuffer;
      ofstream outfile;
      outfile.rdbuf()->pubsetbuf(outputBuffer.data(), outputBuffer.size());
      outfile.open(filename.c_str());
      outfile << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge \n";
      outfile << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e \n";
      outfile << "# event " <<  n_event << " out\n";
      while (particle_number > 0)
      {
        Particle particle_to_dump = particles->front();
        double t = particle_to_dump.getT();
        double x = particle_to_dump.getX();
        double y = particle_to_dump.getY();
        double z = particle_to_dump.getZ();
        double mass = particle_to_dump.getM();
        double p0 = particle_to_dump.getE();
        double px = particle_to_dump.getPx();
        double py = particle_to_dump.getPy();
        double pz = particle_to_dump.getPz();
        int pdg = particle_to_dump.getPdg();
        int id = 0;
        int charge = particle_to_dump.getQ();
        while (particle_to_dump.getEventNo() > n_event)
        {
          outfile << "# event " <<  n_event << " end\n";
          n_event++;
          outfile << "# event " <<  n_event << " out\n"; 
        }
        outfile << t << " " << x << " " << y << " " << z << " " << mass << " " << 
                    p0 << " " << px << " " << py << " " << pz << " " << pdg << 
                    " " << id << " " << charge << "\n";
        n_part += 1;
        particles->pop_front();
        particle_number = particles->size();
      }
      outfile << "# event " <<  n_event << " end\n";
    }
    // check-in about the procedure
    std::cout << n_part << " particles were left after hydro evolution.\n" <<
                 "At the end of the run, particle queue is empty: " <<
                 particles->empty() << "\n";
}
