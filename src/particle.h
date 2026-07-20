class Fluid;

class Particle {
    private:
     int B, Q, S;            // baryon number, charge, strangeness
     double t, x, y, z;      // particle position
     double e, px, py, pz;   // particle momentum
     int pdg;                // pdg code
     int ixc, iyc, izc;      // cell coordinates of the particle
     int eventNo;            // event number: important for averaged IC
     double R;               // sigma of the gaussian
     int nsmoothx, nsmoothy,
      nsmoothz;              // smoothly distribute to +- this many cells
     double gauss_norm;      // normalization of the smearing kernel
     double scale;           // scaling factor for the case Nevents>1
     double calculateNorm(Fluid *f, double R);
                             // function to calculate the norm of the kernel
         
    public:
     Particle();
     Particle(Fluid *f, double _R, int _B, int _Q, int _S, double _t, double _x, 
      double _y, double _z, double _e, double _px, double _py, double _pz, int _pdg, int _eventNo);
     ~Particle(){};
     void setScale(double _scale) { scale = _scale; }
     double getScale() { return scale; }
     double getWeight(int ix, int iy, int iz, Fluid *f, double R);
     double getWeight(double _xdiff, double _ydiff, double _zdiff);

     inline int getIxc(void) { return ixc; }
     inline int getIyc(void) { return iyc; }
     inline int getIzc(void) { return izc; }

     inline int getNsmoothX(void) { return nsmoothx; }
     inline int getNsmoothY(void) { return nsmoothy; }
     inline int getNsmoothZ(void) { return nsmoothz; }

     inline int getB(void) { return B; }
     inline int getQ(void) { return Q; }
     inline int getS(void) { return S; }
     inline int getPdg(void) { return pdg; }

     inline double getT(void) const { return t; }
     inline double getX(void) const { return x; }
     inline double getY(void) const { return y; }
     inline double getZ(void) const { return z; }
     
     inline double getE(void) { return e; }
     inline double getPx(void) { return px; }
     inline double getPy(void) { return py; }
     inline double getPz(void) { return pz; }

     inline double getM(void) { return sqrt(e*e - px*px - py*py - pz*pz); }
     inline int getEventNo(void) const { return eventNo; }
};