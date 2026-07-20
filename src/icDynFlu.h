#include <vector>
#include <deque>

class Fluid;
class EoS;
class Particle;

class IcDynFlu {
// works only in CARTESIAN coordinates
private:
 int nx, ny, nz, nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***T00, ***T0x, ***T0y, ***T0z, ***QB, ***QE, ***QS;
 // auxiliary particle values for reading from file
 double T_val, X_val, Y_val, Z_val, E_val, Px_val, Py_val, Pz_val, M_val;
 int Id_val, Baryon_val, Charge_val, Strangeness_val;
 // auxiliary particle arrays
 std::vector<double> T, X, Y, Z, E, Px, Py, Pz;
 std::vector<int> Id, Charge, Baryon, Strangeness;

 double R;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;

public:
 IcDynFlu(Fluid *f, const char *filename, double gaussian_sigma, std::deque<Particle>* particles);
 ~IcDynFlu();
 void setIC(Fluid *f, EoS *eos, std::deque<Particle>* particles, double &ctime, int min_particles_FO, double &ftime);
};

void outputCoronaParticles(std::deque<Particle>* particles, std::string outputDir);