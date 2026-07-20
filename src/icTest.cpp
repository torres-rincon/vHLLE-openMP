#include <fstream>
#include <iomanip>
#include <iostream>

#include "fld.h"
#include "eos.h"
#include "icTest.h"
#include "inc.h"
#include "s95p.h"

using namespace std;

ICTest::ICTest() {}

ICTest::~ICTest(void) {}

// This is to provide standard testing scenarios
// for the cartesian version of vHLLE:
// 1D shock-tube test - numerical against analytic solution,
// radial expansion - check if solution is symmetric


void ICTest::setIC(Fluid *f, EoS *eos, int test_option) {
 double e, nb, nq, ns, vx, vy, vz;
 Cell *c;
 double r, x, y, z;
 double e_high, e_low;
 
 // Test 1: 1D Shock-tube test
 if (test_option == 1) {
   
   e_high = 16.0;
   e_low = 1.0;
   double x_half = f->getX(f->getNX()) / 2.; 
   std::cout << "x_half = " << x_half << std::endl;

   for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
     for (int iz = 0; iz < f->getNZ(); iz++) {
      
      c = f->getCell(ix, iy, iz);
      x = f->getX(ix);
      if (x < x_half) e = e_high; 
      else e = e_low; 
      nb = nq = ns = 0.0;
      vx = vy = vz = 0.0;
    
      // cartesian coordinates: tau = 1.0
      c->setPrimVar(eos, 1.0, e, nb, nq, ns, vx, vy, vz);
      
      if (e > 0.) c->setAllM(1.); 
     }
   }



 // Test 2: Radial expansion
 if (test_option == 2) {

  double R = 2.0;  // sphere radius
  e_high = 8.0;
  e_low = 1.0;
    
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++)
    for (int iz = 0; iz < f->getNZ(); iz++) {
    c = f->getCell(ix, iy, iz);
    x = f->getX(ix);
    y = f->getY(iy);
    z = f->getZ(iz);
    r = x*x + y*y + z*z;
    if (r < R*R) e = e_high; 
    else e = e_low;
    nb = nq = ns = 0.0;
    vx = vy = vz = 0.0;
    
    // cartesian coordinates: tau = 1.0
    c->setPrimVar(eos, 1.0, e, nb, nq, ns, vx, vy, vz);
    //c->saveQprev();
    
    if (e > 0.) c->setAllM(1.);     
  }
 } 

}





















