
class EoS;
class Fluid;

// This is to provide standard testing scenarios
// for the cartesian version of vHLLE:
// 1D shock-tube test, radial expansion (check symmetry)

class ICTest {
public:
 ICTest(void);
 ~ICTest(void);
 // setIC: initializes entire hydro grid for a given test
 // option 1: shocktube test, option 2: radial expansion
 void setIC(Fluid *f, EoS *eos, int test_option);
};
