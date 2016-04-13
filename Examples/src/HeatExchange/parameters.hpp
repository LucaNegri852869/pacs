#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
struct parameters
{
  //! max number of iteration for Gauss-Siedel
  int   itermax;
  //! Tolerance for stopping criterion
  double  toler;
  //! Bar length
   double L;
  //! First longitudinal dimension
  double a1;
 //! Second longitudinal dimension
  double a2;
  //! Dirichlet condition
  double To;
  //! External temperature 
  double Te;
  //! Conductivity
  double k;
  //! Convection coefficient
  double hc;
  //! Number of elements
  int M;
  //! Filename of the output file
  std::string out_file;
  //! Parameter that specify the norm to use
  int norm;
  //! Constructor takes default values
  parameters():
    itermax(1000000),
    toler(1.e-8),
    L(40.),
    a1(4.),
    a2(50.),
    To(46.),
    Te(20.),
    k(0.164),
    hc(1.e-6*200.),
    M(100),
    out_file("result.dat"),
    norm(0)
  {}
};
//! Prints parameters
std::ostream & operator << (std::ostream &,const parameters &);
#endif
