#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  const double& toler=param.toler;   // Tolerance for stopping criterion
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  const int& norm=param.norm;
  std::string out_file(param.out_file); //Name of the result file
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);
  
  // Gauss Siedel is initialised with a linear variation
  // of T
  
  for(unsigned int m=0;m <= M;++m)
     theta[m]=(1.-m*h)*(To-Te)/Te;
  
  // Gauss-Seidel
  // epsilon=||x^{k+1}-x^{k}||
  // Stopping criteria epsilon<=toler
  
  int iter=0;
  double xnew, epsilon;
  switch(norm)
  {
  case 0:    //I use the Rn norm as stopping criterion
  {
  cout<<"Norm used: Rn"<<endl;
     do
       { epsilon=0.;

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   epsilon += (xnew-theta[m])*(xnew-theta[m]);
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 epsilon += (xnew-theta[M])*(xnew-theta[M]);
	 theta[M]=  xnew; 

	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );
      break;
      }
      
      case 1:    //I use norm L2 as stopping criterion 
      {
      cout<<"Norm used: L2"<<endl;
      double theta_back;
      do
       { epsilon=0.;
       	theta_back=theta[0];

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   epsilon += (((xnew-theta[m])*(xnew-theta[m]))+((theta[m-1]-theta_back)*(theta[m-1]-theta_back)))*0.5*h;  //I used the trapezoidal rule for computing the integral
	   theta_back=theta[m];
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 epsilon += ((xnew-theta[M])*(xnew-theta[M])+(theta[M-1]-theta_back)*(theta[M-1]-theta_back))*0.5*h;    //I used the trapezoidal rule for computing the integral
	 theta[M]=  xnew; 

	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );
       break;
       }
       
       
       case 2:    //I use norm H1 as stopping criterion
       {
       cout<<"Norm used: H1"<<endl;
       double theta_back;
       double epsilonL2;
       double epsilon_grad_L2;
       double xnew_forw;
       double theta_back2;
       double xn_der;
       double xo_der;
         do
       { epsilon=0.;
       epsilonL2=0.;
       epsilon_grad_L2=0.;
       theta_back=theta[0];
       theta_back2=theta[0];
       xnew_forw=(theta[0]+theta[2])/(2.+h*h*act);

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = xnew_forw;
	   m < (M-1) ? xnew_forw=(xnew+theta[m+2])/(2.+h*h*act) : xnew_forw=xnew;
	   epsilonL2 += (((xnew-theta[m])*(xnew-theta[m]))+((theta[m-1]-theta_back)*(theta[m-1]-theta_back)))*0.5*h;   //I used the trapezoidal rule for computing the integral
	   xn_der=((xnew_forw-theta[m-1])*0.5/h)-((theta[m+1]-theta_back)*0.5/h);
	   xo_der=((xnew-theta[m-2])*0.5/h)-((theta[m]-theta_back2)*0.5/h);
	   m > 1 ? epsilon_grad_L2 += (xn_der*xn_der+xo_der*xo_der)*0.5*h : epsilon_grad_L2 += (xn_der*xn_der +(((xnew-theta[0])*0.5/h)-((theta[m]-theta_back2)*0.5/h))*(((xnew-theta[0])*0.5/h)-((theta[m]-theta_back2)*0.5/h)))*0.5*h ; 	//I used the trapezoidal rule for computing the integral and the second order central approximation for the derivatives
	   theta_back2=theta_back;
	   theta_back=theta[m];
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 epsilonL2 += ((xnew-theta[M])*(xnew-theta[M])+(theta[M-1]-theta_back)*(theta[M-1]-theta_back))*0.5*h;     //I used the trapezoidal rule for computing the integral
	 epsilon_grad_L2 += ((((xnew-theta[M-1])/h)-(theta[M]-theta_back)/h)*(((xnew-theta[M-1])/h)-(theta[M]-theta_back)/h)+(((xnew-theta[M-2])*0.5/h)-((theta[M]-theta_back2)*0.5/h))*(((xnew-theta[M-2])*0.5/h)-((theta[M]-theta_back2)*0.5/h)))*0.5*h;	//I used the trapezoidal rule for computing the integral and both the second order central approximation and the first order left approximation for the derivatives	
	 theta[M]=  xnew; 

	 iter=iter+1;
	 epsilon=sqrt(epsilonL2)+sqrt(epsilon_grad_L2);     
       }while((epsilon > toler) && (iter < itermax) );
       break;
       }
    }
    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }
      

 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 

     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);

    cout<<"Result file: "<<out_file<<endl;
      ofstream f(out_file);
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
       }
     // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<std::endl;
     f.close();
     return status;
}
