#include"thomas.hpp"

void thomas(std::vector<double> &A, std::vector<double> &sol, std::vector<double> &F, const unsigned int nnodes){

double* alfa=new double [nnodes];
double* beta=new double [nnodes-1];
double* c=new double [nnodes-1];	

for(unsigned int i=0; i<nnodes-1; i++){	
	c[i]=A[i*nnodes+i+1];
}
alfa[0]=A[0];
for(unsigned int i=0; i<nnodes-1; i++){
	beta[i]=A[(i+1)*nnodes+i]/alfa[i];
	alfa[i+1]=A[(i+1)*nnodes+(i+1)]-beta[i]*c[i];
}
double* y=new double [nnodes];

y[0]=F[0];
for(unsigned int i=1; i<nnodes; i++)
	y[i]=F[i]-beta[i-1]*y[i-1];

sol[nnodes-1]=y[nnodes-1]/alfa[nnodes-1];
for(int i=nnodes-2; i>=0; i--)
	sol[i]=(y[i]-c[i]*sol[i+1])/alfa[i];

delete [] alfa;
delete [] beta;
delete [] c;
delete [] y;
}
