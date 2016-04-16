#ifndef HAVE_THOMAS_ALGORITHM
#define HAVE_THOMAS_ALGORITHM
#include <vector>
void thomas(std::vector<double> &A, std::vector<double> &sol, std::vector<double> &F, const unsigned int nnodes);		//A è la matrice da scomporre, sol è il vettore dove mettere la soluzione, F è il vettore dei termini noti, nnodes la dimensione della matrice (quadrata).

#endif
