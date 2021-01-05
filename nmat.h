#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <cstdio>



using namespace std;
using namespace arma;

//this file contains the functions needed for the Baer_Nunziato equation to work


class nmat{

  //This is a class for an N-dimensional Tensor/matrix class, with relevant subfunctions. 
  //element access is slightly faster than the field class, but slower than cube class.

  public:
    vector<int> dims; //vector containing the lengths of each dimension of the matrix
    vector<int> cumdims; //vector containing cumulative lengths, for fast indexing of the NMat
    int i_end; //number of dimensions

  

    vector <double> vectorize; //vector containing the nmat, can be requested by the user in a similar manner to the vectorize
                  //function in armadillo.

    void sizemat(vector<int> uilengths){
      //sets the size of the nmat, by inputting a vector: {N,M,P,Q,R....,Z} for a NxMxPxQxRx...xZ size tensor

        dims = uilengths;
        i_end = uilengths.size();
        cumdims.resize(i_end+1);
        //creating cumulative lengths for indexing.
        for(int i; i<i_end+1; i++){
          if(i==0){
            cumdims[0] = 1;
          }else{
            cumdims[i]=dims[i-1]*cumdims[i-1];
          }
        }
        vectorize.resize((cumdims[i_end]));
    };

    void valfill(double val){
      //fills the nmat with a particular value

      std::fill(vectorize.begin(), vectorize.end(), val);
    }

    void setval(vector<int> ind_vec,double val){
      //sets the value of a particular element within the tensor

      int vect_ind =0;
      for(int i=0; i<i_end;i++){
        vect_ind += ind_vec[i]*cumdims[i];
      }     
      vectorize[vect_ind] = val;
    };
      
    double getval(vector<int> ind_vec){
      //accesses a particular element within the tensor

      int vect_ind = 0;
      for(int i=0; i<i_end;i++){
        vect_ind += ind_vec[i]*cumdims[i];
      }
      return vectorize[vect_ind];
    };

    vec nslice(int dir,vector<int> i0_vect,int iend){
      //takes a slice of the matrix accoss dimension 'dir' starting at element i0_vect,
      // through to the element iend in dimension dir.

      int i0 = i0_vect[dir];
      vec sliced_vec = zeros<vec>(iend-i0+1); //vector class containing the slice
      int j = 0;
      for(int i = i0;i<=iend;i++){
        i0_vect[dir] = i;
        sliced_vec(j) = getval(i0_vect);
        j++;
      }
      return sliced_vec;
    };
};




cube fliplrcube(cube preflip){
  int N = preflip.n_cols;
  int M = preflip.n_rows;
  int L = preflip.n_slices;
  cube postflip = preflip;
  for(int i = 0; i<N; i++){
    postflip(span(0,M-1),span(i,i),span(0,L-1)) = preflip(span(0,M-1),span(N-i-1,N-i-1),span(0,L-1));
  }
  return postflip;
}

cube flipudcube(cube preflip){
  int N = preflip.n_cols;
  int M = preflip.n_rows;
  int L = preflip.n_slices;
  cube postflip = preflip;
  for(int i = 0; i<M; i++){
    postflip(span(i,i),span(0,N-1),span(0,L-1)) = preflip(span(M-i-1,M-i-1),span(0,N-1),span(0,L-1));
  }
  return postflip;
}
       
    

mat flux_eval(mat Q_0,int DIM){
  //matrix for evaluating fluxes, using an input of conserved variables
  //using DIM==1 gives the flux in the x direction and DIM==2 gives the flux in the y direction
  //each column vector within the matrix refers to a conserved variable. 

  mat F_Q0 = Q_0; //flux matrix
  F_Q0.zeros();
  double gamma = 1.4;
  //finding conserved variables
  rowvec rho = Q_0.row(0);
  rowvec u = Q_0.row(1)/rho;
  rowvec v = Q_0.row(2)/rho;
  rowvec p = (Q_0.row(3) - 0.5*rho%(u%u + v%v))*(gamma-1);
  //calculating fluxes (for ideal gas)
  if(DIM==1){
    F_Q0.row(0) = rho%u;
    F_Q0.row(1) = rho%u%u + p;
    F_Q0.row(2) = rho%v%u;
    F_Q0.row(3) = u%(Q_0.row(3)+p);

  }else{
    F_Q0.row(0) = rho%v;
    F_Q0.row(1) = rho%u%v;
    F_Q0.row(2) = rho%v%v + p;
    F_Q0.row(3) = v%(Q_0.row(3)+p);
  }

  return F_Q0;
}

mat revolve(mat A){
  mat B = A; //creating a copy
  double N = A.n_cols;

  for(int i =1; i<N; i++){
    A.col(i) = B.col(i-1);
  }
  A.col(0) = B.col(N-1);
  return A;
}

cx_vec revolve_vec(cx_vec A){
  cx_vec B = A;
  double N = A.n_elem;
  //cout<<N<<endl;
    for(int i =1; i<N; i++){
      A(i) = B(i-1);
    }
    A(0) = B(N-1);
    return A;

}


double polyval(rowvec Poly, double x){
  //this evaluates a polynomial of the form P(x) = P[0] + P[1]*x^1 + P[2]*x^2 +... at a point x.

  double Pval = 0;
  for (int i = 0; i<Poly.n_elem;i++){
    Pval = Pval + pow(x,i)*Poly(i);
  }

  return Pval;
}

cube relaxor(cube U,double mu,double nu){

  double alpha1 = U(0,0,0);
  double rho1 = U(0,0,1);
  double u1 = U(0,0,2);
  double v1 = U(0,0,3);
  double p1 = U(0,0,4);
  double rho2 = U(0,0,5);
  double u2 = U(0,0,6);
  double v2 = U(0,0,7);
  double p2 = U(0,0,8);
  double alpha2 = 1-alpha1;
  double err = 1;
  int counter = 0;

  while(err>0.01){

    alpha1 += nu*(p1 - p2);
    rho1 += -(nu*rho1*(p1 - p2))/alpha1;
    u1 += -(mu*(u1 - u2))/(alpha1*rho1);
    v1 += -(mu*(v1 - v2))/(alpha1*rho1);
    p1 += (mu*(u1*(u1 - u2) + v1*(v1 - v2)) + nu*(p1 - p2))/alpha1 - (nu*p1*(p1 - p2))/alpha1 - (mu*u1*(u1 - u2))/alpha1 - (mu*v1*(v1 - v2))/alpha1;
    rho2 += -(nu*rho2*(p1 - p2))/(alpha1 - 1);
    u2 += -(mu*(u1 - u2))/(rho2*(alpha1 - 1));
    v2 += -(mu*(v1 - v2))/(rho2*(alpha1 - 1));
    p2 += (mu*(u1*(u1 - u2) + v1*(v1 - v2)) + nu*(p1 - p2))/(alpha1 - 1) - (mu*u2*(u1 - u2))/(alpha1 - 1) - (mu*v2*(v1 - v2))/(alpha1 - 1) - (nu*p2*(p1 - p2))/(alpha1 - 1);

    err = abs(p1-p2) + abs(u1-u2) + abs(v1-v2);
    counter++ ;
    if(counter>5000){
      cout<<"relaxation failure detected"<<endl;
      err = 0;
      }
    }

  U(0,0,0) = alpha1;
  U(0,0,1) = rho1;
  U(0,0,2) = u1;
  U(0,0,3) = v1;
  U(0,0,4) = p1;
  U(0,0,5) = rho2;
  U(0,0,6) = u2;
  U(0,0,7) = v2;
  U(0,0,8) = p2;

  return U;


}

mat meshgrid(int& N, int&M, double Lx, double Ly, double dx, mat& xx, mat& yy){
  //creates 2 mesh-grid matrices containing x and y coordinates, similar to the matlab function

  N = floor(Lx/dx)+1; //number of x direction elements
  M = floor(Ly/dx)+1; //number of y direction elements
  xx = zeros<mat>(M,N);
  yy = zeros<mat>(M,N);
  for ( int i = 0; i<N; i++){
    for(int j = 0; j<M; j++){
      xx(j,i) = i*dx;
      yy(j,i) = j*dx;
    }
  }
}

vec incFill(double x0,double xend, double increm){
  //this function creates a vector with incrementally increasing value, increm between x0 and xend,
  //similar to meshgrid but for a 1D vector.

  int j = 0;
  vec incvec = zeros<vec>(int(floor((xend-x0)/increm)+1));

  for (double i = x0; i<=xend;i = i+increm ){
    incvec(j) = i;
    j++;
  }

  return incvec;
}



rowvec polyint(rowvec Poly){
  //This function integrates a polynomial of the form P(x) = P[0] + P[1]*x^1 + P[2]*x^2 +... at a point x.

  rowvec Poly_integ(Poly.n_elem+1);
  Poly_integ(0) = 0;
  for (int i = 1;i<Poly.n_elem+1;i++){
    Poly_integ(i) = Poly(i-1)/i;
  }

  return Poly_integ;
}

rowvec polyder(rowvec Poly){
  //Find the derivative of a polynomial of the form P(x) = P[0] + P[1]*x^1 + P[2]*x^2 +... at a point x.

  rowvec Poly_deriv(Poly.n_elem);
  for (int i = 1; i<Poly.n_elem; i++){
    Poly_deriv(i-1) = Poly(i)*(i);
  }
  Poly_deriv(Poly.n_elem-1) = 0;

  return Poly_deriv;
}

void glquad(rowvec& xC, rowvec& wC, int NGP){
  //gaussian quadrature rules up to 3rd degree legendre polynomial

  xC.resize(NGP);
  wC.resize(NGP);
  if(NGP==1){
    xC(0) = 0.0;
    wC(0) = 2.0;
  }else if(NGP==2){
    xC(0) = sqrt(1.0/3.0);
    wC(0) = 1.0;
    xC(1) = -xC(0);
    wC(1) = wC(0);
  }else if(NGP==3){
    xC(0) = sqrt(3.0/5.0);
    xC(1) = 0.0;
    xC(2) = -sqrt(3.0/5.0);
    wC(0) = 5.0/9.0;
    wC(1) = 8.0/9.0;
    wC(2) = 5.0/9.0;
  }else if(NGP ==4){

  }

  xC = (xC + 1)/2;
  wC = wC/2;

}



