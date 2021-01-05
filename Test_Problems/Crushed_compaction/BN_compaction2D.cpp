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

  //Stiffened equation of state parameters:

  double gamma1 = 5.0;
  double pi1 = 0.16021;
  double q1 = 0.11795;
  double CV1 = 0.625;
  //double gamma2 = 4.4;
  double gamma2 = 1.35;
  double b2 = 0.37778;
  double CV2 = 1;


  //upstream states:


  double alpha1_0 = 0.73;
  double rho1_0 = 5.0293;
  double p1_0 = 3.5681e-4;
  double T1_0 = 0.01277;
  double c1_0 = 0.39952;

  double alpha2_0 = 0.27;
  double rho2_0 = 2.6470e-3;
  double p2_0 = 1.1843e-5;
  double T2_0 = 0.012770;
  double c2_0 = 0.077756;

  //granular parameters:

  double delta = 20;
  double mu_F = 0.05;
  double H = 0.2;

  //double pi2 = 6e8;


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


mat A_calc(mat Q_0,int NGP,mat d_psi,mat psi,rowvec xC,int N_var, int ForG){

  //this calculates the quantity int_0^1 A*dq/dzeta dzeta for use in the Discrete Galerkin predictor.
  mat F_Q0 = Q_0;
  mat A = zeros<mat>(N_var,N_var); //Jacobian matrix
  F_Q0.zeros();
  vec dq_dzeta = zeros<vec>(N_var);

  //summing over for each Q in the Q_0 matrix
  for(int p = 0; p<NGP; p++){
    for(int q = 0; q<NGP; q++){
      for(int r = 0; r<NGP; r++){
        
        double alpha1 = Q_0(p*NGP*NGP +q*NGP + r,0);
        double rho1 = Q_0(p*NGP*NGP + q*NGP + r,1);
        double u1 = Q_0(p*NGP*NGP + q*NGP + r,2);
        double v1 = Q_0(p*NGP*NGP +NGP*q + r,3);
        double p1 = Q_0(p*NGP*NGP +NGP*q + r,4);
        double rho2 = Q_0(p*NGP*NGP +NGP*q + r,5);
        double u2 = Q_0(p*NGP*NGP +NGP*q + r,6);
        double v2 = Q_0(p*NGP*NGP +NGP*q + r,7);
        double p2 = Q_0(p*NGP*NGP +NGP*q + r,8);

        double a1 = sqrt(gamma1*(p1+pi1)/rho1);
        double a2 = sqrt((gamma2*p2/rho2)*(1+b2*rho2 - b2*b2*rho2*rho2/(gamma2*(1+b2*rho2)))         );

          double dp = p2-p1;
          double du = u2-u1;
          double dv = v2-v1;


        if(ForG==1){ //if x direction
          //flux jacobian in the x-direction
          A = { {                          u1,  0,         0,  0,      0,  0,         0,  0,      0},
              {                           0, u1,      rho1,  0,      0,  0,         0,  0,      0},
              {           -dp/(alpha1*rho1),  0,        u1,  0, 1/rho1,  0,         0,  0,      0},
              {                           0,  0,         0, u1,      0,  0,         0,  0,      0},
              {                           0,  0, a1*a1*rho1,  0,     u1,  0,         0,  0,      0},
              {      (du*rho2)/(alpha1 - 1),  0,         0,  0,      0, u2,      rho2,  0,      0},
              {                           0,  0,         0,  0,      0,  0,        u2,  0, 1/rho2},
              {                           0,  0,         0,  0,      0,  0,         0, u2,      0},
              { (a2*a2*du*rho2)/(alpha1 - 1),  0,         0,  0,      0,  0, a2*a2*rho2,  0,     u2} };





          dq_dzeta.zeros();
          //cout<<"Acalc"<<endl;
          //integrated derivative using quadrature scheme:

          for(int var = 0; var<N_var; var++){
            for(int pb = 0; pb<NGP; pb++){
              for(int qb = 0; qb<NGP; qb++){
                for(int rb = 0; rb<NGP; rb++){
                  dq_dzeta(var) += polyval(d_psi.row(pb),xC(p))*polyval(psi.row(qb),xC(q))*polyval(psi.row(rb),xC(r))*Q_0(NGP*NGP*pb + NGP*qb + rb,var);
                }
              }
            }
          }

          F_Q0.row(NGP*NGP*p + NGP*q + r) = trans(A*dq_dzeta);

          
        }else{ //if y direction



            A = { {                          v1,  0,  0,         0,      0,  0,  0,         0,      0},
                {                           0, v1,  0,      rho1,      0,  0,  0,         0,      0},
                {                           0,  0, v1,         0,      0,  0,  0,         0,      0},
                {           -dp/(alpha1*rho1),  0,  0,        v1, 1/rho1,  0,  0,         0,      0},
                {                           0,  0,  0, a1*a1*rho1,     v1,  0,  0,         0,      0},
                {      (dv*rho2)/(alpha1 - 1),  0,  0,         0,      0, v2,  0,      rho2,      0},
                {                           0,  0,  0,         0,      0,  0, v2,         0,      0},
                {                           0,  0,  0,         0,      0,  0,  0,        v2, 1/rho2},
                { (a2*a2*dv*rho2)/(alpha1 - 1),  0,  0,         0,      0,  0,  0, a2*a2*rho2,     v2} };


          dq_dzeta.zeros();



          for(int var = 0; var<N_var; var++){
            for(int pb = 0; pb<NGP; pb++){
              for(int qb = 0; qb<NGP; qb++){
                for(int rb = 0; rb<NGP; rb++){
                  dq_dzeta(var) += polyval(psi.row(pb),xC(p))*polyval(d_psi.row(qb),xC(q))*polyval(psi.row(rb),xC(r))*Q_0(NGP*NGP*pb + NGP*qb + rb,var);
                }
              }
            }
          }

          F_Q0.row(NGP*NGP*p + NGP*q + r) = trans(A*dq_dzeta);
        }
      }
    }
  }

  return F_Q0;

}

mat Afunc(vec q_eval, int dir){
  mat A(9,9);
  double alpha1 = q_eval(0);
  double rho1 = q_eval(1);
  double u1 = q_eval(2);
  double v1 = q_eval(3);
  double p1 = q_eval(4);
  double rho2 = q_eval(5);
  double u2 = q_eval(6);
  double v2 = q_eval(7);
  double p2 = q_eval(8);
  double a1 = sqrt(gamma1*(p1+pi1)/rho1);
  double a2 = sqrt((gamma2*p2/rho2)*(1+b2*rho2 - b2*b2*rho2*rho2/(gamma2*(1+b2*rho2)))         );
    double dp = p2-p1;
  double du = u2-u1;
  double dv = v2 - v1;

  if(dir==1){

              //flux jacobian in the x-direction
        A = { {                          u1,  0,         0,  0,      0,  0,         0,  0,      0},
              {                           0, u1,      rho1,  0,      0,  0,         0,  0,      0},
              {           -dp/(alpha1*rho1),  0,        u1,  0, 1/rho1,  0,         0,  0,      0},
              {                           0,  0,         0, u1,      0,  0,         0,  0,      0},
              {                           0,  0, a1*a1*rho1,  0,     u1,  0,         0,  0,      0},
              {      (du*rho2)/(alpha1 - 1),  0,         0,  0,      0, u2,      rho2,  0,      0},
              {                           0,  0,         0,  0,      0,  0,        u2,  0, 1/rho2},
              {                           0,  0,         0,  0,      0,  0,         0, u2,      0},
              { (a2*a2*du*rho2)/(alpha1 - 1),  0,         0,  0,      0,  0, a2*a2*rho2,  0,     u2} };

               // cout<<"Afunc"<<endl;
 

  }else{

          A = { {                          v1,  0,  0,         0,      0,  0,  0,         0,      0},
                {                           0, v1,  0,      rho1,      0,  0,  0,         0,      0},
                {                           0,  0, v1,         0,      0,  0,  0,         0,      0},
                {           -dp/(alpha1*rho1),  0,  0,        v1, 1/rho1,  0,  0,         0,      0},
                {                           0,  0,  0, a1*a1*rho1,     v1,  0,  0,         0,      0},
                {      (dv*rho2)/(alpha1 - 1),  0,  0,         0,      0, v2,  0,      rho2,      0},
                {                           0,  0,  0,         0,      0,  0, v2,         0,      0},
                {                           0,  0,  0,         0,      0,  0,  0,        v2, 1/rho2},
                { (a2*a2*dv*rho2)/(alpha1 - 1),  0,  0,         0,      0,  0,  0, a2*a2*rho2,     v2} };
    
    }

  return A;         

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

cube compactor(cube U,double dtau){

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
  double beta;
  double F;
  double Mu;
  double Mv;
  double Eps;
  double T1 = 0;
  double T2 = 0;
  double mu = mu_F*100;
  double nu = mu*0.0001;

  //compaction coefficients:

  for(int t = 0; t<10000; t++){

    T1 = (p1+pi1)/((gamma1-1)*rho1*CV1);

    T2 = p2/((gamma2-1)*rho2*(1+b2*rho2)*CV2);

    alpha2 = 1-alpha1;

    beta = -(p2_0 - p1_0) * alpha1 * rho1 * (2-alpha1_0)*(2-alpha1_0) * log(1-alpha1) / ( (alpha1_0*rho1_0) * (2-alpha1) * (2-alpha1) * log(1-alpha1_0)  );

    F = alpha1*alpha2*(p1-p2-beta)/mu_F;

    Mu = delta*(u2-u1);
    Mv = delta*(v2-v1);
    //Mv = 0;
    Eps = Mu*u1 + Mv*v1 + H*(T2-T1);

    U(0,0,0) += dtau*F;
    U(0,0,1) += -dtau*(F*rho1)/alpha1;
    U(0,0,2) += dtau*Mu/(alpha1*rho1);
    U(0,0,3) += dtau*Mv/(alpha1*rho1);

    U(0,0,4) += dtau*(   -(Eps - Eps*gamma1 + F*p1 - F*p2 - Mu*u1 - Mv*v1 + F*gamma1*p2 + F*gamma1*pi1 + Mu*gamma1*u1 + Mv*gamma1*v1)/alpha1   );

    U(0,0,5) += -dtau*(F*rho2)/(alpha1 - 1);
    U(0,0,6) +=  dtau*Mu/(rho2*(alpha1 - 1));
    U(0,0,7) +=   dtau*Mv/(rho2*(alpha1 - 1));


    U(0,0,8) += dtau*(    -(Eps - Eps*gamma2 - Mu*u2 - Mv*v2 + Eps*b2*b2*rho2*rho2 + 2*Eps*b2*rho2 + F*gamma2*p2 + Mu*gamma2*u2 + Mv*gamma2*v2 - Eps*b2*b2*gamma2*rho2*rho2 - F*b2*b2*p2*rho2*rho2 - Mu*b2*b2*rho2*rho2*u2 - Mv*b2*b2*rho2*rho2*v2 - 2*Eps*b2*gamma2*rho2 - 2*Mu*b2*rho2*u2 - 2*Mv*b2*rho2*v2 + F*b2*b2*gamma2*p2*rho2*rho2 + Mu*b2*b2*gamma2*rho2*rho2*u2 + Mv*b2*b2*gamma2*rho2*rho2*v2 + 2*F*b2*gamma2*p2*rho2 + 2*Mu*b2*gamma2*rho2*u2 + 2*Mv*b2*gamma2*rho2*v2)/((b2*rho2 + 1)*(alpha1 - 1))     );
    //U(0,0,4) += dtau*( (mu*(u1*(u1 - u2) + v1*(v1 - v2)) + nu*(p1 - p2))/alpha1 - (nu*p1*(p1 - p2))/alpha1 - (mu*u1*(u1 - u2))/alpha1 - (mu*v1*(v1 - v2))/alpha1  );
    //U(0,0,8) +=  dtau*(  (mu*(u1*(u1 - u2) + v1*(v1 - v2)) + nu*(p1 - p2))/(alpha1 - 1) - (mu*u2*(u1 - u2))/(alpha1 - 1) - (mu*v2*(v1 - v2))/(alpha1 - 1) - (nu*p2*(p1 - p2))/(alpha1 - 1)   );

    alpha1 = U(0,0,0);
    rho1 = U(0,0,1);
    u1 = U(0,0,2);
    v1 = U(0,0,3);
    p1 = U(0,0,4);
    rho2 = U(0,0,5);
    u2 = U(0,0,6);
    v2 = U(0,0,7);
    p2 = U(0,0,8);





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

void BNeigfunc(vec& q_eval,mat& eigvec,mat& eigval,int dir){

  double gamma = 1.4;
  double alpha1 = q_eval(0);
  double rho1 = q_eval(1);
  double u1 = q_eval(2);
  double v1 = q_eval(3);
  double p1 = q_eval(4);
  double rho2 = q_eval(5);
  double u2 = q_eval(6);
  double v2 = q_eval(7);
  double p2 = q_eval(8);
  double a1 = sqrt(gamma1*(p1+pi1)/rho1);
  double a2 = sqrt((gamma2*p2/rho2)*(1+b2*rho2 - b2*b2*rho2*rho2/(gamma2*(1+b2*rho2)))         );
  double dp = p2-p1;
  double du = u2-u1;
  double dv = v2-v1;

  if(dir==1){

         
       eigvec ={{                                        1,          0, 0, 0,         0,             0, 0, 0,         0},
                {                                             0,       rho1, 1, 0,      rho1,             0, 0, 0,         0},
                {                                             0,        -a1, 0, 0,         a1,            0, 0, 0,         0},
                {                                             0,          0, 0, 1,        0,              0, 0, 0,         0},
                {                                     dp/alpha1, a1*a1*rho1, 0, 0, a1*a1*rho1,            0, 0, 0,         0},
                {   (du*du*rho2)/((a2*a2 - du*du)*(alpha1 - 1)),          0, 0, 0,         0,          rho2, 1, 0,      rho2},
                {    -(a2*a2*du)/((a2*a2 - du*du)*(alpha1 - 1)),          0, 0, 0,         0,           -a2, 0, 0,        a2},
                {                                             0,          0, 0, 0,         0,             0, 0, 1,         0},
            { (a2*a2*du*du*rho2)/((a2*a2 - du*du)*(alpha1 - 1)),          0, 0, 0,         0,    a2*a2*rho2, 0, 0, a2*a2*rho2} };


        eigval.zeros();
        eigval(5,5) = u2-a2;
        eigval(6,6) = u2;
        eigval(7,7) = u2;
        eigval(8,8) = u2+a2;
        eigval(1,1) = u1-a1;
        eigval(2,2) = u1;
        eigval(3,3) = u1;
        eigval(4,4) = u1+a1;
        eigval(0,0) = u1;

        //cout<<"Eigen"<<endl;

 
  }else{


    eigvec = {  {                                             1,         0, 0, 0,         0,         0, 0, 0,         0},
                {                                             0,      rho1, 0, 1,      rho1,         0, 0, 0,         0},
                {                                             0,         0, 1, 0,         0,         0, 0, 0,         0},
                {                                             0,       -a1, 0, 0,        a1,         0, 0, 0,         0},
                {                                     dp/alpha1, a1*a1*rho1, 0, 0, a1*a1*rho1,         0, 0, 0,         0},
                {      (dv*dv*rho2)/((a2*a2 - dv*dv)*(alpha1 - 1)),         0, 0, 0,         0,      rho2, 0, 1,      rho2},
                {                                             0,         0, 0, 0,         0,         0, 1, 0,         0},
                {       -(a2*a2*dv)/((a2*a2 - dv*dv)*(alpha1 - 1)),         0, 0, 0,         0,       -a2, 0, 0,        a2},
                { (a2*a2*dv*dv*rho2)/((a2*a2 - dv*dv)*(alpha1 - 1)),         0, 0, 0,         0, a2*a2*rho2, 0, 0, a2*a2*rho2} };

      eigval.zeros();
      eigval(5,5) = v2-a2;
      eigval(6,6) = v2;
      eigval(7,7) = v2;
      eigval(8,8) = v2+a2;
      eigval(1,1) = v1-a1;
      eigval(2,2) = v1;
      eigval(3,3) = v1;
      eigval(4,4) = v1+a1;
      eigval(0,0) = v1;
  }



  //normalising
  
    for(int i = 0; i<9; i++){
    eigvec.col(i) = eigvec.col(i)/sqrt(sum(eigvec.col(i)%eigvec.col(i)));
    }
  //eigvec.elem( find_nonfinite(eigvec) ).zeros();




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

int main(){


  //This function uses the ADER-WENO finite volume method to solver the 2D BAER NUNZIATO equations in cartesian geometries
  //This acts as a base code for an abritrary geometry/ arbitrary equation of state multiphase ADER-WENO solver


  //This uses a primitive variable formulation:




  //Domain properties
  int N_var = 9; //number of conserved variables
  double dx = 0.16; //element size
  double Lx = 20; //x width of domains
  double Ly = 20; //y width of domain
  int NGP = 2; //order of accuracy / 'N'umber of 'G'aussian 'P'oints
  mat xx; //matrix of x coordinates
  mat yy; //matrix of y coordinates
  int N; //number of cells in x direction
  int M; //number of cells in y direction
  meshgrid(N,M,Lx,Ly,dx,xx, yy); //filling x and y matrices
  
  //non-linear weighting parameters
  double eps_WP = 1e-14;
  double r_WP = 8;
  double u_piston = 0.013318;
  double us = u_piston; //x velocity of the solid body
  double vs = 0; //y velocity of the solid body
  
  

  //initial conditions for baer nunziato equations:

  

  double CFL = 0.95/(2*(NGP-1) +1);
  //Riemann problem tests override:
  //Vector given by:
  vec WL(9);
  vec WR(9);
  vec WB(9);
  double t_end =150;
  
  double x0;
  double save_t_inc = 5;
  double save_t_accu =0;
  int save_number = 1;
  //WL = {0.75, 1999.939402, 49.998485, 0, 4999849.5, 1.0, 0, 0, 1};
  WR = {alpha1_0, rho1_0, 0, 0, p1_0, rho2_0, 0 ,0,p2_0 };
  WL = WR;
 // mat alpha1 = WL(0) + (WR(0)-WL(0))*(sign(xx-0.5001)+1)/2;
  mat alpha1 = xx*0 + alpha1_0;
  mat rho1 = WL(1) + (WR(1)-WL(1))*(sign(xx-0.5001)+1)/2;
  mat u1 = WL(2) + (WR(2)-WL(2))*(sign(xx-0.5001)+1)/2;
  mat v1 = WL(3) + (WR(3)-WL(3))*(sign(xx-0.5001)+1)/2;
  mat p1 = WL(4) + (WR(4)-WL(4))*(sign(xx-0.5001)+1)/2;
  mat rho2 = WL(5) + (WR(5)-WL(5))*(sign(xx-0.5001)+1)/2;
  mat u2 = WL(6) + (WR(6)-WL(6))*(sign(xx-0.5001)+1)/2;
  mat v2 = WL(7) + (WR(7)-WL(7))*(sign(xx-0.5001)+1)/2;
  mat p2 = WL(8) + (WR(8)-WL(8))*(sign(xx-0.5001)+1)/2;

  //bubble parameters:
  //alpha1 = alpha1*(sign(0.25*0.25 - pow(xx-1,2) - pow(yy-0.75,2))-1)/2;
  //WL = {alpha1L, rho1L, u1L, v1L, p1L, rho2L, u2L, v2L p2L};
  //test taken from Tokareva paper
 /* W1 = {0.4, 800, 0, 0, 500, 0.3, 1000, 0, 0, 600};
  W2 = {0.6, 1.5, 0,0,2,0.7,1,0,0,1};
  //W2 = W1;
  
  mat alpha1 = W1(0) + (W1(5)-W1(0))*(sign(xx-0.5001)+1)/2;
  mat rho1 = W1(1) + (W1(6)-W1(1))*(sign(xx-0.501)+1)/2;
  mat u1 = W1(2) + (W1(7)-W1(2))*(sign(xx-0.5001)+1)/2;
  mat v1 = W1(3) + (W1(8)-W1(3))*(sign(xx-0.5001)+1)/2;
  mat p1 = W1(4) + (W1(9)-W1(4))*(sign(xx-0.5001)+1)/2;
  mat rho2 = W2(1) + (W2(6)-W2(1))*(sign(xx-0.5001)+1)/2;
  mat u2 = W2(2) + (W2(7)-W2(2))*(sign(xx-0.50001)+1)/2;
  mat v2 = W2(3) + (W2(8)-W2(3))*(sign(xx-0.50001)+1)/2;
  mat p2 = W2(4) + (W2(9)-W2(4))*(sign(xx-0.50001)+1)/2;*/
  mat a1 = zeros<mat>(M,N);
  mat a2 = zeros<mat>(M,N);
  


  //creating cube vector of conserved variables:

  cube U = zeros<cube>(M,N,N_var); //matrix of primative variables

  //packing conserved variables matrix
  U.slice(0) = alpha1;
  U.slice(1) = rho1;
  U.slice(2) = u1;
  U.slice(3) = v1;
  U.slice(4) = p1;
  U.slice(5) = rho2;
  U.slice(6) = u2;
  U.slice(7) = v2;
  U.slice(8) = p2;


  
  //generating Gauss legendre nodes and quadrature weights

  //checking whether the scheme is odd (create 3 stencils) or even (create 4 stencils)

  int isodd = NGP%2;
  int N_sten;
  vec lambda_s; // non-linear weighting paramater for combining stencils
  vec R; //Right hand boundary of each stencil
  vec L; //Left hand boundary of each stencil

  //generating stencil sets

  if(isodd==1){
    lambda_s = { 1, 1e5, 1 };
    R = {0, (double(NGP)-1)/2, double(NGP)-1};
    L = {double(NGP)-1, (double(NGP)-1)/2, 0};
    N_sten = 3;

  }else if(NGP == 1){
    lambda_s = { 1, 1e5, 1 };
    R = {0, (double(NGP)-1)/2, double(NGP)-1};
    L = {double(NGP)-1, (double(NGP)-1)/2, 0};
    N_sten = 3;
  }else{
    lambda_s = { 1, 1e5, 1e5, 1 };
    R = {0, floor((double(NGP)-1)/2), floor((double(NGP)-1)/2)+1, double(NGP)-1};
    L = {double(NGP)-1, floor((double(NGP)-1)/2)+1, floor((double(NGP)-1)/2), 0};
    N_sten = 4;
  }


  //generating lagrange polynomials of the form:

  // psi_int_i = [0]*psi^0 + [1]*psi^1 + [2]...

  rowvec psi_int = zeros<rowvec>(NGP+1); //intestitial vector containing preconvoluted polynomials
  psi_int(0) = 1; //initially set so psi = [1 0 0 ... 0]
  rowvec conv_row(2); //each iteration of the Lagrange polynomial constructor
  rowvec xC; //gauss legendre quadrature nodes
  
  mat psi = zeros<mat>(NGP,NGP+1); //matrix containing rows of lagrange polynomials
  mat psi_h = zeros<mat>(NGP,NGP+2); // matrix containing rows of antiderivatives of lagrange polynomials
  rowvec wC;

  glquad(xC,wC,NGP);

  for ( int j = 0; j<NGP; j++){

      for ( int  i = 0; i<NGP; i++){

        if(i!=j){

          conv_row = {-xC(i)/(xC(j)-xC(i)), 1/(xC(j)-xC(i))}; //the next product to be included in the lagrange polynomial

          //multiplying the two polynomials together (convolution is polynomial multiplication)
          rowvec interstitial_psi_int = conv(psi_int,conv_row,"full"); 
          //taking the non-zero elements
          psi_int = interstitial_psi_int(span(0,NGP));
        }
      }
      //collecting in a data matrix:
    psi.row(j) = psi_int;
    //matrix of lagrange polynomial integrals
    psi_h.row(j) = polyint(psi_int);
    //reseting
    psi_int.zeros();
    psi_int(0) = 1;
  }

  mat lim_S = zeros<mat>(N_sten,NGP); //vector containing integral limits
  cube S_integrals = zeros<cube>(NGP,NGP,N_sten); //each slice contains the integral of matrices in eqn (15) of dumbser
  cube S_inverse = zeros<cube>(NGP,NGP,N_sten);
  for (int S = 0; S<N_sten; S++){

    lim_S.row(S) = trans(incFill(-L(S),R(S),1)); //creating vector containing all increments
    for (int i = 0; i<NGP; i++){
      for(int j = 0; j<NGP; j++){

        //evaluating polynomial integrals
        S_integrals(i,j,S) = polyval(psi_h.row(j),lim_S(S,i)+1) - polyval(psi_h.row(j),lim_S(S,i));
      }
    }
    //inverting each integral matrix for fast solution:
    S_inverse.slice(S) = inv(S_integrals.slice(S));
  }

  //generating the oscillator index matrices, these are used to chose the lowest dispersion stencil

  mat Sigma_pm = zeros<mat>(NGP,NGP); //oscillation index
  rowvec del_psi_p = zeros<rowvec>(NGP+1); //nth derivative of psi WRT zeta
  rowvec del_psi_m = zeros<rowvec>(NGP+1); //nth derivative of psi WRT eta
  rowvec psi_psi_int; //interstitial value of the oscillation index as it is continuosly differentiated and multiplied

  for( int p=0;p<NGP;p++){
    for(int m = 0; m<NGP; m++){
      //initial we are taking the zeroth derivative of psi:
      del_psi_p = psi.row(p);
      del_psi_m = psi.row(m);
      //looping over each derivative
      for( int alpha = 0; alpha<NGP-1; alpha++){
        del_psi_p = polyder(del_psi_p);
        del_psi_m = polyder(del_psi_m);
        //multiplying the two polynomials together using a convolution
        psi_psi_int = polyint(conv(del_psi_p,del_psi_m,"fill"));
        //summation of the oscillation index
        Sigma_pm(p,m) = Sigma_pm(p,m) + polyval(psi_psi_int,1) - polyval(psi_psi_int,0);
      }
    }
  }
  //generating tensors from eqn (36)-(40) from Dumbser et al.
  //due to the complicated underlying integrals, these were difficult to split for readability

  mat K1_mat = zeros<mat>(NGP*NGP*NGP,NGP*NGP*NGP); //tensor multiplied by U^n+1
  mat K_zeta = zeros<mat>(NGP*NGP*NGP,NGP*NGP*NGP); //tensor multiplied by F(U);
  mat K_eta = zeros<mat>(NGP*NGP*NGP,NGP*NGP*NGP); // tensor multiplied by G(U);
  mat F_qm = zeros<mat>(NGP*NGP*NGP,NGP*NGP); // tensor multiplied by the ADER reconstructed U^n;
  mat KA_mat = zeros<mat>(NGP*NGP*NGP,NGP*NGP*NGP); // tensor to be multiplied by Ax(U)du/dx + Ay(U)du/dy;

  //differentiating psi lagrange polynomials to use in Galerkin integrals
  mat d_psi = psi;
  for(int i = 0; i<NGP;i++){
    d_psi.row(i) = polyder(psi.row(i));
  }
  //looping over all tensor indices and allocating within each Tensor, eg T_{(p,q,r),(p',q',r')}
  for (int p = 0; p<NGP; p++){ 
    for (int q = 0; q<NGP; q++){
      for (int r = 0; r<NGP; r++){
        for (int pd = 0; pd<NGP; pd++){
          for (int qd = 0; qd<NGP; qd++){
            for (int rd = 0; rd<NGP; rd++){
              //using the quadrature points to evaluate integrals.
              for (int gpi = 0; gpi <NGP; gpi++){
                for(int gpj = 0; gpj<NGP; gpj++){

                  //indexing tensors

                  //using gauss-legendre quadrature to integrate terms:

                  //integrals of u^n+1
                  K1_mat(p*NGP*NGP + q*NGP + r , pd*NGP*NGP + qd*NGP + rd ) = K1_mat(p*NGP*NGP + q*NGP + r , pd*NGP*NGP + qd*NGP + rd ) + wC(gpi)*wC(gpj)*polyval(psi.row(p),xC(gpi))*polyval(psi.row(q),xC(gpj))*polyval(psi.row(r),1)*polyval(psi.row(pd),xC(gpi))*polyval(psi.row(qd),xC(gpj))*polyval(psi.row(rd),1);

                  for(int gpk = 0; gpk<NGP;gpk++){

                    K1_mat(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd ) = K1_mat(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd ) - wC(gpi)*wC(gpj)*wC(gpk)*polyval(psi.row(p),xC(gpi))*polyval(psi.row(q),xC(gpj))*polyval(d_psi.row(r),xC(gpk))*polyval(psi.row(pd),xC(gpi))*polyval(psi.row(qd),xC(gpj))*polyval(psi.row(rd),xC(gpk));
                    //x flux integral
                    K_zeta(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd ) = K_zeta(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd ) + wC(gpi)*wC(gpj)*wC(gpk)*polyval(psi.row(p),xC(gpi))*polyval(psi.row(q),xC(gpj))*polyval(psi.row(r),xC(gpk))*polyval(d_psi.row(pd),xC(gpi))*polyval(psi.row(qd),xC(gpj))*polyval(psi.row(rd),xC(gpk));
                    //y flux integral
                    K_eta(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd ) = K_eta(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd )   + wC(gpi)*wC(gpj)*wC(gpk)*polyval(psi.row(p),xC(gpi))*polyval(psi.row(q),xC(gpj))*polyval(psi.row(r),xC(gpk))*polyval(psi.row(pd),xC(gpi))*polyval(d_psi.row(qd),xC(gpj))*polyval(psi.row(rd),xC(gpk));

                    KA_mat(p*NGP*NGP + q*NGP + r, pd*NGP*NGP + qd*NGP + rd) += wC(gpi)*wC(gpj)*wC(gpk)*polyval(psi.row(p),xC(gpi))*polyval(psi.row(q),xC(gpj))*polyval(psi.row(r),xC(gpk))*polyval(psi.row(pd),xC(gpi))*polyval(psi.row(qd),xC(gpj))*polyval(psi.row(rd),xC(gpk));
                  }
                }
              }
            }
            for(int gpi = 0; gpi<NGP;gpi++){
              for (int gpj = 0; gpj<NGP; gpj++){

                //initial condition integrals
                F_qm(p*NGP*NGP + q*NGP + r, pd*NGP + qd) = F_qm(p*NGP*NGP + q*NGP + r, pd*NGP + qd) + wC(gpi)*wC(gpj)*polyval(psi.row(p),xC(gpi))*polyval(psi.row(q),xC(gpj))*polyval(psi.row(r),0)*polyval(psi.row(pd),xC(gpi))*polyval(psi.row(qd),xC(gpj));
              }
            }
          }
        }   
      }
    }   
  }

  //Inverting K1_mat to increase speed in non-linear solver. 
  mat K1_inv = inv(K1_mat); 

  double t = 0.0; //time variable
  mat omega_ts = zeros<mat>(N_sten,1); //matrix containing non-linear weighting with oscillator index

  nmat W_p; //creating an N-dimensional matrix containing ADER reconstructed variables in the x-direction
  W_p.sizemat({M,N,N_var,NGP});

  vec W_inter; //interstitial value used for ADER reconstruction for each stencil

  nmat W_sx; //linear reconstruction by each stencil in the x direction, to then by non-linearly recombined
  W_sx.sizemat({M,N,N_sten,N_var,NGP});

  double sigma_s = 0; //weighting function
  mat omega_s = omega_ts; //matrix containing non-linear weights

  nmat W_h; //complete non-linear ADER reconstruction in Lagrange polynomial form
  W_h.sizemat({M,N,N_var,NGP,NGP});

  nmat W_s; //linear reconstruction of x-reconstructed ADER variables, one for each stencil
  W_s.sizemat({M,N,N_sten,N_var,NGP,NGP});

  nmat w_m; 
  w_m.sizemat({NGP*NGP,N_var,M,N});

  mat Q_0 = zeros<mat>(NGP*NGP*NGP,N_var); //initial guess for WENO non-linear solver

  nmat Q_pqr; //tensor containing evolution of lagrange polynomials over space and time for each element
  Q_pqr.sizemat({M,N,N_var,NGP,NGP,NGP});

  mat Q_old = Q_0;  //older version of ADER estimate for estimating error
  mat F_Qvect = Q_0; //x-direction flux vector of Q_0;
  mat G_Qvect = Q_0; //y-direction flux vector of Q_0;
  double err =1; //error estimate
  double Smax1 =0; //fastest wave speed estimate
  double Smax2 = 0;
  double dt; //time step
  ; //Courant Friedrich Lewis number
  mat a; //matrix containing speeds of sound.

  double S_max, rhoL, rhoR, uL, uR, vL, vR, pL, pR, aL, aR; //left and right hand state properties for Riemann rpoblem.

  double rhobar, abar, ppvrs, pstar, qL, qR, SL, SR, EL, ER, Sstar;


  vec psi_L = zeros<vec>(NGP);
  vec psi_R = zeros<vec>(NGP);


  //changed to vecs from rowvec
  vec q_minus = zeros<vec>(N_var);
  vec q_plus = zeros<vec>(N_var);
  vec qs = zeros<vec>(N_var);
  vec Dplus = zeros<vec>(N_var);
  vec Dminus = zeros<vec>(N_var);
  
  mat As = zeros<mat>(N_var,N_var);
  mat AsT = zeros<mat>(N_var,N_var);
  //cx_vec eigval_cx = zeros<cx_vec>(N_var);
  //cx_mat eigvec_cx = zeros<cx_mat>(N_var,N_var);
  mat eigval = zeros<mat>(N_var,N_var);
  mat eigvec = zeros<mat>(N_var,N_var);
  mat abs_eigval = zeros<mat>(N_var,N_var);
  mat abs_A = zeros<mat>(N_var,N_var);
  mat eigval_plus = zeros<mat>(N_var,N_var);
  mat eigval_minus = zeros<mat>(N_var,N_var);
  mat eigvecINV = zeros<mat>(N_var,N_var);
  rowvec WS;
  rowvec XS;
  glquad(XS,WS,3);
  //WS << 0.5 << 1 << 1 << 1 <<  0.5 << endr;
  vec AQx_vect = zeros<vec>(N_var);
  vec q_eval = zeros<vec>(N_var);

  vec dqdeta = zeros<vec>(N_var);
  vec Qx_vect = zeros<vec>(N_var);

  cube Dxmi = zeros<cube>(M,N,N_var);
  cube Dxpi = zeros<cube>(M,N,N_var);
  cube Dymi = zeros<cube>(M,N,N_var);
  cube Dypi = zeros<cube>(M,N,N_var);
  cube AQx = zeros<cube>(M,N,N_var);
  cube AQy = zeros<cube>(M,N,N_var);


  rowvec f_tild = zeros<rowvec>(N_var);
  rowvec UL = zeros<rowvec>(N_var);
  rowvec UR = zeros<rowvec>(N_var);
  rowvec FL = zeros<rowvec>(N_var);
  rowvec FR = zeros<rowvec>(N_var);
  rowvec Fhllc = zeros<rowvec>(N_var);
  rowvec D = zeros<rowvec>(N_var);
  rowvec FstarL = zeros<rowvec>(N_var);
  rowvec FstarR = zeros<rowvec>(N_var);
  mat psi_QD = zeros<mat>(NGP,NGP);

  nmat qLx; //evaluation of Q_pqr tensor and left hand state in the x-direction
  qLx.sizemat({M,N,N_var,NGP,NGP});

  nmat qRx; //evaluation of Q_pqr tensor at right hand state in the x-direction
  qRx.sizemat({M,N,N_var,NGP,NGP});

  nmat qLy;
  qLy.sizemat({M,N,N_var,NGP,NGP});

  nmat qRy;
  qRy.sizemat({M,N,N_var,NGP,NGP});

  //Ghost fluid parameters

  //creating the signed distance function phi:
  //this is specifically for having two circles

  //defining extra phi parameter to stitch together:

  mat phi = 0*rho1-1;


  //finding its derivatives:
  mat phi_x = zeros<mat>(M,N); //x derivative
  mat phi_y = zeros<mat>(M,N); //y derivative
  //using central differences:
  double dx2 = 2*dx; //this is used to increase speed in central difference operations.

  //finding normal vectors from equation XX
  mat nx = zeros<mat>(M,N); //x component of the normal vector
  mat ny = zeros<mat>(M,N); //y component of the normal vector

  //setting up a cylinder at (0.5,0.5)


  //cout<<ny<<endl;
  //sleep(5);

  double dtau = 0; //step size of imaginary time for ghost fluid


  //spatial derivatives for ghost fluid method:
  double drho1dx = 0; //density
  double dp1dx = 0; //pressure
  double du1dx = 0; //x velocity 
  double dv1dx = 0; // y velocity
  double dalpha1dx = 0;

  double drho1dy = 0;
  double dp1dy = 0;
  double du1dy = 0;
  double dv1dy = 0;
  double dalpha1dy = 0;


  double drho2dx = 0; //density
  double dp2dx = 0; //pressure
  double du2dx = 0; //x velocity 
  double dv2dx = 0; // y velocity

  double drho2dy = 0;
  double dp2dy = 0;
  double du2dy = 0;
  double dv2dy = 0;


  //relative velocities
  double urel = 0;
  double vrel = 0;
  //reflected velocities
  double uref = 0;
  double vref = 0;
  //extrapolated scalar fields: (for use in GFM)
  mat alpha1_ext = zeros<mat>(M,N);
  mat rho1_ext = zeros<mat>(M,N); 
  mat p1_ext = zeros<mat>(M,N);
  mat u1_ext = zeros<mat>(M,N);
  mat v1_ext = zeros<mat>(M,N);
  mat rho2_ext = zeros<mat>(M,N); 
  mat p2_ext = zeros<mat>(M,N);
  mat u2_ext = zeros<mat>(M,N);
  mat v2_ext = zeros<mat>(M,N);
  mat psi_eval = zeros<mat>(NGP,NGP);
  mat dpsi_eval = zeros<mat>(NGP,NGP);

  //generating evaluation matrices for faster polynomial evaluation:
  for(int pd = 0; pd<NGP; pd++){
    for(int p = 0; p<NGP; p++){
      psi_eval(pd,p) = polyval(psi.row(pd),xC(p));
      dpsi_eval(pd,p) = polyval(d_psi.row(pd),xC(p));
    }
  }

  //polyval(psi.row(pd),xC(p))



  while( t< t_end){
    //time-stepping variables
    //reseting


        phi = sign(4*4- pow(xx-14 +(t_end*us - us*t),2) - pow(yy- 10,2))%sqrt(abs(2*2- pow(xx-14+(t_end*us - us*t),2) - pow(yy-10,2)));
    //phi = 0.1 - abs(yy-0.2);

    //phi = (2-xx);



    for( int i = 1; i<N-1; i++ ) { //looping from the second element to the second last
      for(int j = 1; j<M-1; j++){
        //finding central difference derivatives
        phi_x(j,i) = (phi(j,i+1)-phi(j,i-1))/(dx2);
        phi_y(j,i) = (phi(j+1,i)-phi(j-1,i))/(dx2);

      
        //creating the normal vector, whilst avoiding production of NaN values:
        if (phi_x(j,i) == 0){
          nx(j,i) = 0;
        }else{
          nx(j,i) = phi_x(j,i)/sqrt(phi_x(j,i)*phi_x(j,i) + phi_y(j,i)*phi_y(j,i));
      
           }
           
        if (phi_y(j,i) == 0){
          ny(j,i) = 0;
        }else{
          ny(j,i) = phi_y(j,i)/sqrt(phi_x(j,i)*phi_x(j,i) + phi_y(j,i)*phi_y(j,i));
        }

      }
    }




    std::clock_t start1;
    double duration1;
    start1 = std::clock();
    omega_ts.zeros();
    W_p.valfill(0.0);
    //initially reconstructing in the x-direction

    //looping over elements and stencils
    for(int i = L(0); i<N-R(R.n_elem-1);i++){
      for(int j = L(0);j<M-R(R.n_elem-1);j++){
        for(int var=0;var<N_var;var++){
          for(int S=0;S<N_sten;S++){
            //solving SOE using each stencil to fit a lagrange polynomial over each element 
            W_inter = S_inverse.slice(S)*vectorise((U(span(j,j),span(i+lim_S(S,0),i+lim_S(S,lim_S.n_cols-1)),span(var,var))));
                
            for(int p = 0; p<NGP;p++){
              //placing these lagrange polynomials into a data-structure
              W_sx.setval({j,i,S,var,p},W_inter(p));
            }
            //resetting
            sigma_s = 0;
            for(int p=0; p<NGP;p++){
              for(int m =0; m<NGP;m++){
                //multiplying each stencil by an oscillation index
                sigma_s+= Sigma_pm(p,m)*W_sx.getval({j,i,S,var,p})*W_sx.getval({j,i,S,var,m});
              }
            }
            //non-linear combination of stencils
            omega_ts(S) = lambda_s(S)/pow((sigma_s + eps_WP),r_WP);
          }
          //summation of combinations
          omega_s = omega_ts/accu(omega_ts);

          for(int p =0; p<NGP;p++){
            for(int S=0;S<N_sten;S++){
              //placing these non-linear combination of polynomial coefficients into the x reconstruction Tensor
              W_p.setval({j,i,var,p},W_p.getval({j,i,var,p})+ omega_s(S)*W_sx.getval({j,i,S,var,p}) );
            }
          }
        }
      }
    }
    

    //reconstruction of x reconstruction in the y direction:
    //resetting
    W_h.valfill(0);
    W_s.valfill(0);
    int ki; //indexing ariable
    vec intervec = zeros<vec>(NGP); //vector containing x-reconstruction to then be reconstructed along each stencil

    //looping over elements, stencils, variables and lagrange coefficients of x-direction polynomials
    for(int i=L(0);i<N-L(0);i++){
      for(int j = L(0); j<M-L(0);j++){
        for(int var = 0; var<N_var; var++){
          for(int p = 0; p<NGP; p++){
            for(int S=0; S<N_sten;S++){
              
              //creating vector for W_inter:
              for(int k = lim_S(S,0); k<=lim_S(S,lim_S.n_cols-1);k++){
                intervec(k-lim_S(S,0)) = W_p.getval({j+k,i,var,p});
              }
              W_inter = S_inverse.slice(S)*intervec; //linear reconstruction for each stencil
              
              for(int q=0; q<NGP;q++){
                W_s.setval({j,i,S,var,p,q},W_inter(q)); //placing into a data structure
              }

              //calculating oscillation index
              sigma_s = 0;
              for(int q= 0; q<NGP;q++){
                for(int m =0; m<NGP; m++){
                  
                  sigma_s += Sigma_pm(q,m)*W_s.getval({j,i,S,var,p,q})*W_s.getval({j,i,S,var,p,m});
                }
              }
              //non-linear combination of oscillation limited polynomials
              omega_ts(S) = lambda_s(S)/pow(sigma_s + eps_WP,r_WP);
            }
            //summation of stencils
            omega_s = omega_ts/accu(omega_ts);

            //placing into the reconstruction Tensor:
            for(int q=0; q<NGP;q++){
              for(int S=0; S<N_sten; S++){
                W_h.setval({j,i,var,p,q},W_h.getval({j,i,var,p,q}) + W_s.getval({j,i,S,var,p,q})*omega_s(S) );
              }
            }
          }
        }
      }
    }
    //variables are now reconstruction at time t = n*dt, now move on to solution using a discrete Galerkin predictor
    //finding raw variables:

    rho1 = U.slice(1);
    u1 = U.slice(2);
    v1 = U.slice(3);
    p1 = U.slice(4);
    rho2 = U.slice(5);
    u2 = U.slice(6);
    v2 = U.slice(7);
    p2 = U.slice(8);
    //speed of sound for ideal gas
    
    a1 = sqrt(gamma1*(p1+pi1)/(rho1));
    a2 = - (b2*b2*rho2%rho2);
    a2 = a2/(gamma2*(1+b2*rho2));
    a2 += 1+b2*rho2; 
    cout<<"test"<<endl;
    a2 = sqrt((gamma2*p2/rho2)%a2);
    //estimating maximum wave speed for finding Courant number
    Smax1 = max(max(a1+sqrt(u1%u1 + v1%v1) +us));
    Smax2 = max(max(a2 + sqrt(u2%u2 + v2%v2) +us));
    //finding timestep
    dt = CFL*dx/std::max(Smax1,Smax2);
    
    //resetting solutions
    w_m.valfill(0);
    Q_0.zeros();


    for(int i= L(0); i<N-L(0); i++){
      for(int j=L(0); j<M-L(0);j++){
        for(int var=0; var<N_var; var++){
          for(int pd=0; pd<NGP; pd++){
            for(int qd=0; qd<NGP; qd++){
              w_m.setval({pd*NGP + qd,var,j,i},W_h.getval({j,i,var,pd,qd}) );

            }
          }
        }
      }
    }
    //resetting
    Q_pqr.valfill(0);

    std::clock_t start;
    double duration;
    start = std::clock();
    
    //creating an initial guess for vector containing polynomial coefficients over time and space for each variable:
    for(int i= L(0); i<N-L(0); i++){
      for(int j=L(0); j<M-L(0);j++){
        for(int var=0; var<N_var; var++){
          for(int pd=0; pd<NGP; pd++){
            for(int qd=0; qd<NGP; qd++){
              for(int rd=0; rd<NGP; rd++){

                Q_0(pd*NGP*NGP+qd*NGP+rd,var) = W_h.getval({j,i,var,pd,qd});
              }
            }
          }
        }

        err = 1; //creating initial large error estimate

        //iterating initial guess until error is bellow threshold
        //cout<<Q_0<<endl;
        int loop_count = 0;
        while(err>3e-14){
          loop_count++;
          //evaluating fluxs for each lagrange coefficient over space and time

          Q_old = Q_0;
          F_Qvect = dt*A_calc(Q_0,NGP,d_psi,psi,xC,N_var,1)/dx; 
          //sleep(1);

          G_Qvect = dt*A_calc(Q_0,NGP,d_psi,psi,xC,N_var,2)/dx;

          //generating update for matrix equation
          for(int var = 0; var<N_var; var++){
            Q_0.col(var) = K1_inv*(F_qm*w_m.nslice(0,{0,var,j,i},NGP*NGP-1) - KA_mat*F_Qvect.col(var) - KA_mat*G_Qvect.col(var));
            //Q_0.col(var) = K1_inv*(F_qm*w_m.nslice(0,{0,var,j,i},NGP*NGP-1) - KA_mat*F_Qvect.col(var)  );
          }
          //error estimate
          err = accu(abs(Q_0-Q_old)/(abs(Q_0) + 1));
          //cout<<err<<endl;
          //cout<<Q_0<<endl; 
          if(loop_count>1000){
            cout<<"error is " <<err<<endl;
            //cout<<Q_old<<endl;
            //cout<<i<<"   "<<j<<endl;
            err = 0;
          }
        }
        //cout<<"done"<<endl;
        //storing within an index for each element
        for(int var = 0; var<N_var; var++){
          for(int p = 0; p<NGP; p++){
            for(int q = 0; q<NGP; q++){
              for(int r = 0; r<NGP; r++){
                Q_pqr.setval({j,i,var,p,q,r}, Q_0(NGP*NGP*p + NGP*q + r,var));

              }
            }
          }
        }
      }
    } 
    
    //cout<<"Run time is" duration<<endl;
    //sleep(100)
    //sleep(10);
    //finding the spatial average over time for use in the finite volume update scheme

    for(int i = 0; i<NGP; i++){
      //evaluating legengdre polynomials at left and right hand states
      psi_L(i) = polyval(psi.row(i),0);
      psi_R(i) = polyval(psi.row(i),1);
      for(int j = 0; j<NGP; j++){
        psi_QD(i,j) = polyval(psi.row(i),xC(j));
      }

    }

    //std::clock_t start;


    start = std::clock();

    for(int i= L(0) ; i<N-L(0); i++){
      for(int j=L(0) ; j<M-L(0);j++){
        for(int qd = 0; qd<NGP; qd++){
          for(int rd = 0; rd<NGP; rd++){
            q_minus.zeros();
            q_plus.zeros();

            //looping over each gauss-legendre point in the y - t direction to solve the riemann problem integral
            for(int var = 0; var<N_var; var++){
              for(int p = 0; p<NGP; p++){
                for(int q = 0; q<NGP; q++){
                  for(int r = 0; r<NGP; r++){

                    //evaluating left and right hand predictions of variables at the cell interface of i and i+1

                    q_minus(var) += Q_pqr.getval({j,i,var,p,q,r})*psi_R(p)*psi_QD(q,qd)*psi_QD(r,rd);
                    q_plus(var) += Q_pqr.getval({j,i+1,var,p,q,r})*psi_L(p)*psi_QD(q,qd)*psi_QD(r,rd);
                  }

                }
              }
            }

            //finding left and right hand side variables
            //reseting
            //DS.zeros();
            
            Dplus.zeros();
            Dminus.zeros();
            
            for(int s = 0; s<3; s++){
              //linear parametric path:
              qs = q_minus + XS(s)*(q_plus-q_minus);

              
              //finding flux jacobian
              As = Afunc(qs,1);

              BNeigfunc(qs,eigvec,eigval,1);

              if(det(eigvec)==0|isnan(det(eigvec))==1){
                //cout<<eigvec<<endl;
                //cout<<As<<endl;
              }
              
              abs_A = eigvec*abs(eigval)*inv(eigvec);
             
              
              Dplus += WS(s)*(As + abs_A)*(q_plus - q_minus)/2;
              Dminus += WS(s)*(As - abs_A)*(q_plus - q_minus)/2;
              
              //Dplus += WS(s)*(As)*(q_plus - q_minus)/12;
              
            }
                     
              
            for(int var = 0; var<N_var; var++){
              //cubes of .riemann solution of the flux jacobian at each cell boundary in x direction
              Dxpi(j,i,var) += wC(qd)*wC(rd)*Dplus(var);
              Dxmi(j,i,var) += wC(qd)*wC(rd)*Dminus(var);
            }



            Dplus.zeros();
            Dminus.zeros();
            

          }
        }
        
        AQx_vect.zeros();
        Qx_vect.zeros();
        for(int p = 0; p<NGP; p++){
          for(int q = 0; q<NGP; q++){
            for(int r = 0; r<NGP; r++){
              q_eval.zeros();
              dqdeta.zeros();
              for(int var = 0; var<N_var; var++){
                for(int pd = 0; pd<NGP; pd++){
                  for(int qd = 0; qd<NGP; qd++){
                    for(int rd = 0; rd<NGP; rd++){
                     //q_eval(var) += polyval(psi.row(pd),xC(p))*polyval(psi.row(qd),xC(q))*polyval(psi.row(rd),xC(r))*Q_pqr.getval({j,i,var,pd,qd,rd});
                     // dqdeta(var) += polyval(d_psi.row(pd),xC(p))*polyval(psi.row(qd),xC(q))*polyval(psi.row(rd),xC(r))*Q_pqr.getval({j,i,var,pd,qd,rd});

                      q_eval(var) += psi_eval(pd,p)*psi_eval(qd,q)*psi_eval(rd,r)*Q_pqr.getval({j,i,var,pd,qd,rd});
                      dqdeta(var) += dpsi_eval(pd,p)*psi_eval(qd,q)*psi_eval(rd,r)*Q_pqr.getval({j,i,var,pd,qd,rd});

                    }
                  }
                }
              }

              Qx_vect += wC(p)*wC(q)*wC(r)*Afunc(q_eval,1)*dqdeta;
              //cout<<Qx_vect<<endl;


            }
          }
        }

        
        for(int var = 0; var<N_var;var++){
          AQx(j,i,var) = Qx_vect(var);
          //cout<<AQx(j,i,var)<<endl;
        }
        




      }
    }
    cout<< "Time for GPU is " <<(std::clock() - start)/ (double) CLOCKS_PER_SEC <<endl;
    sleep(5);

    for(int i= L(0); i<N-L(0) ; i++){
      for(int j=L(0); j<M-L(0);j++){

        for(int pd = 0; pd<NGP; pd++){
          for(int rd = 0; rd<NGP; rd++){
            q_minus.zeros();
            q_plus.zeros();


            for(int var = 0; var<N_var; var++){
              for(int p = 0; p<NGP; p++){
                for(int q = 0; q<NGP; q++){
                  for(int r = 0; r<NGP; r++){

                    //evaluating left and right hand states for flux calculations in the y direction
                    //q_minus(var) += Q_pqr.getval({j,i,var,p,q,r});//*polyval(psi.row(p),xC(pd))*polyval(psi.row(q),1)*polyval(psi.row(r),xC(rd));
                    //q_plus(var) += Q_pqr.getval({j+1,i,var,p,q,r});//*polyval(psi.row(p),xC(pd))*polyval(psi.row(q),0)*polyval(psi.row(r),xC(rd));
                    q_minus(var) += Q_pqr.getval({j,i,var,p,q,r})*psi_R(q)*psi_QD(p,pd)*psi_QD(r,rd);
                    q_plus(var) += Q_pqr.getval({j+1,i,var,p,q,r})*psi_L(q)*psi_QD(p,pd)*psi_QD(r,rd);

                  }

                }
              }
            }

            Dplus.zeros();
            Dminus.zeros();
            //cout<<q_minus<<endl;

            //integrating through phase space
            
            for(int s = 0; s<3; s++){
              //linear parametric path:
              qs = q_minus + XS(s)*(q_plus-q_minus);
              //qs = floor(qs*1000)/1000;
              //qs = round(qs*10000)/10000;
              //finding flux jacobian
              As = Afunc(qs,2);

              BNeigfunc(qs,eigvec,eigval,2);
              if(det(eigvec)==0){
                //cout<<eigvec<<endl;
                //cout<<As<<endl;
              }
              //cout<<eigvec<<endl;
              abs_A = eigvec*abs(eigval)*inv(eigvec);


              Dplus += WS(s)*(As + abs_A)*(q_plus - q_minus)/2;
              Dminus += WS(s)*(As - abs_A)*(q_plus - q_minus)/2;
            }
                     
              
            for(int var = 0; var<N_var; var++){
              //cubes of riemann solution of the flux jacobian at each cell boundary in x direction
              Dypi(j,i,var) += wC(pd)*wC(rd)*Dplus(var);
              Dymi(j,i,var) += wC(pd)*wC(rd)*Dminus(var);
            }



            Dplus.zeros();
            Dminus.zeros();
            //cout<<q_minus<<endl;

            //integrating through phase space

          }

    
          
        }

        AQx_vect.zeros();
        Qx_vect.zeros();
        for(int p = 0; p<NGP; p++){
          for(int q = 0; q<NGP; q++){
            for(int r = 0; r<NGP; r++){
              q_eval.zeros();
              dqdeta.zeros();
              for(int var = 0; var<N_var; var++){
                for(int pd = 0; pd<NGP; pd++){
                  for(int qd = 0; qd<NGP; qd++){
                    for(int rd = 0; rd<NGP; rd++){
                     // q_eval(var) += polyval(psi.row(pd),xC(p))*polyval(psi.row(qd),xC(q))*polyval(psi.row(rd),xC(r))*Q_pqr.getval({j,i,var,pd,qd,rd});
                    //  dqdeta(var) += polyval(psi.row(pd),xC(p))*polyval(d_psi.row(qd),xC(q))*polyval(psi.row(rd),xC(r))*Q_pqr.getval({j,i,var,pd,qd,rd});

                      q_eval(var) += psi_eval(pd,p)*psi_eval(qd,q)*psi_eval(rd,r)*Q_pqr.getval({j,i,var,pd,qd,rd});
                      dqdeta(var) += dpsi_eval(qd,q)*psi_eval(pd,p)*psi_eval(rd,r)*Q_pqr.getval({j,i,var,pd,qd,rd});

                    }
                  }
                }
              }

              Qx_vect += wC(p)*wC(q)*wC(r)*Afunc(q_eval,2)*dqdeta;
              //cout<<Qx_vect<<endl;


            }
          }
        }


        for(int var = 0; var<N_var;var++){
          AQy(j,i,var) = Qx_vect(var);
         // cout<<AQx(j,i,var)<<endl;
        }



      }
    }
    
    duration = (std::clock() - start)/ (double) CLOCKS_PER_SEC;
    //updating conservative scheme
    //cout<<trans(Dxmi.slice(0))<<endl;
    //sleep(10);
    //cout<<"poo"<<endl;


    U(span(1,M-2),span(1,N-2),span(0,N_var-1)) -= (dt/dx)*(AQx(span(1,M-2),span(1,N-2),span(0,N_var-1)) + Dxmi(span(1,M-2),span(1,N-2),span(0,N_var-1)) + Dxpi(span(1,M-2),span(0,N-3),span(0,N_var-1)));
    for(int var = 0; var<N_var; var++){
      //cout<<trans(U.slice(var))<<endl;
    }
    //cout<<U<<endl;

    U(span(1,M-2),span(1,N-2),span(0,N_var-1)) -= (dt/dx)*(AQy(span(1,M-2),span(1,N-2),span(0,N_var-1)) + Dymi(span(1,M-2),span(1,N-2),span(0,N_var-1)) + Dypi(span(0,M-3),span(1,N-2),span(0,N_var-1)));

    //U(span(1,M-2),span(1,N-2),span(0,N_var-1)) -= (dt/dx)*(AQx(span(1,M-2),span(1,N-2),span(0,N_var-1)) + Dxmi(span(1,M-2),span(1,N-2),span(0,N_var-1)) + Dxpi(span(1,M-2),span(0,N-3),span(0,N_var-1)));
    //adding gravity:
    double g = 9.8;
   // U.slice(3) = U.slice(3) - dt*g;
   // U.slice(7) = U.slice(7) - dt*g;

/*
    muR = dt*0.01;
    //applying relaxation parameters:
    for(int i= L(0)+1; i<N-R(NGP-1)-1; i++){
      for(int j=L(0)+1; j<M-R(NGP-1)-1;j++){
        U.slice(2) += -(muR*(u1 - u2))/(alpha1*rho1);
        U.slice(3) 


*/
        for(int i = 0; i<N; i++){
          for(int j=0; j<M; j++){
            if(phi(j,i)<=0){
              U(span(j,j),span(i,i),span(0,N_var-1)) = compactor(U(span(j,j),span(i,i),span(0,N_var-1)),dt/10000);
            }
            }
        }

  

    //defining boundary conditions:
    
    
    cube BCLx = U(span(0,M-1),span(L(0)+1,2*L(0)+1),span(0,N_var-1));
    cube BCRx = U(span(0,M-1),span(N-2*L(0)-2,N-L(0)-2),span(0,N_var-1));
    //cout<<AQx<<endl;
    
    BCLx = fliplrcube(BCLx);
    BCRx = fliplrcube(BCRx);

    U(span(0,M-1),span(0,L(0)),span(0,N_var-1)) = BCLx;
    U(span(0,M-1),span(N-L(0)-1,N-1),span(0,N_var-1)) = BCRx;
    U(span(0,M-1),span(N-L(0)-1,N-1),span(2,2)) = -BCRx.slice(2);
    U(span(0,M-1),span(N-L(0)-1,N-1),span(6,6)) = -BCRx.slice(6);

    cube BCLy = U(span(L(0)+1,2*L(0)+1),span(0,N-1),span(0,N_var-1));
    cube BCRy = U(span(M-2*L(0)-2,M-L(0)-2),span(0,N-1),span(0,N_var-1));

    BCLy = flipudcube(BCLy);
    BCRy = flipudcube(BCRy);
    

    U(span(0,L(0)),span(0,N-1),span(0,N_var-1)) = BCLy;
    U(span(M-L(0)-1,M-1),span(0,N-1),span(0,N_var-1)) = BCRy;

     //U(span(0,M-1),span(0,L(0)),span(2,2)) = -BCLx(span(0,M-1),span(0,L(0)),span(2,2));
    // U(span(0,M-1),span(0,L(0)),span(6,6)) = -BCLx(span(0,M-1),span(0,L(0)),span(6,6));
    //zero velocity conditions:
    /*
    for(int i = 0; i<N_var; i++){
      U(span(0,L(0)),span(0,N-1),span(i,i)) = U(span(0,L(0)),span(0,N-1),span(i,i))*0 + WL(i);
    }

    U(span(0,L(0)),span(0,N-1),span(3,3)) = -U(span(0,L(0)),span(0,N-1),span(3,3))*0;
    U(span(0,L(0)),span(0,N-1),span(7,7)) = -U(span(0,L(0)),span(0,N-1),span(7,7))*0;
    
    */

    //applying pressure and velocity relaxation:

    //assuming no surface tension:
    //double nu= dt*5;
    //double mu = dt*5;






    Dxpi.zeros();
    Dxmi.zeros();
    Dymi.zeros();
    Dypi.zeros();
    AQx.zeros();
    AQy.zeros();

    //finding primitive variables for interpolation
    alpha1 = U.slice(0);
    if(min(min(alpha1))<0){
      cout<<min(min(alpha1))<<endl;
    }
    rho1 = U.slice(1);
    u1 = U.slice(2);
    v1 = U.slice(3);
    p1 = U.slice(4);
    rho2 = U.slice(5);
    u2 = U.slice(6);
    v2 = U.slice(7);
    p2 = U.slice(8);
    
    dtau = dt/10;
    
    for(int tau=0; tau<500;tau++){
      

      alpha1_ext = alpha1;
      u1_ext = u1;
      v1_ext = v1;
      p1_ext = p1;
      rho1_ext = rho1;
      u2_ext = u2;
      v2_ext = v2;
      p2_ext = p2;
      rho2_ext = rho2;
    
        for(int i=1;i<N-1;i++){
            for(int j=1; j<M-1; j++){
                if(phi(j,i)>=0){ //ie cell ij is a Ghost cell
                    if(nx(j,i)>0){//calculating which direction derivative to use
                      //calculating derivatives:
                        drho1dx = (rho1(j,i) - rho1(j,i-1))/dx;
                        dp1dx = (p1(j,i) - p1(j,i-1))/dx;

                        drho2dx = (rho2(j,i) - rho2(j,i-1))/dx;
                        dp2dx = (p2(j,i) - p2(j,i-1))/dx;
                        dalpha1dx = (alpha1(j,i) - alpha1(j,i-1))/dx;
              //if the velocity derivative uses a value from outside of the Ghost region
              //the relative tangential component of velocity (compared to the solid body)
              //must be reflected:
                  if(phi(j,i-1)>=0){ //checking if both cells are ghost cells
                    du1dx = (u1(j,i) - u1(j,i-1))/dx;
                    dv1dx = (v1(j,i) - v1(j,i-1))/dx;

                    du2dx = (u2(j,i) - u2(j,i-1))/dx;
                    dv2dx = (v2(j,i) - v2(j,i-1))/dx;
                  }else{
                    //The cell being used for interpolation is not a ghost cell:
                    //Changing to solid body's frame of reference
                    urel = u1(j,i-1) - us;
                    vrel = v1(j,i-1) - vs;
                    //calculating reflected velocity relative to the body's frame of referece
                    uref = urel - 2*(urel*nx(j,i-1) + vrel*ny(j,i-1))*nx(j,i-1);
                    vref = vrel - 2*(urel*nx(j,i-1) + vrel*ny(j,i-1))*ny(j,i-1);
                    
                    //changing back into the domains frame of reference and finding velocity derivatives:
                    du1dx = (u1(j,i) - (uref+us))/dx;
                    dv1dx = (v1(j,i) - (vref+vs))/dx;


                    urel = u2(j,i-1) - us;
                    vrel = v2(j,i-1) - vs;
                    //calculating reflected velocity relative to the body's frame of referece
                    uref = urel - 2*(urel*nx(j,i-1) + vrel*ny(j,i-1))*nx(j,i-1);
                    vref = vrel - 2*(urel*nx(j,i-1) + vrel*ny(j,i-1))*ny(j,i-1);
                    
                    //changing back into the domains frame of reference and finding velocity derivatives:
                    du2dx = (u2(j,i) - (uref+us))/dx;
                    dv2dx = (v2(j,i) - (vref+vs))/dx;


                  }
                  
              }else if(nx(j,i)<0){//using the other derivative direction
                  drho1dx = (rho1(j,i+1) - rho1(j,i))/dx;
                  dp1dx = (p1(j,i+1) - p1(j,i))/dx;


                  drho2dx = (rho2(j,i+1) - rho2(j,i))/dx;
                  dp2dx = (p2(j,i+1) - p2(j,i))/dx;
                  dalpha1dx = (alpha1(j,i+1) - alpha1(j,i))/dx;
                //checking whether cell i+1,j is a ghost cell
                  if(phi(j,i+1)>=0){
                    
                    du1dx = (u1(j,i+1) - u1(j,i))/dx;
                    dv1dx = (v1(j,i+1) - v1(j,i))/dx;


                    du2dx = (u2(j,i+1) - u2(j,i))/dx;
                    dv2dx = (v2(j,i+1) - v2(j,i))/dx;
                  }else{
                    //tangential velocty is reflected:
                    //calculating velocity components:
                    
                    urel = u1(j,i+1) - us;
                    vrel = v1(j,i+1) - vs;
                    
                    uref = urel - 2*(urel*nx(j,i+1) + vrel*ny(j,i+1))*nx(j,i+1);
                    vref = vrel - 2*(urel*nx(j,i+1) + vrel*ny(j,i+1))*ny(j,i+1);
                    
                    //calculating derivatives in the fluid frame of reference:
                    
                    du1dx = (uref+us - u1(j,i))/dx;
                    dv1dx = (vref+vs - v1(j,i))/dx;
                    

                    urel = u2(j,i+1) - us;
                    vrel = v2(j,i+1) - vs;
                    
                    uref = urel - 2*(urel*nx(j,i+1) + vrel*ny(j,i+1))*nx(j,i+1);
                    vref = vrel - 2*(urel*nx(j,i+1) + vrel*ny(j,i+1))*ny(j,i+1);
                    
                    //calculating derivatives in the fluid frame of reference:
                    
                    du2dx = (uref+us - u2(j,i))/dx;
                    dv2dx = (vref+vs - v2(j,i))/dx;




                  }
              }else{
                  drho1dx = 0;
                  dp1dx = 0;
                  du1dx = 0;
                  dv1dx = 0;
                  dalpha1dx = 0;
                  du2dx = 0;
                  dv2dx = 0;
                  dp2dx = 0;
              }
              
              //interpolation is finished in the x direction, now interpolating in the y direction
              //more detailed comments are available in the previous section, we are now doing exactly as previously, but
              //in the y direction.
              if(ny(j,i)>0){ //calculating derivative direction
                  //scalar derivatives:
                  drho1dy = (rho1(j,i) - rho1(j-1,i))/dx;
                  dp1dy = (p1(j,i) - p1(j-1,i))/dx;

                  drho2dy = (rho2(j,i) - rho2(j-1,i))/dx;
                  dp2dy = (p2(j,i) - p2(j-1,i))/dx;
                  dalpha1dy = (alpha1(j,i) - alpha1(j-1,i))/dx;
                  
                  if(phi(j-1,i)>=0){//ie cell i,j-1 is a ghost cell
                    du1dy = (u1(j,i) - u1(j-1,i))/dx;
                    dv1dy = (v1(j,i) - v1(j-1,i))/dx;

                    du2dy = (u2(j,i) - u2(j-1,i))/dx;
                    dv2dy = (v2(j,i) - v2(j-1,i))/dx;


                  }else{ //not a ghost cell
                    //relative velocities:
                    urel = u1(j-1,i) - us;
                    vrel = v1(j-1,i) - vs;
                    
                    //reflected against the normal <nx,ny>
                    
                    uref = urel - 2*(urel*nx(j-1,i) + vrel*ny(j-1,i))*nx(j-1,i);
                    vref = vrel - 2*(urel*nx(j-1,i) + vrel*ny(j-1,i))*ny(j-1,i);
                    
                    //derivatives
                    
                    du1dy = (u1(j,i) - uref - us)/dx;
                    dv1dy = (v1(j,i) - vref - vs)/dx;

                    urel = u2(j-1,i) - us;
                    vrel = v2(j-1,i) - vs;
                    
                    //reflected against the normal <nx,ny>
                    
                    uref = urel - 2*(urel*nx(j-1,i) + vrel*ny(j-1,i))*nx(j-1,i);
                    vref = vrel - 2*(urel*nx(j-1,i) + vrel*ny(j-1,i))*ny(j-1,i);
                    
                    //derivatives
                    
                    du2dy = (u2(j,i) - uref - us)/dx;
                    dv2dy = (v2(j,i) - vref - vs)/dx;
                  }
                  
              }else if(ny(j,i)<0){
                  //scalar derivatives:
                  drho1dy = (rho1(j+1,i) - rho1(j,i))/dx;
                  dp1dy = (p1(j+1,i) - p1(j,i))/dx;

                  drho2dy = (rho2(j+1,i) - rho2(j,i))/dx;
                  dp2dy = (p2(j+1,i) - p2(j,i))/dx;
                  dalpha1dy = (alpha1(j+1,i) - alpha1(j,i))/dx;


                  if(phi(j+1,i)>=0){ //i,j+1 is in the ghost region
                    
                    du1dy = (u1(j+1,i) - u1(j,i))/dx;
                    dv1dy = (v1(j+1,i) - v1(j,i))/dx;

                    du2dy = (u2(j+1,i) - u2(j,i))/dx;
                    dv2dy = (v2(j+1,i) - v2(j,i))/dx;
                    
                  }else{ //i,j+1 is in the fluid
                    //relative velocities:
                    urel = u1(j+1,i) - us;
                    vrel = v1(j+1,i) - vs;
                    
                    //refleced relative velocities
                    
                    uref = urel - 2*(urel*nx(j+1,i) + vrel*ny(j+1,i))*nx(j+1,i);
                    vref = vrel - 2*(urel*nx(j+1,i) + vrel*ny(j+1,i))*ny(j+1,i);
                    
                    //changing frame of reference to fluid
                    
                    du1dy = (uref + us - u1(j,i))/dx;
                    dv1dy = (vref + vs - v1(j,i))/dx;


                    urel = u2(j+1,i) - us;
                    vrel = v2(j+1,i) - vs;
                    
                    //refleced relative velocities
                    
                    uref = (urel - 2*(urel*nx(j+1,i) + vrel*ny(j+1,i))*nx(j+1,i));
                    vref = vrel - 2*(urel*nx(j+1,i) + vrel*ny(j+1,i))*ny(j+1,i);
                    
                    //changing frame of reference to fluid
                    
                    du2dy = (uref + us - u2(j,i))/dx;
                    dv2dy = (vref + vs - v2(j,i))/dx;
                  }
                  
              }else{
                  //derivatives do not need to be interpolated
                  du1dy = 0;
                  dv1dy = 0;
                  dp1dy = 0;
                  drho1dy = 0;

                  du2dy = 0;
                  dv2dy = 0;
                  dp2dy = 0;
                  drho2dy = 0;

                  dalpha1dy = 0;
                  
              }
              
              //stepping the interpolation equation forward one imaginary time step XX
              
              u1_ext(j,i) = u1(j,i) - dtau*(du1dy*ny(j,i) + du1dx*nx(j,i));
              v1_ext(j,i) = v1(j,i) - dtau*(dv1dy*ny(j,i) + dv1dx*nx(j,i));
              p1_ext(j,i) = p1(j,i) - dtau*(dp1dy*ny(j,i) + dp1dx*nx(j,i));
              rho1_ext(j,i) = rho1(j,i) - dtau*(drho1dy*ny(j,i) + drho1dx*nx(j,i));

              u2_ext(j,i) = u2(j,i) - dtau*(du2dy*ny(j,i) + du2dx*nx(j,i));
              v2_ext(j,i) = v2(j,i) - dtau*(dv2dy*ny(j,i) + dv2dx*nx(j,i));
              p2_ext(j,i) = p2(j,i) - dtau*(dp2dy*ny(j,i) + dp2dx*nx(j,i));
              rho2_ext(j,i) = rho2(j,i) - dtau*(drho2dy*ny(j,i) + drho2dx*nx(j,i));

              alpha1_ext(j,i) = alpha1(j,i) - dtau*(dalpha1dy*ny(j,i) + dalpha1dx*nx(j,i));
              //cout<<rho_ext(i,j)<<endl;
                } //if phi ...
                
            }//for j = ...
            
        } //for i = ...
        
        
        u1 = u1_ext;
        v1 = v1_ext;
        p1 = p1_ext;
        rho1 = rho1_ext;

        u2 = u2_ext;
        v2 = v2_ext;
        p2 = p2_ext;
        rho2 = rho2_ext;

        //alpha1 = alpha1_ext;
  
    }   //for tau = ...   
    
  
    //placing primitives in vector of conserved variables:
    U.slice(0) = alpha1;
    U.slice(1) = rho1;
    U.slice(2) = u1;
    U.slice(3) = v1;
    U.slice(4) = p1;
    U.slice(5) = rho2;
    U.slice(6) = u2;
    U.slice(7) = v2;
    U.slice(8) = p2;
    

    //sleep(1);
    cout<<t+dt<<endl;
    cout<<std::max(abs(p1).max(),abs(p2).max())<<endl;
    //sleep(2);
    t=t+dt;
    save_t_accu += dt;
    
    if(save_t_accu>save_t_inc){
      save_t_accu -= save_t_inc;
      std::string s;
      std::ostringstream SAVENAME1;
      std::ostringstream SAVENAME2;
      std::ostringstream SAVENAME3;
      std::ostringstream SAVENAME4;
      std::ostringstream SAVENAME5;
      std::ostringstream SAVENAME6;
      std::ostringstream SAVENAME7;
      std::ostringstream SAVENAME8;
      std::ostringstream SAVENAME9;

        SAVENAME1 << "alpha1_" << save_number << ".dat";
        s = SAVENAME1.str();

        alpha1.save(s,raw_ascii);

        //SAVENAME.clear();

        SAVENAME2 << "rho1_" << save_number << ".dat";
        s = SAVENAME2.str();

        rho1.save(s,raw_ascii);
        //SAVENAME.clear();


        SAVENAME3 << "u1_" << save_number << ".dat";
        s = SAVENAME3.str();

        u1.save(s,raw_ascii);

        //SAVENAME.clear();

        SAVENAME4 << "v1_" << save_number << ".dat";
        s = SAVENAME4.str();

        v1.save(s,raw_ascii);
        //SAVENAME.clear();


        SAVENAME5 << "p1_" << save_number << ".dat";
        s = SAVENAME5.str();

        p1.save(s,raw_ascii);


        //SAVENAME.clear();
        SAVENAME6 << "rho2_" << save_number << ".dat";
        s = SAVENAME6.str();

        rho2.save(s,raw_ascii);
        //SAVENAME.clear();
        SAVENAME7 << "u2_" << save_number << ".dat";
        s = SAVENAME7.str();

        u2.save(s,raw_ascii);
       // SAVENAME.clear();
        SAVENAME8 << "v2_" << save_number << ".dat";
        s = SAVENAME8.str();

        v2.save(s,raw_ascii);
        //SAVENAME9.clear();
        SAVENAME9 << "p2_" << save_number << ".dat";
        s = SAVENAME9.str();

        p2.save(s,raw_ascii);

        cout<<"saving"<<endl;
        save_number++;
    
    }
    
    duration1 = (std::clock() - start1)/(double) CLOCKS_PER_SEC;

    //cout<<duration1<<endl<<duration<<endl;
    //cout<<trans(p)<<endl;
    //sleep(10);


  }
  //cout<<phi<<endl;
  //sleep(5);

  alpha1 = U.slice(0);
  rho1 = U.slice(1);
  u1 = U.slice(2);
  v1 = U.slice(3);
  p1 = U.slice(4);
  rho2 = U.slice(5);
  u2 = U.slice(6);
  v2 = U.slice(7);
  p2 = U.slice(8);
  mat p = alpha1%p1 + (1-alpha1)%p2;
  mat rho = alpha1%rho1 + (1-alpha1)%rho2;
/*  vec rho_si = trans(rho1.row(5));
  vec p_si = trans(p1.row(5));
  vec u_si = trans(v1.row(5));
  vec e_si = trans(p1.row(5)/(rho1.row(5)*(gamma-1)));
  vec phi_si = phi.col(5);
  //cout<< rho_si<<endl;
  //p = p;//%(1-sign(phi))/2;
  rho_si.save("rho2.dat",raw_ascii);
  p1.save("p2.dat",raw_ascii);
  p_si.save("prow.dat",raw_ascii);
  u_si.save("u2.dat",raw_ascii);
  e_si.save("e2.dat",raw_ascii);
  phi_si.save("phi_si.dat",raw_ascii);*/
  phi.save("phi.dat",raw_ascii);
  alpha1.save("alpha1.dat",raw_ascii);
  rho1.save("rho1.dat",raw_ascii);
  u1.save("u1.dat",raw_ascii);
  v1.save("v1.dat",raw_ascii);
  p1.save("p1.dat",raw_ascii);
  rho2.save("rho2.dat",raw_ascii);
  u2.save("u2.dat",raw_ascii);
  v2.save("v2.dat",raw_ascii);
  p2.save("p2.dat",raw_ascii);
  return 0;
}