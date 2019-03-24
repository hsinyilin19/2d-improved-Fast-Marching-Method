//travel_time_marcher.cpp

#include "travel_time_marcher.h"
#include "math.h"
#include "heap.h"
#include <stdexcept>

void travelTimeMarcher::initalizeFrozen()
{
  distanceMarcher::initalizeFrozen();
  for (int i=0; i<size_; i++)
  {
    if (flag_[i]==Frozen)
    {
      // convert distance to time
      distance_[i]=fabs(distance_[i]/speed_[i]);
    }
  }
}


double travelTimeMarcher::updatePointOrderOne(int i)
{
  double a,b,c;
  double tmp;
  tmp = 0;
  a=b=c=0;
  double u_small=0;
  int naddr=0;
  if (dim_ == 2)
  {
    double *u_vec = new double[2];
    for (int dim=0; dim<dim_; dim++)
    {
      double value = maxDouble;
      for (int j=-1; j<2; j+=2) // each direction
      {
        naddr = _getN(i,dim,j,Mask);
        if (naddr!=-1 && flag_[naddr]==Frozen)
        {
          if (fabs(distance_[naddr])<fabs(value))
          {
            value = distance_[naddr];
          }
        }
      }

      if (value<maxDouble)
      {
        u_vec[dim] = value;
        a+=idx2_[dim];
        b-=idx2_[dim]*2*value;
        c+=idx2_[dim]*pow(value,2);
      }
      else
      {
        u_vec[dim] = maxDouble;
      }
    
    }

    double u_small=u_vec[0], u_big=u_vec[1];
    if (u_small > u_big)
    {
      double tem = u_small;
      double u_small = u_big;
      double u_big = tem;
    }

    
    tmp = solveQuadratic(i,a,b,c,u_small);
     //cout<<"dim=2 and positive discriminant!!!\n"; 
    delete[] u_vec;  
    
  }
  else{
    for (int dim=0; dim<dim_; dim++)
    {
      double value = maxDouble;
      for (int j=-1; j<2; j+=2) // each direction
      {
        naddr = _getN(i,dim,j,Mask);
        if (naddr!=-1 && flag_[naddr]==Frozen)
        {
          if (fabs(distance_[naddr])<fabs(value))
          {
            value = distance_[naddr];
          }
        }
      }
      if (value<maxDouble)
      {
        a+=idx2_[dim];
        b-=idx2_[dim]*2*value;
        c+=idx2_[dim]*pow(value,2);
      }
    }
    tmp = solveQuadratic(i,a,b,c,u_small);
    //cout<<"dim !=2 and positive discriminant!!!\n"; 
  }
  //double tmp = solveQuadratic(i,a,b,c);
  return tmp;
}





// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double travelTimeMarcher::updatePointOrderTwo(int i)
{
  double a,b,c;
  a=b=c=0;
  int naddr=0;
  double u_small=0;
  for (int dim=0; dim<dim_; dim++)
  {
    double value1 = maxDouble;
    double value2 = maxDouble;
    for (int j=-1; j<2; j+=2) // each direction
    {
      naddr = _getN(i,dim,j,Mask);
      if (naddr!=-1 && flag_[naddr]==Frozen)
      {
        if (fabs(distance_[naddr])<fabs(value1))
        {
          value1 = distance_[naddr];
          int naddr2 = _getN(i,dim,j*2,Mask);
          if (naddr2!=-1 &&
              flag_[naddr2]==Frozen &&
              ((distance_[naddr2]<=value1 && value1 >=0) ||
               (distance_[naddr2]>=value1 && value1 <=0)))
          {
            value2=distance_[naddr2];
            if (phi_[naddr2] * phi_[naddr] < 0  || phi_[naddr2] * phi_[i] < 0)
              value2 *= -1;
          }
        }
      }
    }
    if (value2<maxDouble)
    {
      double tp = oneThird*(4*value1-value2);
      a+=idx2_[dim]*aa;
      b-=idx2_[dim]*2*aa*tp;
      c+=idx2_[dim]*aa*pow(tp,2);
    }
    else if (value1<maxDouble)
    {
      a+=idx2_[dim];
      b-=idx2_[dim]*2*value1;
      c+=idx2_[dim]*pow(value1,2);
    }
  }
  return solveQuadratic(i,a,b,c,u_small);
}


double travelTimeMarcher::solveQuadratic(int i, const double &a,
                                         const double &b,
                                         double &c, double &u_small )
{
  c -= 1/pow(speed_[i],2);
  double r0=0;
  double det = pow(b,2)-4*a*c;
  if (det>=0)
  {
    r0 = (-b+sqrt(det))/2.0/a;
  }
  else
  {
    //throw std::runtime_error("Negative discriminant in time marcher quadratic.");
    r0 = u_small + 1/sqrt(idx2_[0]);    
  }
  return r0;
}
