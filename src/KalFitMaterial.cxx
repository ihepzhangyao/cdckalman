//-----------------------------------------------------------------------
// File from KalFit module
//
// Filename : KalFitMaterial.cc
//------------------------------------------------------------------------
// Description : 
// Material is a class which describes the properties of a given material, 
// for instance atomic number, atomic weight and so on.
//------------------------------------------------------------------------
// Modif :
//------------------------------------------------------------------------
#include <math.h>
#include <iostream>
#include "KalFitAlg/KalFitMaterial.h"

using namespace std;

KalFitMaterial::KalFitMaterial(double z, double a, double i,
			   double rho, double x0)
  : x0_(x0), z_(z)  // rho is the density, z is the atomic number, a is the weight
                    // i is mean excitation potention, x0 is the radiation length
{
  rza_ = rho * z / a;
  isq_  = i * i * 1e-18;
}

KalFitMaterial::KalFitMaterial(const KalFitMaterial& mat)
  : rza_(mat.rza_), isq_(mat.isq_),
    x0_(mat.x0_), z_(mat.z_)
{
   
}

// calculate dE
double KalFitMaterial::dE(double mass, double path,
			double p) const
{
  if (!(p>0))
    return 0;
 
  //cout<<"this material:x0 "<< x0_ << " Z  " << z_ << endl
  //    <<"  rho*Z/A  "<< rza_ << " I^2 "<< isq_ << endl;

  
  const double Me = 0.000510999;
  double psq = p * p;
  double bsq = psq / (psq + mass * mass);
  double esq = psq / (mass * mass);

  double s = Me / mass;
  double w = (4 * Me * esq
	      / (1 + 2 * s * sqrt(1 + esq)
		 + s * s));

  // Density correction :
  double cc, x0;
  cc = 1+2*log(sqrt(isq_)/(28.8E-09*sqrt(rza_)));
  if (cc < 5.215) 
    x0 = 0.2;
  else
    x0 = 0.326*cc-1.5;
  double x1(3), xa(cc/4.606), aa;
  aa = 4.606*(xa-x0)/((x1-x0)*(x1-x0)*(x1-x0));
  double delta(0);
  double x(log10(sqrt(esq)));
  if (x > x0){    
    delta = 4.606*x-cc;
    if (x < x1) delta=delta+aa*(x1-x)*(x1-x)*(x1-x);
  }

 // Shell correction :  
  float f1, f2, f3, f4, f5, ce;
  f1 = 1/esq;
  f2 = f1*f1;
  f3 = f1*f2;
  f4 = (f1*0.42237+f2*0.0304-f3*0.00038)*1E12;
  f5 = (f1*3.858-f2*0.1668+f3*0.00158)*1E18;
  ce = f4*isq_+f5*isq_*sqrt(isq_);

  return (0.0001535 * rza_ / bsq
	  * (log(Me * esq * w / isq_)
	     - 2 * bsq-delta-2.0*ce/z_)) * path;
}


// calculate Multiple Scattering angle
double KalFitMaterial::mcs_angle(double mass, double path,
			       double p) const
{
  //cout<<"this material:x0 "<< x0_ << " Z  " << z_ << endl
  //    <<"  rho*Z/A  "<< rza_ << " I^2 "<< isq_ << endl;
  double t = path / x0_;
  double psq = p*p;
  return 0.0136 * sqrt(t * (mass*mass + psq)) / psq
    * (1 + 0.038 * log(t));
}

// calculate the RMS of straggling for the energy loss (del_E) :
double KalFitMaterial::del_E(double mass, double path,
                           double p) const
{
  double sigma0_2 = 0.1569*rza_*path;

  if (sigma0_2<0) return 0;

  double psq = p * p;
  double bsq = psq / (psq + mass * mass);

  // Correction for relativistic particles :
  double sigma_2 = sigma0_2*(1-0.5*bsq)/(1-bsq);

  if (sigma_2<0) return 0;

  // Return sigma in GeV !!
  return sqrt(sigma_2)*0.001;

}
