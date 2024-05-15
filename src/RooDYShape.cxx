/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   TS, Thomas Schietinger, SLAC,           schieti@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
 
/** \class RooDYShape
    \ingroup Roofit
 
    RooDYShape is an efficient implementation of the convolution of a
    Breit-Wigner with a Gaussian, making use of the complex error function.
    RooFitCore provides two algorithms for the evaluation of the complex error
    function (the default CERNlib C335 algorithm, and a faster, look-up-table
    based method). By default, RooDYShape employs the default (CERNlib)
    algorithm. Select the faster algorithm either in the constructor, or with
    the selectFastAlgorithm() method.
**/
 
#include <cmath>
#include <complex>
 
#include "RooFit.h"
 
#include "Riostream.h"

#include "TMath.h"

#include "RooDYShape.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "RooArgList.h"
 
#include "TError.h"
using namespace std;
 
ClassImp(RooDYShape);
 
////////////////////////////////////////////////////////////////////////////////
 
RooDYShape::RooDYShape(const char *name, const char *title,
           RooAbsReal& x, RooAbsReal& mean,
           RooAbsReal& width, RooAbsReal& power,
           const RooArgList& coefList,
           int BWType,
           int EPFType,
           int polyType) :
  RooAbsPdf(name,title),
  _x("x","Dependent",this,x),
  _mean("mean","Mean",this,mean),
  _width("width","Breit-Wigner Width",this,width),
  _power("power","Power in CMS function",this,power),
  _coefList("coefficients","List of coefficients",this),
  _BWType(BWType),
  _EPFType(EPFType),
  _polyType(polyType)
{
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << __PRETTY_FUNCTION__<<" ERROR: coefficient " << coef->GetName()
     << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;
  _invRootPi= 1./sqrt(atan2(0.,-1.));
}
 
////////////////////////////////////////////////////////////////////////////////
 
RooDYShape::RooDYShape(const RooDYShape& other, const char* name) :
  RooAbsPdf(other,name),
  _x("x",this,other._x),
  _mean("mean",this,other._mean),
  _width("width",this,other._width),
  _power("power",this,other._power),
  _coefList("coefList",this,other._coefList),
  _BWType(other._BWType),
  _EPFType(other._EPFType),
  _polyType(other._polyType)
{
  _invRootPi= 1./sqrt(atan2(0.,-1.));
}
 
////////////////////////////////////////////////////////////////////////////////
 
Double_t RooDYShape::evaluate() const
{
  double poly=0;
  if(_polyType==0) poly=bernstein();
  else if(_polyType==1) poly=chebychev();
  else poly=polynomial();

  double add=1;
  if(_EPFType==0) add=exp(poly);
  else if(_EPFType==1) add=poly;
  else if(_EPFType==2) add=pow(_x,poly);
  add=std::max(add, 1e-32);

  if(_BWType<0){
    Double_t w = (_width>0) ? _width : -_width ;
    Double_t p = _power;
    Double_t arg = _x - _mean;
    return (add/(pow(arg,p)+pow(0.5*w,p)));   // Simple Breit-Wigner times function
  }

  // Physics constants
  double alphaEM=7.297352e-3, Nc=3, Qd=-1./3., Qu=2./3., Qs=-1./3., sin2W=0.23126, GF=1.1663787e-5;
  double gAmu=-0.5, gAu=0.5, gAd=-0.5, gAs=-0.5, gVmu=-0.5+2*sin2W, gVu=0.5-4./3.*sin2W, gVd=-0.5+2./3.*sin2W, gVs=-0.5+2./3.*sin2W;

  // Useful variables
  double s=_x*_x, alphaEM2=alphaEM*alphaEM, kappa=sqrt(2)*GF*_mean*_mean/(4*TMath::Pi()*alphaEM);
  double S_Zgam=kappa*s*(s-_mean*_mean)*BW(), S_Z=kappa*kappa*s*s*BW();
  
  double ME_ddbar=4*TMath::Pi()*alphaEM2/(3*Nc*s)*((Qd*Qd)-2*Qd*gVmu*gVd*S_Zgam+(gAmu*gAmu+gVmu*gVmu)*(gAd*gAd+gVd*gVd)*S_Z);
  double ME_uubar=4*TMath::Pi()*alphaEM2/(3*Nc*s)*((Qu*Qu)-2*Qu*gVmu*gVu*S_Zgam+(gAmu*gAmu+gVmu*gVmu)*(gAu*gAu+gVu*gVu)*S_Z);
  double ME_ssbar=4*TMath::Pi()*alphaEM2/(3*Nc*s)*((Qs*Qs)-2*Qs*gVmu*gVs*S_Zgam+(gAmu*gAmu+gVmu*gVmu)*(gAs*gAs+gVs*gVs)*S_Z);

  
  return add*(ME_ddbar*partonLumi_ddbar()+ME_uubar*partonLumi_uubar()+ME_ssbar*partonLumi_ssbar())*_x;
}

Double_t RooDYShape::bernstein() const {
  Double_t xmin = _x.min(), xmax = _x.max();
  Double_t x_norm = (_x - xmin) / (xmax - xmin); // rescale to [0,1]
  Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n
  RooFIter iter = _coefList.fwdIterator();
 
  if(degree == 0) {
 
    return ((RooAbsReal *)iter.next())->getVal();
 
  } else if(degree == 1) {
 
    Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); // c0
    Double_t a1 = ((RooAbsReal *)iter.next())->getVal() - a0; // c1 - c0
    return a1 * x_norm + a0;
 
  } else if(degree == 2) {
 
    Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); // c0
    Double_t a1 = 2 * (((RooAbsReal *)iter.next())->getVal() - a0); // 2 * (c1 - c0)
    Double_t a2 = ((RooAbsReal *)iter.next())->getVal() - a1 - a0; // c0 - 2 * c1 + c2
    return (a2 * x_norm + a1) * x_norm + a0;
 
  } else if(degree > 2) {
 
    Double_t t = x_norm;
    Double_t s = 1 - x_norm;
    
    Double_t result = ((RooAbsReal *)iter.next())->getVal() * s;
    for(Int_t i = 1; i < degree; i++) {
      result = (result + t * TMath::Binomial(degree, i) * ((RooAbsReal *)iter.next())->getVal()) * s;
      t *= x_norm;
    }
    result += t * ((RooAbsReal *)iter.next())->getVal();
 
    return result;
  }
 
  // in case list of arguments passed is empty
  return TMath::SignalingNaN();
}


Double_t RooDYShape::polynomial() const{
  // Calculate and return value of polynomial
  const unsigned sz = _coefList.getSize();
  const int lowestOrder = 0;  // Lowest order contains a constant
  if (!sz) return lowestOrder ? 1. : 0.;
  _wksp.clear();
  _wksp.reserve(sz);
  {
    const RooArgSet* nset = _coefList.nset();
    RooFIter it = _coefList.fwdIterator();
    RooAbsReal* c;
    while ((c = (RooAbsReal*) it.next())) _wksp.push_back(c->getVal(nset));
  }
  const Double_t x = _x;
  Double_t retVal = _wksp[sz - 1];
  for (unsigned i = sz - 1; i--; ) retVal = _wksp[i] + x * retVal;
  return retVal * std::pow(x, lowestOrder) + (lowestOrder ? 1.0 : 0.0);
}

Double_t RooDYShape::PL() const {
  // Parameterization of qqbar parton luminosity
  double x_norm=_x/13000.;  // Normalize mgg by sqrt(s)
  double L = pow(1-pow(x_norm,1./3),8.0919)*pow(x_norm,-2.26565)*(8.43489e-8);
  if(!std::isfinite(L)) L=0;
  return L;
}

Double_t RooDYShape::BW() const {
  double s=_x*_x;
  if(_BWType==0) return 1/((s-_mean*_mean)*(s-_mean*_mean)+(_mean*_mean*_width*_width));
  if(_BWType==1) return 1/((s-_mean*_mean)*(s-_mean*_mean)+(s*s/(_mean*_mean)*_width*_width));
  else return 1/((_x-_mean)*(_x-_mean)+(_width/2.)*(_width/2.));
}

double RooDYShape::partonLumi_uubar() const {
  double p0                        =   0.00703528;
  double p1                        = -0.000174816;
  double p2                        =  1.89616e-06;
  double p3                        = -1.12038e-08;
  double p4                        =   3.7537e-11;
  double p5                        = -6.71198e-14;  
  double p6                        =  4.98316e-17;
  
  return p0+p1*_x+p2*_x*_x+p3*_x*_x*_x+p4*_x*_x*_x*_x+p5*_x*_x*_x*_x*_x+p6*_x*_x*_x*_x*_x*_x;
}

double RooDYShape::partonLumi_ddbar() const{
  double p0                        =    0.0052163;
  double p1                        = -0.000130677;
  double p2                        =  1.42682e-06;
  double p3                        = -8.48218e-09;
  double p4                        =  2.85847e-11;
  double p5                        = -5.14001e-14;
  double p6                        =  3.83665e-17;
  
  return p0+p1*_x+p2*_x*_x+p3*_x*_x*_x+p4*_x*_x*_x*_x+p5*_x*_x*_x*_x*_x+p6*_x*_x*_x*_x*_x*_x;
}

double RooDYShape::partonLumi_ssbar() const{
  double p0                        =   0.00227443;
  double p1                        = -5.80923e-05;
  double p2                        =  6.40022e-07;
  double p3                        = -3.81913e-09;
  double p4                        =  1.28814e-11;
  double p5                        =  -2.3144e-14;
  double p6                        =  1.72444e-17;
  
  return p0+p1*_x+p2*_x*_x+p3*_x*_x*_x+p4*_x*_x*_x*_x+p5*_x*_x*_x*_x*_x+p6*_x*_x*_x*_x*_x*_x;
}

Double_t RooDYShape::chebychev() const
{
  Double_t xmin = _x.min() ; Double_t xmax = _x.max();
  Double_t x(-1+2*(_x-xmin)/(xmax-xmin));
  Double_t x2(x*x);
  Double_t sum(0) ;
  switch (_coefList.getSize()) {
  case  7: sum+=((RooAbsReal&)_coefList[6]).getVal()*x*p3(x2,64,-112,56,-7);
  case  6: sum+=((RooAbsReal&)_coefList[5]).getVal()*p3(x2,32,-48,18,-1);
  case  5: sum+=((RooAbsReal&)_coefList[4]).getVal()*x*p2(x2,16,-20,5);
  case  4: sum+=((RooAbsReal&)_coefList[3]).getVal()*p2(x2,8,-8,1);
  case  3: sum+=((RooAbsReal&)_coefList[2]).getVal()*x*p1(x2,4,-3);
  case  2: sum+=((RooAbsReal&)_coefList[1]).getVal()*p1(x2,2,-1);
  case  1: sum+=((RooAbsReal&)_coefList[0]).getVal()*x;
  case  0: sum+=1; break;
  default: std::cerr << "In " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ", line " <<
      __LINE__ << "): Higher order Chebychev polynomials currently "
      "unimplemented." << std::endl;
    R__ASSERT(false);
  }
  return sum;
}
