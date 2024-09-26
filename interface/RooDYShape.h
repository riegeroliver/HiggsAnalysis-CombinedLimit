/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooVoigtian.h,v 1.7 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   TS, Thomas Schietinger, SLAC,          schieti@slac.stanford.edu        *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_DYSHAPE
#define ROO_DYSHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;

class RooDYShape : public RooAbsPdf {
public:
  RooDYShape() {} ;
  RooDYShape(const char *name, const char *title,
     RooAbsReal& x, RooAbsReal& mean,
     RooAbsReal& width, RooAbsReal& power,
     const RooArgList& coefList,
     int BWType,  // 0 for relativistic BW with running scale, 1 for relativistic BW with fixed scale, 2 for non-relativistic, negative number for simple BWxfunc
     int EPFType, // 0 for exponentiated polynomial, 1 for poly nominal, 2 for power law
     int polyType // 0 for bernstein, 1 for chebychev, 2 for polynomial
       );
  RooDYShape(const RooDYShape& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDYShape(*this,newname); }
  inline virtual ~RooDYShape() { }

protected:

  RooRealProxy _x ;
  RooRealProxy _mean ;
  RooRealProxy _width ;
  RooRealProxy _power ;
  RooListProxy _coefList ;
  int _BWType;
  int _EPFType;
  int _polyType;
  
  Double_t evaluate() const ;
private:

  Double_t _invRootPi;

  mutable std::vector<Double_t> _wksp; //! do not persist
  
  Double_t bernstein() const ;
  Double_t polynomial() const ;
  Double_t chebychev() const ;
  Double_t asymptote() const ;
  
  inline static double p1(double t,double a,double b) { return a*t+b; }
  inline static double p2(double t,double a,double b,double c) { return p1(t,p1(t,a,b),c); }
  inline static double p3(double t,double a,double b,double c,double d) { return p2(t,p1(t,a,b),c,d); }
  
  Double_t BW() const ;
  Double_t PL() const ;
  Double_t partonLumi_uubar() const ;
  Double_t partonLumi_ddbar() const ;
  Double_t partonLumi_ssbar() const ;
  ClassDef(RooDYShape,1) // VoigtianExt PDF (Gauss (x) BreitWigner)
};

#endif

