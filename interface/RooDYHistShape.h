#ifndef ROO_HISTSHAPE
#define ROO_HISTSHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAICRegistry.h"

// class RooRealProxy;
// class RooListProxy;

class RooDYHistShape : public RooAbsPdf {
public:
  RooDYHistShape();
  RooDYHistShape(const char *name, const char *title,
	       RooAbsReal& x,
	       const RooArgList& coefList,
	       int histType,	// 0: smeared DY in 0.1 GeV binning. 1: smeared DY in 0.2 GeV binning. 2: smeared DY with alternative smearing function from Giacomo, 0.1 GeV binning. 6: smeared BW in 0.5 GeV binning. 7: smeared DY in 0.5 GeV binning.
	       int EPFType, // 0 for exponentiated polynomial, 1 for poly nominal, 2 for power law
	       int polyType	// 0 for bernstein, 1 for chebychev, 2 for polynomial
	     );
  RooDYHistShape(const RooDYHistShape& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDYHistShape(*this,newname); }
  inline virtual ~RooDYHistShape() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  virtual Bool_t isBinnedDistribution(const RooArgSet&) const {return true;}
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ;
  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;
protected:

  RooRealProxy _x ;
  RooListProxy _coefList ;
  int _histType;
  int _EPFType;
  int _polyType;
  
  Double_t evaluate() const ;
private:

  Double_t _xMin;
  Double_t _xMax;
  mutable TNamed* _refRangeName ; 
  mutable std::vector<Double_t> _wksp; //! do not persist
  
  Double_t bernstein(Double_t x) const ;
  Double_t polynomial(Double_t x) const ;
  Double_t chebychev(Double_t x) const ;
  Double_t empiricalFunc(Double_t x) const ;
  
  void setupHist();
  inline static double p1(double t,double a,double b) { return a*t+b; }
  inline static double p2(double t,double a,double b,double c) { return p1(t,p1(t,a,b),c); }
  inline static double p3(double t,double a,double b,double c,double d) { return p2(t,p1(t,a,b),c,d); }

  std::vector<double> histContent[10];
  double binWidth[10];
  ClassDef(RooDYHistShape,4) // VoigtianExt PDF (Gauss (x) BreitWigner)
};

#endif

