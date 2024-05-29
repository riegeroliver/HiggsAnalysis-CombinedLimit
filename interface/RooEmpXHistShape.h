#ifndef ROO_EMPXHISTSHAPE
#define ROO_EMPXHISTSHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAICRegistry.h"

#include "TDOMParser.h"
#include "TXMLNode.h"
#include "TXMLAttr.h"

// class RooRealProxy;
class RooListProxy;

class RooEmpXHistShape : public RooAbsPdf {

public:
  enum EPFType{Epoly, Poly, Power};
  enum PolyType{Bern, Cheby, Plain};

  RooEmpXHistShape();
  RooEmpXHistShape(const char *name, const char *title,
		   RooAbsReal& x,
		   const RooArgList& coefList,
		   const char* inputXMLFileName,
		   const char* histType,	// explained in input XML file
		   int EPFType, // 0 for exponentiated polynomial, 1 for poly nominal, 2 for power law
		   int polyType	// 0 for bernstein, 1 for chebychev, 2 for polynomial
		   );
  RooEmpXHistShape(const RooEmpXHistShape& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooEmpXHistShape(*this,newname); }
  inline virtual ~RooEmpXHistShape() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  virtual Bool_t isBinnedDistribution(const RooArgSet&) const {return true;}
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ;
  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;

  void DumpHist() const;			// Dump hist info out
protected:

  RooRealProxy _x ;
  RooListProxy _coefList ;
  TString _histType;
  int _EPFType;
  int _polyType;
  mutable TNamed* _refRangeName ;
  
  Double_t evaluate() const ;

private:

  Double_t _xMin;
  Double_t _xMax;
  std::vector<Double_t> _histContent;
  Double_t _binWidth;
  
  mutable std::vector<Double_t> _wksp; //! do not persist
  
  Double_t bernstein(Double_t x) const ;
  Double_t polynomial(Double_t x) const ;
  Double_t chebychev(Double_t x) const ;
  Double_t empiricalFunc(Double_t x) const ;
  
  void setupHist(const char* inputXMLFileName);
  inline static double p1(double t,double a,double b) { return a*t+b; }
  inline static double p2(double t,double a,double b,double c) { return p1(t,p1(t,a,b),c); }
  inline static double p3(double t,double a,double b,double c,double d) { return p2(t,p1(t,a,b),c,d); }

  TString getAttributeValue(TXMLNode* rootNode, TString attributeKey) const;

  static constexpr double epsilon = 1e-8;
  ClassDef(RooEmpXHistShape,4) 
};

#endif

