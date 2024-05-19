#ifndef ROO_MULTIPDF
#define ROO_MULTIPDF

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooAbsCategory.h"
#include "RooCategory.h"
#include "RooCategoryProxy.h"
#include "RooArgProxy.h"
#include "RooAbsProxy.h"
#include "RooFormulaVar.h"
#include "RooLinkedList.h"
#include "RooConstVar.h"
#include "RooListProxy.h"

#include <iostream>
#include <vector>

using namespace std;

class RooMultiPdf : public RooAbsPdf {
public:
  enum PenaltyScheme { PVAL, AIC };

  RooMultiPdf(const char *name, const char *title, RooCategory &, const RooArgList &_c);
  RooMultiPdf(const RooMultiPdf &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const { return new RooMultiPdf(*this, newname); }
  inline virtual ~RooMultiPdf() {}

  bool checkIndexDirty() const;
  double getCorrection() const;
  RooAbsPdf *getCurrentPdf() const;
  int getNumPdfs() const { return nPdfs; }
  void setCorrectionFactor(PenaltyScheme penal);
  void setCorrectionFactor(double penal);
  int getCurrentIndex() const;
  RooAbsPdf *getPdf(int index) const;
  virtual Double_t getValV(const RooArgSet *nset) const;
  virtual void printValue(ostream &os) const { getCurrentPdf()->printValue(os); }

protected:
  RooListProxy c;
  RooListProxy corr;
  RooCategoryProxy x;

  int fIndex;
  int nPdfs;
  mutable Int_t _oldIndex;

  Double_t evaluate() const;
  Double_t getLogVal(const RooArgSet *set = 0) const;
  double cFactor;

private:
  ClassDef(RooMultiPdf, 1) // Multi PDF
};
#endif
