#include "RooFit.h"
#include "Riostream.h"
#include <math.h>

#include "../interface/RooMultiPdf.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include <stdexcept>

ClassImp(RooMultiPdf)

//_____________________________________________________________________________
RooMultiPdf::RooMultiPdf(const char *name, const char *title, RooCategory &_x, const RooArgList &_c) 
    : RooAbsPdf(name, title), 
      c("c", "The list of pdfs", this), 
      corr("corr", "The correction factors", this), 
      x("x", "the pdf index", this, _x) {

    int count = 0;
    for (RooAbsArg *fPdf : _c) {
        c.add(*fPdf);
        _x.defineType(Form("_pdf%d", count), count);
        std::unique_ptr<RooArgSet> variables(fPdf->getVariables());
        std::unique_ptr<RooAbsCollection> nonConstVariables(variables->selectByAttrib("Constant", false));
        RooConstVar *tmp = new RooConstVar(Form("const%s", fPdf->GetName()), "", nonConstVariables->getSize());
        corr.add(*tmp);
        count++;
    }
    nPdfs = c.getSize();
    cFactor = 0.5;
    _oldIndex = fIndex;
}

//_____________________________________________________________________________
RooMultiPdf::RooMultiPdf(const RooMultiPdf &other, const char *name)
    : RooAbsPdf(other, name), 
      c("c", this, other.c), 
      corr("corr", this, other.corr), 
      x("x", this, other.x) {

    fIndex = other.fIndex;
    nPdfs = other.nPdfs;
    cFactor = other.cFactor;
    _oldIndex = fIndex;
}

bool RooMultiPdf::checkIndexDirty() const {
    return _oldIndex != x;
}

//_____________________________________________________________________________
void RooMultiPdf::setCorrectionFactor(PenaltyScheme penal) {
    if (penal == AIC) {
        cFactor = 1.0;
    } else if (penal == PVAL) {
        cFactor = 0.5;
    }
}

//_____________________________________________________________________________
void RooMultiPdf::setCorrectionFactor(double penal) {
    cFactor = penal;
}

//_____________________________________________________________________________
double RooMultiPdf::getCorrection() const {
    double val = ((RooAbsReal *)corr.at(x))->getVal();
    return cFactor * val;
}

//_____________________________________________________________________________
RooAbsPdf *RooMultiPdf::getCurrentPdf() const {
    RooAbsPdf *cPdf = ((RooAbsPdf *)c.at(x));
    return cPdf;
}

RooAbsPdf *RooMultiPdf::getPdf(int index) const {
    RooAbsPdf *cPdf = ((RooAbsPdf *)c.at(index));
    return cPdf;
}

int RooMultiPdf::getCurrentIndex() const {
    Int_t index = x;
    return index;
}

//_____________________________________________________________________________
Double_t RooMultiPdf::getValV(const RooArgSet *nset) const {
    RooAbsPdf *cPdf = ((RooAbsPdf *)c.at(x));
    double val = cPdf->getVal(nset);
    _oldIndex = x;
    return val;
}

//_____________________________________________________________________________
Double_t RooMultiPdf::evaluate() const {
    throw std::invalid_argument("RooMultiPdf::evaluate() called\n");
}

//_____________________________________________________________________________
Double_t RooMultiPdf::getLogVal(const RooArgSet *nset) const {
    RooAbsPdf *cPdf = ((RooAbsPdf *)c.at(x));
    double logval = cPdf->getLogVal(nset);
    _oldIndex = x;
    return logval;
}
