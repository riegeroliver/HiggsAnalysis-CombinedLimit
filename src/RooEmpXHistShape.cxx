#include <cmath>
#include <complex>
 
#include "RooFit.h"
 
#include "Riostream.h"

#include "TMath.h"

#include "RooEmpXHistShape.h"
#include "RooMsgService.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "RooArgList.h"
 
#include "TError.h"

ClassImp(RooEmpXHistShape);

RooEmpXHistShape::RooEmpXHistShape() : _refRangeName(0)
{
}

////////////////////////////////////////////////////////////////////////////////
 
RooEmpXHistShape::RooEmpXHistShape(const char *name, const char *title,
           RooAbsReal& x,
           const RooArgList& coefList,
           const char* inputXMLFileName,
           const char* histType,
           int EPFType,
           int polyType) :
  RooAbsPdf(name,title),
  _x("x","Dependent",this,x),
  _coefList("coefficients","List of coefficients",this),
  _histType(histType),
  _EPFType(EPFType),
  _polyType(polyType),
  _refRangeName(0)
{
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cout << __PRETTY_FUNCTION__<<" ERROR: coefficient " << coef->GetName()
    << " is not of type RooAbsReal" << std::endl ;
      R__ASSERT(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;
  setupHist(inputXMLFileName);
}
 
////////////////////////////////////////////////////////////////////////////////
 
RooEmpXHistShape::RooEmpXHistShape(const RooEmpXHistShape& other, const char* name) :
  RooAbsPdf(other,name),
  _x("x",this,other._x),
  _coefList("coefList",this,other._coefList),
  _histType(other._histType),
  _EPFType(other._EPFType),
  _polyType(other._polyType),
  _refRangeName(other._refRangeName),
  _xMin(other._xMin),
  _xMax(other._xMax),
  _histContent(other._histContent),
  _binWidth(other._binWidth)
{
}
 
////////////////////////////////////////////////////////////////////////////////
 
Double_t RooEmpXHistShape::evaluate() const
{
  double x=_x; //, xmin=_x.min(), xmax=_x.max();
  x=std::max(x, _xMin);
  x=std::min(x, _xMax);
  unsigned int bin = int((x-_xMin)/_binWidth); // Always round down

  return empiricalFunc(_xMin+(bin+0.5)*_binWidth)*_histContent[bin];
}

Double_t RooEmpXHistShape::bernstein(Double_t x_) const {
  const Double_t xmax = _x.max(_refRangeName?_refRangeName->GetName():0);
  const Double_t xmin = _x.min(_refRangeName?_refRangeName->GetName():0);

  Double_t x_norm = (x_ - xmin) / (xmax - xmin); // rescale to [0,1]
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


Double_t RooEmpXHistShape::polynomial(Double_t x_) const{
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
  const Double_t x = x_;
  Double_t retVal = _wksp[sz - 1];
  for (unsigned i = sz - 1; i--; ) retVal = _wksp[i] + x * retVal;
  return retVal * std::pow(x, lowestOrder) + (lowestOrder ? 1.0 : 0.0);
}

Double_t RooEmpXHistShape::chebychev(Double_t x_) const
{
  const Double_t xmax = _x.max(_refRangeName?_refRangeName->GetName():0);
  const Double_t xmin = _x.min(_refRangeName?_refRangeName->GetName():0);

  Double_t x(-1+2*(x_-xmin)/(xmax-xmin));
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

Double_t RooEmpXHistShape::empiricalFunc(Double_t x) const{
  double poly=0;
  if(_polyType==Bern) poly=bernstein(x);
  else if(_polyType==Cheby) poly=chebychev(x);
  else if(_polyType==Plain) poly=polynomial(x);
  else{
    std::cerr << "Unknown polynomial type " << _polyType << std::endl;
    R__ASSERT(false);
  }
  double ret=1;
  if(_EPFType==Epoly) ret=exp(poly);
  else if(_EPFType==Poly) ret=poly;
  else if(_EPFType==Power) ret=pow(x,poly);
  else{
    std::cerr << "Unknown empirical function type " << _EPFType << std::endl;
    R__ASSERT(false);
  }
  
  ret=std::max(ret, std::numeric_limits<Double_t>::min());

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine integration scenario. If no interpolation is used,
/// RooEmpXHistShape can perform all integrals over its dependents
/// analytically via partial or complete summation of the input
/// histogram. If interpolation is used on the integral over
/// all histogram observables is supported
 
Int_t RooEmpXHistShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if(matchArgs(allVars, analVars, _x)) return 1;
  else return 0;
}
 
 
 
////////////////////////////////////////////////////////////////////////////////
/// Return integral identified by 'code'. The actual integration
/// is deferred to RooDataHist::sum() which implements partial
/// or complete summation over the histograms contents
 
Double_t RooEmpXHistShape::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1);

  double xmin=std::max(_x.min(rangeName), _xMin), xmax=std::min(_x.max(rangeName), _xMax);
  Double_t total(0);
  // Check whether xmin and xmax are at boundaries
  unsigned int offset = int( (xmin - _xMin) / _binWidth );  
  if(fabs( (xmin-_xMin)/_binWidth - offset ) > epsilon){ // if xmin is in the middle of one bin
    total += (_xMin + (offset+1) * _binWidth - xmin) * empiricalFunc(_xMin+(offset+0.5)*_binWidth) *_histContent[offset];
    offset += 1;               // Shift up the offset by 1
    xmin = _xMin + offset * _binWidth;           // Shift also xmin to match bin boundary
  }
  
  unsigned int nbin = int( (xmax - xmin) / _binWidth );
  if(fabs( nbin - (xmax-xmin)/_binWidth ) > epsilon){ // if xmax is in the middle of one bin
    total += ( xmax - (_xMin + (nbin + offset) * _binWidth) ) * empiricalFunc(_xMin + (nbin + offset + 0.5) * _binWidth) *_histContent[nbin + offset];
    xmax = _xMin + (nbin + offset) * _binWidth;
  }
  
  for(unsigned int ibin=offset; ibin<nbin+offset; ibin++){
    total+=empiricalFunc(_xMin+(ibin+0.5)*_binWidth)*_histContent[ibin]*_binWidth;
  }
  return total ;
}

////////////////////////////////////////////////////////////////////////////////
/// Return sampling hint for making curves of (projections) of this function
/// as the recursive division strategy of RooCurve cannot deal efficiently
/// with the vertical lines that occur in a non-interpolated histogram
 
std::list<Double_t>* RooEmpXHistShape::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  // Check that observable is in dataset, if not no hint is generated
  if (TString(obs.GetName())!=TString(_x.arg().GetName())) {
    std::cout<<"Does not match!"<<std::endl;
    std::cout<<obs.GetName()<<std::endl;
    std::cout<<_x.GetName()<<std::endl;
    return 0 ;
  }
  const RooAbsLValue* lva=dynamic_cast<const RooAbsLValue*>(&_x.arg());
  if(!lva){
    std::cout<<"Not RooAbsLValue!"<<std::endl;
    return 0;
  }
  
  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lva->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;
 
  std::list<Double_t>* hint = new std::list<Double_t> ;
 
  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      hint->push_back(boundaries[i]) ;
    }
  }
 
  return hint ;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Return sampling hint for making curves of (projections) of this function
/// as the recursive division strategy of RooCurve cannot deal efficiently
/// with the vertical lines that occur in a non-interpolated histogram
 
std::list<Double_t>* RooEmpXHistShape::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{ 
  // Check that observable is in dataset, if not no hint is generated
  if (TString(obs.GetName())!=TString(_x.arg().GetName())) {
    std::cout<<"Does not match!"<<std::endl;
    std::cout<<obs.GetName()<<std::endl;
    std::cout<<_x.GetName()<<std::endl;
    return 0 ;
  }
  const RooAbsLValue* lva=dynamic_cast<const RooAbsLValue*>(&_x.arg());
  if(!lva){
    std::cout<<"Not RooAbsLValue!"<<std::endl;
    return 0;
  }
  // Retrieve position of all bin boundaries
   
  std::list<Double_t>* hint = new std::list<Double_t> ;
 
  // Widen range slighty
  xlo = xlo - 0.01*(xhi-xlo) ;
  xhi = xhi + 0.01*(xhi-xlo) ;
 
  Double_t delta = (xhi-xlo)*epsilon ;
  
  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (double boundary=_xMin; boundary<_xMax; boundary+=_binWidth) {
    if (boundary>=xlo && boundary<=xhi) {
      hint->push_back(boundary-delta) ;
      hint->push_back(boundary+delta) ;
    }
  }
 
  return hint ;
}

////////////////////////////////////////////////////////////////////////////////
 
void RooEmpXHistShape::selectNormalizationRange(const char* rangeName, Bool_t force)
{
  if (rangeName && (force || !_refRangeName)) {
    _refRangeName = (TNamed*) RooNameReg::instance().constPtr(rangeName) ;
  }
  if (!rangeName) {
    _refRangeName = 0 ;
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooEmpXHistShape::DumpHist() const{
  std::cout << "xMin " << _xMin
      << ", xMax " << _xMax
      << ", bin width " << _binWidth
      << std::endl;
  for(auto bin : _histContent) std::cout << bin << ", ";
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

TString RooEmpXHistShape::getAttributeValue(TXMLNode* rootNode, TString attributeKey) const{
  TListIter attribIt = rootNode->GetAttributes();
  TXMLAttr* curAttr = NULL;

  while (( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 )
    if ( curAttr->GetName() == attributeKey ) break;

  if(!curAttr){
    std::cerr << "\tERROR: cannot find attribute \"" << attributeKey
        << "\"" << std::endl;
    R__ASSERT(false);
  }
  return curAttr->GetValue();
}

////////////////////////////////////////////////////////////////////////////////

void RooEmpXHistShape::setupHist(const char* inputXMLFileName){
  // Clear vectors
  _histContent.clear();

  // Fill in histogram info from XML file
  TDOMParser xmlparser;
  Int_t parseError = xmlparser.ParseFile( inputXMLFileName );

  if ( parseError ){
    std::cerr << "\tERROR: Loading of xml document \"" << inputXMLFileName
        << "\" failed" << std::endl;
    R__ASSERT(false);
  }

  TXMLDocument* xmldoc=xmlparser.GetXMLDocument();
  TXMLNode* rootNode=xmldoc->GetRootNode();
  TXMLNode* node=rootNode->GetChildren();

  while ( node != 0 ){
    TString nodeName=node->GetNodeName();
    if(nodeName=="Data"){
      TString dataName = getAttributeValue(node, "Name");
      if(dataName == _histType){
  _xMin = getAttributeValue(node, "Min").Atof();
  _xMax = getAttributeValue(node, "Max").Atof();
  _binWidth = getAttributeValue(node, "BinWidth").Atof();
  
  std::vector<std::pair<double, double> > dataList;
  TXMLNode* subnode=node->GetChildren();
  while ( subnode != 0 ){
    TString subnodeName=subnode->GetNodeName();
    if(subnodeName=="Bin")
      dataList.push_back(std::make_pair(getAttributeValue(subnode, "X").Atof(), getAttributeValue(subnode, "Y").Atof()));
    subnode=subnode->GetNextNode();
  }
  std::sort(dataList.begin(), dataList.end());
  for(auto data : dataList) _histContent.push_back(data.second);
  break;
      }
    } 
    node=node->GetNextNode();
  }

  // Sanity check
  if(!_histContent.size()){
    std::cerr << "\tERROR: empty histogram. Please check inputs" << std::endl;
    R__ASSERT(false);
  }

  if(nearbyint((_xMax-_xMin)/_binWidth) != _histContent.size()){
    std::cerr << "\tERROR: inconsistent bin width provided in XML file. Please check inputs" << std::endl;
    R__ASSERT(false);
  }

  // Duplicate the last element to have contingence
  _histContent.push_back(_histContent.back());
}  
