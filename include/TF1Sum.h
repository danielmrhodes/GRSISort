#ifndef TF1Sum_H
#define TF1Sum_H

/** \addtogroup Fitting Fitting & Analysis
 *  @{
 */

#include <vector>
#include <string>

#include <TF1.h>

class TF1Sum : public TNamed {

public:
  TF1Sum():TNamed("TF1Sum","TF1Sum"),fFit(0),npars(0) { }
  TF1Sum(const TF1Sum& other) 
    : TNamed(other), fFit(0), npars(other.npars), xlow(other.xlow), xhigh(other.xhigh),
      fTF1s(other.fTF1s) { }

  ~TF1Sum() { if(fFit) delete fFit; } 

  void AddTF1(TF1 *f);
  void AddExclusion(double x1, double x2);
  
  void Print(Option_t *opt="") const override;

  Double_t EvalPar(const Double_t *x,const Double_t *params=0);

  double operator()(double *x, double *p) { return EvalPar(x,p); }
  
  void SetRange(double l,double h) { xlow =l; xhigh=h; }
  int  GetNpar() const  { return npars; }

  operator TF1*() { return fFit;}
        
  TF1 *GetFunc() { return fFit; }

  virtual void Draw(Option_t *opt="all") override;


private:
    
  TF1 *fFit;
  int npars;
  double xlow;
  double xhigh;

  std::vector<double> fParam;
  std::vector<double> fParErr;
  std::vector<double> fParMax;
  std::vector<double> fParMin;
  std::vector<std::string> fParName;

  //std::vector<std::array<double,2>> fExclude;
  std::vector<TF1*> fTF1s;

  /// \cond CLASSIMP
  ClassDefOverride(TF1Sum, 0)
  /// \endcond
};
/*! @} */
#endif
