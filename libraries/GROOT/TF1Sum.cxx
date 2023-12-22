#include "TF1Sum.h"

void TF1Sum::AddTF1(TF1 *f) {
  
  if(!f)
    return;
  
  if(fTF1s.empty()) {
    double min,max;
    f->GetRange(min,max);
    SetRange(min,max);
  }

  fTF1s.push_back(f);
  int num_old_pars = GetNpar(); 
  int num_new_pars = f->GetNpar();
  npars= num_new_pars +num_old_pars;
	
  double lmin,lmax;

  for(int i=0; i<num_new_pars; i++) {
    fParam.push_back(f->GetParameter(i));  
    fParErr.push_back(f->GetParError(i));    
    f->GetParLimits(i,lmin,lmax);
    fParMin.push_back(lmin);          
    fParMax.push_back(lmax);
    fParName.push_back(Form("%s_par%i_%s",f->GetName(),i,f->GetParName(i)));
  }

  if(fFit)
    delete fFit;
    
  fFit = new TF1(Form("%s_obj",GetName()),new TF1Sum(*this),xlow, xhigh, npars);
  for(int i=0;i<npars;i++) {
    fFit->SetParameter(i,fParam.at(i));
    fFit->SetParError(i,fParErr.at(i));
    fFit->SetParLimits(i,fParMin.at(i),fParMax.at(i));
    fFit->SetParName(i,fParName.at(i).c_str());
  }
  

  return;
}

std::vector<std::array<double,2>> exclude;
void TF1Sum::AddExclusion(double x1, double x2) {
  
  exclude.push_back({x1,x2});
  
  return;
}

void TF1Sum::Print(Option_t *opt) const { 

  if(fFit)
    fFit->Print(opt);

  return;
}

Double_t TF1Sum::EvalPar(const Double_t *x,const Double_t *params) {
  int parnum = 0;
  double sum = 0.0;
  
  for(auto fit : fTF1s) {
    
    for(auto exc : exclude)
      if((x[0] > exc[0]) && (x[0] < exc[1]))
	fit->RejectPoint();
    
    if(params==0)
      sum += fit->EvalPar(x,params);
    else
      sum += fit->EvalPar(x, params+parnum);


    parnum += fit->GetNpar();
  }
  
  return sum;
}


void TF1Sum::Draw(Option_t *opt) {
  
  if(!fFit)
    return;

  TString sopt = opt;
  int poffset=0;
  int color=2;
  
  for(int i=0;i<fFit->GetNpar();i++)
    fParam[i] = fFit->GetParameter(i);
  
  if(sopt.Contains("all")) { 
    fFit->Draw(); 
  }
  else {
    for(unsigned int j=0;j<fTF1s.size();j++) {
      fTF1s[j]->SetParameters(fParam.data()+poffset);
      poffset += fTF1s[j]->GetNpar();
      fTF1s[j]->SetLineColor(color++);
      fTF1s[j]->Draw("same");
    }
  }

  return;
}
