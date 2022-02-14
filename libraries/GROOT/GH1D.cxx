#include "GH1D.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include "TVirtualPad.h"
#include "TString.h"
#include "TF1.h"
#include "TFrame.h"
//#include "TROOT.h"
//#include "TSystem.h"

#include "GCanvas.h"
#include "GH2I.h"
#include "GH2D.h"

GH1D::GH1D(const TH1& source) : parent(nullptr), projection_axis(-1)
{
   source.Copy(*this);
}

GH1D::GH1D(const TF1& function, Int_t nbinsx, Double_t xlow, Double_t xup)
   : TH1D(Form("%s_hist", function.GetName()), Form("%s_hist", function.GetName()), nbinsx, xlow, xup), parent(nullptr),
     projection_axis(-1)
{

   // TF1 *f = (TF1*)function.Clone();
   // f->SetRange(xlow,xup);

   for(int i = 0; i < nbinsx; i++) {
      double x = GetBinCenter(i);
      Fill(x, function.Eval(x));
   }
   // f->Delete();
}

GH1D* GH1D::Calibrate(double m, double b) {

  GH1D* htmp = (GH1D*)Clone(Form("%s_%s",GetName(),"cal"));
  htmp->SetTitle(Form("%s_%s",GetTitle(),"cal"));
  htmp->GetXaxis()->SetLimits(GetXaxis()->GetXmin()*m + b,GetXaxis()->GetXmax()*m + b);
  
  new GCanvas();
  htmp->Draw("hist");

  return htmp;
  
}

GH1D* GH1D::Calibrate(std::vector<double> raw, std::vector<double> cal) {

  int size = raw.size();
  if(cal.size() != (unsigned int)size) {
    std::cout << "Must be same number of raw values and calibrated values" << std::endl;
    return (GH1D*)NULL;
  }

  double xmean = 0.0;
  for(auto x : raw)
    xmean += x;
  xmean /= double(size);

  double ymean = 0.0;
  for(auto y : cal)
    ymean += y;
  ymean /= double(size);

  double nume = 0.0;
  double denom = 0.0;
  for(int i=0;i<size;i++) {
    nume += (raw.at(i) - xmean)*(cal.at(i) - ymean);
    denom += TMath::Power(raw.at(i) - xmean,2.0);
  }

  double m = nume/denom;
  double b = ymean - m*xmean;
  
  GH1D* htmp = (GH1D*)Clone(Form("%s_%s",GetName(),"cal1"));
  htmp->SetTitle(Form("%s_%s",GetTitle(),"cal1"));
  htmp->GetXaxis()->SetLimits(GetXaxis()->GetXmin()*m + b,GetXaxis()->GetXmax()*m + b);
  
  new GCanvas();
  htmp->Draw("hist");

  return htmp;
  
}

bool GH1D::WriteDatFile(const char* outFile)
{
   if(strlen(outFile) < 1) {
      return false;
   }

   std::ofstream out;
   out.open(outFile);

   if(!(out.is_open())) {
      return false;
   }

   out << GetNbinsX() << "\t" << GetXaxis()->GetXmin() << "\t" << GetXaxis()->GetXmax() << std::endl;
   for(int i = 1; i < GetNbinsX() + 1; i++) {
      out << GetBinContent(i) << std::endl;
   }
   out<<std::endl;
   out.close();

   return true;
}

bool GH1D::WriteDatFileError(const char* outFile)
{
   if(strlen(outFile) < 1) {
      return false;
   }

   std::ofstream out;
   out.open(outFile);

   if(!(out.is_open())) {
      return false;
   }

   out << GetNbinsX() << "\t" << GetXaxis()->GetXmin() << "\t" << GetXaxis()->GetXmax() << std::endl;
   for(int i = 1; i < GetNbinsX() + 1; i++) {
      out << GetBinContent(i) << "\t" << GetBinError(i) << std::endl;
   }
   out<<std::endl;
   out.close();

   return true;
}

GH1D* GH1D::ReadDatFile(const char* inFile, int nSkip) {

  GH1D* htmp = (GH1D*)NULL;
  if(strlen(inFile) < 1) {
      return htmp;
   }

   std::ifstream in;
   in.open(inFile);

   if(!(in.is_open())) {
      return htmp;
   }

   std::string line;
   if(nSkip > -1) { //Go to where the column of numbers starts, can be the first line (nSkip = 0)

     for(int i=0;i<nSkip;i++) {
       std::getline(in,line);
     }

     std::vector<double> vals;
     while(std::getline(in,line)) {

       double val;
       std::stringstream ss(line);
       ss >> val;
       
       vals.push_back(val);
     }

     int size = vals.size();
     htmp = new GH1D(inFile,inFile,size,0,size);
     
     for(int i=0;i<size;i++) {
       htmp->SetBinContent(i+1,vals.at(i));
       htmp->SetBinError(i+1,TMath::Sqrt(TMath::Abs(vals.at(i))));
     }
   }
   
   else { //The first line is: nbins xlow xhigh. Next line is the column of numbers

     std::getline(in,line);
     std::stringstream ss1(line);

     int nbins;
     ss1 >> nbins;

     double xlow;
     ss1 >> xlow;

     double xhigh;
     ss1 >> xhigh;

     htmp = new GH1D(inFile,inFile,nbins,xlow,xhigh);

     int linenum = 1;
     while(std::getline(in,line)) {

       double val;
       std::stringstream ss2(line);
       ss2 >> val;
       
       htmp->SetBinContent(linenum,val);
       htmp->SetBinError(linenum,TMath::Sqrt(TMath::Abs(val)));

       linenum++;
     }
     
   }
   
   return htmp;
  
}

GH1D* GH1D::ReadDatFileError(const char* inFile, int nSkip) {

  GH1D* htmp = (GH1D*)NULL;
  if(strlen(inFile) < 1) {
      return htmp;
   }

   std::ifstream in;
   in.open(inFile);

   if(!(in.is_open())) {
      return htmp;
   }

   std::string line;
   if(nSkip > -1) { //Go to where the columns of numbers and errors start, can be the first line (nSkip = 0)

     for(int i=0;i<nSkip;i++) {
       std::getline(in,line);
     }

     std::vector<double> vals;
     std::vector<double> errs;
     while(std::getline(in,line)) {

       std::stringstream ss1(line);
       
       double val;
       ss1 >> val;
       vals.push_back(val);
       
       ss1 >> val;
       errs.push_back(val);
       
     }

     int size = vals.size();
     htmp = new GH1D(inFile,inFile,size,0,size);
     
     for(int i=0;i<size;i++) {
       htmp->SetBinContent(i+1,vals.at(i));
       htmp->SetBinError(i+1,errs.at(i));
     }
   }
   
   else { //The first line is: nbins xlow xhigh. Next line is the column of numbers and errors

     std::getline(in,line);
     std::stringstream ss1(line);

     int nbins;
     ss1 >> nbins;

     double xlow;
     ss1 >> xlow;

     double xhigh;
     ss1 >> xhigh;
    
     htmp = new GH1D(inFile,inFile,nbins,xlow,xhigh);

     int linenum = 1;
     while(std::getline(in,line)) {

       std::stringstream ss4(line);
       
       double val; 
       ss4 >> val;
       htmp->SetBinContent(linenum,val);

       ss4 >> val;
       htmp->SetBinError(linenum,val);

       linenum++;
     }
     
   }
   
   return htmp;
  
}

/*
GH1D::GH1D(const TH1 *source)
  : parent(nullptr), projection_axis(-1) {
  if(source->GetDiminsion()>1) {
    return;
  }

  // Can copy from any 1-d TH1, not just a TH1D
  source->Copy(*this);

  // Force a refresh of any parameters stored in the option string.
  SetOption(GetOption());
}

void GH1D::SetOption(Option_t* opt) {
  fOption = opt;

  TString sopt = opt;
  if(sopt.Index("axis:")) {
    projection_axis = 0;// TODO
  }
}
*/

void GH1D::Clear(Option_t* opt)
{
   TH1D::Clear(opt);
   parent = nullptr;
}

void GH1D::Print(Option_t* opt) const
{
   TH1D::Print(opt);
   std::cout<<"\tParent: "<<parent.GetObject()<<std::endl;
}

void GH1D::Copy(TObject& obj) const
{
   TH1D::Copy(obj);

   static_cast<GH1D&>(obj).parent = parent;
}

void GH1D::Draw(Option_t* opt)
{
   TString option(opt);
   if(option.Contains("new", TString::kIgnoreCase)) {
      option.ReplaceAll("new", "");
      new GCanvas;
   }
   TH1D::Draw(option.Data());
   if(gPad) {
      gPad->Update();
      gPad->GetFrame()->SetBit(TBox::kCannotMove);
   }
}

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
TH1* GH1D::DrawCopy(Option_t* opt) const
{
   TH1* h = TH1D::DrawCopy(opt);
#else
TH1* GH1D::DrawCopy(Option_t* opt, const char* name_postfix) const
{
   TH1* h = TH1D::DrawCopy(opt, name_postfix);
#endif
   if(gPad) {
      gPad->Update();
      gPad->GetFrame()->SetBit(TBox::kCannotMove);
   }
   return h;
}

TH1* GH1D::DrawNormalized(Option_t* opt, Double_t norm) const
{
   TH1* h = TH1D::DrawNormalized(opt, norm);
   if(gPad) {
      gPad->Update();
      gPad->GetFrame()->SetBit(TBox::kCannotMove);
   }
   return h;
}

GH1D* GH1D::GetPrevious(bool DrawEmpty) const
{
   if((parent.GetObject() != nullptr) && parent.GetObject()->InheritsFrom(GH2Base::Class())) {
      GH2D* gpar  = static_cast<GH2D*>(parent.GetObject());
      int   first = GetXaxis()->GetFirst();
      int   last  = GetXaxis()->GetLast();
      GH1D* prev  = gpar->GetPrevious(this, DrawEmpty);
      prev->GetXaxis()->SetRange(first, last);
      return prev; // gpar->GetPrevious(this,DrawEmpty);
   }
   return nullptr;
}

GH1D* GH1D::GetNext(bool DrawEmpty) const
{
   if((parent.GetObject() != nullptr) && parent.GetObject()->InheritsFrom(GH2Base::Class())) {
      GH2D* gpar  = static_cast<GH2D*>(parent.GetObject());
      int   first = GetXaxis()->GetFirst();
      int   last  = GetXaxis()->GetLast();
      GH1D* next  = gpar->GetNext(this, DrawEmpty);
      next->GetXaxis()->SetRange(first, last);
      return next; // gpar->GetNext(this,DrawEmpty);
   }
   return nullptr;
}

GH1D* GH1D::Project(int low, int high) const
{

  if((parent.GetObject() == nullptr) || !(parent.GetObject()->InheritsFrom(GH2Base::Class())) || 
     projection_axis == -1)
    return nullptr;

  if(low > high) {
    std::swap(low,high);
  }
  
  GH2D* gpar = static_cast<GH2D*>(parent.GetObject());
  if(projection_axis == 0)
    return gpar->ProjectionY("_py",low,high);
  
  return gpar->ProjectionX("_px",low,high);
   
}

 GH1D* GH1D::Project(int low, int high, int bg_low, int bg_high, double scale) const
{

  if((parent.GetObject() == nullptr) || !(parent.GetObject()->InheritsFrom(GH2Base::Class())) || 
     projection_axis == -1)
    return nullptr;

  if(low > high)
    std::swap(low,high);
  
  GH2D* gpar = static_cast<GH2D*>(parent.GetObject());
  if(projection_axis == 0)
    return gpar->ProjectionY_BG("_py",low,high,bg_low,bg_high,scale,"keep+");
  
  return gpar->ProjectionX_BG("_px",low,high,bg_low,bg_high,scale,"keep+");
   
}

GH1D* GH1D::Project_Background(double value_low, double value_high, double bg_value_low, double bg_value_high,
                               EBackgroundSubtraction mode) const
{
   if((parent.GetObject() != nullptr) && parent.GetObject()->InheritsFrom(GH2Base::Class()) && projection_axis != -1) {
      if(value_low > value_high) {
         std::swap(value_low, value_high);
      }
      if(bg_value_low > bg_value_high) {
         std::swap(bg_value_low, bg_value_high);
      }

      GH2D* gpar = static_cast<GH2D*>(parent.GetObject());
      if(projection_axis == 0) {
         int bin_low     = gpar->GetXaxis()->FindBin(value_low);
         int bin_high    = gpar->GetXaxis()->FindBin(value_high);
         int bg_bin_low  = gpar->GetXaxis()->FindBin(bg_value_low);
         int bg_bin_high = gpar->GetXaxis()->FindBin(bg_value_high);

         return gpar->ProjectionY_Background(bin_low, bin_high, bg_bin_low, bg_bin_high, mode);
      }
      int bin_low     = gpar->GetYaxis()->FindBin(value_low);
      int bin_high    = gpar->GetYaxis()->FindBin(value_high);
      int bg_bin_low  = gpar->GetYaxis()->FindBin(bg_value_low);
      int bg_bin_high = gpar->GetYaxis()->FindBin(bg_value_high);

      return gpar->ProjectionX_Background(bin_low, bin_high, bg_bin_low, bg_bin_high, mode);
   }
   return nullptr;
}

GH1D* GH1D::Project(int bins)
{
   GH1D*  proj = nullptr;
   double ymax = GetMinimum();
   double ymin = GetMaximum();
   if(bins == -1) {
      bins = static_cast<int>(std::abs(ymax - ymin));
      if(bins < 1) {
         bins = 100;
      }
   }
   proj = new GH1D(Form("%s_y_axis_projection", GetName()), Form("%s_y_axis_projection", GetName()), bins, ymin, ymax);
   for(int x = 0; x < GetNbinsX(); x++) {
      if(GetBinContent(x) != 0) {
         proj->Fill(GetBinContent(x));
      }
   }

   return proj;
}
