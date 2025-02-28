#include "GH2D.h"

#include <iostream>

#include <TDirectory.h>

#include "GH1D.h"

ClassImp(GH2D)

   GH2D::GH2D(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy,
              const Double_t* ybins)
   : TH2D(name, title, nbinsx, xbins, nbinsy, ybins), GH2Base()
{
}

GH2D::GH2D(const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins)
   : TH2D(name, title, nbinsx, xbins, nbinsy, ybins), GH2Base()
{
}

GH2D::GH2D(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, Double_t ylow,
           Double_t yup)
   : TH2D(name, title, nbinsx, xbins, nbinsy, ylow, yup), GH2Base()
{
}

GH2D::GH2D(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy,
           Double_t* ybins)
   : TH2D(name, title, nbinsx, xlow, xup, nbinsy, ybins), GH2Base()
{
}

GH2D::GH2D(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow,
           Double_t yup)
   : TH2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup), GH2Base()
{
}

GH2D::GH2D(const TObject& obj)
{
   if(obj.InheritsFrom(TH2::Class())) {
      obj.Copy(*this);
   }
}

GH2D::~GH2D() = default;

void GH2D::Copy(TObject& obj) const
{
   TH2::Copy(obj);
   // fProjections->Copy(*(((GH2D&)obj).fProjections));
   // fSummaryProjections->Copy(*(((GH2D&)obj).fSummaryProjections));
}

TObject* GH2D::Clone(const char* newname) const
{
   std::string name = newname;
   if(name.length() == 0u) {
      name = Form("%s_clone", GetName());
   }
   return TH2::Clone(name.c_str());
}

void GH2D::Clear(Option_t* opt)
{
   TString sopt(opt);
   if(!sopt.Contains("projonly")) {
      TH2D::Clear(opt);
   }
   GH2Clear();
}

void GH2D::Print(Option_t*) const
{
}

void GH2D::Draw(Option_t* opt)
{
   std::string option = opt;
   if(option == "") {
      option = "colz";
   }
   TH2D::Draw(option.c_str());
   if(gPad) {
      gPad->Update();
      gPad->GetFrame()->SetBit(TBox::kCannotMove);
   }
}

void GH2D::Draw(TCutG* cut)
{
   if(cut == nullptr) {
      return;
   }
   std::string option = Form("colz [%s]", cut->GetName());
   TH2D::Draw(option.c_str());
}

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
TH1* GH2D::DrawCopy(Option_t* opt) const
{
   TH1* h = TH2D::DrawCopy(opt);
#else
TH1* GH2D::DrawCopy(Option_t* opt, const char* name_postfix) const
{
   TH1* h = TH2D::DrawCopy(opt, name_postfix);
#endif
   if(gPad) {
      gPad->Update();
      gPad->GetFrame()->SetBit(TBox::kCannotMove);
   }
   return h;
}

TH1* GH2D::DrawNormalized(Option_t* opt, Double_t norm) const
{
   TH1* h = TH2D::DrawNormalized(opt, norm);
   if(gPad) {
      gPad->Update();
      gPad->GetFrame()->SetBit(TBox::kCannotMove);
   }
   return h;
}

GH1D* GH2D::ProjectionX(const char* name, int firstbin, int lastbin, Option_t* option)
{
   return GH2ProjectionX(name, firstbin, lastbin, option);
}

GH1D* GH2D::ProjectionY(const char* name, int firstbin, int lastbin, Option_t* option)
{
   return GH2ProjectionY(name, firstbin, lastbin, option);
}

GH1D *GH2D::ProjectionX_BG(const char *name, int ylowbin, int yhighbin, int ylowbgbin, int yhighbgbin,
			   double scale, Option_t *opt) {
  
  GH1D *add = ProjectionX("_px",ylowbin,yhighbin);
  GH1D *sub = ProjectionX("_px",ylowbgbin,yhighbgbin);

  if(scale>0)
    scale*=-1;
  
  add->Add(sub,scale);

  add->SetName(Form("%s_bg",add->GetName()));
  add->SetTitle(Form("%s - %s",add->GetTitle(),sub->GetTitle()));

  sub->Delete();
  
  return add;
}

GH1D *GH2D::ProjectionY_BG(const char *name, int xlowbin, int xhighbin, int xlowbgbin, int xhighbgbin,
			  double scale, Option_t *opt) {
  
  GH1D *add = ProjectionY("_py",xlowbin,xhighbin);
  GH1D *sub = ProjectionY("_py",xlowbgbin,xhighbgbin);
  
  if(scale>0)
    scale*=-1;
  
  add->Add(sub,scale);
  
  add->SetName(Form("%s_bg",add->GetName()));
  add->SetTitle(Form("%s - %s",add->GetTitle(),sub->GetTitle()));

  sub->Delete();

  return add;
}
