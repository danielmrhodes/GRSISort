#include "Globals.h"
#include "GCanvas.h"

#include "TClass.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TText.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "Buttons.h"
#include "KeySymbols.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPolyMarker.h"
#include "TSpectrum.h"
#include "TPython.h"
#include "TCutG.h"
#include "TGInputDialog.h"

#include "TApplication.h"
#include "TContextMenu.h"
#include "TGButton.h"

#include "GPopup.h"

#include "GRootCommands.h"
#include "GH2I.h"
#include "GH2D.h"
#include "GH1D.h"
#include "GCutG.h"

#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"

#include "TGRSIint.h"

#ifndef kArrowKeyPress
#define kArrowKeyPress 25
#define kArrowKeyRelease 26
#endif

enum MyArrowPress { kMyArrowLeft = 0x1012, 
		    kMyArrowUp = 0x1013, 
		    kMyArrowRight = 0x1014, 
		    kMyArrowDown = 0x1015 };

/// \cond CLASSIMP
ClassImp(GMarker)
/// \endcond

GMarker::GMarker(int x, int y, TH1* hist)
: fHist(hist)
{
  if(fHist->GetDimension() == 1) {
    double localX = gPad->AbsPixeltoX(x);

    fLineX = new TLine(localX, gPad->GetUymin(), localX, gPad->GetUymax());
    fLineY = nullptr;
    SetColor(kRed);
    Draw();
  } else if(fHist->GetDimension() == 2) {
    double localX     = gPad->AbsPixeltoX(x);
    double localY     = gPad->AbsPixeltoY(y);

    fLineX = new TLine(localX, gPad->GetUymin(), localX, gPad->GetUymax());
    fLineY = new TLine(gPad->GetUxmin(), localY, gPad->GetUxmax(), localY);

    SetColor(kRed);
    Draw();
  }
}

void GMarker::Copy(TObject& object) const
{
  TObject::Copy(object);
  (static_cast<GMarker&>(object)).fLineX  = nullptr;
  (static_cast<GMarker&>(object)).fLineY  = nullptr;
  (static_cast<GMarker&>(object)).fHist   = fHist;
}

int GCanvas::lastx = 0;
int GCanvas::lasty = 0;

GCanvas::GCanvas(Bool_t build) : TCanvas(build)
{
  GCanvasInit();
}

GCanvas::GCanvas(const char* name, const char* title, Int_t form) : TCanvas(name, title, form)
{
  GCanvasInit();
}

GCanvas::GCanvas(const char* name, const char* title, Int_t ww, Int_t wh) : TCanvas(name, title, ww, wh)
{
  GCanvasInit();
}

GCanvas::GCanvas(const char* name, Int_t ww, Int_t wh, Int_t winid) : TCanvas(name, ww, wh, winid)
{
  // this constructor is used to create an embedded canvas
  // I see no reason for us to support this here.  pcb.
  GCanvasInit();
  fGuiEnabled = true;
}

GCanvas::GCanvas(const char* name, const char* title, Int_t wtopx, Int_t wtopy, Int_t ww, Int_t wh, bool gui)
  : TCanvas(name, title, wtopx, wtopy, ww, wh)
{
  GCanvasInit();
  fGuiEnabled = gui;
}

GCanvas::~GCanvas()
{
  // TCanvas::~TCanvas();
  delete[] fCutName;
}

void GCanvas::GCanvasInit()
{
  // ok, to interact with the default TGWindow
  // stuff from the root gui we need our own GRootCanvas.
  // We make this using GROOTGuiFactory, which replaces the
  // TRootGuiFactory used in the creation of some of the
  // default gui's (canvas,browser,etc).
  // fStatsDisplayed = true;
  fMarkerMode     = true;
  fGuiEnabled     = false;
  fBackgroundMode = EBackgroundSubtraction::kNoBackground;
  fCutName = new char[256];
  // if(gVirtualX->InheritsFrom("TGX11")) {
  //    printf("\tusing x11-like graphical interface.\n");
  //}
  // SetCrosshair(true);
  SetBit(kNotDeleted, false); // root voodoo.
}

void GCanvas::AddMarker(int x, int y, TH1* hist)
{
  auto* mark = new GMarker(x, y, hist);
  unsigned int max_number_of_markers = (hist->GetDimension() == 1) ? 4 : 2;

  fMarkers.push_back(mark);

  if(fMarkers.size() > max_number_of_markers) {
    delete fMarkers.at(0);
    fMarkers.erase(fMarkers.begin());
  }
  return;
}

void GCanvas::RemoveMarker(Option_t* opt)
{
  TString options(opt);

  if(options.Contains("all")) {
    for(auto marker : fMarkers) {
      delete marker;
    }
    for(auto marker : fBackgroundMarkers) {
      delete marker;
    }
    fMarkers.clear();
    fBackgroundMarkers.clear();
  } else {
    if(fMarkers.empty()) {
      return;
    }
    delete fMarkers.back();
    // printf("Marker %i Removed\n");
    fMarkers.pop_back();
  }
}

void GCanvas::OrderMarkers()
{
  std::sort(fMarkers.begin(), fMarkers.end());
}

void GCanvas::RedrawMarkers()
{
  gPad->Update();
  for(auto marker : fMarkers) {
    marker->Update(GetUxmin(), GetUxmax(), GetUymin(), GetUymax());
    marker->Draw();
  }

  for(auto marker : fBackgroundMarkers) {
    marker->Update(GetUxmin(), GetUxmax(), GetUymin(), GetUymax());
    marker->Draw();
  }
}

bool GCanvas::SetBackgroundMarkers()
{
  if(GetNMarkers() < 2) {
    return false;
  }

  // Delete previous background, if any.
  for(auto marker : fBackgroundMarkers) {
    delete marker;
  }
  fBackgroundMarkers.clear();

  // Push last two markers into the background.
  fBackgroundMarkers.push_back(fMarkers.back());
  fMarkers.pop_back();
  fBackgroundMarkers.push_back(fMarkers.back());
  fMarkers.pop_back();

  // Change background marker color.
  for(auto marker : fBackgroundMarkers) {
    marker->SetColor(kBlue);
  }

  fBackgroundMode = EBackgroundSubtraction::kRegionBackground;

  return true;
}

bool GCanvas::CycleBackgroundSubtraction()
{
  if(fBackgroundMarkers.size() < 2) {
    return false;
  }

  Color_t color = 0;

  switch(fBackgroundMode) {
  case EBackgroundSubtraction::kNoBackground:
    fBackgroundMode = EBackgroundSubtraction::kRegionBackground;
    printf("hello??\n");
    Prompt();
    color = kBlue;
    break;
  case EBackgroundSubtraction::kRegionBackground:
    fBackgroundMode = EBackgroundSubtraction::kTotalFraction;
    color           = kGreen;
    break;
  case EBackgroundSubtraction::kTotalFraction:
    fBackgroundMode = EBackgroundSubtraction::kMatchedLowerMarker;
    color           = kOrange;
    break;
  case EBackgroundSubtraction::kMatchedLowerMarker:
    fBackgroundMode = EBackgroundSubtraction::kSplitTwoMarker;
    color           = kMagenta;
    break;
  case EBackgroundSubtraction::kSplitTwoMarker:
    fBackgroundMode = EBackgroundSubtraction::kNoBackground;
    color           = 0;
    break;
  };

  for(auto marker : fBackgroundMarkers) {
    marker->SetColor(color);
  }

  return true;
}

GCanvas* GCanvas::MakeDefCanvas()
{
  // Static function to build a default canvas.

  const char* defcanvas = gROOT->GetDefCanvasName();
  char*       cdef;
  TList*      lc = static_cast<TList*>(gROOT->GetListOfCanvases());
  if(lc->FindObject(defcanvas) != nullptr) {
    Int_t n = lc->GetSize() + 1;
    cdef    = new char[strlen(defcanvas) + 15];
    do {
      strlcpy(cdef, Form("%s_n%d", defcanvas, n++), strlen(defcanvas) + 15);
    } while(lc->FindObject(cdef) != nullptr);
  } else {
    cdef = StrDup(Form("%s", defcanvas));
  }
  auto* c = new GCanvas(cdef, cdef, 1);
  delete[] cdef;
  return c;
}

void GCanvas::HandleInput(int event, Int_t x, Int_t y)
{
  // If the below switch breaks. You need to upgrade your version of ROOT
  // Version 5.34.24 works. //older version should work now too pcb (8/2015)
  bool used = false;
  switch(event) {
  case kButton1Down:   // single click
  case kButton1Double: // double click
    used = HandleMousePress(event, x, y);
    break;
  case kButton1Shift: // shift-click
    used = HandleMouseShiftPress(event, x, y);
    break;
  case 9: // control-click
    used = HandleMouseControlPress(event, x, y);
    break;
  };
  if(!used) {
    TCanvas::HandleInput(static_cast<EEventType>(event), x, y);
  }
  return;
}

void GCanvas::Draw(Option_t* opt)
{
  printf("GCanvas Draw was called.\n");
  TCanvas::Draw(opt);
  if(FindObject("TFrame") != nullptr) {
    FindObject("TFrame")->SetBit(TBox::kCannotMove);
  }
}

std::vector<TH1*> GCanvas::FindHists(int dim)
{
  std::vector<TH1*> tempvec;
  TIter             iter(gPad->GetListOfPrimitives());
  while(TObject* obj = iter.Next()) {
    if(obj->InheritsFrom(TH1::Class())) {
      TH1* hist = static_cast<TH1*>(obj);
      if(hist->GetDimension() == dim) {
	tempvec.push_back(hist);
      }
    }
  }
  return tempvec;
}

std::vector<TH1*> GCanvas::FindAllHists()
{
  std::vector<TH1*> tempvec;
  TIter             iter(gPad->GetListOfPrimitives());
  while(TObject* obj = iter.Next()) {
    if(obj->InheritsFrom("TH1")) {
      tempvec.push_back(static_cast<TH1*>(obj));
    }
  }
  return tempvec;
}

bool GCanvas::HandleArrowKeyPress(Event_t* event, UInt_t* keysym)
{

  bool edited = Process1DArrowKeyPress(event, keysym);
  if(!edited) {
    edited = Process2DArrowKeyPress(event, keysym);
  }

  if(edited) {
    gPad->Modified();
    gPad->Update();
  }
  return true;
}

bool GCanvas::HandleKeyboardPress(Event_t* event, UInt_t* keysym)
{
  bool edited = false;

  edited = ProcessNonHistKeyboardPress(event, keysym);

  if(!edited) {
    edited = Process1DKeyboardPress(event, keysym);
  }
  if(!edited) {
    edited = Process2DKeyboardPress(event, keysym);
  }

  if(edited) {
    gPad->Modified();
    gPad->Update();
  }
  return true;
}

bool GCanvas::HandleMousePress(Int_t event, Int_t x, Int_t y)
{
  if(GetSelected() == nullptr) {
    return false;
  }

  TH1* hist = nullptr;
  if(GetSelected()->InheritsFrom(TH1::Class())) {
    hist = static_cast<TH1*>(GetSelected());
  } else if(GetSelected()->IsA() == TFrame::Class()) {
    std::vector<TH1*> hists = FindAllHists();
    if(!hists.empty()) {
      hist = hists.front();

      // Let everybody know that the histogram is selected
      SetSelected(hist);
      SetClickSelected(hist);
      Selected(GetSelectedPad(), hist, event);
    }
  }

  if((hist == nullptr) || hist->GetDimension() > 2) {
    return false;
  }

  bool used = false;

  if(fMarkerMode) {
    AddMarker(x, y, hist);
    used = true;
  }

  if(used) {
    gPad->Modified();
    gPad->Update();
  }

  return used;
}

bool GCanvas::HandleMouseShiftPress(Int_t, Int_t, Int_t)
{
  TH1*  hist = nullptr;
  TIter iter(gPad->GetListOfPrimitives());
  while(TObject* obj = iter.Next()) {
    if(obj->InheritsFrom(TH1::Class())) {
      hist = static_cast<TH1*>(obj);
    }
  }
  if(hist == nullptr) {
    return false;
  }

  TString options;
  switch(hist->GetDimension()) {
  case 1: {
    if(hist->InheritsFrom(GH1D::Class())) {
      new GCanvas();
      (static_cast<GH1D*>(hist))->GetParent()->Draw("colz");
      return true;
    }
    std::vector<TH1*> hists = FindHists();
    new GCanvas();
    // options.Append("HIST");
    hists.at(0)->DrawCopy(options.Data());
    for(unsigned int j = 1; j < hists.size(); j++) {
      hists.at(j)->DrawCopy("same");
    }
  }
    return true;
  case 2:
    options.Append("colz");
    auto* ghist = new GH2D(*(static_cast<TH2*>(hist)));
    new GCanvas();
    ghist->Draw();
    return true;
  };
  return false;
}

bool GCanvas::HandleMouseControlPress(Int_t, Int_t, Int_t)
{
  // printf("GetSelected() = 0x%08x\n",GetSelected());
  if(GetSelected() == nullptr) {
    return false;
  }
  // printf("GetSelected()->GetName() = %s\n",GetSelected()->GetName());
  if(GetSelected()->InheritsFrom(TCutG::Class())) {
    // TODO: Bring this back, once we have brought over more from GRUTinizer
    // if(TRuntimeObjects::Get())
    //   TRuntimeObjects::Get()->GetGates().Add(GetSelected());
  }
  return true;
}

TF1* GCanvas::GetLastFit()
{
  TH1*  hist = nullptr;
  TIter iter(gPad->GetListOfPrimitives());
  while(TObject* obj = iter.Next()) {
    if(obj->InheritsFrom("TH1") && !obj->InheritsFrom("TH2") && !obj->InheritsFrom("TH3")) {
      hist = static_cast<TH1*>(obj);
    }
  }
  if(hist == nullptr) {
    return nullptr;
  }
  if(hist->GetListOfFunctions()->GetSize() > 0) {
    TF1* tmpfit = static_cast<TF1*>(hist->GetListOfFunctions()->Last());
    return tmpfit;
  }
  return nullptr;
}

bool GCanvas::Process1DArrowKeyPress(Event_t*, UInt_t* keysym)
{
  bool              edited = false;
  std::vector<TH1*> hists  = FindHists();
  if(hists.empty()) {
    return edited;
  }

  int first = hists.at(0)->GetXaxis()->GetFirst();
  int last  = hists.at(0)->GetXaxis()->GetLast();

  int min = std::min(first, 0);
  int max = std::max(last, hists.at(0)->GetXaxis()->GetNbins() + 1);
  // int max = std::max(last,axis->GetNbins()+1);

  int xdiff = last - first;
  int mdiff = max - min - 2;

  switch(*keysym) {
  case kMyArrowLeft: {
    if(mdiff > xdiff) {
      if(first == (min + 1)) {
	//
      } else if((first - (xdiff / 2)) < min) {
	first = min + 1;
	last  = min + (xdiff) + 1;
      } else {
	first = first - (xdiff / 2);
	last  = last - (xdiff / 2);
      }
    }
    for(auto& hist : hists) {
      hist->GetXaxis()->SetRange(first, last);
    }

    edited = true;
  } break;
  case kMyArrowRight: {
    if(mdiff > xdiff) {
      if(last == (max - 1)) {
	//
      } else if((last + (xdiff / 2)) > max) {
	first = max - 1 - (xdiff);
	last  = max - 1;
      } else {
	last  = last + (xdiff / 2);
	first = first + (xdiff / 2);
      }
    }
    for(auto& hist : hists) {
      hist->GetXaxis()->SetRange(first, last);
    }

    edited = true;
  } break;

  case kMyArrowUp: {
    GH1D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH1D::Class())) {
	ghist = static_cast<GH1D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      TH1* prev = ghist->GetNext();
      if(prev != nullptr) {
	prev->GetXaxis()->SetRange(first, last);
	prev->Draw("");
	RedrawMarkers();
	edited = true;
      }
    }
  } break;

  case kMyArrowDown: {
    GH1D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH1D::Class())) {
	ghist = static_cast<GH1D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      TH1* prev = ghist->GetPrevious();
      if(prev != nullptr) {
	prev->GetXaxis()->SetRange(first, last);
	prev->Draw("");
	RedrawMarkers();
	edited = true;
      }
    }
  } break;
  default: printf("keysym = %i\n", *keysym); break;
  }
  return edited;
}

bool GCanvas::ProcessNonHistKeyboardPress(Event_t*, UInt_t* keysym)
{
  bool edited = false;

  switch(*keysym) {
  case kKey_F2:
    GetCanvasImp()->ShowEditor(!GetCanvasImp()->HasEditor());
    edited = true;
    break;
  case kKey_F9: 
    SetCrosshair(static_cast<Int_t>(!HasCrosshair()));
    edited = true;
    break;
  }

  return edited;
}

bool GCanvas::Process1DKeyboardPress(Event_t*, UInt_t* keysym)
{
  bool edited = false;
  std::vector<TH1*> hists  = FindHists();
  if(hists.empty()) {
    return edited;
  }

  switch(*keysym) {

  case kKey_a:

    {
      GetContextMenu()->Action(hists.at(0)->GetXaxis(),hists.at(0)->GetXaxis()->Class()->GetMethodAny("SetRangeUser"));
      
      double x1 = hists.at(0)->GetXaxis()->GetBinCenter(hists.at(0)->GetXaxis()->GetFirst());
      double x2 = hists.at(0)->GetXaxis()->GetBinCenter(hists.at(0)->GetXaxis()->GetLast());

      //printf("(x1,x2)=(%f,%f)\n",x1,x2);
    
      for(unsigned int i=1;i<hists.size();i++) {
	hists.at(i)->GetXaxis()->SetRangeUser(x1,x2);
      }
      
    }
    edited = true;
    break;

  case kKey_A:
    GetContextMenu()->Action(hists.back()->GetXaxis(),hists.back()->GetXaxis()->Class()->GetMethodAny("SetRangeUser"));
    {
      double x1 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetFirst());
      double x2 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetLast());
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->GetXaxis()->SetRangeUser(x1,x2);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited = true;
    break;
     
  case kKey_b: edited = SetBackgroundMarkers(); break;

  case kKey_B: edited = CycleBackgroundSubtraction(); break;

  case kKey_c:

    {
      int xbins_before = hists.back()->GetNbinsX();
      GetContextMenu()->Action(hists.at(0),hists.at(0)->Class()->GetMethodAny("Rebin"));
      int xbins_after  = hists.at(0)->GetNbinsX();
      int factor = xbins_before/xbins_after;
      
      for(unsigned int i=1;i<hists.size();i++) {
	hists.at(i)->Rebin(factor);
      }
    }
    edited = true;
    break;
   
  case kKey_C: 
    {
      int xbins_before = hists.back()->GetNbinsX();
      GetContextMenu()->Action(hists.back(),hists.back()->Class()->GetMethodAny("Rebin"));
      int xbins_after  = hists.back()->GetNbinsX();
      int factor = xbins_before/xbins_after;
   
      TH1 *start = hists.back();

      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      if(hist==start) continue;
	      hist->Rebin(factor);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited = true;
    break;

  case kKey_d: {
    //new GPopup(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(), 500, 200);
    hists.at(0)->DrawPanel();
    edited = false;
  } break;

  case kKey_e:
    if(GetNMarkers() < 2) {
      break;
    }
    {
      if(fMarkers.at(fMarkers.size() - 1)->GetLocalX() < fMarkers.at(fMarkers.size() - 2)->GetLocalX()) {
	for(auto& hist : hists) {
	  hist->GetXaxis()->SetRangeUser(fMarkers.at(fMarkers.size() - 1)->GetLocalX(),
					 fMarkers.at(fMarkers.size() - 2)->GetLocalX());
	}
      } else {
	for(auto& hist : hists) {
	  hist->GetXaxis()->SetRangeUser(fMarkers.at(fMarkers.size() - 2)->GetLocalX(),
					 fMarkers.at(fMarkers.size() - 1)->GetLocalX());
	}
      }
    }
    edited = true;
    RemoveMarker("all");
    break;
  case kKey_E:
    // GetListOfPrimitives()->Print();
    GetContextMenu()->Action(hists.back()->GetXaxis(),
			     hists.back()->GetXaxis()->Class()->GetMethodAny("SetRangeUser"));
    {
      double x1 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetFirst());
      double x2 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetLast());
      TIter  iter(GetListOfPrimitives());
      while(TObject* obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad* pad = static_cast<TPad*>(obj);
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject* obj2 = iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = static_cast<TH1*>(obj2);
	      hist->GetXaxis()->SetRangeUser(x1, x2);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }

    }
    edited = true;
    break;
  case kKey_f:
    if(!hists.empty() && GetNMarkers() > 1) {
      printf("x low = %.1f\t\txhigh = %.1f\n",fMarkers.at(fMarkers.size()-2)->GetLocalX(),fMarkers.back()->GetLocalX());
      if(PhotoPeakFit(hists.back(), fMarkers.at(fMarkers.size() - 2)->GetLocalX(), fMarkers.back()->GetLocalX()) != nullptr) {
	edited = true;
      }
    }
    break;

  case kKey_F:
    if(!hists.empty() && GetNMarkers() > 1) {
      printf("x low = %.1f\t\txhigh = %.1f\n",fMarkers.at(fMarkers.size()-2)->GetLocalX(),fMarkers.back()->GetLocalX());
      if(AltPhotoPeakFit(hists.back(), fMarkers.at(fMarkers.size() - 2)->GetLocalX(), fMarkers.back()->GetLocalX(), "q+") !=
	 nullptr) {
	edited = true;
      }
    }
    break;

  case kKey_g:
    hists.back()->GetListOfFunctions()->Clear();
    if(GausFit(hists.back(),fMarkers.at(fMarkers.size() - 2)->GetLocalX(),fMarkers.back()->GetLocalX()) != nullptr) {
      hists.back()->Draw();
      edited = true;
    }
    break;

  case kKey_G:
    hists.back()->GetListOfFunctions()->Clear();
    if(!hists.back() || !(fMarkers.size()==4)) {
      printf( CYAN "must have a a1 hist with 4 markers drawn" RESET_COLOR "\n");
    } else  {
      std::vector<double> xvalues;
      xvalues.push_back(fMarkers.at(0)->GetLocalX());
      xvalues.push_back(fMarkers.at(1)->GetLocalX());
      xvalues.push_back(fMarkers.at(2)->GetLocalX());
      xvalues.push_back(fMarkers.at(3)->GetLocalX());
      RemoveMarker("all");
      std::sort(xvalues.begin(),xvalues.end());
      //std::cout << xvalues.at(1)<<"\t"<<xvalues.at(2)<<"\t"<<xvalues.at(0)<<"\t"<<xvalues.at(3)<<std::endl;
      edited = DoubleGausFit(hists.back(),xvalues.at(1),xvalues.at(2),xvalues.at(0),xvalues.at(3));
      hists.back()->Draw();
    }
    break;

  case kKey_h: {
    this->Clear();
    RemoveMarker("all");
    RemovePeaks(hists.data(),hists.size());
    hists.at(0)->Draw("hist");
    for(unsigned int i=1;i<hists.size();i++) {
      hists.at(i)->Draw("hist same");
    }
    for(unsigned int i=0;i<hists.size();i++) {
      for(const auto&& func : *(hists.at(i)->GetListOfFunctions())) {
	func->Draw("same");
      }
    }
    edited=true;
  }
    break;

  case kKey_H: {
    hists.back()->SetDrawOption("hist");

    TIter iter(this->GetListOfPrimitives());
    while(TObject *obj = iter.Next()) {
      if(obj->InheritsFrom(TPad::Class())) {
	TPad *pad = (TPad*)obj;
	pad->cd();
	TIter iter2(pad->GetListOfPrimitives());
	while(TObject *obj2=iter2.Next()) {
	  if(obj2->InheritsFrom(TH1::Class())) {
	    
	    TH1* hist = (TH1*)obj2;
	    hist->SetDrawOption("hist");

	    for(const auto&& func : *(hist->GetListOfFunctions()))
	      func->Draw("same");
	
	    pad->Modified();
	    pad->Update();
	  }
	}
      }
    }
  }
    edited=true;
    break;

  case kKey_i: {

    if(GetNMarkers() < 1)
      break;

    int bin1 = fMarkers.at(fMarkers.size() - 1)->GetBinX();
    int bin2 = fMarkers.at(fMarkers.size() - 2)->GetBinX();
    if(bin1 > bin2)
      std::swap(bin1,bin2);
    
    std::cout << BLUE "\nSum between bins [" << bin1 << "," << bin2 << "]";
    for(TH1* hist : hists) {
      double sum = hist->Integral(bin1,bin2);
      //printf(BLUE "\n  %s: Sum between bins [%i : %i] = %.01f" RESET_COLOR "\n",
      //     hist->GetName(),bin1,bin2,sum);

      std::cout << "\n  " << hist->GetName() << ": " << sum;
    }
    std::cout << RESET_COLOR << std::endl;
  }
    break;

  case kKey_I:
    break;

  case kKey_l:
    this->Clear();
    hists.at(0)->Draw("hist");
    for(unsigned int i=1;i<hists.size();i++) {
      hists.at(i)->Draw("histsame");
    }
    SetLogy(!GetLogy());
    edited=true;
    break;

  case kKey_L:
    {
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  pad->SetLogy(!pad->GetLogy());
	}
      }
    }
    edited=true;
    break; 

  case kKey_m: SetMarkerMode(true); break;
  case kKey_M: SetMarkerMode(false); break;

  case kKey_n:
    this->Clear();
    RemoveMarker("all");
    RemovePeaks(hists.data(),hists.size());
    hists.at(0)->Draw();
    for(unsigned int i=1;i<hists.size();i++) {
      hists.at(i)->Draw("same");
    }
    edited=true;
    break;
   
  case kKey_N:
    this->Clear();
    RemoveMarker("all");
    RemovePeaks(hists.data(),hists.size());
    hists.at(0)->Draw("hist");
    for(unsigned int i=1;i<hists.size();i++) {
      hists.at(i)->Draw("hist same");
    }
    edited=true;
    break;

  case kKey_o:
    for(auto& hist : hists) {
      hist->GetXaxis()->UnZoom();
      hist->GetYaxis()->UnZoom();
    }
    //RemoveMarker("all");
    edited = true;
    break;

  case kKey_O:
    {
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->GetXaxis()->UnZoom();
	      hist->GetYaxis()->UnZoom();
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited=true;
    break; 

  case kKey_p: {
    if(GetNMarkers() < 2)
      break;
    
    GH1D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH1D::Class())) {
	ghist = static_cast<GH1D*>(hist);
	break;
      }
    }
    
    if(!ghist)
      break;
    
    TObject* parent = ghist->GetParent();
    if(!parent)
      break;

    if(!(parent->InheritsFrom(GH2Base::Class())))
      break;

    int nbins;
    if(ghist->GetProjectionAxis() == 0)
      nbins = ((GH2D*)parent)->GetNbinsX();
    else
      nbins = ((GH2D*)parent)->GetNbinsY();
    
    if(nbins != ghist->GetNbinsX()) {
      std::cout << RED "Cannot rebin child hist or zoom parent hist before projecting... sorry " 
	RESET_COLOR << std::endl;
      break;
    }

    int binlow  = fMarkers.at(fMarkers.size() - 1)->GetBinX();
    int binhigh = fMarkers.at(fMarkers.size() - 2)->GetBinX();
    if(binlow > binhigh)
      std::swap(binlow,binhigh);
    
    GH1D* proj = nullptr;
    if(fBackgroundMarkers.size() >= 2) {

      int bg_binlow  = fBackgroundMarkers.at(0)->GetBinX();
      int bg_binhigh = fBackgroundMarkers.at(1)->GetBinX();
      if(bg_binlow > bg_binhigh)
	std::swap(bg_binlow,bg_binhigh);

      double scale = -1*(double(binhigh) - double(binlow))/(double(bg_binhigh) - double(bg_binlow));
      proj = ghist->Project(binlow,binhigh,bg_binlow,bg_binhigh,scale);
    }
    else
      proj = ghist->Project(binlow,binhigh);

    if(proj) {
      proj->Draw("");
      edited=true;
    }
    
    /*
    //  ok, i found a bug.  if someone tries to gate on a histogram
    //  that is already zoomed, bad things will happen; namely the bins
    //  in the zoomed histogram will not map correctly to the parent. To get
    //  around this we need the bin value, not the bin!   pcb.
    //
    if(ghist != nullptr) {
      GH1D* proj    = nullptr;
      int   binlow  = fMarkers.at(fMarkers.size() - 1)->GetBinX();
      int   binhigh = fMarkers.at(fMarkers.size() - 2)->GetBinX();
      if(binlow > binhigh) {
	std::swap(binlow, binhigh);
      }
      double value_low  = ghist->GetXaxis()->GetBinLowEdge(binlow);
      double value_high = ghist->GetXaxis()->GetBinLowEdge(binhigh);

      {
	double epsilon = 16 * (std::nextafter(value_low, INFINITY) - value_low);
	value_low += epsilon;
      }

      {
	double epsilon = 16 * (value_high - std::nextafter(value_high, -INFINITY));
	value_high -= epsilon;
      }

      if(fBackgroundMarkers.size() >= 2 && fBackgroundMode != EBackgroundSubtraction::kNoBackground) {
	int bg_binlow  = fBackgroundMarkers.at(0)->GetBinX();
	int bg_binhigh = fBackgroundMarkers.at(1)->GetBinX();
	if(bg_binlow > bg_binhigh) {
	  std::swap(bg_binlow, bg_binhigh);
	}
	double bg_value_low  = ghist->GetXaxis()->GetBinCenter(bg_binlow);
	double bg_value_high = ghist->GetXaxis()->GetBinCenter(bg_binhigh);
	// Using binhigh-1 instead of binhigh,
	//  because the ProjectionX/Y functions from ROOT use inclusive bin numbers,
	//  rather than exclusive.
	//
	proj = ghist->Project_Background(value_low, value_high, bg_value_low, bg_value_high, fBackgroundMode);
      } else {
	proj = ghist->Project(value_low, value_high);
      }
      if(proj != nullptr) {
	proj->Draw("");
	edited = true;
      }
      }
    */
  } 
    edited = true;
    break;

  case kKey_P: {
    GH1D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH1D::Class())) {
	ghist = static_cast<GH1D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      ghist->GetParent()->Draw();
      edited = true;
    }
  } break;
  case kKey_q: {
    TH1* ghist = hists.at(0);
    if(GetNMarkers() > 1) {
      edited = (PhotoPeakFit(ghist, fMarkers.at(fMarkers.size() - 2)->GetLocalX(), fMarkers.back()->GetLocalX()) != nullptr);
    }
    if(edited) {
      ghist->Draw("hist");

      TIter iter(ghist->GetListOfFunctions());
      while(TObject* o = iter.Next()) {
	if(o->InheritsFrom(TF1::Class())) {
	  (static_cast<TF1*>(o))->Draw("same");
	}
      }
    }
  }

    break;

  case kKey_r:
    {
      GetContextMenu()->Action(hists.at(0)->GetYaxis(),hists.at(0)->GetYaxis()->Class()->GetMethodAny("SetRangeUser"));
      gPad->Modified();
      gPad->Update();
      
      double y1=gPad->GetUymin();
      double y2=gPad->GetUymax();
      
      for(unsigned int i=1;i<hists.size();i++) {
	hists.at(i)->GetYaxis()->SetRangeUser(y1,y2);
      }
    }
    edited = true;
    break;

  case kKey_R:
    {
      GetContextMenu()->Action(hists.at(0)->GetYaxis(),hists.at(0)->GetYaxis()->Class()->GetMethodAny("SetRangeUser"));
      gPad->Modified();
      gPad->Update();
      
      double y1=gPad->GetUymin();
      double y2=gPad->GetUymax();
      
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->GetYaxis()->SetRangeUser(y1,y2);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited = true;
    break;

  case kKey_s: {

    if(hists.size() == 1) {
      edited = false;
      break;
    }

    int bin1 = 1;
    int bin2 = 0;
    if(GetNMarkers() > 1) {
      bin1 = fMarkers.at(GetNMarkers() - 1)->GetBinX();
      bin2 = fMarkers.at(GetNMarkers() - 2)->GetBinX();
      
      if(bin2 < bin1) {
	std::swap(bin1,bin2);
      }
    }
    
    double smallest = std::numeric_limits<double>::infinity();
    for(TH1* hist : hists) {
      
      double sum = hist->Integral(bin1,bin2);
      if(sum < smallest)
	smallest = sum;
    }
    
    if(smallest < 0) {
	std::cout << "Negative integral found! Won't rescale." << std::endl;
	edited = false;
	break;
    }
      
    for(TH1* hist : hists) {
      
      double sum = hist->Integral(bin1,bin2);	
      if(!(hist->GetSumw2N()))
	hist->Sumw2(true);
      
      hist->Scale(smallest/sum);
      
    }
    
  }
    edited = true;
    break;

  case kKey_S: {

    TIter iter(this->GetListOfPrimitives());
    while(TObject *obj = iter.Next()) {
      
      if(!(obj->InheritsFrom(TPad::Class())))
	continue;
	
      TPad *pad = (TPad*)obj;
      TIter iter2(pad->GetListOfPrimitives());

      double smallest = std::numeric_limits<double>::infinity();
      while(TObject *obj2=iter2.Next()) {
	  
	if(!(obj2->InheritsFrom(TH1::Class())))
	  continue;
	    
	TH1* hist = (TH1*)obj2;
	double sum = hist->Integral();
	if(sum < smallest)
	  smallest = sum;
	  
      }

      if(smallest < 0) {
	std::cout << "Negative integral found!" << std::endl;
	edited = false;
	break;
      }

      TIter iter3(pad->GetListOfPrimitives());
      while(TObject *obj2=iter3.Next()) {
	
	if(!(obj2->InheritsFrom(TH1::Class())))
	  continue;
	    
	TH1* hist = (TH1*)obj2;
	double sum = hist->Integral();
	
	if(!(hist->GetSumw2N()))
	  hist->Sumw2(true);
	
	hist->Scale(smallest/sum);
	    
      }
    }	
      
  }
    edited = true;
    break;

  case kKey_x: {
    this->Clear();
    RemoveMarker("all");
    RemovePeaks(hists.data(),hists.size());
    hists.at(0)->GetListOfFunctions()->Clear();
    hists.at(0)->Draw("hist");
    for(unsigned int i=1;i<hists.size();i++) {
      hists.at(i)->GetListOfFunctions()->Clear();
      hists.at(i)->Draw("histsame");
    }
    edited=true;
  }
    break;

  case kKey_X: {
    
    RemoveMarker("all");
    RemovePeaks(hists.data(),hists.size());
    
    TIter iter(this->GetListOfPrimitives());
    while(TObject *obj = iter.Next()) {
      if(obj->InheritsFrom(TPad::Class())) {
	TPad *pad = (TPad*)obj;
	TIter iter2(pad->GetListOfPrimitives());
	while(TObject *obj2=iter2.Next()) {
	  if(obj2->InheritsFrom(TH1::Class())) {
	    GH1D* hist = (GH1D*)obj2;
	    hist->GetListOfFunctions()->Clear();
	    pad->Modified();
	    pad->Update();
	  }
	}
      }
    } 
  }
    edited=true;
    break;
    
  case kKey_z: {
    if(hists.size() == 1) {
       
      int color = hists.at(0)->GetLineColor() + 1;
      if(color == 9)
	color += 30;
       
      hists.at(0)->SetLineColor(color);
    }
    else {

      std::vector<std::string> colors = {"White","Black","Red","Green","Blue","Yellow","Pink",
					 "Cyan","Dark Green","Purple"};

      unsigned int color = 1;
      for(TH1* h : hists) {
	
	if(color == 10)
	  color += 20;
	if(color == 50)
	  color = 1;

	h->SetLineColor(color);
	
	std::cout << "\nHist " << h->GetName() << ": ";
	if(color < 10) 
	  std::cout << colors.at(color) << std::endl;
	else 
	  std::cout << "Some weird color" << std::endl;
	
	color++;
      }
    }
    edited = true;

  } break;
     
  };
  return edited;
}

bool GCanvas::Process1DMousePress(Int_t, Int_t, Int_t)
{
  bool edited = false;
  return edited;
}

bool GCanvas::Process2DArrowKeyPress(Event_t*, UInt_t* keysym)
{
  /// Moves displayed 2D histograms by 50% of the visible range left, right, up, or down
	
  bool              edited = false;
  std::vector<TH1*> hists  = FindHists(2);
  if(hists.empty()) {
    return edited;
  }

  int firstX = hists.at(0)->GetXaxis()->GetFirst();
  int lastX  = hists.at(0)->GetXaxis()->GetLast();
  int firstY = hists.at(0)->GetYaxis()->GetFirst();
  int lastY  = hists.at(0)->GetYaxis()->GetLast();

  int minX = std::min(firstX, 0);
  int maxX = std::max(lastX, hists.at(0)->GetXaxis()->GetNbins() + 1);
  int minY = std::min(firstY, 0);
  int maxY = std::max(lastY, hists.at(0)->GetYaxis()->GetNbins() + 1);

  int xdiff = lastX - firstX;
  int mxdiff = maxX - minX - 2;
  int ydiff = lastY - firstY;
  int mydiff = maxY - minY - 2;

  switch(*keysym) {
  case kMyArrowLeft: {
    if(mxdiff > xdiff) {
      if(firstX == (minX + 1)) {
	//
      } else if((firstX - (xdiff / 2)) < minX) {
	firstX = minX + 1;
	lastX  = minX + (xdiff) + 1;
      } else {
	firstX = firstX - (xdiff / 2);
	lastX  = lastX  - (xdiff / 2);
      }
    }
    for(auto& hist : hists) {
      hist->GetXaxis()->SetRange(firstX, lastX);
    }

    edited = true;
  } break;
  case kMyArrowRight: {
    if(mxdiff > xdiff) {
      if(lastX == (maxX - 1)) {
	//
      } else if((lastX + (xdiff / 2)) > maxX) {
	firstX = maxX - 1 - (xdiff);
	lastX  = maxX - 1;
      } else {
	lastX  = lastX  + (xdiff / 2);
	firstX = firstX + (xdiff / 2);
      }
    }
    for(auto& hist : hists) {
      hist->GetXaxis()->SetRange(firstX, lastX);
    }

    edited = true;
  } break;

  case kMyArrowUp: {
    if(mydiff > ydiff) {
      if(lastY == (maxY - 1)) {
	//
      } else if((lastY + (ydiff / 2)) > maxY) {
	firstY = maxY - 1 - ydiff;
	lastY  = maxY - 1;
      } else {
	firstY = firstY + (ydiff / 2);
	lastY  = lastY  + (ydiff / 2);
      }
    }
    for(auto& hist : hists) {
      hist->GetYaxis()->SetRange(firstY, lastY);
    }

    edited = true;
  } break;

  case kMyArrowDown: {
    if(mydiff > ydiff) {
      if(firstY == (minY + 1)) {
	//
      } else if((firstY - (ydiff / 2)) < minY) {
	firstY = minY + 1;
	lastY  = minY + (ydiff) + 1;
      } else {
	firstY = firstY - (ydiff / 2);
	lastY  = lastY  - (ydiff / 2);
      }
    }
    for(auto& hist : hists) {
      hist->GetYaxis()->SetRange(firstY, lastY);
    }

    edited = true;
  } break;
  default: printf("keysym = %i\n", *keysym); break;
  }
  return edited;
}

bool GCanvas::Process2DKeyboardPress(Event_t*, UInt_t* keysym)
{
  bool edited = false;
  // printf("2d hist key pressed.\n");
  std::vector<TH1*> hists = FindHists(2);
  if(hists.empty()) {
    return edited;
  }
  switch(*keysym) {

  case kKey_1:
    {
      hists.at(0)->SetMinimum(1);
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->SetMinimum(1);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited=true;
    break;

  case kKey_a:
    GetContextMenu()->Action(hists.back()->GetXaxis(),hists.back()->GetXaxis()->Class()->GetMethodAny("SetRangeUser"));
    {
      double x1 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetFirst());
      double x2 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetLast());
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->GetXaxis()->SetRangeUser(x1,x2);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited = true;
    break;

  case kKey_c:
    { 
      int xbins_before = hists.back()->GetNbinsX();
      GetContextMenu()->Action(hists.back(),hists.back()->Class()->GetMethodAny("RebinX"));
      int xbins_after  = hists.back()->GetNbinsX();
      int bin_factor = xbins_before/xbins_after;
      
      TH1 *start = hists.back();

      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      if(hist==start) continue;
	      hist->RebinX(bin_factor);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited = true;
    break;

  case kKey_C:
    {
      int xbins_before = hists.back()->GetNbinsY();
      GetContextMenu()->Action((GH2D*)hists.back(),((GH2D*)hists.back())->Class()->GetMethodAny("RebinY"));
      int xbins_after  = hists.back()->GetNbinsY();
      int bin_factor = xbins_before/xbins_after;
   
      TH1 *start = hists.back();

      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH2::Class())) {
	      GH2D* hist = (GH2D*)obj2;
	      if(hist==start) continue;
	      hist->RebinY(bin_factor);
	      pad->Modified();
	      pad->Update();
	    }
	    
	    else if(obj2->InheritsFrom(GH2D::Class())) {
	      GH2D* hist = (GH2D*)obj2;
	      if(hist==start) continue;
	      hist->RebinY(bin_factor);
	      pad->Modified();
	      pad->Update();
	    }
	    
	  }
	}
      }
       
    }
    edited = true;
    break;

  case kKey_e:
    if(GetNMarkers() < 2) {
      break;
    }
    {
      double x1 = fMarkers.at(fMarkers.size() - 1)->GetLocalX();
      double y1 = fMarkers.at(fMarkers.size() - 1)->GetLocalY();
      double x2 = fMarkers.at(fMarkers.size() - 2)->GetLocalX();
      double y2 = fMarkers.at(fMarkers.size() - 2)->GetLocalY();
      if(x1 > x2) {
	std::swap(x1, x2);
      }
      if(y1 > y2) {
	std::swap(y1, y2);
      }
      for(auto& hist : hists) {
	hist->GetXaxis()->SetRangeUser(x1, x2);
	hist->GetYaxis()->SetRangeUser(y1, y2);
      }
    }
    edited = true;
    RemoveMarker("all");
    break;

  case kKey_E:
    // GetListOfPrimitives()->Print();
    GetContextMenu()->Action(hists.back()->GetXaxis(),
			     hists.back()->GetXaxis()->Class()->GetMethodAny("SetRangeUser"));
    {
      double x1 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetFirst());
      double x2 = hists.back()->GetXaxis()->GetBinCenter(hists.back()->GetXaxis()->GetLast());
      TIter  iter(GetListOfPrimitives());
      while(TObject* obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad* pad = static_cast<TPad*>(obj);
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject* obj2 = iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = static_cast<TH1*>(obj2);
	      hist->GetXaxis()->SetRangeUser(x1, x2);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }

    edited = true;
    break;

  case kKey_g:
    if(GetNMarkers() < 2) {
      break;
    }
    {
      static int cutcounter = 0;
      GCutG*      cut        = new GCutG(Form("_cut%i", cutcounter++), 9);
      // cut->SetVarX("");
      // cut->SetVarY("");
      //
      double x1 = fMarkers.at(fMarkers.size() - 1)->GetLocalX();
      double y1 = fMarkers.at(fMarkers.size() - 1)->GetLocalY();
      double x2 = fMarkers.at(fMarkers.size() - 2)->GetLocalX();
      double y2 = fMarkers.at(fMarkers.size() - 2)->GetLocalY();
      if(x1 > x2) {
	std::swap(x1, x2);
      }
      if(y1 > y2) {
	std::swap(y1, y2);
      }
      double xdist = (x2 - x1) / 2.0;
      double ydist = (y2 - y1) / 2.0;
      //
      //
      cut->SetPoint(0, x1, y1);
      cut->SetPoint(1, x1, y1 + ydist);
      cut->SetPoint(2, x1, y2);
      cut->SetPoint(3, x1 + xdist, y2);
      cut->SetPoint(4, x2, y2);
      cut->SetPoint(5, x2, y2 - ydist);
      cut->SetPoint(6, x2, y1);
      cut->SetPoint(7, x2 - xdist, y1);
      cut->SetPoint(8, x1, y1);
      cut->SetLineColor(kBlack);
      hists.at(0)->GetListOfFunctions()->Add(cut);

      TGRSIint::instance()->LoadGCutG(cut);
    }
    edited = true;
    RemoveMarker("all");
    break;

  case kKey_i:
    {
      std::vector<std::string> names; 
      TIter iter(gROOT->GetListOfSpecials());
      while(TObject *obj =iter.Next()) {
	if(!(obj->InheritsFrom(GCutG::Class())))
	  continue;
	   
	std::string name = obj->GetName();

	bool good = true;
	for(std::string nm : names) {
	  if(!strcmp(nm.c_str(),name.c_str())) {
	    good = false;
	    break;
	  }
	}
	   
	if(!good)
	  continue;

	names.push_back(name);

	double sum = ((GCutG*)obj)->IntegralHist((TH2*)hists.at(0));
	std::string tag = ((GCutG*)obj)->GetTag();
	 
	printf(CYAN "\n  %s %s: %.02f" RESET_COLOR "\n",tag.c_str(),name.c_str(),sum);
	 
      }
    }
    break;

  case kKey_l:
    {
      this->SetLogz(!this->GetLogz());
      gPad->Update();
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  pad->SetLogz(!pad->GetLogz());
	}
      }
    }
    edited = true;
    break;

  case kKey_m:
    {
    
      GetContextMenu()->Action(hists.back(),hists.back()->Class()->GetMethodAny("SetMinimum"));
      double min = hists.back()->GetMinimum();
      gPad->Update();
      
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->SetMinimum(min);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited=true;
    break;

  case kKey_M:
    {
      GetContextMenu()->Action(hists.back(),hists.back()->Class()->GetMethodAny("SetMaximum"));
      double max = hists.back()->GetMaximum();
      gPad->Update();

      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->SetMaximum(max);
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited=true;
    break;

  case kKey_n:
    RemoveMarker("all");
    RemovePeaks(hists.data(),hists.size());
    this->Clear();
    hists.at(0)->Draw("colz");
    for(unsigned int i=1;i<hists.size();i++) {
      hists.at(i)->Draw("colz same");
    }
    edited=true;
    break;
   
  case kKey_o:

    {

      hists.at(0)->GetXaxis()->UnZoom();
      hists.at(0)->GetYaxis()->UnZoom();
      
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->GetXaxis()->UnZoom();
	      hist->GetYaxis()->UnZoom();
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
    }
    edited=true;
    break;
   
  case kKey_P: {
    GH2D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH2Base::Class())) {
	ghist = static_cast<GH2D*>(hist);
	break;
      }
    }

    if((ghist != nullptr) && (ghist->GetProjections()->GetSize() != 0)) {
      ghist->GetProjections()->At(0)->Draw("");
      edited = true;
    }
  } break;
  case kKey_r:
    GetContextMenu()->Action(hists.back()->GetYaxis(),hists.back()->GetYaxis()->Class()->GetMethodAny("SetRangeUser"));
    {
      double y1 = hists.back()->GetYaxis()->GetBinCenter(hists.back()->GetYaxis()->GetFirst());
      double y2 = hists.back()->GetYaxis()->GetBinCenter(hists.back()->GetYaxis()->GetLast());
      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      if(hist->GetDimension()>1) {
		hist->GetYaxis()->SetRangeUser(y1,y2);
		pad->Modified();
		pad->Update();
	      }
	    }
	  }
	}
      }
    }
    edited = true;
    break;
    /*
      case kKey_s: 
      {
      TDirectory* oldDir = gDirectory;
      TString defaultName = "CutFile.cuts";
      char* fileName = new char[256];
      new TGInputDialog(nullptr, static_cast<TRootCanvas*>(GetCanvasImp()), "Enter file name to save cuts to", defaultName, fileName);
      if(strlen(fileName) == 0) {
      break;
      }
      TFile f(fileName, "update");
      if(!f.IsOpen()) {
      std::cout<<RESET_COLOR<<"Failed to open file '"<<fileName<<"', not saving cuts!"<<std::endl;
      break;
      }
      std::cout<<RESET_COLOR<<"Writing the following cuts to '"<<fileName<<"':"<<std::endl;
      for(auto cut : fCuts) {
      std::cout<<cut->GetName()<<std::endl;
      cut->Write();
      }
      f.Close();
      delete[] fileName;
      oldDir->cd();
      }
      break;
    */

  case kKey_s:
    {

      hists.back()->SetDrawOption("scat");

      TIter iter(this->GetListOfPrimitives());
      while(TObject *obj = iter.Next()) {
	if(obj->InheritsFrom(TPad::Class())) {
	  TPad *pad = (TPad*)obj;
	  pad->cd();
	  TIter iter2(pad->GetListOfPrimitives());
	  while(TObject *obj2=iter2.Next()) {
	    if(obj2->InheritsFrom(TH1::Class())) {
	      TH1* hist = (TH1*)obj2;
	      hist->SetDrawOption("scat");
	      pad->Modified();
	      pad->Update();
	    }
	  }
	}
      }
      
    }
    edited=true;
    break;
    
  case kKey_x: {
    GH2D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH2Base::Class())) {
	ghist = static_cast<GH2D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      ghist->SetSummary(false);
      TH1* phist = ghist->ProjectionX(); //->Draw();
      if(phist != nullptr) {
	new GCanvas();
	phist->Draw("");
      }
      edited = true;
    }
  } break;

  case kKey_X: {
    GH2D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH2Base::Class())) {
	ghist = static_cast<GH2D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      ghist->SetSummary(true);
      ghist->SetSummaryDirection(EDirection::kYDirection);
      TH1* phist = ghist->GetNextSummary(nullptr, false);
      if(phist != nullptr) {
	new GCanvas();
	phist->Draw("");
      }
      edited = true;
    }
  } break;

  case kKey_y: {
    GH2D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH2Base::Class())) {
	ghist = static_cast<GH2D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      ghist->SetSummary(false);
      // printf("ghist = 0x%08x\n",ghist);
      TH1* phist = ghist->ProjectionY(); //->Draw();
      // printf("phist = 0x%08x\n",phist);
      // printf("phist->GetName() = %s\n",phist->GetName());
      if(phist != nullptr) {
	new GCanvas();
	phist->Draw("");
      }
      edited = true;
    }
  } break;

  case kKey_Y: {
    GH2D* ghist = nullptr;
    for(auto hist : hists) {
      if(hist->InheritsFrom(GH2Base::Class())) {
	ghist = static_cast<GH2D*>(hist);
	break;
      }
    }

    if(ghist != nullptr) {
      ghist->SetSummary(true);
      ghist->SetSummaryDirection(EDirection::kXDirection);
      // TH1* phist = ghist->SummaryProject(1);
      TH1* phist = ghist->GetNextSummary(nullptr, false);
      if(phist != nullptr) {
	new GCanvas();
	phist->Draw("");
      }
      edited = true;
    }
  } break;
   
  case kKey_z: {
    hists.back()->SetDrawOption("colz");

    TIter iter(this->GetListOfPrimitives());
    while(TObject *obj = iter.Next()) {
      if(obj->InheritsFrom(TPad::Class())) {
	TPad *pad = (TPad*)obj;
	pad->cd();
	TIter iter2(pad->GetListOfPrimitives());
	while(TObject *obj2=iter2.Next()) {
	  if(obj2->InheritsFrom(TH1::Class())) {
	    TH1* hist = (TH1*)obj2;
	    hist->SetDrawOption("colz");
	    pad->Modified();
	    pad->Update();
	  }
	}
      }
    }
  }
    edited=true;
    break;

  };
  return edited;
}

bool GCanvas::Process2DMousePress(Int_t, Int_t, Int_t)
{
  bool edited = false;
  return edited;
}
