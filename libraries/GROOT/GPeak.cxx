
#include <GPeak.h>
#include <TGraph.h>
#include <TVirtualFitter.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "Globals.h"
#include "GRootFunctions.h"
#include "TGRSIFunctions.h"
#include "GCanvas.h"

/// \cond CLASSIMP
ClassImp(GPeak)
/// \endcond

GPeak* GPeak::fLastFit = nullptr;

GPeak::GPeak(Double_t cent, Double_t xlow, Double_t xhigh, Option_t*)
   : TF1("photopeakbg", GRootFunctions::PhotoPeakBG, xlow, xhigh, 7),
     fBGFit("background", GRootFunctions::StepBG, xlow, xhigh, 6)
{
   Clear("");
   if(cent > xhigh || cent < xlow) {
      // out of range...
      if(xlow > cent) {
         std::swap(xlow, cent);
      }
      if(xlow > xhigh) {
         std::swap(xlow, xhigh);
      }
      if(cent > xhigh) {
         std::swap(cent, xhigh);
      }
   }

   TF1::SetRange(xlow, xhigh);
   TF1::SetNpx(10000);
   
   fBGFit.SetNpx(10000);
   fBGFit.SetLineStyle(2);
   fBGFit.SetLineColor(kBlack);

   SetName(Form("Chan%d_%d_to_%d", static_cast<Int_t>(cent), static_cast<Int_t>(xlow), static_cast<Int_t>(xhigh)));
   InitNames();
   TF1::SetParameter("centroid", cent);

   SetParent(nullptr);
   // TF1::SetDirectory(0);
   fBGFit.SetParent(nullptr);
   fBGFit.SetBit(TObject::kCanDelete, false);
   // fBGFit.SetDirectory(0);
}

GPeak::GPeak(Double_t cent, Double_t xlow, Double_t xhigh, TF1* bg, Option_t*)
   : TF1("photopeakbg", GRootFunctions::PhotoPeakBG, xlow, xhigh, 7)
{
   Clear("");
   if(cent > xhigh || cent < xlow) {
      // out of range...
      if(xlow > cent) {
         std::swap(xlow, cent);
      }
      if(xlow > xhigh) {
         std::swap(xlow, xhigh);
      }
      if(cent > xhigh) {
         std::swap(cent, xhigh);
      }
   }
   TF1::SetRange(xlow, xhigh);
   SetName(Form("Chan%d_%d_to_%d", static_cast<Int_t>(cent), static_cast<Int_t>(xlow), static_cast<Int_t>(xhigh)));
   InitNames();
   TF1::SetParameter("centroid", cent);

   if(bg != nullptr) {
      fBGFit.Clear();
      fBGFit.Copy(*bg);
   } else {
      fBGFit = TF1("BGFit", GRootFunctions::StepBG, xlow, xhigh, 10);
   }

   fBGFit.SetNpx(1000);
   fBGFit.SetLineStyle(2);
   fBGFit.SetLineColor(kBlack);

   SetParent(nullptr);
   // SetDirectory(0);
   fBGFit.SetParent(nullptr);
   // fBGFit.SetDirectory(0);
}

GPeak::GPeak()
   : TF1("photopeakbg", GRootFunctions::PhotoPeakBG, 0, 1000, 10),
     fBGFit("background", GRootFunctions::StepBG, 0, 1000, 10)
{

   Clear();
   InitNames();
   fBGFit.SetNpx(1000);
   fBGFit.SetLineStyle(2);
   fBGFit.SetLineColor(kBlack);

   SetParent(nullptr);
   // SetDirectory(0);
   fBGFit.SetParent(nullptr);
   // fBGFit.SetDirectory(0);
}

GPeak::GPeak(const GPeak& peak) : TF1(peak)
{

   SetParent(nullptr);
   // SetDirectory(0);
   fBGFit.SetParent(nullptr);
   // fBGFit.SetDirectory(0);
   peak.Copy(*this);
}

GPeak::~GPeak()
{
   // gROOT->RecursiveRemove(&fBGFit);
   // gROOT->RecursiveRemove(this);
   // if(background)
   //  delete background;
}

// void GPeak::Fcn(Int_t &npar,Double_t *gin,Double_T &f,Double_t *par,Int_t iflag) {
// chisquared calculator
//
//  int i=0;
//  double chisq = 0;
//  double delta = 0;
//  for(i=0;i<nbins;i++) {
//    delta = (data[i] - GRootFunctions::PhotoPeakBG((x+i),par))/error[i];
//    chisq += delta*delta;
//  }
//  f=chisq;
//}

void GPeak::InitNames()
{
   TF1::SetParName(0, "Height");
   TF1::SetParName(1, "centroid");
   TF1::SetParName(2, "sigma");
   TF1::SetParName(3, "R");
   TF1::SetParName(4, "beta");
   // TF1::SetParName(5,"step");
   // TF1::SetParName(6,"A");
   // TF1::SetParName(7,"B");
   // TF1::SetParName(8,"C");
   TF1::SetParName(5, "step");
   TF1::SetParName(6, "bg_offset");
   // TF1::SetParName(7,"bg_slope");
}

void GPeak::Copy(TObject& obj) const
{
   // printf("0x%08x\n",&obj);
   // fflush(stdout);
   // printf("%s\n",obj.GetName());
   // fflush(stdout);

   TF1::Copy(obj);
   (static_cast<GPeak&>(obj)).init_flag = init_flag;
   (static_cast<GPeak&>(obj)).fArea     = fArea;
   (static_cast<GPeak&>(obj)).fDArea    = fDArea;
   (static_cast<GPeak&>(obj)).fSum      = fSum;
   (static_cast<GPeak&>(obj)).fDSum     = fDSum;
   (static_cast<GPeak&>(obj)).fChi2     = fChi2;
   (static_cast<GPeak&>(obj)).fNdf      = fNdf;

   fBGFit.Copy(((static_cast<GPeak&>(obj)).fBGFit));
}

bool GPeak::InitParams(TH1* fithist)
{
   if(fithist == nullptr) {
      printf("No histogram is associated yet, no initial guesses made\n");
      return false;
   }
   // printf("%s called.\n",__PRETTY_FUNCTION__); fflush(stdout);
   // Makes initial guesses at parameters for the fit. Uses the histogram to
   Double_t xlow, xhigh;
   GetRange(xlow, xhigh);

   // Int_t bin = fithist->GetXaxis()->FindBin(GetParameter("centroid"));
   Int_t binlow  = fithist->GetXaxis()->FindBin(xlow);
   Int_t binhigh = fithist->GetXaxis()->FindBin(xhigh);

   Double_t highy = fithist->GetBinContent(binlow);
   Double_t lowy  = fithist->GetBinContent(binhigh);
   for(int x = 1; x < 5; x++) {
      highy += fithist->GetBinContent(binlow - x);
      lowy += fithist->GetBinContent(binhigh + x);
   }
   highy = highy / 5.0;
   lowy  = lowy / 5.0;

   //  Double_t yhigh  = fithist->GetBinContent(binhigh);
   //  Double_t ylow   = fithist->GetBinContent(binlow);
   if(lowy > highy) {
      std::swap(lowy, highy);
   }

   double largestx = 0.0;
   double largesty = 0.0;
   int    i        = binlow;
   for(; i <= binhigh; i++) {
      if(fithist->GetBinContent(i) > largesty) {
         largesty = fithist->GetBinContent(i);
         largestx = fithist->GetXaxis()->GetBinCenter(i);
      }
   }

   // - par[0]: height of peak
   // - par[1]: cent of peak
   // - par[2]: sigma
   // - par[3]: R:    relative height of skewed gaus to gaus
   // - par[4]: beta: "skewedness" of the skewed gaussin
   // - par[5]: step: size of stepfunction step.

   // - par[6]: base bg height.

   // limits.
   TF1::SetParLimits(0, 0, largesty * 2);
   TF1::SetParLimits(1, xlow, xhigh);
   TF1::SetParLimits(2, 0.1, xhigh - xlow);
   TF1::SetParLimits(3, 0.0, 40);
   TF1::SetParLimits(4, 0.01, 5);
   double step = ((highy - lowy) / largesty) * 50;

   // TF1::SetParLimits(5,step-step*.1,step+.1*step);
   TF1::SetParLimits(5, 0.0, step + step);

   // double slope  = (yhigh-ylow)/(xhigh-xlow);
   // double offset = yhigh-slope*xhigh;
   double offset = lowy;
   TF1::SetParLimits(6, offset - 0.5 * offset, offset + offset);
   // TF1::SetParLimits(7,-2*slope,2*slope);

   // Make initial guesses
   TF1::SetParameter(0, largesty);                // fithist->GetBinContent(bin));
   TF1::SetParameter(1, largestx);                // GetParameter("centroid"));
   TF1::SetParameter(2, (largestx * .01) / 2.35); // 2,(xhigh-xlow));     //2.0/binWidth); //
   TF1::SetParameter(3, 5.);
   TF1::SetParameter(4, 1.);
   TF1::SetParameter(5, step);
   TF1::SetParameter(6, offset);
   // TF1::SetParameter(7,slope);

   TF1::SetParError(0, 0.10 * largesty);
   TF1::SetParError(1, 0.25);
   TF1::SetParError(2, 0.10 * ((largestx * .01) / 2.35));
   TF1::SetParError(3, 5);
   TF1::SetParError(4, 0.5);
   TF1::SetParError(5, 0.10 * step);
   TF1::SetParError(6, 0.10 * offset);

   // TF1::Print();

   SetInitialized();
   return true;
}

Bool_t GPeak::Fit(TH1* fithist, Option_t* opt)
{
   if(fithist == nullptr) {
      return false;
   }
   TString options = opt;
   options.ToLower();
   if(!IsInitialized()) {
      InitParams(fithist);
   }
   TVirtualFitter::SetMaxIterations(100000);

   bool verbose = !options.Contains("q");
   bool noprint = options.Contains("no-print");
   if(noprint) {
      options.ReplaceAll("no-print", "");
      verbose = false;
   }
   
   bool retryFit = options.Contains("retryfit");
   options.ReplaceAll("retryfit", "");

   if(fithist->GetSumw2()->fN != fithist->GetNbinsX() + 2) {
      fithist->Sumw2();
   }

   TFitResultPtr fitres = fithist->Fit(this, Form("%sQLRSE", options.Data()));

   if(verbose && !noprint)
     printf("chi^2/NDF = %.02f\n", GetChisquare() / static_cast<double>(GetNDF()));

   if(!fitres.Get()->IsValid()) {
     if(!noprint)
       printf(RED "fit has failed, trying refit... " RESET_COLOR);
     
     fithist->GetListOfFunctions()->Last()->Delete();
     fitres = fithist->Fit(this, Form("%sQLRSME", options.Data()));
     
      if(fitres.Get()->IsValid()) {
	if(!noprint)
	  printf(DGREEN " refit passed!" RESET_COLOR "\n");
      }
      else {
	if(!noprint)
	  printf(DRED " refit also failed :( " RESET_COLOR "\n");
      }
   }

   // check parameter errors and re-try using minos instead of minuit
   if(!TGRSIFunctions::CheckParameterErrors(fitres)) {
     if(!noprint)
       std::cout<<YELLOW
		<<"Re-fitting with \"E\" option to get better error estimation using Minos technique."
		<<RESET_COLOR<<std::endl;
      
      fitres = fithist->Fit(this, Form("%sQLRSME", options.Data()));
   }
   
   // check the parameter errors again and see if we need to fit again with all parameters released
   if(!TGRSIFunctions::CheckParameterErrors(fitres, "q") && retryFit) {
     
     if(!noprint)
       std::cout<<GREEN<<"Re-fitting with released parameters (without any limits):"
		<<RESET_COLOR<<std::endl;
     
      for(int i = 0; i < GetNpar(); ++i) {
         ReleaseParameter(i);
      }
      
      fitres = fithist->Fit(this, Form("%sQLRSM", options.Data()));
      TGRSIFunctions::CheckParameterErrors(fitres);
   }

   Double_t xlow, xhigh;
   TF1::GetRange(xlow, xhigh);

   double bgpars[5];
   bgpars[0] = TF1::GetParameters()[0];
   bgpars[1] = TF1::GetParameters()[1];
   bgpars[2] = TF1::GetParameters()[2];
   bgpars[3] = TF1::GetParameters()[5];
   bgpars[4] = TF1::GetParameters()[6];

   fBGFit.SetParameters(bgpars);

   fChi2 = GetChisquare();
   fNdf  = GetNDF();

   fArea         = Integral(xlow, xhigh) / fithist->GetBinWidth(1);
   double bgArea = fBGFit.Integral(xlow, xhigh) / fithist->GetBinWidth(1);
   fArea -= bgArea;

   if(xlow > xhigh) {
      std::swap(xlow, xhigh);
   }
   fSum = fithist->Integral(fithist->GetXaxis()->FindBin(xlow),
                            fithist->GetXaxis()->FindBin(xhigh)); //* fithist->GetBinWidth(1);
   if(verbose) printf("sum between markers: %02f\n", fSum);
   fDSum = TMath::Sqrt(fSum);
   fSum -= bgArea;
   if(verbose) printf("sum after subtraction: %02f\n", fSum);

   // Make a function that does not include the background
   // Intgrate the background.
   // TPeak* tmppeak = new TPeak(*this);

   Double_t range_low, range_high;
   GetRange(range_low,range_high);

   GPeak* tmppeak = new GPeak;
   Copy(*tmppeak);
   tmppeak->SetParameter("step", 0.0);
   //tmppeak->SetParameter("A", 0.0);
   //tmppeak->SetParameter("B", 0.0);
   //tmppeak->SetParameter("C", 0.0);
   tmppeak->SetParameter("bg_offset", 0.0);
   tmppeak->SetRange(range_low, range_high); // This will help get the true area of the gaussian 200 ~ infinity in a gaus
   tmppeak->SetName("tmppeak");

   // This is where we will do integrals and stuff.
   TMatrixDSym CovMat = fitres->GetCovarianceMatrix();
   CovMat(5, 5) = 0.0;
   CovMat(6, 6) = 0.0;
   fDArea = (tmppeak->IntegralError(xlow, xhigh, tmppeak->GetParameters(), CovMat.GetMatrixArray())) / fithist->GetBinWidth(1);

   delete tmppeak;

   // print the results of the fit even if not verbose
   if(!noprint)
     Print();

   Copy(*fithist->GetListOfFunctions()->FindObject(GetName()));
   fithist->GetListOfFunctions()->Add(fBGFit.Clone());
   fLastFit = this;

   SetParent(nullptr);

   return true;
}

void GPeak::Clear(Option_t* opt)
{
   TString options = opt;
   // Clear the GPeak including functions and histogram
   if(options.Contains("all")) {
      TF1::Clear();
   }
   init_flag = false;
   fArea     = 0.0;
   fDArea    = 0.0;
   fSum      = 0.0;
   fDSum     = 0.0;
   fChi2     = 0.0;
   fNdf      = 0.0;
}

void GPeak::Print(Option_t* opt) const
{
   TString options = opt;
   printf(GREEN);
   printf("Name: %s \n", GetName());
   printf("Centroid:  %1f +/- %1f \n", GetParameter("centroid"), GetParError(GetParNumber("centroid")));
   printf("Area:      %1f +/- %1f \n", fArea, fDArea);
   printf("Sum:       %1f +/- %1f \n", fSum, fDSum);
   printf("FWHM:      %1f +/- %1f \n", GetFWHM(), GetFWHMErr());
   printf("Reso:      %1f%%  \n", GetFWHM() / GetParameter("centroid") * 100.);
   printf("Chi^2/NDF: %1f\n", fChi2 / fNdf);
   if(options.Contains("all")) {
      TF1::Print(opt);
   }
   printf(RESET_COLOR);
}

void GPeak::DrawResiduals(TH1* hist) const
{
   if(hist == nullptr) {
      return;
   }
   if(fChi2 < 0.000000001) {
      printf("No fit performed\n");
      return;
   }
   Double_t xlow, xhigh;
   GetRange(xlow, xhigh);
   Int_t nbins  = hist->GetXaxis()->GetNbins();
   auto* res    = new Double_t[nbins];
   auto* bin    = new Double_t[nbins];
   Int_t points = 0;
   for(int i = 1; i <= nbins; i++) {
      if(hist->GetBinCenter(i) <= xlow || hist->GetBinCenter(i) >= xhigh) {
         continue;
      }
      res[points] = (hist->GetBinContent(i) - Eval(hist->GetBinCenter(i))) + GetParameter("Height") / 2;
      bin[points] = hist->GetBinCenter(i);
      points++;
   }
   new GCanvas();
   auto* residuals = new TGraph(points, bin, res);
   residuals->Draw("*AC");
   delete[] res;
   delete[] bin;
}
