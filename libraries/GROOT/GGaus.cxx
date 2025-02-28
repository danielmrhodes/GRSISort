#include "GGaus.h"
#include "TGraph.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "Globals.h"
#include "GRootFunctions.h"
#include "GCanvas.h"

ClassImp(GGaus)

GGaus::GGaus(Double_t xlow, Double_t xhigh, Option_t*)
   : TF1("gausbg", "gaus(0)+pol1(3)", xlow, xhigh), fBGFit("background", "pol1", xlow, xhigh)
{
   Clear("");
   if(xlow > xhigh) {
      std::swap(xlow, xhigh);
   }

   TF1::SetRange(xlow, xhigh);
   TF1::SetNpx(10000);

   fBGFit.SetNpx(1000);
   fBGFit.SetLineStyle(2);
   fBGFit.SetLineColor(kBlack);

   // Changing the name here causes an infinite loop when starting the FitEditor
   // SetName(Form("gaus_%d_to_%d",(Int_t)(xlow),(Int_t)(xhigh)));
   InitNames();
   // TF1::SetParameter("centroid",cent);
}

GGaus::GGaus(Double_t xlow, Double_t xhigh, TF1* bg, Option_t*) : TF1("gausbg", "gaus(0)+pol1(3)", xlow, xhigh)
{
   Clear("");
   if(xlow > xhigh) {
      std::swap(xlow, xhigh);
   }
   TF1::SetRange(xlow, xhigh);
   // Changing the name here causes an infinite loop when starting the FitEditor
   // SetName(Form("gaus_%d_to_%d",(Int_t)(xlow),(Int_t)(xhigh)));
   InitNames();

   if(bg != nullptr) {
      fBGFit.Clear();
      fBGFit.Copy(*bg);
   } else {
      fBGFit = TF1("BGFit", "pol1", xlow, xhigh);
   }

   fBGFit.SetNpx(1000);
   fBGFit.SetLineStyle(2);
   fBGFit.SetLineColor(kBlack);
}

GGaus::GGaus() : TF1("gausbg", "gaus(0)+pol1(3)", 0, 1000), fBGFit("background", "pol1", 0, 1000)
{
   Clear();
   InitNames();
   fBGFit.SetNpx(1000);
   fBGFit.SetLineStyle(2);
   fBGFit.SetLineColor(kBlack);
}

GGaus::GGaus(const GGaus& peak) : TF1(peak)
{
   peak.Copy(*this);
}

GGaus::~GGaus()
{
   // if(background)
   //  delete background;
}

// void GGaus::Fcn(Int_t &npar,Double_t *gin,Double_T &f,Double_t *par,Int_t iflag) {
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

void GGaus::InitNames()
{
   TF1::SetParName(0, "Height");
   TF1::SetParName(1, "centroid");
   TF1::SetParName(2, "sigma");
   TF1::SetParName(3, "bg_offset");
   TF1::SetParName(4, "bg_slope");
}

void GGaus::Copy(TObject& obj) const
{
   // printf("0x%08x\n",&obj);
   // fflush(stdout);
   // printf("%s\n",obj.GetName());
   // fflush(stdout);

   TF1::Copy(obj);
   (static_cast<GGaus&>(obj)).init_flag = init_flag;
   (static_cast<GGaus&>(obj)).fArea     = fArea;
   (static_cast<GGaus&>(obj)).fDArea    = fDArea;
   (static_cast<GGaus&>(obj)).fSum      = fSum;
   (static_cast<GGaus&>(obj)).fDSum     = fDSum;
   (static_cast<GGaus&>(obj)).fChi2     = fChi2;
   (static_cast<GGaus&>(obj)).fNdf      = fNdf;

   fBGFit.Copy(((static_cast<GGaus&>(obj)).fBGFit));
}

bool GGaus::InitParams(TH1* fithist)
{
   if(fithist == nullptr) {
      printf("No histogram is associated yet, no initial guesses made\n");
      return false;
   }
   // printf("%s called.\n",__PRETTY_FUNCTION__); fflush(stdout);
   // Makes initial guesses at parameters for the fit. Uses the histogram to
   Double_t xlow, xhigh;
   GetRange(xlow, xhigh);

   Int_t binlow  = fithist->GetXaxis()->FindBin(xlow);
   Int_t binhigh = fithist->GetXaxis()->FindBin(xhigh);
   fNdf = binhigh - binlow - 5;
   
   double x1 = fithist->GetBinLowEdge(binlow);
   double x2 = fithist->GetBinLowEdge(binhigh+1);
   SetRange(x1,x2);
   
   /*
   Double_t y1 = fithist->GetBinContent(binlow);
   Double_t y2  = fithist->GetBinContent(binhigh);
   for(int x = 1; x < 5; x++) {
      y1 += fithist->GetBinContent(binlow - x);
      y2 += fithist->GetBinContent(binhigh + x);
   }
   y1 = y1 / 5.0;
   y2  = y2 / 5.0;

   double slope_guess = (y2-y1)/(x2-x1);
   double offset_guess = y1 - slope_guess*x1;
   */
   
   double largestx = 0.0;
   double largesty = 0.0;
   for(int i = binlow; i <= binhigh; i++) {
      if(fithist->GetBinContent(i) > largesty) {
         largesty = fithist->GetBinContent(i);
         largestx = fithist->GetXaxis()->GetBinCenter(i);
      }
   } 
   
   // - par[0]: height of peak
   // - par[1]: cent of peak
   // - par[2]: sigma
   // - par[3]: bg constant
   // - par[4]: bg slope

   // limits.
   TF1::SetParLimits(0, 0, largesty * 2);
   TF1::SetParLimits(1, xlow+2., xhigh-2.);
   TF1::SetParLimits(2, 0, xhigh - xlow);
   // TF1::SetParLimits(3,0.0,40);
   // TF1::SetParLimits(4,0.01,5);

   // Make initial guesses
   TF1::SetParameter(0, largesty);                // fithist->GetBinContent(bin));
   TF1::SetParameter(1, largestx);                // GetParameter("centroid"));
   TF1::SetParameter(2, (largestx * .01) / 2.35); // 2,(xhigh-xlow));     //2.0/binWidth); //
   //TF1::SetParameter(3, offset_guess);
   //TF1::SetParameter(4, slope_guess);

   //TF1::SetParError(0, 0.10 * largesty);
   //TF1::SetParError(1, 0.25);
   //TF1::SetParError(2, 0.10 * ((largestx * .01) / 2.35));
   // TF1::SetParError(3,5);
   // TF1::SetParError(4,0.5);

   // TF1::Print();

   SetInitialized();
   return true;
}

Bool_t GGaus::Fit(TH1* fithist, Option_t* opt)
{
   if(fithist == nullptr) {
      return false;
   }
   TString options = opt;
   if(!IsInitialized()) {
      InitParams(fithist);
   }
   TVirtualFitter::SetMaxIterations(100000);

   bool verbose = !options.Contains("Q");
   bool noprint = options.Contains("no-print");
   bool no_area = options.Contains("no-area");
   if(noprint) {
      options.ReplaceAll("no-print","");
      verbose = false;
   }
   if(no_area)
     options.ReplaceAll("no-area","");

   if(fithist->GetSumw2()->fN != fithist->GetNbinsX() + 2)
      fithist->Sumw2();
   
   TFitResultPtr fitres = fithist->Fit(this, Form("%sQRS+", options.Data()));
   if(!fitres.Get()->IsValid()) {

     if(verbose)
       printf(RED "fit has failed, trying refit... " RESET_COLOR);

     fithist->GetListOfFunctions()->Last()->Delete();
     fitres = fithist->Fit(this, Form("%sQRSME+", options.Data()));
     
     if(fitres.Get()->IsValid()) {
       if(verbose && !noprint)
	 printf(DGREEN " refit passed!" RESET_COLOR "\n");
     }
     else {
       if(verbose && !noprint) 
	 printf(DRED " refit also failed :( " RESET_COLOR "\n");
     }
   }
   
   fChi2 = fithist->Chisquare(this,"R");
   
   // if(fitres->ParError(2) != fitres->ParError(2)) { // checks if nan.
   //  if(fitres->Parameter(3)<1) {
   //    FixParameter(4,0);
   //    FixParameter(3,1);
   // printf("Beta may have broken the fit, retrying with R=0);
   //    fithist->GetListOfFunctions()->Last()->Delete();
   //    fitres = fithist->Fit(this,Form("%sRSM",options.Data()));
   //  }
   //}

   // TF1::Print();

   Double_t xlow, xhigh;
   TF1::GetRange(xlow, xhigh);

   // Make a function that does not include the background
   // Intgrate the background.
   GGaus *tmppeak = new GGaus;
   Copy(*tmppeak);

   tmppeak->SetParameter("bg_offset",0.0);
   tmppeak->SetParameter("bg_slope",0.0);
   tmppeak->SetRange(xlow,xhigh);//This will help get the true area of the gaussian 200 ~ infinity in a gaus
   tmppeak->SetName("tmppeak");

   //fArea = (tmppeak->Integral(xlow,xhigh))/fithist->GetBinWidth(1);
   TMatrixDSym CovMat = fitres->GetCovarianceMatrix();
   CovMat(3,3) = 0.0;
   CovMat(4,4) = 0.0;

   if(no_area)
     fDArea = 0.0;
   else
     fDArea = tmppeak->IntegralError(xlow,xhigh,tmppeak->GetParameters(),
				     CovMat.GetMatrixArray())/fithist->GetBinWidth(1);

   double bgpars[2];
   bgpars[0] = TF1::GetParameters()[3];
   bgpars[1] = TF1::GetParameters()[4];
   // bgpars[5] = TF1::GetParameters()[7];

   fBGFit.SetParameters(bgpars);
   // fithist->GetListOfFunctions()->Print();

   double bgArea = 0.0;
   if(no_area)
     fArea = 0.0;
   else {
     fArea         = Integral(xlow, xhigh) / fithist->GetBinWidth(1);
     bgArea = fBGFit.Integral(xlow, xhigh) / fithist->GetBinWidth(1);
     fArea -= bgArea;
   }
   
   if(xlow > xhigh)
     std::swap(xlow, xhigh);
   
   fSum = fithist->Integral(fithist->GetXaxis()->FindBin(xlow),
                            fithist->GetXaxis()->FindBin(xhigh));
   if(!noprint)
     printf(RESET_COLOR "\nsum between markers: %02f\n", fSum);
   
   fDSum = TMath::Sqrt(fSum);
   fSum -= bgArea;
   if(!noprint && !no_area)
     printf("sum after subtraction: %02f\n", fSum);
   
   if(verbose)
     Print();
   
   fithist->GetListOfFunctions()->Add(Background());
   
   return true;
}

void GGaus::Clear(Option_t* opt)
{
   TString options = opt;
   // Clear the GGaus including functions and histogram
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

void GGaus::Print(Option_t* opt) const
{
   TString options = opt;
   printf(GREEN);
   printf("Name:      %s \n", GetName());
   printf("Centroid:  %1f +/- %1f \n", GetParameter("centroid"), GetParError(GetParNumber("centroid")));
   printf("Sigma:     %1f +/- %1f \n", GetParameter("sigma"), GetParError(GetParNumber("sigma")));
   printf("FWHM:      %1f +/- %1f \n", GetFWHM(), GetFWHMErr());
   printf("Reso:      %1f%%  \n", GetFWHM() / GetParameter("centroid") * 100.);
   printf("Area:      %1f +/- %1f \n", fArea, fDArea);
   printf("Sum:       %1f +/- %1f \n", fSum, fDSum);
   printf("Chi^2/NDF: %1f/%1f = %1f\n", fChi2, fNdf, fChi2/double(fNdf));
   if(options.Contains("all")) {
      TF1::Print(opt);
   }
   printf(RESET_COLOR);
   printf("\n");
}

void GGaus::DrawResiduals(TH1* hist) const
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
