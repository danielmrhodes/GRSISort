#include "GH1D.h"
#include "TF1Sum.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TNtuple.h"
#include "TMath.h"

#include <fstream>

int REBIN_FACTOR;

const std::string gate_name = "78Kr";
const std::string targ_name = "208Pb";

std::string data_file_name = "/home/drhodes/Analysis/S1937/78Kr/208PbTarget/HistFiles/sum.root";
//std::string data_file_name = "/home/drhodes/Analysis/S1937/78Kr/HistFiles/sum.root";
//std::string data_file_name = "/home/drhodes/Analysis/S1937/78Kr/HistFiles/good/sum.root";

std::string fit_dir = "/home/drhodes/Analysis/S1937/78Kr/208PbTarget/Simulation/Yields/";
//std::string fit_dir = "/home/drhodes/Analysis/S1937/78Kr/Simulation/Yields/";

TF1* fEff = NULL;
double getEff(double energy) {
  
  return fEff->Eval(energy);
  
}

void MakeEffCurve() {

  std::string func = "exp([0]";
  int Npars = 6;
  for(int k=1;k<Npars;k++)
    func += "+([" + std::to_string(k) + "]*(log(x/200))**" + std::to_string(k) + ")";
  func += ")";

  double xlow = 40.0;
  double xhigh = 3600.0;

  if(fEff)
    delete fEff;
  
  fEff = new TF1("fEff",func.c_str(),xlow,xhigh);
  fEff->SetParameters(-1.28339,-0.483539,-0.218486,0.111893,-0.0139730,-0.00271216);
  
  return;
}

int getEnergy(const std::string &name) {
  return std::stoi(name.substr(4,4));
}

void printResults(std::vector<double> fep_counts, std::vector<double> fep_counts_unc,
		  std::vector<int> energies) {
  
  std::vector<double> intensities;
  std::vector<double> intensities_unc;
  
  std::cout << "\n========== COUNTS ================\n";
  for (unsigned int i=0;i<energies.size();i++){
    
    std::cout << energies.at(i) << " keV peak:\t" << fep_counts.at(i) << "\t+- "
	      << fep_counts_unc.at(i) << "\n";
    
    double eff = getEff(energies.at(i));
    intensities.push_back(fep_counts.at(i)/eff);

    double rel_err = sqrt(pow(fep_counts_unc.at(i)/fep_counts.at(i),2.) + 0.015*0.015);
    intensities_unc.push_back(intensities.at(i)*rel_err);
    
  }
  std::cout << "==================================\n";
  
  std::cout << "\n========== INTENSITIES ===========\n";
  for (unsigned int i = 0; i < intensities.size(); i++){
    
    std::cout << energies.at(i) << " keV peak:\t" << intensities.at(i) << "\t+- "
	      << intensities_unc.at(i) << "\n";
    
  }
  std::cout << "==================================\n\n";

}

/*
TF1* ExpBackground(double p0, double p1, bool fix) {

  TF1* bg = new TF1("Single_Exp","[0]*TMath::Exp([1]*x)",0,4000);
  //bg->SetParameter(0,p0);
  //bg->SetParameter(1,p1);
  bg->SetParameters(p0,p1);
  
  if(fix) {
    bg->FixParameter(0,p0);
    bg->FixParameter(1,p1);
  }

  bg->SetRange(0,4000);
  bg->SetNpx(4000/REBIN_FACTOR);

  return bg;
}
*/

TF1* ExpBackground(double p0, double p1, double p2, double p3, int num, bool fix) {

  TF1* bg = NULL;
  if(num == 1) {
    bg = new TF1("Single_Exp","[0]*TMath::Exp([1]*x)",0,4000);
    bg->SetParameters(p0,p1);
  }
  else if(num == 2) {
    bg = new TF1("Double_Exp","[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)",0,4000);
    bg->SetParameters(p0,p1,p2,p3);
    bg->SetParLimits(0,0.0,1000.0);
    bg->SetParLimits(1,-1.0,0.0);
    bg->SetParLimits(2,0.0,1000.0);
    bg->SetParLimits(3,-1.0,0.0);
  }

  double pars[4] = {p0,p1,p2,p3};
  if(fix) 
    for(int i=0;i<bg->GetNpar();i++)
      bg->FixParameter(0,pars[i]);

  bg->SetRange(0,4000);
  bg->SetNpx(4000/REBIN_FACTOR);

  return bg;
}

TF1* TimeRandomBackground(std::string hname1, bool fix) {

  std::string fn1 = "/home/drhodes/Analysis/S1937/HistFiles/InvertedTE/sum.root";
  TFile* file1 = new TFile(fn1.c_str(),"read");
  if(file1->IsZombie()) {
    std::cout << "Time-random background ROOT file not found!" << std::endl;
  }

  std::string fn2 = "/home/drhodes/Analysis/S1937/HistFiles/Off_Prompt/sum.root";
  TFile* file2 = new TFile(fn2.c_str(),"read");
  if(file2->IsZombie()) {
    std::cout << "Off-prompt ROOT file not found!" << std::endl;
  }

  GH1D* hist1 = (GH1D*)file1->Get(hname1.c_str());
  if(!hist1) {
    std::cout << "Histogram not properly retrieved from time-random ROOT file!" << std::endl;
  }
  hist1->Sumw2();

  size_t pos = hname1.find("TE") + 2;
  std::string  hname2 = hname1.substr(0,pos);
  hname2 += "_OP" + hname1.substr(pos);
  
  GH1D* hist2 = (GH1D*)file2->Get(hname2.c_str());
  if(!hist2) {
    std::cout << "Histogram not properly retrieved from off-prompt ROOT file!" << std::endl;
  }
  
  hist1->SetName("time_random");
  hist1->SetTitle("Time Random");
  hist1->Rebin(REBIN_FACTOR);

  double scale = hist2->Integral()/hist1->Integral();
  
  TF1* func = hist1->ConstructTF1();
  func->SetParameter(0,scale);

  if(fix)
    func->FixParameter(0,scale);
  
  return func;
  
}

TF1* Se78Background(std::string hname1, int gate, int range, int nbins, double scale, bool fix) {

  std::string fn1 = "/home/drhodes/Analysis/S1937/78Kr/Simulation/78Se/";
  
  if(range == 0) {
    fn1 += "full.root";
  }
  else {
    
    if(gate == 1)
      fn1 += "PDS/";
    else if(gate == 2)
      fn1 += "REC/";
    
    fn1 += Form("Y%d/full.root",range);
  }
  
  TFile* file1 = new TFile(fn1.c_str(),"read");
  if(file1->IsZombie())
    std::cout << "78Se background ROOT file not found!" << std::endl;

  GH1D* hist1 = (GH1D*)file1->Get(hname1.c_str());
  if(!hist1)
    std::cout << "Histogram not properly retrieved from 78Se ROOT file!" << std::endl;

  int rebin = hist1->GetNbinsX() / nbins;
  if(rebin > 1)
    hist1->Rebin(rebin);
  
  hist1->Sumw2();

  TF1* func = hist1->ConstructTF1();
  func->SetName("78Se");
  func->SetTitle("78Se");
  
  func->SetParameter(0,scale);
  if(fix)
    func->FixParameter(0,scale);
  
  return func;
}

TF1Sum fitAllPeaks(GH1D* data_hist, const std::vector<TF1*>& fit_funcs,
		   std::vector<std::array<double,2>> exclude, double low, double high) {
  
  TF1Sum fullSum;
  for(unsigned int i=0;i<fit_funcs.size();i++) {
    fullSum.AddTF1(fit_funcs.at(i)); 
  }
  fullSum.GetFunc()->SetRange(0,4000);

  for(auto exc : exclude) 
    fullSum.AddExclusion(exc[0],exc[1]);

  int counter = 0;
  TFitResultPtr r;
  while(true) {
    
    r = data_hist->Fit(fullSum.GetFunc(),"MESQL0","",low,high);
    if(r->IsValid() || counter > 4)
      break;
    
    counter++;
  }

  r->Print();
  
  std::cout << "\nReduced Chi2 = " << r->Chi2()/r->Ndf() << "\nr->Status() = " << r->Status()
	    << ", r->IsValid() = " <<  r->IsValid()
	    << "\n****************************************\n";

  std::cout << "\n*************MINOS ERRORS*************\n";
  for(int i=0;i<fullSum.GetNpar();i++) {

    std::cout << r->ParName(i) << ": " << r->Parameter(i) << " + " << r->UpperError(i) << " - "
	      << std::abs(r->LowerError(i)) << "\n" ; 
    
  }
  std::cout << "**************************************\n";
  
  return fullSum;
}

void fitCoulEx(std::string output_fn, std::string peak_input, std::string bg_params,
	       std::string exclusion, int gate, int dop, int range, int fit_low_x, int fit_high_x) {
  
  TFile *data_file = new TFile(data_file_name.c_str(),"read");
  if(data_file->IsZombie()) {
    std::cout << "Data File " << data_file_name << " does not exist!" << std::endl;
    return;
  }
  
  std::string data_hist_name = "S3_Tigress";
  if(range > 0)
    data_hist_name = "Angle_Ranges";
  
  std::string fit_hist_name = "Coincidence/";
  std::string fit_tag;
  
  std::string data_tag;
  std::string data_title;
  
  if(gate == 1) {

    data_tag = Form("_%s",gate_name.c_str());
    data_title = Form("%s DS Gate",gate_name.c_str());

    if(range > 0)
      fit_dir += "PDS/";
    else
      fit_dir += "all/";
    
    fit_hist_name += "ProjectileDS/";
    fit_tag = "DS";
    
  }
  else if(gate == 2) {

    data_tag = Form("_%s",targ_name.c_str());
    data_title = Form("%s Gate",targ_name.c_str());
    
    if(range > 0)
      fit_dir += "REC/";
    else
      fit_dir += "all/";
    
    fit_hist_name += "Recoil/";
    fit_tag = "Rec";
  }
  else if(gate == 3) {

    data_tag = Form("_%sUS",gate_name.c_str());
    data_title = Form("%s US Gate",gate_name.c_str());

    if(range > 0)
      fit_dir += "PUS/";
    else
      fit_dir += "all/";
    
    fit_hist_name += "ProjectileUS/";
    fit_tag = "US";
  }
  else {
    std::cout << "Gate input = " << gate << " It should be 1, 2, or 3" << std::endl;
    return;
  }

  data_hist_name += data_tag;
  if(dop == 0) {
    
    data_title += " Tigress Energy";
    if(range == 0)
      data_hist_name += "/Tigress_Energy";
    else
      data_hist_name += "/tigEn";      
    
    fit_hist_name += "Core_Energy";
  }
  else if(dop == 1) {

    data_title += " Doppler Energy";
    if(range == 0)
      data_hist_name += "/Doppler_Energy";
    else
      data_hist_name += "/dopEn";

    fit_hist_name += "Doppler/Dop_Energy";
  }
  else if(dop == 2) {

    data_title += " Recon Energy";
    if(range == 0)
      data_hist_name += "/Recon_Energy";
    else
      data_hist_name += "/recEn";
    
    fit_hist_name += "Recon/Recon_Energy";
  }
  else {
    std::cout << "Dop input = " << dop << " It should be 0, 1, or 2" << std::endl;
    return;
  }

  data_hist_name += data_tag;
  fit_hist_name += fit_tag;

  if(range > 0) {
    data_hist_name += Form("_Y%d",range);
    data_title += Form(" Range %d",range);
    
    fit_dir += Form("Y%01d/",range);
  }
  
  GH1D* data_hist((GH1D*)data_file->Get(data_hist_name.c_str()));
  if(!data_hist) {
    std::cout << "Histogram " << data_hist_name << " not properly retrieved from data file."
	      << std::endl;
    return;
  }
  data_hist->SetName("hData");
  data_hist->SetTitle(data_title.c_str());

  if(REBIN_FACTOR > 1)
    data_hist->Rebin(REBIN_FACTOR);
  int nbins = data_hist->GetNbinsX();
  
  TNtuple* tpl = new TNtuple("input","input","energy:shift:parameter:ip:ic:fix");
  tpl->ReadFile(peak_input.c_str());

  std::vector<int> energies;
  std::map<int,double> shifts;
  std::map<int,double> params;
  std::map<int,bool> peaks;
  std::map<int,bool> comps;
  std::map<int,bool> fix;

  float en; //energy
  float sh; //shift
  float pa; //parameter
  float ip; //include fep
  float ic; //include compt
  float fx; //fix parameter

  tpl->SetBranchAddress("energy",&en);
  tpl->SetBranchAddress("shift",&sh);
  tpl->SetBranchAddress("parameter",&pa);
  tpl->SetBranchAddress("ip",&ip);
  tpl->SetBranchAddress("ic",&ic);
  tpl->SetBranchAddress("fix",&fx);

  std::cout << "\nPeak Info: " << "\nEnergy\tShift\tParam\t\tPeak\tComp\tFix\n";
  for(int i=0;i<tpl->GetEntries();i++) {

    tpl->GetEntry(i);
    if(!(bool)ip && !(bool)ic) {
      continue;
    }
      
    energies.push_back((int)en);
    shifts[energies.back()] = (int)sh;
    params[energies.back()] = (double)pa; 
    peaks[energies.back()] = (bool)ip;
    comps[energies.back()] = (bool)ic;
    fix[energies.back()] = (bool)fx;
    
    std::cout << energies.back() << "\t" << shifts[energies.back()] << "\t" << params[energies.back()]
	      << "\t" << peaks[energies.back()] << "\t" << comps[energies.back()] << "\t"
	      << fix[energies.back()];
      
    if(en < fit_low_x || en > fit_high_x) {
      std::cout << "   Outside fit range (" << fit_low_x << "," << fit_high_x << ")";
    }
    std::cout << "\n";

    if(!fix[energies.back()] && (!peaks[energies.back()] || !comps[energies.back()])) {
      std::cout << "\t!!!This peak has a free scaling parameter but is missing a component!!!\n\n";
    }
    
  } //end loop over input NTuple (peak input)

  std::string fep_hist_name = fit_hist_name + "_fep";
  std::string com_hist_name = fit_hist_name + "_nfep";
  std::vector<GH1D*> fit_hists;
  std::vector<GH1D*> fep_hists;
  std::vector<GH1D*> com_hists;
  
  for(unsigned int i=0;i<energies.size();i++) {

    std::string fit_file_name = fit_dir + Form("hist%04d.root",energies.at(i));
    
    TFile* f = new TFile(fit_file_name.c_str(),"read");
    if(f->IsZombie()) {
      std::cout << "The ROOT file for the " << energies.at(i) << " keV peak was not found!"
		<< std::endl;
      
      continue;
    }
    
    if(peaks[energies.at(i)] && comps[energies.at(i)]) {
      fit_hists.push_back((GH1D*)f->Get(fit_hist_name.c_str()));
    }
    else if(peaks[energies.at(i)] && !comps[energies.at(i)]) {
      fit_hists.push_back((GH1D*)f->Get(fep_hist_name.c_str()));
    }
    else if(!peaks[energies.at(i)] && comps[energies.at(i)]) {
      fit_hists.push_back((GH1D*)f->Get(com_hist_name.c_str()));
    }

    if(!fit_hists.at(i)) {

      std::cout << "Fit template hist " << fit_hist_name << " for " << energies.at(i)
		<< " keV peak was " << " not found!" << std::endl;
      
      fit_hists.pop_back();
      continue;
      
    }

    fep_hists.push_back((GH1D*)(f->Get(fep_hist_name.c_str())->Clone()));
    com_hists.push_back((GH1D*)(f->Get(com_hist_name.c_str())->Clone()));
      
    fit_hists.back()->Sumw2();
    fep_hists.back()->Sumw2();
    com_hists.back()->Sumw2();
      
    //Do not change this naming convention. It will break the next loop.
    fit_hists.back()->SetName(Form("hist%04i",energies.at(i)));
    fep_hists.back()->SetName(Form("fep_hist%04i",energies.at(i)));
    com_hists.back()->SetName(Form("com_hist%04i",energies.at(i)));

    fit_hists.back()->GetXaxis()->SetLimits(fit_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
					    fit_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]
					    );

    fep_hists.back()->GetXaxis()->SetLimits(fep_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
					    fep_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]
					    );

    com_hists.back()->GetXaxis()->SetLimits(com_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
					    com_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]
					    );

    

    int rebin = fep_hists.back()->GetNbinsX() / nbins;
    if(rebin > 1) {
      fit_hists.back()->Rebin(rebin);
      fep_hists.back()->Rebin(rebin);
      com_hists.back()->Rebin(rebin); 
    }
    
  } //end loop over peaks
  
  std::vector<TF1*> fit_funcs;
  for(unsigned int i=0;i<fit_hists.size();i++) {

    fit_funcs.push_back(fit_hists.at(i)->ConstructTF1());
    fit_funcs.back()->SetParameter(0,params[getEnergy(fit_hists.at(i)->GetName())]);
    
    if(fix[getEnergy(fit_hists.at(i)->GetName())]) {
      fit_funcs.back()->FixParameter(0,params[getEnergy(fit_hists.at(i)->GetName())]);
    }
    
    fit_funcs.back()->SetName(Form("func%04d",getEnergy(fit_hists.at(i)->GetName()))); 
  }
  
  TNtuple* t = new TNtuple("t","t","parameters");
  t->ReadFile(bg_params.c_str());
  t->Draw("parameters","","goff");

  if(t->GetSelectedRows() != 11) {
    std::cout << "Invalid background parameter list" << std::endl;
    return;
  }

  bool inc_tr = bool(t->GetV1()[0]);

  //Fix binning and ranges
  if(false)
    fit_funcs.push_back(TimeRandomBackground(data_hist_name,bool(t->GetV1()[1])));
  
  int inc_bg = int(t->GetV1()[2]);
  if(inc_bg)
    fit_funcs.push_back(ExpBackground(t->GetV1()[3],t->GetV1()[4],t->GetV1()[5],t->GetV1()[6],inc_bg,
				      bool(t->GetV1()[7])));

  bool inc_se = bool(t->GetV1()[8]);
  if(inc_se)
    fit_funcs.push_back(Se78Background(fit_hist_name,gate,range,nbins,t->GetV1()[9],bool(t->GetV1()[10])));

  std::ifstream file;
  file.open(exclusion.c_str(),std::ios::in);

  std::string line;
  std::getline(file,line);

  std::stringstream ss(line);
  int add_exclude;
  ss >> add_exclude;

  std::vector<std::array<double,2>> exclude;
  if(add_exclude) {
    while(std::getline(file,line)) {

      std::stringstream ss1(line);

      double x1, x2;
      ss1 >> x1 >> x2;

      exclude.push_back({x1,x2});
      
    }
  }
  file.close();
  
  TF1Sum fSum(fitAllPeaks(data_hist,fit_funcs,exclude,fit_low_x,fit_high_x));
  fSum.GetFunc()->SetNpx(data_hist->GetNbinsX());
  
  std::vector<double> fep_counts;
  std::vector<double> fep_counts_unc;
  std::vector<int> ens;
  for(unsigned int i=0;i<fep_hists.size();i++){
    
    fep_counts.push_back(fep_hists.at(i)->Integral()*fSum.GetFunc()->GetParameter(i));
    fep_counts_unc.push_back(fep_counts.at(i)*
			     (fSum.GetFunc()->GetParError(i)/fSum.GetFunc()->GetParameter(i)));
    
    ens.push_back(getEnergy(fit_hists.at(i)->GetName()));
    
  }
  printResults(fep_counts,fep_counts_unc,ens);
  
  GH1D* hf = (GH1D*)fSum.GetFunc()->GetHistogram();
  //data_hist->GetListOfFunctions()->Clear();
  
  GH1D* htr = (GH1D*)NULL;
  if(inc_tr) {

    TF1* ftr = fit_funcs.at(fit_hists.size());
    ftr->SetParameter(0,fSum.GetFunc()->GetParameter(fit_hists.size()));
    ftr->SetNpx(4000/REBIN_FACTOR);

    htr = (GH1D*)ftr->GetHistogram();
    htr->SetName("time_random");
    htr->SetTitle("Time Random");
    
  }
  
  GH1D* hbg = (GH1D*)NULL;
  if(inc_bg) {

    int off = 0;
    if(inc_tr)
      off += 1;
    /*
    TF1* bg = new TF1("exp_bg","[0]*TMath::Exp([1]*x)",0,4000);
    bg->SetParameter(0,fSum.GetFunc()->GetParameter(fit_hists.size()+off));
    bg->SetParameter(1,fSum.GetFunc()->GetParameter(fit_hists.size()+off+1));
    bg->SetNpx(4000/REBIN_FACTOR);
  
    hbg = (GH1D*)bg->GetHistogram();
    hbg->SetName("Exp_Bg");
    hbg->SetTitle("Exp_Bg");
    */

    TF1* bg = fit_funcs.at(fit_hists.size()+off);
    bg->SetNpx(4000/REBIN_FACTOR);
    
    for(int i=0;i<bg->GetNpar();i++) {
      //std::cout << "Par " << i << ": " >> ;
      bg->SetParameter(i,fSum.GetFunc()->GetParameter(fit_hists.size()+i+off));
    }
    
    hbg = (GH1D*)bg->GetHistogram();
    hbg->SetName("Exp_Bg");
    hbg->SetTitle("Exp_Bg");
    
  }

  GH1D* hse = (GH1D*)NULL;
  if(inc_se) {

    int off = 0;
    if(inc_tr)
      off += 1;
    
    if(inc_bg)
      off += 2*inc_bg;

    TF1* fse = fit_funcs.at(fit_funcs.size()-1);
    fse->SetParameter(0,fSum.GetFunc()->GetParameter(fit_hists.size()+off));
    fse->SetNpx(4000/REBIN_FACTOR);

    hse = (GH1D*)fse->GetHistogram();
    hse->SetName("Se78");
    hse->SetTitle("Se78");
    
  }
  
  GH1D* hr = new GH1D("Residuum","Residuals",data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),
		      data_hist->GetXaxis()->GetXmax()
		      );

  GH1D* he = new GH1D("Error","Error",data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),
		      data_hist->GetXaxis()->GetXmax()
		      );
  
  GH1D* he1 = new GH1D("Error1","Error1",data_hist->GetNbinsX(),
		       data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		       );
  
  GH1D* he2 = new GH1D("Error2","Error2",data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),
		       data_hist->GetXaxis()->GetXmax()
		       );

  GH1D* he21 = new GH1D("Error21","Error21",data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),
			data_hist->GetXaxis()->GetXmax()
			);
  
  for(int i=1;i<data_hist->GetNbinsX()+1;i++) {
    
    hr->SetBinContent(i,data_hist->GetBinContent(i) - hf->GetBinContent(i));
    hr->SetBinError(i,data_hist->GetBinError(i));

    he->SetBinContent(i,data_hist->GetBinError(i));
    he->SetBinError(i,0);

    he1->SetBinContent(i,-data_hist->GetBinError(i));
    he1->SetBinError(i,0);

    he2->SetBinContent(i,2*data_hist->GetBinError(i));
    he2->SetBinError(i,0);

    he21->SetBinContent(i,-2*data_hist->GetBinError(i));
    he21->SetBinError(i,0);
    
  }
  
  TFile *outfile = new TFile(output_fn.c_str(),"recreate");
  data_hist->Write();
  hf->Write();

  if(htr)
    htr->Write();
  
  if(hbg)
    hbg->Write();

  if(hse)
    hse->Write();
  
  hr->Write();
  he->Write();
  he1->Write();
  he2->Write();
  he21->Write();

  for(unsigned int i=0;i<fit_hists.size();i++) {
    
    fit_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    fit_hists.at(i)->SetName(Form("Peak%02i",i));
    fit_hists.at(i)->SetTitle(Form("hist%04i",energies.at(i)));
    fit_hists.at(i)->Write();

    fep_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    fep_hists.at(i)->SetName(Form("FEP%02i",i));
    fep_hists.at(i)->SetTitle(Form("FEP_hist%04i",energies.at(i)));
    fep_hists.at(i)->Write();

    com_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    com_hists.at(i)->SetName(Form("COM%02i",i));
    com_hists.at(i)->SetTitle(Form("COM_hist%04i",energies.at(i)));
    com_hists.at(i)->Write();  
    
  }

  data_file->Close();
  outfile->Close();
  
  return;
}


int main(int argc, char** argv) { 

  if(argc != 11) {
    
    std::cout << "USAGE: fitCoulEx output_file_name peak_input bg_model_params exclusion_ranges "
	      << "rebin_factor gate dop range x_low x_high" << std::endl;
    
    return 1;
  }

  REBIN_FACTOR = std::stoi(argv[5]);
  
  int gate = std::stoi(argv[6]);
  int dop = std::stoi(argv[7]);
  int range = std::stoi(argv[8]);
  
  int x_low = std::stoi(argv[9]);
  int x_high = std::stoi(argv[10]);
  
  MakeEffCurve(); 
  fitCoulEx(argv[1],argv[2],argv[3],argv[4],gate,dop,range,x_low,x_high);
  
  return 0;
  
}
