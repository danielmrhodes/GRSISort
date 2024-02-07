#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "GH1D.h"
#include "GH2D.h"
#include "TSRIM.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMath.h"

#include "/opt/G4TCX/include/Data_Format.hh"

////These should match the parameters defined in the simulation input////
//Masses in MeV/c^2

const int beamZ = 36;
//const double beam_mass = 74442.6; //80Kr
const double beam_mass = 72582.36; //78Kr
//const double beam_mass = 74441.6; //80Ge
//const double beam_mass = 74449.2; //80Sr
//const double beam_mass = 63274.6; //68Ge
//const double beam_mass = 98626.9; //106Cd
//const double beam_mass = 117279.1; //126Xe

//MeV
//const double beam_en = 278.8; //68Ge
const double beam_en = 331.5; //78Kr Aug 2023
//const double beam_en = 329.6; //78Kr Aug 2023 Pb Target
//const double beam_en = 429.0; //78Kr HE
//const double beam_en = 296.4; //78Kr LE
//const double beam_en = 340.6; //80Sr
//const double beam_en = 318.0; //106Cd

const int targZ = 82;
//const double targ_mass = 182540.0; //196Pt
//const double targ_mass = 183473.2; //197Au
//const double targ_mass = 102376.0; //110Pd
const double targ_mass = 193688.0; //208Pb
//const double targ_mass = 44652.0; //48Ti

//Silicon detector Z-offsets (downstream and upstream) (cm)
//const double DS_Offset = 2.835; //78Kr
const double DS_Offset = 3.3;
const double US_Offset = 3.0;

//Tigress configuration. 0 for detectors forward, 1 for detectors back
const int tigConfig = 0;

//Tigress Z-offset (cm)
const double Tigress_Offset = 0.0;

//Beam spot position (cm)
//const double beam_X = 0.0323; //78Kr
//const double beam_Y = 0.0946; //78Kr
//const double beam_X = -0.03; //78Kr Aug 2023 196Pt
//const double beam_Y = 0.06; //78Kr Aug 2023 196Pt
const double beam_X = -0.12; //78Kr Aug 2023 208Pb
const double beam_Y = 0.05; //78Kr Aug 2023 208Pb
//const double beam_X = 0.01; //78Kr Aug 2023 194Pt
//const double beam_Y = 0.06; //78Kr Aug 2023 194Pt
//const double beam_X = -0.16; //84Kr Aug 2023 196Pt
//const double beam_Y = 0.03; //84Kr Aug 2023 196Pt
//const double beam_X = 0.0;
//const double beam_Y = 0.0;
const double beam_Z = 0.00;

//Total linear target thickness (um)
//const double target_width = 1.0;
//const double target_width = 2.20; //48Ti
//const double target_width = 0.738; //196Pt
//const double target_width = 0.984; //197Au
//const double target_width = 0.805; //110Pd
//const double target_width = 0.882; //208Pb 1.0 mg/cm2
const double target_width = 0.970; //208Pb 1.1 mg/cm2
//const double target_width = 1.36; //196Pt 2.93 mg/cm2

/////////////////////////////////////////////////////////////////////////

//SRIM file names
//const std::string beam_srim_name = "ge68_in_pt196";
//const std::string beam_srim_name = "ge68_in_pb208";
//const std::string beam_srim_name = "kr78_in_au197";
//const std::string beam_srim_name = "kr78_in_pd110";
//const std::string beam_srim_name =  "kr78_in_pt196";
const std::string beam_srim_name =  "kr78_in_pb208";
//const std::string beam_srim_name = "cd106_in_ti48";
//const std::string targ_srim_name = "ti48_in_ti48";
const std::string targ_srim_name = "pb208_in_pb208";
//const std::string targ_srim_name = "196Pt_on_196Pt";
//const std::string targ_srim_name = "pd110_in_pd110";
//const std::string targ_srim_name = "au197_in_au197";

std::map<int,int> secMap = {{0,1},{1,0},{2,31},{3,30},{4,29},{5,28},{6,27},{7,26},{8,25},{9,24},{10,23},{11,22},{12,21},{13,20},{14,19},{15,18},{16,17},{17,16},{18,15},{19,14},{20,13},{21,12},{22,11},{23,10},{24,9},{25,8},{26,7},{27,6},{28,5},{29,4},{30,3},{31,2}};

//Gate on S3 rings (inclusive)
const int Ri = 1;
const int Rf = 24;

////Tigress Intrinsic Energy Resolution////
double Sigma(double en) {
  return 0.7088 + en*0.00034535; //S1893
}
////////////////////////////////////////

////Kinematics////
double Theta_CM_FP(double ThetaLAB, double Ep, bool sol2=false, double Ex=0.0) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  if(std::sin(ThetaLAB) > 1.0/tau) {
    
    ThetaLAB = std::asin(1.0/tau);
    if(ThetaLAB < 0)
      ThetaLAB += TMath::Pi();

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2)
    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;

  return std::asin(tau*std::sin(-ThetaLAB)) + ThetaLAB + TMath::Pi();

}

double Theta_CM_FR(double ThetaLAB, double Ep, bool sol2=false, double Ex=0.0) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  if(std::sin(ThetaLAB) > 1.0/tau) {
    
    ThetaLAB = std::asin(1.0/tau);
    if(ThetaLAB < 0)
      ThetaLAB += TMath::Pi();

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2)
    return TMath::Pi() - (std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB);

  return -std::asin(tau*std::sin(-ThetaLAB)) - ThetaLAB;

}

double Theta_LAB_Max(double Ep, double Ex=0.0) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  if(tau < 1.0)
    return TMath::Pi();
  
  return std::asin(1.0/tau);
  
}

double Theta_LAB(double thetaCM, double Ep, double Ex=0.0) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(thetaCM)/(std::cos(thetaCM) + tau);

  if(tanTheta > 0)
    return std::atan(tanTheta);

  return std::atan(tanTheta) + TMath::Pi();
  
}

double Recoil_Theta_LAB(double thetaCM, double Ep, double Ex=0.0) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(TMath::Pi() - thetaCM)/(std::cos(TMath::Pi() - thetaCM) + tau);
  
  return std::atan(tanTheta);
  
}

double KE_LAB(double thetaCM, double Ep, double Ex=0.0) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = std::pow(targ_mass/(beam_mass + targ_mass),2);
  double term2 = 1 + tau*tau + 2*tau*std::cos(thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

double Recoil_KE_LAB(double thetaCM, double Ep, double Ex=0.0) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = beam_mass*targ_mass/std::pow(beam_mass + targ_mass,2);
  double term2 = 1 + tau*tau + 2*tau*std::cos(TMath::Pi() - thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}
//////////////////

////Positions////
//S3 segment position
TVector3 GetPos(const int det, const int ring, const int sec) {

   if(det < 0 || det > 1 || ring < 1 || ring > 24 || sec < 1 || sec > 32) {
     std::cout << "Bad det, ring, sec (" << det << "," << ring << "," << sec << ")"<< std::endl;
     return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  const double PI = TMath::Pi();
  
  double phi_offset = 0.5*PI; // Phi of sector 1 of downstream detector 
  bool clockwise; // Winding direction of sectors.
  if(det==0)
    clockwise = false;
  else
    clockwise = true;

  double janus_outer_radius = 3.5;
  double janus_inner_radius = 1.1;

  double rad_slope = (janus_outer_radius - janus_inner_radius)/24.;
  double rad_offset = janus_inner_radius;

  TVector3 pos(1.,0,0);
  pos.SetPerp((ring - 0.5)*rad_slope + rad_offset);
  
  double phi = phi_offset + (clockwise ? -1 : 1) * 2.*PI/32. * (sec - 1);
  pos.SetPhi(phi);

  double zoff;
  if(det == 1)
    zoff = DS_Offset;
  else
    zoff = US_Offset;
  
  pos.SetZ(zoff);
  if(det == 0)
    pos.RotateY(PI);

  return pos;
  
}

//Tigress positions
//forward (0) or back (1), det, seg
TVector3 tigPosVec[2][64][9];
void FillPositonVectors() {

  std::ifstream file;
  file.open("/opt/G4TCX/positions.txt",std::ios::in);

  std::string line;
  double x,y,z;
  
  for(int i=0;i<64;i++) {
    for(int j=0;j<9;j++) {

      std::getline(file,line);
      std::stringstream ss(line);
      ss >> x >> y >> z;
	
      TVector3 vec(x/10.,y/10.,z/10.);
      tigPosVec[0][i][j] = vec;

    }
  }

  std::ifstream file1;
  file1.open("/opt/G4TCX/positionsBack.txt",std::ios::in);
  
  for(int i=0;i<64;i++) {
    for(int j=0;j<9;j++) {

      std::getline(file1,line);
      std::stringstream ss(line);
      ss >> x >> y >> z;
	
      TVector3 vec(x/10.,y/10.,z/10.);
      tigPosVec[1][i][j] = vec;

    }
  }

  return;
  
}

//Tigress segment position
TVector3 GetPos(const int det, const int seg) {

  if(det < 1 || det > 64 || seg < 0 || seg > 8)
    return TVector3(0,0,0);
  
  return tigPosVec[tigConfig][det-1][seg];

}
/////////////////

////////////Build data format////////////
struct S3 {

  //realistic info
  int det, ring, sector;
  double rEn, sEn;

  //perfect info
  TVector3 rPos, sPos;
  bool rP, rR, sP, sR;
  
};

struct TIGRESS {

  //realistic info
  int det, nsegs, segs[8];
  double cEn, sEn[8];
  bool sup;
  
  //perfect info
  double x[8], y[8], z[8];
  bool fep, pfep;

  int MainSeg() {

    if(!nsegs)
      return 0;
    
    int index=0;
    for(int i=1;i<nsegs;i++) {
      if(sEn[i] > sEn[index]) {
	index=i;
      }
    }
    
    return segs[index];
  }

  int LastSeg() {

     if(!nsegs)
      return 0;
    
    int index=0;
    for(int i=1;i<nsegs;i++) {
      if(sEn[i] < sEn[index]) {
	index=i;
      }
    }
    
    return segs[index];
  }

  TVector3 MainPos() {
    return GetPos(det,MainSeg());
  }

  TVector3 LastPos() {
    return GetPos(det,LastSeg());
  }

  TVector3 ExactMainPos() {

    if(!nsegs)
      return GetPos(det,0);
    
    int index=0;
    for(int i=1;i<nsegs;i++) {
      if(sEn[i] > sEn[index]) {
	index=i;
      }
    }

    return TVector3(x[index],y[index],z[index]);
  }

  void AddHit(const TIGRESS& hit) {

    cEn += hit.cEn;
    if(!sup)
      sup = hit.sup;
    
    return;
  }

  bool operator>(const TIGRESS& hit) const {
    return cEn > hit.cEn;
  }
  
};

bool AddbackCriterion(TIGRESS& hit1,TIGRESS& hit2) {

  double res = (hit1.LastPos() - hit2.MainPos()).Mag();

  // In clover core separation 54.2564, 76.7367
  // Between clovers core separation 74.2400 91.9550 (high-eff mode)
  double seperation_limit = 93.0/10.0;

  int one_seg = hit1.LastSeg();
  int two_seg = hit2.MainSeg();

  // front segment to front segment OR back segment to back segment
  if((one_seg < 5 && two_seg < 5) || (one_seg > 4 && two_seg > 4))
    seperation_limit = 54.0/10.0;

  // front to back
  else if((one_seg < 5 && two_seg > 4) || (one_seg > 4 && two_seg < 5))
    seperation_limit = 105.0/10.0;
  
  if(res < seperation_limit)
    return true;

  return false;
}

bool BadSeg(TIGRESS tig) {
  
  //int num = hit->GetArrayNumber();
  int num = tig.det;
  if(num != 5 && num != 10)
    return false;

  //int seg_mult = hit->GetNSegments();
  int seg_mult = tig.nsegs;
  for(int j=0;j<seg_mult;j++) {
    
    //TDetectorHit seg_hit = hit->GetSegmentHit(j);
    //int seg = seg_hit.GetSegment();
    int seg = tig.segs[j];
   
    if((seg > 4 && num == 5) || (seg == 1 && num == 10))
      return true;
    
  }

  return false;

}

struct BuiltData {

  int evt;
  int nS3;
  int nTi;
  
  S3 s3[5];
  std::vector<TIGRESS> tigress;
  std::vector<TIGRESS> tigressAB;

  void MakeS3Hit(const S3Data ringHit, const S3Data secHit) {

    s3[nS3].det = ringHit.det; //either hit works
    s3[nS3].ring = ringHit.ring;
    s3[nS3].sector = secHit.sector;
    
    s3[nS3].rEn = ringHit.en;
    s3[nS3].sEn = secHit.en;

    s3[nS3].rPos = TVector3(ringHit.x,ringHit.y,ringHit.z);
    s3[nS3].sPos = TVector3(secHit.x,secHit.y,secHit.z);
    
    s3[nS3].rP = ringHit.proj;
    s3[nS3].rR = ringHit.rec;

    s3[nS3].sP = secHit.proj;
    s3[nS3].sR = secHit.rec;

    nS3++;
  }

  void BuildAddbackHits() {
    
    if(nTi == 0)
      return;

    tigressAB.push_back(tigress[0]);
    if(nTi == 1)
      return;
      
    int i;
    unsigned int j;
    for(i=1;i<nTi;i++) {
      
      for(j=0;j<tigressAB.size();j++) {
	if(AddbackCriterion(tigressAB[j],tigress[i])) {

	  tigressAB[j].AddHit(tigress[i]);
	  break;
	}
      }
      
      // if hit[i] was not added to a higher energy hit, create its own addback hit
      if(j == tigressAB.size())
	tigressAB.push_back(tigress[i]);
      
    }

    return;
  }
  
};
/////////////////////////////////////////

//Unpack raw data into correlated data
BuiltData BuildData(const int nE, const int nSch, const int nTch, const RawData& raw_dat) {
  
  BuiltData data;
  data.evt = nE;
  data.nS3 = 0;
  data.nTi = 0;

  //Correlate S3 Data
  std::vector<S3Data> rings;
  std::vector<S3Data> sectors;
  for(int i=0;i<nSch;i++) {

    S3Data sChan(raw_dat.sData[i]);
    
    if(sChan.IsRing())
      rings.push_back(sChan);
    else
      sectors.push_back(sChan);
    
  }

  std::sort(rings.begin(),rings.end(),std::greater<S3Data>());
  std::sort(sectors.begin(),sectors.end(),std::greater<S3Data>());

  std::vector<bool> used_rings;
  std::vector<bool> used_sectors;
  
  used_rings.resize(rings.size());
  std::fill(used_rings.begin(),used_rings.end(),false);

  used_sectors.resize(sectors.size());
  std::fill(used_sectors.begin(),used_sectors.end(),false);

  for(unsigned int i=0;i<sectors.size();i++) {
    if(used_sectors.at(i))
      continue;
    
    for(unsigned int j=0;j<rings.size();j++) {
      if(used_rings.at(j))
        continue;

      //Same detector
      if((sectors.at(i).det == rings.at(j).det) &&

	 //Same energy
	 (TMath::Abs(sectors.at(i).en - rings.at(j).en) < 1.)) {

	data.MakeS3Hit(rings.at(j),sectors.at(i));

	used_sectors.at(i) = true;
	used_rings.at(j) = true;
	break;

      }
    }
  }

  bool broken = false;
  for(unsigned int i=0;i<sectors.size();i++) {
    broken=false;
    if(used_sectors.at(i))
      continue;
    
    for(unsigned int j=0;j<rings.size();j++) {
      if(used_rings.at(j))
        continue;

      for(unsigned int k=0;k<sectors.size();k++) {
	if(used_sectors.at(k))
	  continue;

	//Same detector
	if((sectors.at(i).det == rings.at(j).det) && (sectors.at(k).det == rings.at(j).det) &&

	   //Sector energies add to ring energy
	   (TMath::Abs(sectors.at(i).en + sectors.at(k).en - rings.at(j).en) < 1.)) {

	  data.MakeS3Hit(rings.at(j),sectors.at(i));
	  data.MakeS3Hit(rings.at(j),sectors.at(k));
	  
	  used_sectors.at(i) = true;
	  used_rings.at(j) = true;
	  used_sectors.at(k) = true;
	  broken = true;
	  break;

	}
      }
      if(broken)
	break;
    }
  }

  //Organize Tigress data
  std::vector<bool> exists;
  exists.resize(64);
  std::fill(exists.begin(),exists.end(),false);
  
  int nT = 0;
  for(int i=0;i<nTch;i++) {

    int detect = raw_dat.tData[i].det;
    int segment = raw_dat.tData[i].seg;
    
    double energy = raw_dat.tData[i].en;
    double x = raw_dat.tData[i].x;
    double y = raw_dat.tData[i].y;
    double z = raw_dat.tData[i].z;
    
    bool FEP = raw_dat.tData[i].fep;
    bool PFEP = raw_dat.tData[i].pfep;
    bool SUP = raw_dat.tData[i].sup;
    
    if(!exists.at(detect-1)) {
      data.tigress.emplace_back();
      data.tigress[nT].det = detect;
      data.tigress[nT].sup = SUP;
      
      if((bool)segment) {
        data.tigress[nT].nsegs = 1;
	data.tigress[nT].segs[0] = segment;
	data.tigress[nT].sEn[0] = energy;

	data.tigress[nT].x[0] = x;
        data.tigress[nT].y[0] = y;
        data.tigress[nT].z[0] = z;
      }
      else {
	data.tigress[nT].nsegs = 0;
	data.tigress[nT].cEn = energy;
	data.tigress[nT].fep = FEP;
	data.tigress[nT].pfep = PFEP;
      }

      nT++;
      exists.at(detect-1) = true;
    }
    else {

      int index = 0;
      for(int j=0;j<nT;j++) {
	if(data.tigress[j].det == detect) {
	  index = j;
	  break;
	}
      }
      data.tigress[index].sup = SUP;
      
      int Nsegs = data.tigress[index].nsegs;
      
      if((bool)segment) {
	data.tigress[index].segs[Nsegs] = segment;
	data.tigress[index].sEn[Nsegs] = energy;

	data.tigress[index].x[Nsegs] = x;
        data.tigress[index].y[Nsegs] = y;
        data.tigress[index].z[Nsegs] = z;
	
	data.tigress[index].nsegs++;
      }
      else {
	data.tigress[index].cEn = energy;
	data.tigress[index].fep = FEP;
	data.tigress[index].pfep = PFEP;
      }
      
    }
    
  }
  data.nTi = nT;
  
  std::sort(data.tigress.begin(),data.tigress.end(),std::greater<TIGRESS>());
  data.BuildAddbackHits();
  
  return data;
}
