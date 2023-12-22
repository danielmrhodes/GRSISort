#include <iostream>
#include <cstdio>
#include <vector>

#include "TFile.h"
//#include "GH1.h"
//#include "GH2.h"
#include "GH1D.h"
#include "GH2D.h"
#include "TSRIM.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMath.h"

#include "/opt/JANUS/include/Data_Format.hh"

////These should match the parameters defined in the simulation input////
//Masses in MeV/c^2

const int beamZ = 48;
//const double beam_mass = 74442.6; //80Kr
//const double beam_mass = 72582.36; //78Kr
//const double beam_mass = 74441.6; //80Ge
const double beam_mass = 98626.9; //106Cd
//const double beam_mass = 117279.1; //126Xe

//MeV
//const double beam_en = 281.6; //80Ge, 80Kr
//const double beam_en = 429.0; //78Kr
//const double beam_en = 462.2; //106Cd HE
//const double beam_en = 471.24; //126Xe
const double beam_en = 318.0; //106Cd-48Ti

const int targZ = 22;
//const double targ_mass = 182540.0; //196Pt
//const double targ_mass = 183473.2; //197Au
//const double targ_mass = 193688.0; //208Pb
const double targ_mass = 44652.0; //48Ti

//Silicon detector Z-offsets (downstream and upstream) (cm)
const double DS_Offset = 2.6;
const double US_Offset = 3.4;

//SeGA Z-offset (cm)
const double SeGA_Offset = 3.1;
//const double SeGA_Offset = 6.5;
//const double SeGA_Offset = 0.0;

//Beam spot position (cm)
//const double beam_X = -0.051; //80Ge
//const double beam_Y = 0.104; //80Ge
//const double beam_X = -0.033; //80Kr
//const double beam_Y = 0.036; //80Kr
//const double beam_X = 0.0323; //78Kr
//const double beam_Y = 0.0946; //78Kr
//const double beam_X = 0.01;
//const double beam_Y = -0.044;
const double beam_X = 0.00;
const double beam_Y = 0.00;
const double beam_Z = 0.00;

//Total linear target thickness (um)
//const double target_width = 1.0;
const double target_width = 2.20; //48Ti
//const double target_width = 0.738; //196Pt
//const double target_width = 0.9834; //197Au
//const double target_width = 0.882; //208Pb

/////////////////////////////////////////////////////////////////////////

//SRIM file names
//const std::string beam_srim_name = "80Ge_on_196Pt";
//const std::string beam_srim_name = "80Kr_on_196Pt";
//const std::string targ_srim_name = "196Pt_on_196Pt";
//const std::string beam_srim_name = "cd106_in_pb208";
const std::string beam_srim_name = "cd106_in_ti48";
const std::string targ_srim_name = "ti48_in_ti48";
//const std::string targ_srim_name = "pb208_in_pb208";
//const std::string beam_srim_name = "Xe126_in_Pt196";

const int Ri = 1; //inclusive
const int Rf = 24; //inclusive

////SeGA Intrinsic Energy Resolution////
double Sigma(double en) {
  //return 1.03753 + en*0.000274797; //2018
  return 1.0*(1.2127 + en*0.000275177); //2020
  //return 0.0;
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
//Bambino2 segment position
TVector3 GetPos(const int det, const int ring, const int sec) {

   if(det < 0 || det > 1 || ring < 1 || ring > 24 || sec < 1 || sec > 32) {
     std::cout << "Bad det, ring, sec (" << det << "," << ring << "," << sec << ")"<< std::endl;
     return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  const double PI = TMath::Pi();
  
  double phi_offset = 0.5*PI; // Phi of sector 1 of downstream detector 
  bool clockwise; // Winding direction of sectors.
  if(det==0) {
    clockwise = false;
  }
  else {
    clockwise = true;
  }

  double janus_outer_radius = 3.5;
  double janus_inner_radius = 1.1;

  double rad_slope = (janus_outer_radius - janus_inner_radius)/24.;
  double rad_offset = janus_inner_radius;

  TVector3 pos(1.,0,0);
  pos.SetPerp((ring - 0.5)*rad_slope + rad_offset);
  
  double phi = phi_offset + (clockwise ? -1 : 1) * 2.*PI/32. * (sec - 1);
  pos.SetPhi(phi);

  double zoff;
  if(det == 1) {
    zoff = DS_Offset;
  }
  else {
    zoff = US_Offset;
  }
  pos.SetZ(zoff);

  if(det == 0) {
    pos.RotateY(PI);
  }

  return pos;
  
}

//SeGA segment position
TVector3 GetPos(const int det, const int seg) {

  if(det < 0 || det > 16 || seg < 1 || seg > 32) {
    std::cout << "Bad det, seg (" << det << "," << seg << ")" << std::endl;
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  const double PI = TMath::Pi();

  int quad = int(seg-1)/(int)8;

  int slice;
  if(seg < 9) {
    slice = seg-1;
  }
  else if(seg < 17) {
    slice = seg-9;
  }
  else if(seg < 25) {
    slice = seg-17;
  }
  else {
    slice = seg-25;
  }

  double length = 4.025;
  double outerRadius = 3.165;
  double innerRadius;
  
  if(slice == 7) {
    innerRadius = 0.0;
  }
  else {
    innerRadius = 0.5 + 0.03; //fingerRadius + DL thickness
  }

  TVector3 pos(1.0,1.0,1.0);
  pos.SetPerp((outerRadius + innerRadius)/2.0);
  pos.SetPhi((quad+0.5)*2.0*PI/4.0);
  pos.SetZ((length/8.0)*(2.0*slice - 7.0));

  double rd = 12.975;
  double phid = (det-1)*(2.0*PI/8.0) + PI/8.0;
  //double phid = (det-1)*(2.0*PI/8.0);
  double zd = length + 2*0.05 + 0.6;
  if(det > 8) {
    zd*=-1;
  }
  
  TVector3 origin(rd*TMath::Cos(phid),rd*TMath::Sin(phid),zd+SeGA_Offset);

  return origin+pos;

}
/////////////////

////////////Build data format////////////
struct BAM2 {

  //realistic info
  int det, ring, sector;
  double rEn, sEn;

  //perfect info
  TVector3 rPos, sPos;
  bool rP, rR, sP, sR;
  
};

struct SEGA {

  //realistic info
  int det, nsegs, segs[32];
  double cEn, sEn[32];

  //perfect info
  double x[32], y[32], z[32];
  bool fep, pfep;

  int MainSeg() {
    
    int index=0;
    for(int i=0;i<nsegs;i++) {
      if(sEn[i] > sEn[index]) {
	index=i;
      }
    }
    
    return segs[index];
  }
  
};

struct BuiltData {

  int evt;
  int nBa;
  int nSe;
  
  BAM2 bam2[50];
  SEGA sega[16];

  void MakeHit(const Bambino2Data ringHit, const Bambino2Data secHit) {

    bam2[nBa].det = ringHit.det; //either one works
    bam2[nBa].ring = ringHit.ring;
    bam2[nBa].sector = secHit.sector;
    
    bam2[nBa].rEn = ringHit.en;
    bam2[nBa].sEn = secHit.en;

    bam2[nBa].rPos = TVector3(ringHit.x,ringHit.y,ringHit.z);
    bam2[nBa].sPos = TVector3(secHit.x,secHit.y,secHit.z);
    
    bam2[nBa].rP = ringHit.proj;
    bam2[nBa].rR = ringHit.rec;

    bam2[nBa].sP = secHit.proj;
    bam2[nBa].sR = secHit.rec;

    nBa++;
  }
  
};
/////////////////////////////////////////

//Unpack raw data into correlated data
BuiltData BuildData(const Header& head, const JANUSData& dat) {

  const int nBch = head.nBdata;
  const int nSch = head.nSdata;
  const int nE = head.evtNum;
  
  BuiltData data;
  data.evt = nE;
  data.nBa = 0;
  data.nSe = 0;

  //Correlate Bambino2 Data
  std::vector<Bambino2Data> rings;
  std::vector<Bambino2Data> sectors;
  for(int i=0;i<nBch;i++) {

    Bambino2Data bChan(dat.bData[i]);
    
    if(bChan.IsRing()) {
      rings.push_back(bChan);
    }
    else {
      sectors.push_back(bChan);
    }
  }

  std::sort(rings.begin(),rings.end());
  std::sort(sectors.begin(),sectors.end());

  std::vector<bool> used_rings;
  std::vector<bool> used_sectors;
  
  used_rings.resize(rings.size());
  std::fill(used_rings.begin(),used_rings.end(),false);

  used_sectors.resize(sectors.size());
  std::fill(used_sectors.begin(),used_sectors.end(),false);

  for(unsigned int i=0;i<sectors.size();i++) {
    if(used_sectors.at(i)) {
      continue;
    }
    
    for(unsigned int j=0;j<rings.size();j++) {
      if(used_rings.at(j)) {
        continue;
      }

         //Same detector
      if((sectors.at(i).det == rings.at(j).det) &&

	 //Same energy
	 (TMath::Abs(sectors.at(i).en - rings.at(j).en) < 1.)) {

	data.MakeHit(rings.at(j),sectors.at(i));

	used_sectors.at(i) = true;
	used_rings.at(j) = true;
	break;

      }
    }
  }

  bool broken = false;
  for(unsigned int i=0;i<sectors.size();i++) {
    broken=false;
    if(used_sectors.at(i)) {
      continue;
    }
    
    for(unsigned int j=0;j<rings.size();j++) {
      if(used_rings.at(j)) {
        continue;
      }

      for(unsigned int k=0;k<sectors.size();k++) {
	if(used_sectors.at(k)) {
	  continue;
	}

	  //Same detector
	if((sectors.at(i).det == rings.at(j).det) && (sectors.at(k).det == rings.at(j).det) &&

	   //Sector energies add to ring energy
	   (TMath::Abs(sectors.at(i).en + sectors.at(k).en - rings.at(j).en) < 1.)) {

	  data.MakeHit(rings.at(j),sectors.at(i));
	  data.MakeHit(rings.at(j),sectors.at(k));
	  
	  used_sectors.at(i) = true;
	  used_rings.at(j) = true;
	  used_sectors.at(k) = true;
	  broken = true;
	  break;

	}
      }
      if(broken) {
	break;
      }
    }
  }

  //Organize SeGA data
  std::vector<bool> exists;
  exists.resize(16);
  std::fill(exists.begin(),exists.end(),false);

  int nS = 0;
  for(int i=0;i<nSch;i++) {

    int detect = dat.sData[i].det;
    int segment = dat.sData[i].seg;
    
    double energy = dat.sData[i].en;
    double x = dat.sData[i].x;
    double y = dat.sData[i].y;
    double z = dat.sData[i].z;
    
    bool FEP = dat.sData[i].fep;
    bool PFEP = dat.sData[i].pfep;
    
    if(!exists.at(detect-1)) {
      data.sega[nS].det = detect;

      if((bool)segment) {
        data.sega[nS].nsegs = 1;
	data.sega[nS].segs[0] = segment;
	data.sega[nS].sEn[0] = energy;

	data.sega[nS].x[0] = x;
	data.sega[nS].y[0] = y;
	data.sega[nS].z[0] = z;
      }
      else {
	data.sega[nS].nsegs = 0;
	data.sega[nS].cEn = energy;
	data.sega[nS].fep = FEP;
	data.sega[nS].pfep = PFEP;
      }

      nS++;
      exists.at(detect-1) = true;
    }
    else {

      int index = 0;
      for(int j=0;j<nS;j++) {
	if(data.sega[j].det == detect) {
	  index = j;
	  break;
	}
      }
      
      int Nsegs = data.sega[index].nsegs;
      
      if((bool)segment) {
	data.sega[index].segs[Nsegs] = segment;
	data.sega[index].sEn[Nsegs] = energy;

	data.sega[index].x[Nsegs] = x;
	data.sega[index].y[Nsegs] = y;
	data.sega[index].z[Nsegs] = z;
	
	data.sega[index].nsegs++;
      }
      else {
	data.sega[index].cEn = energy;
	data.sega[index].fep = FEP;
	data.sega[index].pfep = PFEP;
      }
      
    }
    
  }
  data.nSe = nS;

  return data;
}

int main(int argc, char** argv) {
  
  if(argc != 3) {
    std::cerr << "Usage: correlator INPUT_FILE OUTPUT_FILE" << std::endl;
    return 1;
  }

  const char* input_filename = argv[1];
  const char* output_filename = argv[2];
  if(!strcmp(input_filename,output_filename)) {
    std::cout << "Give your input and output files different names" << std::endl;
    return 1;
  }
  
  //Bambino2 singles
  GH2D* bSum = new GH2D("Summary","Janus Summary",120,1,121,1000,0,1000);
  GH2D* pSum = new GH2D("pSummary","Janus Projectile Summary",120,1,121,1000,0,1000);
  GH2D* rSum = new GH2D("rSummary","Janus Recoil Summary",120,1,121,1000,0,1000);

  /*
  double shift = 0.5*TMath::TwoPi()/32.0;
  GH2D* pPvPDS = new GH2D("pPerpvPhiDS","Projectile Perp_v_Phi",32,-TMath::Pi()-shift,TMath::Pi()-shift,
			 24,1.1,3.5);
  GH2D* rPvP = new GH2D("rPerpvPhi","Recoil Perp_v_Phi",32,-TMath::Pi()-shift,TMath::Pi()-shift,
		       24,1.1,3.5);
  */

  GH2D* pSvRDS = new GH2D("pSvRDS","Projectile Sector vs Ring",24,1,25,32,1,33);
  GH2D* pThvPhDS = new GH2D("pThvPhDS","Projectile #phi-#theta surface",1000,0,90,1000,-200,200);
  GH2D* pThvPhDS1 = new GH2D("pThvPhDS1","Projectile #phi-#theta surface",1000,0,90,1000,-200,200);
  GH2D* pThvPhDS2 = new GH2D("pThvPhDS2","Projectile #phi-#theta surface",1000,0,90,1000,-200,200);

  GH2D* rThvPh = new GH2D("rThvPh","Recoil #phi-#theta surface",1000,0,90,1000,-200,200);
  GH2D* pThvPhUS = new GH2D("pThvPhUS","Projectile US #phi-#theta surface",1000,90,180,1000,-200,200);
  
  GH1D* secD0 = new GH1D("Sectors_Det0","US Sectors",32,1,33);
  GH1D* secD1 = new GH1D("Sectors_Det1","DS Sectors",32,1,33);

  GH1D* pSecDS = new GH1D("pSecDS","DS Projectile Sectors",32,1,33);
  GH1D* rSec = new GH1D("rSec","Recoil Sectors",32,1,33);

  /*
  GH1D* pSecDS_m1 = new GH1D("pSecDS_m1","DS Projectile Sectors Mult1",32,1,33);
  GH1D* rSec_m1 = new GH1D("rSec_m1","Recoil Sectors Mult1",32,1,33);
  
  GH1D* pSecDS_m2 = new GH1D("pSecDS_m2","DS Projectile Sectors Mult2",32,1,33);
  GH1D* rSec_m2 = new GH1D("rSec_m2","Recoil Sectors Mult2",32,1,33);

  GH1D* pPhiDS = new GH1D("pPhiDS","DS Projectile Phi",34,-191.25,191.25);
  GH1D* rPhi = new GH1D("rPhi","Recoil Phi",34,-191.25,191.25);
  */
  
  GH2D* rPid0 = new GH2D("RingPID_Det0","US RingEn PID",24,1,25,1000,0,500);
  GH2D* rPid1 = new GH2D("RingPID_Det1","DS RingEn PID",24,1,25,1500,0,1000);
  GH2D* rPid1m1 = new GH2D("RingPID_Det1m1","DS RingEn PID Mult1",24,1,25,1000,0,1000);
  GH2D* rPid1m2 = new GH2D("RingPID_Det1m2","DS RingEn PID Mult2",24,1,25,1000,0,1000);

  GH2D* sPid0 = new GH2D("SecPID_Det0","US SectorEn PID",24,1,25,1000,0,500);
  GH2D* sPid1 = new GH2D("SecPID_Det1","DS SectorEn PID",24,1,25,1500,0,1000);
  GH2D* sPid1m1 = new GH2D("SecPID_Det1m1","DS SectorEn PID Mult1",24,1,25,1000,0,1000);
  GH2D* sPid1m2 = new GH2D("SecPID_Det1m2","DS SectorEn PID Mult2",24,1,25,1000,0,1000);

  GH2D* sPid1_p = new GH2D("SecPID_Det1_proj","DS SectorEn Projectile PID",24,1,25,1500,0,1000);
  GH2D* sPid1_r = new GH2D("SecPID_Det1_rec","DS SectorEn Recoil PID",24,1,25,1500,0,1000);
  
  //SeGA singles
  GH1D* coreEnergy = new GH1D("Core_Energy","SeGA Core Energy",3000,0,3000);
  GH1D* segEnergy = new GH1D("Seg_Energy","SeGA Segment Energy",3000,0,3000);

  GH2D* coreSum = new GH2D("Core_Summary","Core Energy Summary",16,1,17,3000,0,3000);
  GH2D* segSum = new GH2D("Seg_Summary","Segment Energy Summary",512,1,513,3000,0,3000);

  GH1D* coreEn_Fep = new GH1D("FEP","SeGA FEP",3000,0,3000);
  GH1D* coreEn_NotFep = new GH1D("nFEP","SeGA Not FEP",3000,0,3000);

  GH2D* FEPsum = new GH2D("FEP_Summary","Core Energy Summary FEP",16,1,17,3000,0,3000);

  //GH1D* gamGate = new GH1D("gamGate","CoreEn with 1332 keV Coinc",3000,0,3000);
  //GH1D* gammaAng = new GH1D("gamAng","Gamma Gamma Angle",100,-1.1,1.1);
  //GH1D* gammaAngFEP = new GH1D("gamAngFEP","Gamma Gamma Angle",100,-1.1,1.1);
  
  GH2D* sPosTP = new GH2D("sPosPT","SeGA Phi-Theta Surface",360,0,180,760,-10,370);
  
  //Coincidences
  GH2D* sPidDS_gam = new GH2D("SecPID_DS_gam","SectorEn PID Gamma Gate",24,1,25,1500,0,1000);
  GH2D* psPidDS_gam = new GH2D("pSecPID_DS_gam","Projectile SectorEn PID Gamma Gate",24,1,25,1500,0,1000);
  GH2D* rsPidDS_gam = new GH2D("rSecPID_DS_gam","Recoil SectorEn PID Gamma Gate",24,1,25,1500,0,1000);
   
  //Projectile DS
  GH2D* sPidDS = new GH2D("SecPID_DS","DS SectorEn PID",24,1,25,1500,0,1000);
  GH2D* rPidDS = new GH2D("RingPID_DS","DS RingEn PID",24,1,25,1500,0,1000);
  
  GH1D* pCoreEnergyDS = new GH1D("Core_EnergyDS","SeGA Core Energy",12000,0,4000);
  GH2D* pCoreSumDS = new GH2D("Core_SummaryDS","Core Energy Summary",16,1,17,3000,0,3000);

  GH1D* pDopEnergyDS = new GH1D("Dop_EnergyDS","Doppler Energy",12000,0,4000);
  GH2D* pDopSumDS = new GH2D("Dop_SummaryDS","Doppler Energy Summary",16,1,17,6000,0,3000);
  GH2D* pDopASSumDS = new GH2D("Dop_sSliceSum_DS","Doppler Energy Array Slice Summary",16,1,17,6000,0,3000);

  GH1D* pCoreEnergyDS_fep = new GH1D("Core_EnergyDS_fep","SeGA Core Energy FEP",12000,0,4000);
  GH1D* pCoreEnergyDS_nfep = new GH1D("Core_EnergyDS_nfep","SeGA Core Energy Not FEP",12000,0,4000);
  GH1D* pCoreEnergyDS_pfep = new GH1D("Core_EnergyDS_pfep","SeGA Core Energy Projectile FEP",12000,0,4000);
  GH1D* pCoreEnergyDS_rfep = new GH1D("Core_EnergyDS_rfep","SeGA Core Energy Recoil FEP",12000,0,4000);

  GH1D* pDopEnergyDS_fep = new GH1D("Dop_EnergyDS_fep","Doppler Energy FEP",12000,0,4000);
  GH1D* pDopEnergyDS_pfep = new GH1D("Dop_EnergyDS_pfep","Doppler Energy Projectile FEP",12000,0,4000);
  GH1D* pDopEnergyDS_rfep = new GH1D("Dop_EnergyDS_rfep","Doppler Energy Recoil FEP",12000,0,4000);
  GH1D* pDopEnergyDS_nfep = new GH1D("Dop_EnergyDS_nfep","Doppler Energy Not FEP",12000,0,4000);

  GH2D* pDopvPartDS = new GH2D("DopEn_v_PartEn_DS","Doppler Energy vs Particle Energy",
			      3000,0,3000,1000,0,1000);

  GH2D* pCEnvRingDS = new GH2D("CoreEn_v_Ring_DS","Core Energy vs Ring",
			      24,1,25,4000,0,4000);
  GH2D* pDopvRingDS = new GH2D("DopEn_v_Ring_DS","Doppler Energy vs Ring",
			      24,1,25,4000,0,4000);
  GH2D* pRecvRingDS = new GH2D("ReconEn_v_Ring_DS","Recon Energy vs Ring",
			      24,1,25,4000,0,4000);
  
  GH2D* pThCorDS = new GH2D("Theta_CorrDS","Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pThCrtDS = new GH2D("Theta_CrctDS","Theta Correction",6000,0,3000,90,0,180);

  GH2D* pcThCrtDS = new GH2D("cosTheta_CrctDS","CosTheta Correction",6000,0,3000,200,-1.1,1.1);
  GH2D* pcThCrtDS_fep = new GH2D("cosTheta_CrctDS_fep","FEP CosTheta Correction",6000,0,3000,200,-1.1,1.1);

  double thing1 = 65.*180./32.0;
  GH2D* pPhCorDS = new GH2D("Phi_CorrDS","Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pPhCrtDS = new GH2D("Phi_CrctDS","Phi Correction",6000,0,3000,32,0,thing1);

  GH1D* pReconEnergyDS = new GH1D("Recon_EnergyDS","Recon Energy",12000,0,4000);
  GH2D* pReconSumDS = new GH2D("Recon_SummaryDS","Recon Energy Summary",16,1,17,6000,0,3000);

  GH1D* pReconEnergyDS_fep = new GH1D("Recon_EnergyDS_fep","Recon Energy FEP",12000,0,4000);
  GH1D* pReconEnergyDS_pfep = new GH1D("Recon_EnergyDS_pfep","Recon Energy Projectile FEP",12000,0,4000);
  GH1D* pReconEnergyDS_rfep = new GH1D("Recon_EnergyDS_rfep","Recon Energy Recoil FEP",12000,0,4000);
  GH1D* pReconEnergyDS_nfep = new GH1D("Recon_EnergyDS_nfep","Recon Energy Not FEP",12000,0,4000);

  GH2D* pReconvPartDS = new GH2D("ReconEn_v_partEn_DS","Recon Energy vs Particle Energy",
				3000,0,3000,1000,0,1000);

  GH2D* pReconThCorDS = new GH2D("ReconTheta_CorrDS","Recon Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pReconThCrtDS = new GH2D("ReconTheta_CrctDS","Recon Theta Correction",6000,0,3000,90,0,180);

  GH2D* pReconPhCorDS = new GH2D("ReconPhi_CorrDS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pReconPhCrtDS = new GH2D("ReconPhi_CrctDS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Projectile US
  GH2D* sPidUS = new GH2D("SecPID_US","DS SectorEn PID",24,1,25,1000,0,500);
  GH2D* rPidUS = new GH2D("RingPID_US","DS RingEn PID",24,1,25,1000,0,500);
  
  GH1D* pCoreEnergyUS = new GH1D("Core_EnergyUS","SeGA Core Energy",12000,0,4000);
  GH2D* pCoreSumUS = new GH2D("Core_SummaryUS","Core Energy Summary",16,1,17,3000,0,3000);
  
  GH1D* pDopEnergyUS = new GH1D("Dop_EnergyUS","Doppler Energy",12000,0,4000);
  GH2D* pDopSumUS = new GH2D("Dop_SummaryUS","Doppler Energy Summary",16,1,17,6000,0,3000);

  GH1D* pCoreEnergyUS_fep = new GH1D("Core_EnergyUS_fep","SeGA Core Energy FEP",12000,0,4000);
  GH1D* pCoreEnergyUS_pfep = new GH1D("Core_EnergyUS_pfep","SeGA Core Energy Projectile FEP",12000,0,4000);
  GH1D* pCoreEnergyUS_rfep = new GH1D("Core_EnergyUS_rfep","SeGA Core Energy Recoil FEP",12000,0,4000);
  GH1D* pCoreEnergyUS_nfep = new GH1D("Core_EnergyUS_nfep","SeGA Core Energy Not FEP",12000,0,4000);
  
  GH1D* pDopEnergyUS_fep = new GH1D("Dop_EnergyUS_fep","Doppler Energy FEP",12000,0,4000);
  GH1D* pDopEnergyUS_pfep = new GH1D("Dop_EnergyUS_pfep","Doppler Energy Projectile FEP",12000,0,4000);
  GH1D* pDopEnergyUS_rfep = new GH1D("Dop_EnergyUS_rfep","Doppler Energy Recoil FEP",12000,0,4000);
  GH1D* pDopEnergyUS_nfep = new GH1D("Dop_EnergyUS_nfep","Doppler Energy Not FEP",12000,0,4000);

  GH2D* pDopvPartUS = new GH2D("Dop_PartEn_US","Doppler Energy vs Particle Energy",3000,0,3000,1000,0,1000);

  GH2D* pThCorUS = new GH2D("Theta_CorrUS","Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pThCrtUS = new GH2D("Theta_CrctUS","Theta Correction",6000,0,3000,90,0,180);
  
  GH2D* pPhCorUS = new GH2D("Phi_CorrUS","Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pPhCrtUS = new GH2D("Phi_CrctUS","Phi Correction",6000,0,3000,32,0,thing1);

  GH1D* pReconEnergyUS = new GH1D("Recon_EnergyUS","Recon Energy",12000,0,4000);
  GH2D* pReconSumUS = new GH2D("Recon_SummaryUS","Recon Energy Summary",16,1,17,6000,0,3000);

  GH1D* pReconEnergyUS_fep = new GH1D("Recon_EnergyUS_fep","Recon Energy FEP",12000,0,4000);
  GH1D* pReconEnergyUS_pfep = new GH1D("Recon_EnergyUS_pfep","Recon Energy Projectile FEP",12000,0,4000);
  GH1D* pReconEnergyUS_rfep = new GH1D("Recon_EnergyUS_rfep","Recon Energy Recoil FEP",12000,0,4000);
  GH1D* pReconEnergyUS_nfep = new GH1D("Recon_EnergyUS_nfep","Recon Energy Not FEP",12000,0,4000);

  GH2D* pReconvPartUS = new GH2D("ReconEn_v_partEn_US","Recon Energy vs Particle Energy",
				3000,0,3000,1000,0,1000);

  GH2D* pReconThCorUS = new GH2D("ReconTheta_CorrUS","Recon Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pReconThCrtUS = new GH2D("ReconTheta_CrctUS","Recon Theta Correction",6000,0,3000,90,0,180);

  GH2D* pReconPhCorUS = new GH2D("ReconPhi_CorrUS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pReconPhCrtUS = new GH2D("ReconPhi_CrctUS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Recoil
  GH2D* sPidRec = new GH2D("SecPID_Rec","DS SectorEn PID",24,1,25,1000,0,1000);
  GH2D* rPidRec = new GH2D("RingPID_Rec","DS RingEn PID",24,1,25,1000,0,1000);
  
  GH1D* rCoreEnergy = new GH1D("Core_EnergyRec","SeGA Core Energy",12000,0,4000);
  GH2D* rCoreSum = new GH2D("Core_SummaryRec","Core Energy Summary",16,1,17,3000,0,3000);
  
  GH1D* rDopEnergy = new GH1D("Dop_EnergyRec","Doppler Energy",12000,0,4000);
  GH2D* rDopSum = new GH2D("Dop_SummaryRec","Doppler Energy Summary",16,1,17,6000,0,3000);

  GH1D* rCoreEnergy_fep = new GH1D("Core_EnergyRec_fep","SeGA Core Energy FEP",12000,0,4000);
  GH1D* rCoreEnergy_pfep = new GH1D("Core_EnergyRec_pfep","SeGA Core Energy Projectile FEP",12000,0,4000);
  GH1D* rCoreEnergy_rfep = new GH1D("Core_EnergyRec_rfep","SeGA Core Energy Recoil FEP",12000,0,4000);
  GH1D* rCoreEnergy_nfep = new GH1D("Core_EnergyRec_nfep","SeGA Core Energy Not FEP",12000,0,4000);

  GH1D* rDopEnergy_fep = new GH1D("Dop_EnergyRec_fep","Doppler Energy FEP",12000,0,4000);
  GH1D* rDopEnergy_pfep = new GH1D("Dop_EnergyRec_pfep","Doppler Energy Projectile FEP",12000,0,4000);
  GH1D* rDopEnergy_rfep = new GH1D("Dop_EnergyRec_rfep","Doppler Energy Recoil FEP",12000,0,4000);
  GH1D* rDopEnergy_nfep = new GH1D("Dop_EnergyRec_nfep","Doppler Energy Not FEP",12000,0,4000);

  GH2D* rDopvPart = new GH2D("DopEv_v_PartEn_Rec","Doppler Energy vs Particle Energy",3000,0,3000,1000,0,1000);

  GH2D* rCEnvRing = new GH2D("CoreEn_v_Ring_Rec","Core Energy vs Ring",
			      24,1,25,4000,0,4000);
  GH2D* rDopvRing = new GH2D("DopEn_v_Ring_Rec","Doppler Energy vs Ring",
			      24,1,25,4000,0,4000);
  GH2D* rRecvRing = new GH2D("ReconEn_v_Ring_Rec","Recon Energy vs Ring",
			    24,1,25,4000,0,4000);
  
  GH2D* rThCor = new GH2D("Theta_CorrRec","Theta Correlation",3000,0,3000,90,0,180);
  GH2D* rThCrt = new GH2D("Theta_CrctRec","Theta Correction",6000,0,3000,90,0,180);;
  GH2D* rPhCor = new GH2D("Phi_CorrRec","Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* rPhCrt = new GH2D("Phi_CrctRec","Phi Correction",6000,0,3000,32,0,thing1);

  GH2D* rcThCrt = new GH2D("cosTheta_CrctRec","CosTheta Correction",6000,0,3000,200,-1.1,1.1);
  GH2D* rcReconThCrt = new GH2D("cosRecTheta_CrctRec","CosTheta Correction",6000,0,3000,200,-1.1,1.1);

  GH1D* rReconEnergy = new GH1D("Recon_EnergyRec","Recon Energy",12000,0,4000);
  GH2D* rReconSum = new GH2D("Recon_SummaryRec","Recon Energy Summary",16,1,17,6000,0,3000);

  GH1D* rReconEnergy_fep = new GH1D("Recon_EnergyRec_fep","Recon Energy FEP",12000,0,4000);
  GH1D* rReconEnergy_pfep = new GH1D("Recon_EnergyRec_pfep","Recon Energy Projectile FEP",12000,0,4000);
  GH1D* rReconEnergy_rfep = new GH1D("Recon_EnergyRec_rfep","Recon Energy Recoil FEP",12000,0,4000);
  GH1D* rReconEnergy_nfep = new GH1D("Recon_EnergyRec_nfep","Recon Energy Not FEP",12000,0,4000);

  GH2D* rReconvPart = new GH2D("ReconEn_v_PartEn_Rec","Recon Energy vs Particle Energy",3000,0,3000,1000,0,1000);
  //GH2D* rReconvPartNS2 = new GH2D("ReconEv_v_PartEn_NoS2_Rec","Recon Energy vs Particle Energy (No Sol2)",3000,0,3000,1000,0,1000);

  GH2D* rReconThCor = new GH2D("ReconTheta_CorrRec","Recon Theta Correlation",3000,0,3000,90,0,180);
  GH2D* rReconThCrt = new GH2D("ReconTheta_CrctRec","Recon Theta Correction",6000,0,3000,90,0,180);

  GH2D* rReconPhCor = new GH2D("ReconPhi_CorrRec","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* rReconPhCrt = new GH2D("ReconPhi_CrctRec","Recon Phi Correction",6000,0,3000,32,0,thing1);

  /*
  //SeGA dets
  std::vector<GH1D*> pDetDopEnDS;
  std::vector<GH2D*> pDetDopSumDS;
  std::vector<GH2D*> pDetPhCorDS;
  std::vector<GH2D*> pDetPhCrtDS;
  std::vector<GH2D*> pDetThCorDS;
  std::vector<GH2D*> pDetThCrtDS;

  std::vector<GH1D*> pDetDopEnUS;
  std::vector<GH2D*> pDetDopSumUS;
  std::vector<GH2D*> pDetPhCorUS;
  std::vector<GH2D*> pDetPhCrtUS;
  std::vector<GH2D*> pDetThCorUS;
  std::vector<GH2D*> pDetThCrtUS;

  std::vector<GH1D*> rDetDopEn;
  std::vector<GH2D*> rDetDopSum;
  std::vector<GH2D*> rDetPhCor;
  std::vector<GH2D*> rDetPhCrt;
  std::vector<GH2D*> rDetThCor;
  std::vector<GH2D*> rDetThCrt;
  */

  std::vector<GH2D*> detPosTP;
  for(int i=0;i<16;i++) {

    detPosTP.push_back(new GH2D(Form("posPT_det%02d",i+1),Form("Det%02d Phi-Theta Surface",i+1),
				360,0,180,760,-10,370));
    
    /*
    //Downstream
    pDetDopEnDS.push_back(new GH1D(Form("Dop_EnergyDS_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				6000,0,3000));

    pDetDopSumDS.push_back(new GH2D(Form("Dop_SegSumDS_Det%02i",i+1),
                                    Form("Det%02i Doppler Seg Summary",i+1),
				    32,1,33,6000,0,3000));

    pDetPhCorDS.push_back(new GH2D(Form("PhiCorDS_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				3000,0,3000,32,0,thing1));

    pDetPhCrtDS.push_back(new GH2D(Form("PhiCrtDS_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				6000,0,3000,32,0,thing1));

    pDetThCorDS.push_back(new GH2D(Form("ThetaCorDS_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				3000,0,3000,90,0,180));

    pDetThCrtDS.push_back(new GH2D(Form("ThetaCrtDS_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				6000,0,3000,90,0,180));

    //Upstream
    pDetDopEnUS.push_back(new GH1D(Form("Dop_EnergyUS_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				   6000,0,3000));

    pDetDopSumUS.push_back(new GH2D(Form("Dop_SegSumUS_Det%02i",i+1),
                                    Form("Det%02i Doppler Seg Summary",i+1),
				    32,1,33,6000,0,3000));

    pDetPhCorUS.push_back(new GH2D(Form("PhiCorUS_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				   3000,0,3000,32,0,thing1));

    pDetPhCrtUS.push_back(new GH2D(Form("PhiCrtUS_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				   6000,0,3000,32,0,thing1));

    pDetThCorUS.push_back(new GH2D(Form("ThetaCorUS_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				   3000,0,3000,90,0,180));

    pDetThCrtUS.push_back(new GH2D(Form("ThetaCrtUS_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				   6000,0,3000,90,0,180));

    //Recoil
    rDetDopEn.push_back(new GH1D(Form("Dop_EnergyRec_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				   6000,0,3000));

    rDetDopSum.push_back(new GH2D(Form("Dop_SegSumRec_Det%02i",i+1),Form("Det%02i Doppler Seg Summary",i+1),
				    32,1,33,6000,0,3000));

    rDetPhCor.push_back(new GH2D(Form("PhiCorRec_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				   3000,0,3000,32,0,thing1));

    rDetPhCrt.push_back(new GH2D(Form("PhiCrtRec_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				   6000,0,3000,32,0,thing1));

    rDetThCor.push_back(new GH2D(Form("ThetaCorRec_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				   3000,0,3000,90,0,180));

    rDetThCrt.push_back(new GH2D(Form("ThetaCrtRec_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				   6000,0,3000,90,0,180));
    */
    
  }

  //SeGA Segments
  std::vector<GH2D*> segPosTP02;
  std::vector<GH2D*> segPosTP10;
  
  for(int i=0;i<32;i++) {
    segPosTP02.push_back(new GH2D(Form("posPT_d02_s%02d",i+1),
				  Form("Det02 Seg%02d Phi-Theta Surface",i+1),
				  360,0,180,760,-10,370));

    segPosTP10.push_back(new GH2D(Form("posPT_d10_s%02d",i+1),
				  Form("Det10 Seg%02d Phi-Theta Surface",i+1),
				  360,0,180,760,-10,370));
  } 

  //Bambino2 rings
  //std::vector<GH1D*> pRingSecDS;
  std::vector<GH1D*> pRingCoreEnDS;
  std::vector<GH1D*> pRingDopEnDS;
  std::vector<GH1D*> pRingRecEnDS;
  std::vector<GH1D*> pRingThDopDS;
  std::vector<GH1D*> pRingRecThDopDS;
  std::vector<GH2D*> pRingThvPhDS;
  std::vector<GH2D*> pRingPhCorDS;
  std::vector<GH2D*> pRingPhCrtDS;
  std::vector<GH2D*> pRingThCorDS;
  std::vector<GH2D*> pRingThCrtDS;

  //std::vector<GH1D*> pRingCoreEnUS;
  //std::vector<GH1D*> pRingDopEnUS;
  //std::vector<GH2D*> pRingPhCorUS;
  //std::vector<GH2D*> pRingPhCrtUS;
  //std::vector<GH2D*> pRingThCorUS;
  //std::vector<GH2D*> pRingThCrtUS;

  //std::vector<GH1D*> rRingSec;
  std::vector<GH1D*> rRingCoreEn;
  std::vector<GH1D*> rRingDopEn;
  std::vector<GH1D*> rRingRecEn;
  std::vector<GH1D*> rRingThDop;
  std::vector<GH1D*> rRingRecThDop;
  //std::vector<GH2D*> rRingThvPh;
  //std::vector<GH2D*> rRingPhCor;
  //std::vector<GH2D*> rRingPhCrt;
  //std::vector<GH2D*> rRingThCor;
  //std::vector<GH2D*> rRingThCrt;
  //std::vector<GH2D*> rRingRecThCor;
  //std::vector<GH2D*> rRingRecThCrt;
  
  for(int i=0;i<24;i++) {

    //Downstream
    
    //pRingSecDS.push_back(new GH1D(Form("pSecDS_R%02i",i+1),Form("Ring%02i Sectors",i+1),32,1,33));
    
    pRingCoreEnDS.push_back(new GH1D(Form("Core_EnergyDS_R%02i",i+1),Form("Ring%02i Core Energy",i+1),
				     3000,0,3000));

    pRingDopEnDS.push_back(new GH1D(Form("Dop_EnergyDS_R%02i",i+1),Form("Ring%02i Doppler Energy",i+1),
				    6000,0,3000));

    pRingRecEnDS.push_back(new GH1D(Form("Rec_EnergyDS_R%02i",i+1),Form("Ring%02i Recon Energy",i+1),
				    6000,0,3000));

    pRingThDopDS.push_back(new GH1D(Form("ThetaDS_R%02i",i+1),Form("Ring%02i Doppler Theta",i+1),
				    180,0,180));

    pRingRecThDopDS.push_back(new GH1D(Form("Recon_ThetaDS_R%02i",i+1),
				       Form("Ring%02i Recon Theta",i+1),
				       180,0,180));
    
    pRingThvPhDS.push_back(new GH2D(Form("pThvPhDS_R%02i",i+1),
    				    Form("Projectile Ring%02i #phi-#theta surface",i+1),
    				    1000,0,90,1000,-200,200));

    /*
    pRingPhCorDS.push_back(new GH2D(Form("PhiCorDS_R%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				3000,0,3000,32,0,thing1));

    pRingPhCrtDS.push_back(new GH2D(Form("PhiCrtDS_R%02i",i+1),Form("Ring%02i Phi Correction",i+1),
				6000,0,3000,32,0,thing1));

    pRingThCorDS.push_back(new GH2D(Form("ThetaCorDS_R%02i",i+1),Form("Ring%02i Theta Correlation",i+1),
				3000,0,3000,90,0,180));

    pRingThCrtDS.push_back(new GH2D(Form("ThetaCrtDS_R%02i",i+1),Form("Ring%02i Theta Correction",i+1),
				6000,0,3000,90,0,180));
    */

    //Upstream
    /*
    pRingCoreEnUS.push_back(new GH1D(Form("Core_EnergyUS_R%02i",i+1),Form("Ring%02i Core Energy",i+1),
				     3000,0,3000));

    pRingDopEnUS.push_back(new GH1D(Form("Dop_EnergyUS_R%02i",i+1),Form("Ring%02i Doppler Energy",i+1),
				    6000,0,3000));

    
    pRingPhCorUS.push_back(new GH2D(Form("PhiCorUS_R%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				    3000,0,3000,32,0,thing1));

    pRingPhCrtUS.push_back(new GH2D(Form("PhiCrtUS_R%02i",i+1),Form("Ring%02i Phi Correction",i+1),
				    6000,0,3000,32,0,thing1));

    pRingThCorUS.push_back(new GH2D(Form("ThetaCorUS_R%02i",i+1),Form("Ring%02i Theta Correlation",i+1),
				    3000,0,3000,90,0,180));

    pRingThCrtUS.push_back(new GH2D(Form("ThetaCrtUS_R%02i",i+1),Form("Ring%02i Theta Correction",i+1),
				    6000,0,3000,90,0,180));
    */

    //Recoil
    
    //rRingSec.push_back(new GH1D(Form("rSec_R%02i",i+1),Form("Ring%02i Sectors",i+1),32,1,33));
    
    rRingCoreEn.push_back(new GH1D(Form("Core_EnergyRec_R%02i",i+1),Form("Ring%02i Core Energy",i+1),
				   3000,0,3000));

    rRingDopEn.push_back(new GH1D(Form("Dop_EnergyRec_R%02i",i+1),Form("Ring%02i Doppler Energy",i+1),
				  6000,0,3000));

    rRingRecEn.push_back(new GH1D(Form("Recon_EnergyRec_R%02i",i+1),Form("Ring%02i Recon Energy",i+1),
				  6000,0,3000));

    rRingThDop.push_back(new GH1D(Form("ThetaRec_R%02i",i+1),Form("Ring%02i Doppler Theta",i+1),
				  180,0,180));

    rRingRecThDop.push_back(new GH1D(Form("Recon_ThetaRec_R%02i",i+1),
				     Form("Ring%02i Recon Theta",i+1),
				     180,0,180));
    /*
    rRingThvPh.push_back(new GH2D(Form("rThvPh_R%02i",i+1),Form("Recoil Ring%02i #phi-#theta surface",i+1),
				  1000,0,90,1000,-200,200));

    rRingPhCor.push_back(new GH2D(Form("PhiCorRec_R%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				    3000,0,3000,32,0,thing1));

    rRingPhCrt.push_back(new GH2D(Form("PhiCrtRec_R%02i",i+1),Form("Ring%02i Phi Correction",i+1),
				    6000,0,3000,32,0,thing1));

    rRingThCor.push_back(new GH2D(Form("ThetaCorRec_R%02i",i+1),Form("Ring%02i Theta Correlation",i+1),
				    3000,0,3000,90,0,180));

    rRingThCrt.push_back(new GH2D(Form("ThetaCrtRec_R%02i",i+1),Form("Ring%02i Theta Correction",i+1),
				    6000,0,3000,90,0,180));
    

    rRingRecThCor.push_back(new GH2D(Form("rRecThetaCor_R%02i",i+1),Form("Ring%02i Recon Theta Correlation",
                                          i+1),3000,0,3000,90,0,180));

    rRingRecThCrt.push_back(new GH2D(Form("rRecThetaCrt_R%02i",i+1),Form("Ring%02i Recon Theta Correction",
                                          i+1),6000,0,3000,90,0,180));
    */

  }

  std::cout << "Correlating and histograming data..." << std::endl;
  FILE* input_file = fopen(input_filename,"rb");
  
  const TVector3 incBeam = TVector3(0.0,0.0,1.0);

  TRandom* rand = new TRandom(50747227);
  const double r2d = TMath::RadToDeg();
  
  TSRIM srimB(beam_srim_name.c_str());
  TSRIM srimT(targ_srim_name.c_str());

  const double reac_en = 0.001*srimB.GetAdjustedEnergy(1000.0*beam_en,0.5*target_width,0.01);
  //const double reac_en = beam_en;
  const double Sol2_En = KE_LAB(Theta_CM_FP(Theta_LAB_Max(reac_en),reac_en),reac_en);
  
  Header header;
  JANUSData jData;
  while(fread(&header,header.bytes(),1,input_file)) {

    const int nB = header.nBdata;
    const int nS = header.nSdata;
    //const int nE = header.evtNum;
  
    fread(&jData.bData,sizeof(Bambino2Data),nB,input_file);
    fread(&jData.sData,sizeof(SegaData),nS,input_file);

    BuiltData data = BuildData(header,jData);

    //Bambino2 singles
    for(int i=0;i<data.nBa;i++) {

      int det = data.bam2[i].det;
      int ring = data.bam2[i].ring;
      int sec = data.bam2[i].sector;

      if(ring < Ri || ring > Rf)
	  continue;

      double ring_en = data.bam2[i].rEn;
      double sec_en = data.bam2[i].sEn;

      TVector3 segPos = GetPos(det,ring,sec);
      TVector3 pos = data.bam2[i].sPos;

      segPos.SetX(segPos.X() - beam_X);
      segPos.SetY(segPos.Y() - beam_Y);
      pos.SetX(pos.X() - beam_X);
      pos.SetY(pos.Y() - beam_Y);
      pos.SetZ(pos.Z() - beam_Z);
      
      if(!det) { //Upstream

	bSum->Fill(sec,sec_en);
	bSum->Fill(ring+32,ring_en);
	
	rPid0->Fill(ring,ring_en);
	sPid0->Fill(ring,sec_en);
	secD0->Fill(sec);

	if(data.bam2[i].rP && data.bam2[i].sP) { //Projectile
	  pSum->Fill(sec,sec_en);
	  pSum->Fill(ring+32,ring_en);

	  pThvPhUS->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	}
	if(data.bam2[i].rR && data.bam2[i].sR) {
	  rSum->Fill(sec,sec_en);
	  rSum->Fill(ring+32,ring_en);
	}
      }
      else { //Downstream
	
	rPid1->Fill(ring,ring_en);
	sPid1->Fill(ring,sec_en);
	secD1->Fill(sec);

	bSum->Fill(sec+64,sec_en);
	bSum->Fill(ring+96,ring_en);

	if(data.bam2[i].rP && data.bam2[i].sP) { //Projectile
	  sPid1_p->Fill(ring,sec_en);

	  pSum->Fill(sec+64,sec_en);
	  pSum->Fill(ring+96,ring_en);

	  //pPvPDS->Fill(segPos.Phi(),segPos.Perp());
	  
	  pThvPhDS->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	  if(ring%2) {
	    if(sec%2) {
	      pThvPhDS1->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	    else {
	      pThvPhDS2->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	  }
	  else {
	    if(sec%2) {
	      pThvPhDS2->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	    else {
	      pThvPhDS1->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	  }

	  pSvRDS->Fill(ring,sec);
	  pSecDS->Fill(sec);
	  //pPhiDS->Fill(segPos.Phi()*r2d);
	  
	  //pRingSecDS.at(ring-1)->Fill(sec);
	  pRingThvPhDS.at(ring-1)->Fill(pos.Theta()*r2d,pos.Phi()*r2d);

	  /*
	  if(data.nBa == 1) {
	    pSecDS_m1->Fill(sec);
	  }
	  else if(data.nBa == 2) {
	    pSecDS_m2->Fill(sec);
	  }
	  */
	    
	}

	if(data.bam2[i].rR && data.bam2[i].sR) { //Recoil
	  sPid1_r->Fill(ring,sec_en);

	  rSum->Fill(sec+64,sec_en);
	  rSum->Fill(ring+96,ring_en);

	  rSec->Fill(sec);
	  //rPhi->Fill(segPos.Phi()*r2d);
	  
	  
	  //rPvP->Fill(segPos.Phi(),segPos.Perp());
	  rThvPh->Fill(pos.Theta()*r2d,pos.Phi()*r2d);

	  /*
	  rRingThvPh.at(ring-1)->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	  rRingSec.at(ring-1)->Fill(sec);

	  if(data.nBa == 1) {
	    rSec_m1->Fill(sec);
	  }
	  else if(data.nBa == 2) {
	    rSec_m2->Fill(sec);
	  }	  
	  */
	  
	}

	if(data.nBa == 1) {
	  rPid1m1->Fill(ring,ring_en);
	  sPid1m1->Fill(ring,sec_en);
	  
	}
	else if(data.nBa == 2) {
	  rPid1m2->Fill(ring,ring_en);
	  sPid1m2->Fill(ring,sec_en);
	}

      }
      
    } //End Bambino2 singles
    
    //SeGA singles
    for(int i=0;i<data.nSe;i++) {

      int det = data.sega[i].det;
      double en = data.sega[i].cEn;
      double core_en = rand->Gaus(en,Sigma(en));
      bool fep = data.sega[i].fep;
      
      coreEnergy->Fill(core_en);
      coreSum->Fill(det,core_en);

      for(int j=0;j<data.sega[i].nsegs;j++) {

	double exactX = data.sega[i].x[j];
	double exactY = data.sega[i].y[j];
	double exactZ = data.sega[i].z[j];

	TVector3 exact_pos(exactX,exactY,exactZ);

	double phi = exact_pos.Phi();
	if(phi < 0) {
	  phi += TMath::TwoPi();
	}
	phi *= r2d;
	double theta = exact_pos.Theta()*r2d;
	
	double seg_en = data.sega[i].sEn[j];
	int seg = data.sega[i].segs[j];
	int num = seg + 32*(det-1);
	
	segEnergy->Fill(seg_en);
	segSum->Fill(num,seg_en);

	sPosTP->Fill(theta,phi);
	detPosTP.at(det-1)->Fill(theta,phi);

	if(det == 2) {
	  segPosTP02.at(seg-1)->Fill(theta,phi);
	}
	else if(det == 10) {
	  segPosTP10.at(seg-1)->Fill(theta,phi);
	}
      }
      
      if(fep) {
	coreEn_Fep->Fill(core_en);
	FEPsum->Fill(det,core_en);
      }
      else {
	coreEn_NotFep->Fill(core_en);
      }
      
    } //End loop over sega

    /*
    //Mult2
    if(data.nSe == 2) {

      double coreEn1 = data.sega[0].cEn;
      double coreEn2 = data.sega[1].cEn;

      double en1 = rand->Gaus(coreEn1,Sigma(coreEn1));
      double en2 = rand->Gaus(coreEn2,Sigma(coreEn2));

      bool good = false;
      if(en1 > 1330.0 && en1 < 1335.0) {
      
	gamGate->Fill(en2);
	if(en2 > 1170.8 && en2 < 1175.8) {
	  good = true;
	}
      
      }
      else if(en2 > 1330.0 && en2 < 1335.0) {
	
	gamGate->Fill(en1);
	if(en1 > 1170.8 && en1 < 1175.8) {
	  good = true;
	}
	
      }
      
      if(good) {

	int det1 = data.sega[0].det;
	int seg1 = data.sega[0].MainSeg();
	TVector3 pos1 = GetPos(det1,seg1);

	int det2 = data.sega[1].det;
	int seg2 = data.sega[1].MainSeg();
	TVector3 pos2 = GetPos(det2,seg2);

	double theta = pos1.Angle(pos2);
	
	gammaAng->Fill(std::cos(theta));

      }

      if(data.sega[0].fep && data.sega[1].fep) {

	int det1 = data.sega[0].det;
	int seg1 = data.sega[0].MainSeg();
	TVector3 pos1 = GetPos(det1,seg1);

	int det2 = data.sega[1].det;
	int seg2 = data.sega[1].MainSeg();
	TVector3 pos2 = GetPos(det2,seg2);

	double theta = pos1.Angle(pos2);
      
	gammaAngFEP->Fill(std::cos(theta));
      }
    } //End SeGA singles
    */
    
    //Coincidences
    if(data.nBa > 0 && data.nSe > 0) {

      bool good =false;
      for(int i=0;i<data.nBa;i++) {

	int bDet = data.bam2[i].det;
	int ring = data.bam2[i].ring;
	int sector = data.bam2[i].sector;

	if(ring < Ri || ring > Rf)
	  continue;

	double ring_en = data.bam2[i].rEn;
        double sec_en = data.bam2[i].sEn;

	TVector3 bPos = GetPos(bDet,ring,sector);
	bPos.SetX(bPos.X() - beam_X);
	bPos.SetY(bPos.Y() - beam_Y);
	
	if(bDet && data.bam2[i].rP && data.bam2[i].sP) { //Projectile DS gate

	  sPidDS->Fill(ring,sec_en);
	  rPidDS->Fill(ring,ring_en);

	  bool sol2 = false;
	  if(sec_en < Sol2_En) {
	    sol2 = true;
	  }
	  
	  double thetaCM = Theta_CM_FP(bPos.Theta(),reac_en,sol2);
	  
	  double energy = KE_LAB(thetaCM,reac_en);
	  double distance = 0.5*target_width/std::abs(std::cos(bPos.Theta()));
	  energy = 0.001*srimB.GetAdjustedEnergy(1000.0*energy,distance,0.01);
      
	  double gam = (energy/beam_mass) + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));
	  
	  TVector3 rPos(0,0,1);
	  rPos.SetTheta(Recoil_Theta_LAB(thetaCM,reac_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());

	  double recon_energy = Recoil_KE_LAB(thetaCM,reac_en);
	  distance = 0.5*target_width/std::abs(std::cos(rPos.Theta()));
	  recon_energy = 0.001*srimT.GetAdjustedEnergy(1000.*recon_energy,distance);
	  
	  double recon_gam = (recon_energy)/targ_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));
	  
	  for(int j=0;j<data.nSe;j++) {

	    int det = data.sega[j].det;
	    int seg = data.sega[j].MainSeg();
	    int map_slice;
	    if(seg < 9) {
	      map_slice = 9-seg;
	    }
	    else if(seg < 17) {
	      map_slice = 17-seg;
	    }
	    else if(seg < 25) {
	      map_slice = 25-seg;
	    }
	    else {
	      map_slice = 33-seg;
	    }

	    if(det > 8) {
	      map_slice += 8;
	    }
	    
	    double en = data.sega[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.sega[j].fep;
	    bool PFEP = data.sega[j].pfep;
	  
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);
	    
	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	    if(dopEn > 625.0 && dopEn < 640.0) {
	      good = true;
	    }
	    
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0) {
	      reac_phi += TMath::TwoPi();
	    }

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0) {
	      det_phi += TMath::TwoPi();
	    }

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0) {
	      planeAng += TMath::TwoPi();
	    }

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0) {
	      recon_phi += TMath::TwoPi();
	    }

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0) {
	      reconAng += TMath::TwoPi();
	    }
	    
	    if(FEP) {
	      
	      pCoreEnergyDS_fep->Fill(coreEn);
	      pDopEnergyDS_fep->Fill(dopEn);
	      pReconEnergyDS_fep->Fill(recon_en);

	      pcThCrtDS_fep->Fill(dopEn,TMath::Cos(theta));

	      if(PFEP) {
		pCoreEnergyDS_pfep->Fill(coreEn);
		pDopEnergyDS_pfep->Fill(dopEn);
		pReconEnergyDS_pfep->Fill(recon_en);
	      }
	      else {
		pCoreEnergyDS_rfep->Fill(coreEn);
		pDopEnergyDS_rfep->Fill(dopEn);
		pReconEnergyDS_rfep->Fill(recon_en);
	      }
	    }
	    else {
	      pCoreEnergyDS_nfep->Fill(coreEn);
	      pDopEnergyDS_nfep->Fill(dopEn);
	      pReconEnergyDS_nfep->Fill(recon_en);
	    }

	    pCoreEnergyDS->Fill(coreEn);
	    pCoreSumDS->Fill(det,coreEn);
	    
	    pDopEnergyDS->Fill(dopEn);
	    pDopSumDS->Fill(det,dopEn);
	    pDopASSumDS->Fill(map_slice,dopEn);
	      
	    pDopvPartDS->Fill(dopEn,sec_en);

	    pCEnvRingDS->Fill(ring,coreEn);
	    pDopvRingDS->Fill(ring,dopEn);
	    pRecvRingDS->Fill(ring,recon_en);
	    
	    pThCorDS->Fill(coreEn,theta*r2d);
	    pThCrtDS->Fill(dopEn,theta*r2d);

	    pcThCrtDS->Fill(dopEn,TMath::Cos(theta));

	    pPhCorDS->Fill(coreEn,planeAng*r2d);
	    pPhCrtDS->Fill(dopEn,planeAng*r2d);

	    /*
	    pDetDopEnDS.at(det-1)->Fill(dopEn);
	    pDetDopSumDS.at(det-1)->Fill(seg,dopEn);

	    pDetThCorDS.at(det-1)->Fill(coreEn,theta*r2d);
	    pDetThCrtDS.at(det-1)->Fill(dopEn,theta*r2d);

	    pDetPhCorDS.at(det-1)->Fill(coreEn,planeAng*r2d);
	    pDetPhCrtDS.at(det-1)->Fill(dopEn,planeAng*r2d);
	    */
	    
	    pRingCoreEnDS.at(ring-1)->Fill(coreEn);
	    pRingDopEnDS.at(ring-1)->Fill(dopEn);

	    pRingThDopDS.at(ring-1)->Fill(theta*r2d);

	    //pRingThCorDS.at(ring-1)->Fill(coreEn,theta*r2d);
	    //pRingThCrtDS.at(ring-1)->Fill(dopEn,theta*r2d);

	    //pRingPhCorDS.at(ring-1)->Fill(coreEn,planeAng*r2d);
	    //pRingPhCrtDS.at(ring-1)->Fill(dopEn,planeAng*r2d);  
	    
	    pReconEnergyDS->Fill(recon_en);
	    pReconSumDS->Fill(det,recon_en);

	    pReconvPartDS->Fill(recon_en,sec_en);

	    pReconThCorDS->Fill(coreEn,recon_theta*r2d);
	    pReconThCrtDS->Fill(recon_en,recon_theta*r2d);

	    pReconPhCorDS->Fill(coreEn,reconAng*r2d);
	    pReconPhCrtDS->Fill(recon_en,reconAng*r2d);

	    pRingRecEnDS.at(ring-1)->Fill(recon_en);

	    pRingRecThDopDS.at(ring-1)->Fill(recon_theta*r2d);
	    
	  } //End SeGA loop
	} //End DS projectile gate

	else if(!bDet && data.bam2[i].rP && data.bam2[i].sP) { //Projectile US gate

	  sPidUS->Fill(ring,sec_en);
	  rPidUS->Fill(ring,ring_en);

	  double thetaCM = Theta_CM_FP(bPos.Theta(),reac_en,false);
	  
	  double energy = KE_LAB(thetaCM,reac_en);
	  double distance = 0.5*target_width/std::abs(std::cos(bPos.Theta()));
	  energy += 0.001*srimB.GetEnergyChange(1000.*energy,distance);
	  
	  double gam = (energy)/beam_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  TVector3 rPos(0,0,1);
	  rPos.SetTheta(Recoil_Theta_LAB(thetaCM,reac_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());

	  double recon_energy = Recoil_KE_LAB(thetaCM,reac_en);
	  distance = 0.5*target_width/std::abs(std::cos(rPos.Theta()));
	  recon_energy += 0.001*srimT.GetEnergyChange(1000.*recon_energy,distance);
	  
	  double recon_gam = (recon_energy)/targ_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));
	  
	  for(int j=0;j<data.nSe;j++) {

	    int det = data.sega[j].det;
	    int seg = data.sega[j].MainSeg();
	    
	    double en = data.sega[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.sega[j].fep;
	    bool PFEP = data.sega[j].pfep;
	  
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	  
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0) {
	      reac_phi += TMath::TwoPi();
	    }

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0) {
	      det_phi += TMath::TwoPi();
	    }

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0) {
	      planeAng += TMath::TwoPi();
	    }

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0) {
	      recon_phi += TMath::TwoPi();
	    }

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0) {
	      reconAng += TMath::TwoPi();
	    }

	    pCoreEnergyUS->Fill(coreEn);
	    pCoreSumUS->Fill(det,coreEn);
	    
	    pDopEnergyUS->Fill(dopEn);
	    pDopSumUS->Fill(det,dopEn);

	    if(FEP) {
	      
	      pCoreEnergyUS_fep->Fill(coreEn);
	      pDopEnergyUS_fep->Fill(dopEn);
	      pReconEnergyUS_fep->Fill(recon_en);
	      
	      if(PFEP) {
		pCoreEnergyUS_pfep->Fill(coreEn);
		pDopEnergyUS_pfep->Fill(dopEn);
		pReconEnergyUS_pfep->Fill(recon_en);
	      }
	      else {
		pCoreEnergyUS_rfep->Fill(coreEn);
		pDopEnergyUS_rfep->Fill(dopEn);
		pReconEnergyUS_rfep->Fill(recon_en);
	      }
	    }
	    else {
	      pCoreEnergyUS_nfep->Fill(coreEn);
	      pDopEnergyUS_nfep->Fill(dopEn);
	      pReconEnergyUS_nfep->Fill(recon_en);
	    }
	    
	    pDopvPartUS->Fill(dopEn,sec_en);

	    pThCorUS->Fill(coreEn,theta*r2d);
	    pThCrtUS->Fill(dopEn,theta*r2d);

	    pPhCorUS->Fill(coreEn,planeAng*r2d);
	    pPhCrtUS->Fill(dopEn,planeAng*r2d);  

	    /*
	    pDetDopEnUS.at(det-1)->Fill(dopEn);
	    pDetDopSumUS.at(det-1)->Fill(seg,dopEn);

	    pDetThCorUS.at(det-1)->Fill(coreEn,theta*r2d);
	    pDetThCrtUS.at(det-1)->Fill(dopEn,theta*r2d);

	    pDetPhCorUS.at(det-1)->Fill(coreEn,planeAng*r2d);
	    pDetPhCrtUS.at(det-1)->Fill(dopEn,planeAng*r2d);

	    pRingCoreEnUS.at(ring-1)->Fill(coreEn);
	    pRingDopEnUS.at(ring-1)->Fill(dopEn);

	    pRingThCorUS.at(ring-1)->Fill(coreEn,theta*r2d);
	    pRingThCrtUS.at(ring-1)->Fill(dopEn,theta*r2d);

	    pRingPhCorUS.at(ring-1)->Fill(coreEn,planeAng*r2d);
	    pRingPhCrtUS.at(ring-1)->Fill(dopEn,planeAng*r2d);
	    */

	    pReconEnergyUS->Fill(recon_en);
	    pReconSumUS->Fill(det,recon_en);

	    pReconvPartUS->Fill(recon_en,sec_en);

	    pReconThCorUS->Fill(coreEn,recon_theta*r2d);
	    pReconThCrtUS->Fill(recon_en,recon_theta*r2d);

	    pReconPhCorUS->Fill(coreEn,reconAng*r2d);
	    pReconPhCrtUS->Fill(recon_en,reconAng*r2d);
	    
	  
	  } //End SeGA loop
	} ////End US projectile gate

	else if(data.bam2[i].rR && data.bam2[i].sR) { //Recoil gate

	  sPidRec->Fill(ring,sec_en);
	  rPidRec->Fill(ring,ring_en);

	  double thetaCM = Theta_CM_FR(bPos.Theta(),reac_en);
	  
	  double energy = Recoil_KE_LAB(thetaCM,reac_en);
	  double distance = 0.5*target_width/std::abs(std::cos(bPos.Theta()));
	  energy += 0.001*srimT.GetEnergyChange(1000.*energy,distance);
	  
	  double gam = (energy)/targ_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  TVector3 rPos(0,0,1);
	  rPos.SetTheta(Theta_LAB(thetaCM,reac_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());

	  double recon_energy = KE_LAB(thetaCM,reac_en);
	  distance = 0.5*target_width/std::abs(std::cos(rPos.Theta()));
	  recon_energy += 0.001*srimB.GetEnergyChange(1000.*recon_energy,distance);
	  
	  double recon_gam = (recon_energy)/beam_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));

	  for(int j=0;j<data.nSe;j++) {

	    int det = data.sega[j].det;
	    int seg = data.sega[j].MainSeg();
	    
	    double en = data.sega[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.sega[j].fep;
	    bool PFEP = data.sega[j].pfep;
	  
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	  
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0) {
	      reac_phi += TMath::TwoPi();
	    }

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0) {
	      det_phi += TMath::TwoPi();
	    }

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0) {
	      planeAng += TMath::TwoPi();
	    }

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;
	    if(recon_en > 625.0 && recon_en < 640.0) {
	      good = true;
	    }

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0) {
	      recon_phi += TMath::TwoPi();
	    }

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0) {
	      reconAng += TMath::TwoPi();
	    }
	    
	    rCoreEnergy->Fill(coreEn);
	    rCoreSum->Fill(det,coreEn);
	    
	    rDopEnergy->Fill(dopEn);
	    rDopSum->Fill(det,dopEn);

	    if(FEP) {
	      rCoreEnergy_fep->Fill(coreEn);
	      rDopEnergy_fep->Fill(dopEn);
	      rReconEnergy_fep->Fill(recon_en);

	      if(PFEP) {
		rCoreEnergy_pfep->Fill(coreEn);
		rDopEnergy_pfep->Fill(dopEn);
		rReconEnergy_pfep->Fill(recon_en);
	      }
	      else {
		rCoreEnergy_rfep->Fill(coreEn);
		rDopEnergy_rfep->Fill(dopEn);
		rReconEnergy_rfep->Fill(recon_en);
	      }
	    }
	    else {
	      rCoreEnergy_nfep->Fill(coreEn);
	      rDopEnergy_nfep->Fill(dopEn);
	      rReconEnergy_nfep->Fill(recon_en);
	    }

	    rDopvPart->Fill(dopEn,sec_en);

	    rCEnvRing->Fill(ring,coreEn);
	    rDopvRing->Fill(ring,dopEn);
	    rRecvRing->Fill(ring,recon_en);

	    rThCor->Fill(coreEn,theta*r2d);
	    rThCrt->Fill(dopEn,theta*r2d);

	    rcThCrt->Fill(dopEn,TMath::Cos(theta));
	    
	    rPhCor->Fill(coreEn,planeAng*r2d);
	    rPhCrt->Fill(dopEn,planeAng*r2d);

	    /*
	    rDetDopEn.at(det-1)->Fill(dopEn);
	    rDetDopSum.at(det-1)->Fill(seg,dopEn);

	    rDetThCor.at(det-1)->Fill(coreEn,theta*r2d);
	    rDetThCrt.at(det-1)->Fill(dopEn,theta*r2d);

	    rDetPhCor.at(det-1)->Fill(coreEn,planeAng*r2d);
	    rDetPhCrt.at(det-1)->Fill(dopEn,planeAng*r2d);
	    */
	    
	    rRingCoreEn.at(ring-1)->Fill(coreEn);
	    rRingDopEn.at(ring-1)->Fill(dopEn);

	    rRingThDop.at(ring-1)->Fill(theta*r2d);

	    /*
	    rRingThCor.at(ring-1)->Fill(coreEn,theta*r2d);
	    rRingThCrt.at(ring-1)->Fill(dopEn,theta*r2d);

	    rRingPhCor.at(ring-1)->Fill(coreEn,planeAng*r2d);
	    rRingPhCrt.at(ring-1)->Fill(dopEn,planeAng*r2d);
	    */

	    rReconEnergy->Fill(recon_en);
	    rReconSum->Fill(det,recon_en);

	    rReconvPart->Fill(recon_en,sec_en);

	    rReconThCor->Fill(coreEn,recon_theta*r2d);
	    rReconThCrt->Fill(recon_en,recon_theta*r2d);

	    rcReconThCrt->Fill(recon_en,TMath::Cos(recon_theta));

	    rReconPhCor->Fill(coreEn,reconAng*r2d);
	    rReconPhCrt->Fill(recon_en,reconAng*r2d);

	    rRingRecEn.at(ring-1)->Fill(recon_en);
	    rRingRecThDop.at(ring-1)->Fill(recon_theta*r2d);

	    //rRingRecThCor.at(ring-1)->Fill(coreEn,recon_theta*r2d);
	    //rRingRecThCrt.at(ring-1)->Fill(recon_en,recon_theta*r2d);
	    
	  } //End SeGA loop  
	} //End recoil gate
	
      } //End Bambino2 loop

      if(good) {
	for(int i=0;i<data.nBa;i++) {

	int bDet = data.bam2[i].det;
	int ring = data.bam2[i].ring;

	if(ring < Ri || ring > Rf) {
	  continue;
	}

        double sec_en = data.bam2[i].sEn;

	if(bDet) {

	  sPidDS_gam->Fill(ring,sec_en);
	  if(data.bam2[i].rP && data.bam2[i].sP) {
	    psPidDS_gam->Fill(ring,sec_en);
	  }

	  if(bDet && data.bam2[i].rR && data.bam2[i].sR) {
	    rsPidDS_gam->Fill(ring,sec_en);
	  }
	}
	
	} //End Bambino2 loop
      } // if(good)
      
    } //End coincidences

    
  } //End while loop
  fclose(input_file);

  std::cout << "Writing histograms to file..." << std::endl;
  
  TFile* outFile = new TFile(output_filename,"RECREATE");
  
  outFile->mkdir("SeGA");
  outFile->mkdir("SeGA/Dets");
  outFile->mkdir("SeGA/Dets/Det02");
  outFile->mkdir("SeGA/Dets/Det10");
  
  outFile->mkdir("Bambino2");
  //outFile->mkdir("Bambino2/Rings");
  outFile->mkdir("Bambino2/Rings2D");

  outFile->mkdir("Coincidence");
  outFile->mkdir("Coincidence/ProjectileDS");
  outFile->mkdir("Coincidence/ProjectileDS/Rings");
  outFile->mkdir("Coincidence/ProjectileDS/Doppler");
  outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings");
  //outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings/PhiCorr");
  //outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings/ThetaCorr");

  outFile->mkdir("Coincidence/ProjectileDS/Recon");
  outFile->mkdir("Coincidence/ProjectileDS/Recon/Rings");
  
  /*
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets");
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets/ThetaCorr");
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets/Summaries");
  */

  outFile->mkdir("Coincidence/ProjectileUS");
  outFile->mkdir("Coincidence/ProjectileUS/Doppler");
  outFile->mkdir("Coincidence/ProjectileUS/Recon");
  //outFile->mkdir("Coincidence/ProjectileUS/Doppler/Rings");
  
  /*
  outFile->mkdir("Coincidence/ProjectileUS/Bambino2Rings");
  outFile->mkdir("Coincidence/ProjectileUS/Bambino2Rings/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileUS/Bambino2Rings/ThetaCorr");
  
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets");
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets/ThetaCorr");
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets/Summaries");
  */
  
  outFile->mkdir("Coincidence/Recoil");
  outFile->mkdir("Coincidence/Recoil/Rings");
  outFile->mkdir("Coincidence/Recoil/Doppler");
  outFile->mkdir("Coincidence/Recoil/Doppler/Rings");
  outFile->mkdir("Coincidence/Recoil/Recon");
  outFile->mkdir("Coincidence/Recoil/Recon/Rings");
  
  /*
  outFile->mkdir("Coincidence/Recoil/Bambino2Rings");
  outFile->mkdir("Coincidence/Recoil/Bambino2Rings/PhiCorr");
  outFile->mkdir("Coincidence/Recoil/Bambino2Rings/ThetaCorr");

  outFile->mkdir("Coincidence/Recoil/SegaDets");
  outFile->mkdir("Coincidence/Recoil/SegaDets/PhiCorr");
  outFile->mkdir("Coincidence/Recoil/SegaDets/ThetaCorr");
  outFile->mkdir("Coincidence/Recoil/SegaDets/Summaries");
  */

  outFile->cd("Bambino2");

  bSum->Write();
  pSum->Write();
  rSum->Write();

  //pPvPDS->Write();
  pThvPhDS->Write();
  pThvPhDS1->Write();
  pThvPhDS2->Write();

  //rPvP->Write();
  rThvPh->Write();
  pThvPhUS->Write();
  
  rPid0->Write();
  sPid0->Write();
  secD0->Write();
  
  rPid1->Write();
  sPid1->Write();
  secD1->Write();

  pSvRDS->Write();
  pSecDS->Write();
  //pSecDS_m1->Write();
  //pSecDS_m2->Write();
  
  rSec->Write();
  //rSec_m1->Write();
  //rSec_m2->Write();

  //pPhiDS->Write();
  //rPhi->Write();

  sPid1_p->Write();
  sPid1_r->Write();

  rPid1m1->Write();
  sPid1m1->Write();
  rPid1m2->Write();
  sPid1m2->Write();

  /*
  outFile->cd("Bambino2/Rings");
  for(int i=0;i<24;i++) {
    pRingSecDS.at(i)->Write();
    rRingSec.at(i)->Write();
  }
  */

  outFile->cd("Bambino2/Rings2D");
  for(int i=0;i<24;i++) {
    pRingThvPhDS.at(i)->Write();
    //rRingThvPh.at(i)->Write();
  }
  
  outFile->cd("SeGA");

  coreEnergy->Write();
  coreSum->Write();
  segEnergy->Write();
  segSum->Write();

  coreEn_Fep->Write();
  coreEn_NotFep->Write();
  FEPsum->Write();

  //gamGate->Write();
  //gammaAng->Write();
  //gammaAngFEP->Write();

  sPosTP->Write();

  outFile->cd("SeGA/Dets");
  for(int i=0;i<16;i++) {
    detPosTP.at(i)->Write();
  }

  outFile->cd("SeGA/Dets/Det02");
  for(int i=0;i<32;i++) {
    segPosTP02.at(i)->Write();
  }

  outFile->cd("SeGA/Dets/Det10");
  for(int i=0;i<32;i++) {
    segPosTP10.at(i)->Write();
  }

  outFile->cd("Coincidence");
  sPidDS_gam->Write();
    
  outFile->cd("Coincidence/ProjectileDS");

  sPidDS->Write();
  rPidDS->Write();

  psPidDS_gam->Write();

  pCoreEnergyDS->Write();
  pCoreSumDS->Write();

  pCoreEnergyDS_fep->Write();
  pCoreEnergyDS_pfep->Write();
  pCoreEnergyDS_rfep->Write();
  pCoreEnergyDS_nfep->Write();

  pCEnvRingDS->Write();
  
  outFile->cd("Coincidence/ProjectileDS/Doppler");
  
  pDopEnergyDS->Write();
  pDopSumDS->Write();
  pDopASSumDS->Write();

  pDopEnergyDS_fep->Write();
  pDopEnergyDS_pfep->Write();
  pDopEnergyDS_rfep->Write();
  pDopEnergyDS_nfep->Write();

  pDopvPartDS->Write();
  pDopvRingDS->Write();

  pThCorDS->Write();
  pThCrtDS->Write();

  pcThCrtDS->Write();
  pcThCrtDS_fep->Write();

  pPhCorDS->Write();
  pPhCrtDS->Write();

  
  outFile->cd("Coincidence/ProjectileDS/Rings");
  for(int i=0;i<24;i++) {
    pRingCoreEnDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/Doppler/Rings");
  for(int i=0;i<24;i++) {
    pRingDopEnDS.at(i)->Write();
    pRingThDopDS.at(i)->Write();
  }

  /*
  outFile->cd("Coincidence/ProjectileDS/Doppler/Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    pRingPhCorDS.at(i)->Write();
    pRingPhCrtDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/Doppler/Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    pRingThCorDS.at(i)->Write();
    pRingThCrtDS.at(i)->Write();
  }
  */
  
  outFile->cd("Coincidence/ProjectileDS/Recon");

  pReconEnergyDS->Write();
  pReconSumDS->Write();

  pReconEnergyDS_fep->Write();
  pReconEnergyDS_pfep->Write();
  pReconEnergyDS_rfep->Write();
  pReconEnergyDS_nfep->Write();

  pReconvPartDS->Write();
  pRecvRingDS->Write();

  pReconThCorDS->Write();
  pReconThCrtDS->Write();

  pReconPhCorDS->Write();
  pReconPhCrtDS->Write();

  
  outFile->cd("Coincidence/ProjectileDS/Recon/Rings");
  for(int i=0;i<24;i++) {
    pRingRecEnDS.at(i)->Write();
    pRingRecThDopDS.at(i)->Write();
  }
  
  /*
  outFile->cd("Coincidence/ProjectileDS/SegaDets");
  for(int i=0;i<16;i++) {
    pDetDopEnDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets/Summaries");
  for(int i=0;i<16;i++) {
    pDetDopSumDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets/PhiCorr");
  for(int i=0;i<16;i++) {
    pDetPhCorDS.at(i)->Write();
    pDetPhCrtDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets/ThetaCorr");
  for(int i=0;i<16;i++) {
    pDetThCorDS.at(i)->Write();
    pDetThCrtDS.at(i)->Write();
  }
  */

  outFile->cd("Coincidence/ProjectileUS");

  sPidUS->Write();
  rPidUS->Write();

  pCoreEnergyUS->Write();
  pCoreSumUS->Write();

  pCoreEnergyUS_fep->Write();
  pCoreEnergyUS_pfep->Write();
  pCoreEnergyUS_rfep->Write();
  pCoreEnergyUS_nfep->Write();
  
  outFile->cd("Coincidence/ProjectileUS/Doppler");
  
  pDopEnergyUS->Write();
  pDopSumUS->Write();

  pDopEnergyUS_fep->Write();
  pDopEnergyUS_pfep->Write();
  pDopEnergyUS_rfep->Write();
  pDopEnergyUS_nfep->Write();

  pDopvPartUS->Write();

  pThCorUS->Write();
  pThCrtUS->Write();

  pPhCorUS->Write();
  pPhCrtUS->Write();

  outFile->cd("Coincidence/ProjectileUS/Recon");

  pReconEnergyUS->Write();
  pReconSumUS->Write();

  pReconEnergyUS_fep->Write();
  pReconEnergyUS_pfep->Write();
  pReconEnergyUS_rfep->Write();
  pReconEnergyUS_nfep->Write();

  pReconvPartUS->Write();

  pReconThCorUS->Write();
  pReconThCrtUS->Write();

  pReconPhCorUS->Write();
  pReconPhCrtUS->Write();

  /*
  outFile->cd("Coincidence/ProjectileUS/Doppler/Rings");
  for(int i=0;i<24;i++) {
    pRingCoreEnUS.at(i)->Write();
    pRingDopEnUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/Bambino2Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    pRingPhCorUS.at(i)->Write();
    pRingPhCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/Bambino2Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    pRingThCorUS.at(i)->Write();
    pRingThCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets");
  for(int i=0;i<16;i++) {
    pDetDopEnUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets/Summaries");
  for(int i=0;i<16;i++) {
    pDetDopSumUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets/PhiCorr");
  for(int i=0;i<16;i++) {
    pDetPhCorUS.at(i)->Write();
    pDetPhCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets/ThetaCorr");
  for(int i=0;i<16;i++) {
    pDetThCorUS.at(i)->Write();
    pDetThCrtUS.at(i)->Write();
  }
  */
  
  outFile->cd("Coincidence/Recoil");

  sPidRec->Write();
  rPidRec->Write();

  rsPidDS_gam->Write();
  
  rCoreEnergy->Write();
  rCoreSum->Write();

  rCoreEnergy_fep->Write();
  rCoreEnergy_pfep->Write();
  rCoreEnergy_rfep->Write();
  rCoreEnergy_nfep->Write();

  rCEnvRing->Write();;
  
  outFile->cd("Coincidence/Recoil/Doppler");
  
  rDopEnergy->Write();
  rDopSum->Write();

  rDopEnergy_fep->Write();
  rDopEnergy_pfep->Write();
  rDopEnergy_rfep->Write();
  rDopEnergy_nfep->Write();

  rDopvPart->Write();
  rDopvRing->Write();

  rThCor->Write();
  rThCrt->Write();

  rcThCrt->Write();

  rPhCor->Write();
  rPhCrt->Write();

  outFile->cd("Coincidence/Recoil/Rings");
  for(int i=0;i<24;i++) {
    rRingCoreEn.at(i)->Write();
  }
  
  outFile->cd("Coincidence/Recoil/Doppler/Rings");
  for(int i=0;i<24;i++) {
    rRingDopEn.at(i)->Write();
    rRingThDop.at(i)->Write();
  }
  
  outFile->cd("Coincidence/Recoil/Recon");

  rReconEnergy->Write();
  rReconSum->Write();

  rReconEnergy_fep->Write();
  rReconEnergy_pfep->Write();
  rReconEnergy_rfep->Write();
  rReconEnergy_nfep->Write();

  rReconvPart->Write();
  rRecvRing->Write();

  rReconThCor->Write();
  rReconThCrt->Write();

  rReconPhCor->Write();
  rReconPhCrt->Write();

  rcReconThCrt->Write();

  
  outFile->cd("Coincidence/Recoil/Recon/Rings");
  for(int i=0;i<24;i++) {
    rRingRecEn.at(i)->Write();
    rRingRecThDop.at(i)->Write();

    //rRingRecThCor.at(i)->Write();
    //rRingRecThCrt.at(i)->Write();
  }
  
  /*
  outFile->cd("Coincidence/Recoil/Bambino2Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    rRingPhCor.at(i)->Write();
    rRingPhCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/Bambino2Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    rRingThCor.at(i)->Write();
    rRingThCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets");
  for(int i=0;i<16;i++) {
    rDetDopEn.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets/Summaries");
  for(int i=0;i<16;i++) {
    rDetDopSum.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets/PhiCorr");
  for(int i=0;i<16;i++) {
    rDetPhCor.at(i)->Write();
    rDetPhCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets/ThetaCorr");
  for(int i=0;i<16;i++) {
    rDetThCor.at(i)->Write();
    rDetThCrt.at(i)->Write();
  }
  */

  outFile->Close();

  std::cout << "Done!" << std::endl;

  return 0;
}

