#include "g4tcx_helper.h"
int main(int argc, char** argv) {
  
  if(argc != 3) {
    std::cerr << "Usage: g4tcx_sim_correlator INPUT_FILE OUTPUT_FILE" << std::endl;
    return 1;
  }

  const char* input_filename = argv[1];
  const char* output_filename = argv[2];
  if(!strcmp(input_filename,output_filename)) {
    std::cout << "Give your input and output files different names" << std::endl;
    return 1;
  }

  FILE* input_file = fopen(input_filename,"rb");
  if(input_file == NULL) {
    std::cout << "Could not open file " << input_filename << std::endl;
    return 1;
  }
  
  //S3 singles
  GH2D* bSum = new GH2D("Summary","S3 Summary",120,1,121,1000,0,1000);
  GH2D* pSum = new GH2D("pSummary","S3 Projectile Summary",120,1,121,1000,0,1000);
  GH2D* rSum = new GH2D("rSummary","S3 Recoil Summary",120,1,121,1000,0,1000);

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

  std::vector<GH2D*> pSumRingDS;
  std::vector<GH2D*> rSumRing;

  for(int i=0;i<24;i++) {

    pSumRingDS.push_back(new GH2D(Form("pSumDS_R%02d",i+1),
				  Form("Downstream Projectile Summary Ring %02d",i+1),
				  56,0,56,5000,0,500));

    rSumRing.push_back(new GH2D(Form("rSum_R%02d",i+1),Form("Recoil Summary Ring %02d",i+1),
				56,0,56,5000,0,500));
  }

  std::vector<GH2D*> pSumSecDS;
  std::vector<GH2D*> rSumSec;

  for(int i=0;i<32;i++) {

    pSumSecDS.push_back(new GH2D(Form("pSumDS_S%02d",i+1),
				 Form("Downstream Projectile Summary Sec %02d",i+1),
				 56,0,56,5000,0,500));

    rSumSec.push_back(new GH2D(Form("rSum_S%02d",i+1),Form("Recoil Summary Sec %02d",i+1),
			       56,0,56,5000,0,500));
  }
  
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
  
  //Tigress singles
  GH1D* coreEnergy = new GH1D("Core_Energy","Tigress Core Energy",16000,0,4000);
  GH2D* coreSum = new GH2D("Core_Summary","Core Energy Summary",64,1,65,6000,0,3000);
  
  GH1D* segEnergy = new GH1D("Seg_Energy","Tigress Segment Energy",16000,0,4000);
  GH2D* segSum = new GH2D("Seg_Summary","Segment Energy Summary",512,1,513,3000,0,3000);

  GH1D* addEnergy = new GH1D("Add_Energy","Tigress Addback Energy",16000,0,4000);
  GH2D* addSum = new GH2D("Add_Summary","Addback Energy Summary",64,1,65,6000,0,3000);

  GH1D* supEnergy = new GH1D("Sup_Energy","Tigress Suppressed Energy",16000,0,4000);
  GH2D* supSum = new GH2D("Sup_Summary","Suppressed Energy Summary",64,1,65,6000,0,3000);

  GH1D* addSupEn = new GH1D("AddSup_En","Tigress Suppressed Addback Energy",16000,0,4000);
  GH2D* addSupSum = new GH2D("AddSup_Sum","Suppressed Addback Energy Summary",64,1,65,6000,0,3000);
  
  GH1D* coreEn_Fep = new GH1D("FEP","Tigress FEP",16000,0,4000);
  GH1D* coreEn_NotFep = new GH1D("nFEP","Tigress Not FEP",16000,0,4000);

  GH2D* FEPsum = new GH2D("FEP_Summary","Core Energy Summary FEP",64,1,65,3000,0,6000);

  GH1D* gamGate = new GH1D("gamGate","CoreEn with 1332 keV Coinc",16000,0,4000);
  
  GH1D* gammaAngDetFEP = new GH1D("gamAngDetFEP","Gamma Gamma Angle (Detector pos)",52,-1.1,1.1);
  GH1D* gammaAngFEP = new GH1D("gamAngFEP","Gamma Gamma Angle",110,-1.1,1.1);
  GH1D* gamAngExMsFEP = new GH1D("gamAngExMsFEP","Gamma Gamma Angle (exact main seg pos)",
				 110,-1.1,1.1);
  GH1D* gamAngExFEP = new GH1D("gamAngExFEP","Gamma Gamma Angle (exact pos)",110,-1.1,1.1);
  
  GH2D* sPosTP = new GH2D("sPosPT","Tigress Phi-Theta Surface",360,0,180,760,-10,370);
  
  //Coincidences
  GH2D* sPidDS_gam = new GH2D("SecPID_DS_gam","SectorEn PID Gamma Gate",24,1,25,1500,0,1000);
  GH2D* psPidDS_gam = new GH2D("pSecPID_DS_gam","Projectile SectorEn PID Gamma Gate",
			       24,1,25,1500,0,1000);
  GH2D* rsPidDS_gam = new GH2D("rSecPID_DS_gam","Recoil SectorEn PID Gamma Gate",
			       24,1,25,1500,0,1000);
   
  //Projectile DS
  GH2D* sPidDS = new GH2D("SecPID_DS","DS SectorEn PID",24,1,25,1500,0,1000);
  GH2D* rPidDS = new GH2D("RingPID_DS","DS RingEn PID",24,1,25,1500,0,1000);
  
  GH1D* pCoreEnergyDS = new GH1D("Core_EnergyDS","Tigress Core Energy",16000,0,4000);
  GH2D* pCoreSumDS = new GH2D("Core_SummaryDS","Core Energy Summary",64,1,65,3000,0,3000);

  GH1D* pDopEnergyDS = new GH1D("Dop_EnergyDS","Doppler Energy",16000,0,4000);
  GH2D* pDopSumDS = new GH2D("Dop_SummaryDS","Doppler Energy Summary",64,1,65,6000,0,3000);

  GH1D* pAddDS = new GH1D("Add_EnergyDS","Addback Energy",16000,0,4000);
  GH2D* pAddSumDS = new GH2D("Add_SummaryDS","Addback Energy Summary",64,1,65,6000,0,3000);

  GH1D* pAddDopDS = new GH1D("AddDop_EnergyDS","Addback Doppler Energy",16000,0,4000);
  GH2D* pAddDopSumDS = new GH2D("AddDop_SummaryDS","Addback Doppler Energy Summary",64,1,65,6000,0,3000);

  GH1D* pSupDopDS = new GH1D("SupDop_EnDS","Suppressed Doppler Energy",16000,0,4000);
  GH1D* pAddSupDopDS = new GH1D("AddSupDop_EnDS","Suppressed Addback Doppler Energy",16000,0,4000);
  
  GH1D* pCoreEnergyDS_fep = new GH1D("Core_EnergyDS_fep","Tigress Core Energy FEP",16000,0,4000);
  GH1D* pCoreEnergyDS_nfep = new GH1D("Core_EnergyDS_nfep","Tigress Core Energy Not FEP",16000,0,4000);
  GH1D* pCoreEnergyDS_pfep = new GH1D("Core_EnergyDS_pfep","Tigress Core Energy Projectile FEP",16000,0,4000);
  GH1D* pCoreEnergyDS_rfep = new GH1D("Core_EnergyDS_rfep","Tigress Core Energy Recoil FEP",16000,0,4000);

  GH1D* pDopEnergyDS_fep = new GH1D("Dop_EnergyDS_fep","Doppler Energy FEP",16000,0,4000);
  GH1D* pDopEnergyDS_pfep = new GH1D("Dop_EnergyDS_pfep","Doppler Energy Projectile FEP",16000,0,4000);
  GH1D* pDopEnergyDS_rfep = new GH1D("Dop_EnergyDS_rfep","Doppler Energy Recoil FEP",16000,0,4000);
  GH1D* pDopEnergyDS_nfep = new GH1D("Dop_EnergyDS_nfep","Doppler Energy Not FEP",16000,0,4000);

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

  GH1D* pReconEnergyDS = new GH1D("Recon_EnergyDS","Recon Energy",16000,0,4000);
  GH2D* pReconSumDS = new GH2D("Recon_SummaryDS","Recon Energy Summary",64,1,65,6000,0,3000);

  GH1D* pSupReconEnDS = new GH1D("SupRecon_EnDS","Suppresed Recon Energy",16000,0,4000);
  GH1D* pAddReconEnDS = new GH1D("AddRecon_EnDS","Addback Recon Energy",16000,0,4000);
  GH1D* pAddSupReconEnDS = new GH1D("AddSupRecon_EnDS","Addback Suppresed Recon Energy",16000,0,4000);

  GH1D* pReconEnergyDS_fep = new GH1D("Recon_EnergyDS_fep","Recon Energy FEP",16000,0,4000);
  GH1D* pReconEnergyDS_pfep = new GH1D("Recon_EnergyDS_pfep","Recon Energy Projectile FEP",16000,0,4000);
  GH1D* pReconEnergyDS_rfep = new GH1D("Recon_EnergyDS_rfep","Recon Energy Recoil FEP",16000,0,4000);
  GH1D* pReconEnergyDS_nfep = new GH1D("Recon_EnergyDS_nfep","Recon Energy Not FEP",16000,0,4000);

  GH2D* pReconvPartDS = new GH2D("ReconEn_v_partEn_DS","Recon Energy vs Particle Energy",
				3000,0,3000,1000,0,1000);

  GH2D* pReconThCorDS = new GH2D("ReconTheta_CorrDS","Recon Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pReconThCrtDS = new GH2D("ReconTheta_CrctDS","Recon Theta Correction",6000,0,3000,90,0,180);

  GH2D* pReconPhCorDS = new GH2D("ReconPhi_CorrDS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pReconPhCrtDS = new GH2D("ReconPhi_CrctDS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Projectile US
  GH2D* sPidUS = new GH2D("SecPID_US","DS SectorEn PID",24,1,25,1000,0,500);
  GH2D* rPidUS = new GH2D("RingPID_US","DS RingEn PID",24,1,25,1000,0,500);
  
  GH1D* pCoreEnergyUS = new GH1D("Core_EnergyUS","Tigress Core Energy",16000,0,4000);
  GH2D* pCoreSumUS = new GH2D("Core_SummaryUS","Core Energy Summary",64,1,65,3000,0,3000);
  
  GH1D* pDopEnergyUS = new GH1D("Dop_EnergyUS","Doppler Energy",16000,0,4000);
  GH2D* pDopSumUS = new GH2D("Dop_SummaryUS","Doppler Energy Summary",64,1,65,6000,0,3000);

  GH1D* pAddDopUS = new GH1D("AddDop_EnUS","Addback Doppler Energy",16000,0,4000);
  GH1D* pSupDopUS = new GH1D("SupDop_EnUS","Suppressed Doppler Energy",16000,0,4000);
  GH1D* pAddSupDopUS = new GH1D("AddSupDop_EnUS","Suppressed Addback Doppler Energy",16000,0,4000);

  GH1D* pCoreEnergyUS_fep = new GH1D("Core_EnergyUS_fep","Tigress Core Energy FEP",16000,0,4000);
  GH1D* pCoreEnergyUS_pfep = new GH1D("Core_EnergyUS_pfep","Tigress Core Energy Projectile FEP",16000,0,4000);
  GH1D* pCoreEnergyUS_rfep = new GH1D("Core_EnergyUS_rfep","Tigress Core Energy Recoil FEP",16000,0,4000);
  GH1D* pCoreEnergyUS_nfep = new GH1D("Core_EnergyUS_nfep","Tigress Core Energy Not FEP",16000,0,4000);
  
  GH1D* pDopEnergyUS_fep = new GH1D("Dop_EnergyUS_fep","Doppler Energy FEP",16000,0,4000);
  GH1D* pDopEnergyUS_pfep = new GH1D("Dop_EnergyUS_pfep","Doppler Energy Projectile FEP",16000,0,4000);
  GH1D* pDopEnergyUS_rfep = new GH1D("Dop_EnergyUS_rfep","Doppler Energy Recoil FEP",16000,0,4000);
  GH1D* pDopEnergyUS_nfep = new GH1D("Dop_EnergyUS_nfep","Doppler Energy Not FEP",16000,0,4000);

  GH2D* pDopvPartUS = new GH2D("Dop_PartEn_US","Doppler Energy vs Particle Energy",3000,0,3000,1000,0,1000);

  GH2D* pThCorUS = new GH2D("Theta_CorrUS","Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pThCrtUS = new GH2D("Theta_CrctUS","Theta Correction",6000,0,3000,90,0,180);
  
  GH2D* pPhCorUS = new GH2D("Phi_CorrUS","Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pPhCrtUS = new GH2D("Phi_CrctUS","Phi Correction",6000,0,3000,32,0,thing1);

  GH1D* pReconEnergyUS = new GH1D("Recon_EnergyUS","Recon Energy",16000,0,4000);
  GH2D* pReconSumUS = new GH2D("Recon_SummaryUS","Recon Energy Summary",64,1,65,6000,0,3000);

  GH1D* pAddReconUS = new GH1D("AddRecon_EnUS","Addback Recon Energy",16000,0,4000);
  GH1D* pSupReconUS = new GH1D("SupRecon_EnUS","Suppressed Recon Energy",16000,0,4000);
  GH1D* pAddSupReconUS = new GH1D("AddSupRecon_EnUS","Suppressed Addback Recon Energy",16000,0,4000);

  GH1D* pReconEnergyUS_fep = new GH1D("Recon_EnergyUS_fep","Recon Energy FEP",16000,0,4000);
  GH1D* pReconEnergyUS_pfep = new GH1D("Recon_EnergyUS_pfep","Recon Energy Projectile FEP",16000,0,4000);
  GH1D* pReconEnergyUS_rfep = new GH1D("Recon_EnergyUS_rfep","Recon Energy Recoil FEP",16000,0,4000);
  GH1D* pReconEnergyUS_nfep = new GH1D("Recon_EnergyUS_nfep","Recon Energy Not FEP",16000,0,4000);

  GH2D* pReconvPartUS = new GH2D("ReconEn_v_partEn_US","Recon Energy vs Particle Energy",
				3000,0,3000,1000,0,1000);

  GH2D* pReconThCorUS = new GH2D("ReconTheta_CorrUS","Recon Theta Correlation",3000,0,3000,90,0,180);
  GH2D* pReconThCrtUS = new GH2D("ReconTheta_CrctUS","Recon Theta Correction",6000,0,3000,90,0,180);

  GH2D* pReconPhCorUS = new GH2D("ReconPhi_CorrUS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* pReconPhCrtUS = new GH2D("ReconPhi_CrctUS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Recoil
  GH2D* sPidRec = new GH2D("SecPID_Rec","DS SectorEn PID",24,1,25,1000,0,1000);
  GH2D* rPidRec = new GH2D("RingPID_Rec","DS RingEn PID",24,1,25,1000,0,1000);
  
  GH1D* rCoreEnergy = new GH1D("Core_EnergyRec","Tigress Core Energy",16000,0,4000);
  GH2D* rCoreSum = new GH2D("Core_SummaryRec","Core Energy Summary",64,1,65,3000,0,3000);
  
  GH1D* rDopEnergy = new GH1D("Dop_EnergyRec","Doppler Energy",16000,0,4000);
  GH2D* rDopSum = new GH2D("Dop_SummaryRec","Doppler Energy Summary",64,1,65,6000,0,3000);

  GH1D* rAddDopEn = new GH1D("AddDop_EnRec","Addback Doppler Energy",16000,0,4000);
  GH1D* rSupDopEn = new GH1D("SupDop_EnRec","Suppressed Doppler Energy",16000,0,4000);
  GH1D* rAddSupDopEn = new GH1D("AddSupDop_EnRec","Suppressed Addback Doppler Energy",16000,0,4000);
  
  GH1D* rCoreEnergy_fep = new GH1D("Core_EnergyRec_fep","Tigress Core Energy FEP",16000,0,4000);
  GH1D* rCoreEnergy_pfep = new GH1D("Core_EnergyRec_pfep","Tigress Core Energy Projectile FEP",16000,0,4000);
  GH1D* rCoreEnergy_rfep = new GH1D("Core_EnergyRec_rfep","Tigress Core Energy Recoil FEP",16000,0,4000);
  GH1D* rCoreEnergy_nfep = new GH1D("Core_EnergyRec_nfep","Tigress Core Energy Not FEP",16000,0,4000);

  GH1D* rDopEnergy_fep = new GH1D("Dop_EnergyRec_fep","Doppler Energy FEP",16000,0,4000);
  GH1D* rDopEnergy_pfep = new GH1D("Dop_EnergyRec_pfep","Doppler Energy Projectile FEP",16000,0,4000);
  GH1D* rDopEnergy_rfep = new GH1D("Dop_EnergyRec_rfep","Doppler Energy Recoil FEP",16000,0,4000);
  GH1D* rDopEnergy_nfep = new GH1D("Dop_EnergyRec_nfep","Doppler Energy Not FEP",16000,0,4000);

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

  GH1D* rReconEnergy = new GH1D("Recon_EnergyRec","Recon Energy",16000,0,4000);
  GH2D* rReconSum = new GH2D("Recon_SummaryRec","Recon Energy Summary",64,1,65,6000,0,3000);

  GH1D* rRecAddEn = new GH1D("AddRecon_EnRec","Addback Recon Energy",16000,0,4000);
  GH1D* rRecSupEn = new GH1D("SupRecon_EnRec","Suppressed Recon Energy",16000,0,4000);
  GH1D* rRecAddSupEn = new GH1D("AddSupRecon_EnRec","Suppressed Addback Recon Energy",16000,0,4000);
  
  GH1D* rReconEnergy_fep = new GH1D("Recon_EnergyRec_fep","Recon Energy FEP",16000,0,4000);
  GH1D* rReconEnergy_pfep = new GH1D("Recon_EnergyRec_pfep","Recon Energy Projectile FEP",16000,0,4000);
  GH1D* rReconEnergy_rfep = new GH1D("Recon_EnergyRec_rfep","Recon Energy Recoil FEP",16000,0,4000);
  GH1D* rReconEnergy_nfep = new GH1D("Recon_EnergyRec_nfep","Recon Energy Not FEP",16000,0,4000);

  GH2D* rReconvPart = new GH2D("ReconEn_v_PartEn_Rec","Recon Energy vs Particle Energy",3000,0,3000,1000,0,1000);
  //GH2D* rReconvPartNS2 = new GH2D("ReconEv_v_PartEn_NoS2_Rec","Recon Energy vs Particle Energy (No Sol2)",3000,0,3000,1000,0,1000);

  GH2D* rReconThCor = new GH2D("ReconTheta_CorrRec","Recon Theta Correlation",3000,0,3000,90,0,180);
  GH2D* rReconThCrt = new GH2D("ReconTheta_CrctRec","Recon Theta Correction",6000,0,3000,90,0,180);

  GH2D* rReconPhCor = new GH2D("ReconPhi_CorrRec","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  GH2D* rReconPhCrt = new GH2D("ReconPhi_CrctRec","Recon Phi Correction",6000,0,3000,32,0,thing1);
  
  //Tigress dets
  std::vector<GH1D*> pDetDopEnDS;
  std::vector<GH2D*> pDetDopSumDS;
  std::vector<GH2D*> pDetPhCorDS;
  std::vector<GH2D*> pDetPhCrtDS;
  std::vector<GH2D*> pDetThCorDS;
  std::vector<GH2D*> pDetThCrtDS;

  /*
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
  for(int i=0;i<64;i++) {

    detPosTP.push_back(new GH2D(Form("posPT_det%02d",i+1),Form("Det%02d Phi-Theta Surface",i+1),
				360,0,180,760,-10,370));
    
    
    //Downstream
    pDetDopEnDS.push_back(new GH1D(Form("Dop_EnergyDS_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				6000,0,3000));

    pDetDopSumDS.push_back(new GH2D(Form("Dop_SegSumDS_Det%02i",i+1),
                                    Form("Det%02i Doppler Seg Summary",i+1),
				    8,1,9,6000,0,3000));

    pDetPhCorDS.push_back(new GH2D(Form("PhiCorDS_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				3000,0,3000,32,0,thing1));

    pDetPhCrtDS.push_back(new GH2D(Form("PhiCrtDS_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				6000,0,3000,32,0,thing1));

    pDetThCorDS.push_back(new GH2D(Form("ThetaCorDS_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				3000,0,3000,90,0,180));

    pDetThCrtDS.push_back(new GH2D(Form("ThetaCrtDS_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				6000,0,3000,90,0,180));

    /*
    //Upstream
    pDetDopEnUS.push_back(new GH1D(Form("Dop_EnergyUS_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				   6000,0,3000));

    pDetDopSumUS.push_back(new GH2D(Form("Dop_SegSumUS_Det%02i",i+1),
                                    Form("Det%02i Doppler Seg Summary",i+1),
				    8,1,9,6000,0,3000));

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
				    8,1,9,6000,0,3000));

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

  //Tigress Segments
  std::vector<GH2D*> segPosTP02;
  std::vector<GH2D*> segPosTP10;
  
  for(int i=0;i<8;i++) {
    segPosTP02.push_back(new GH2D(Form("posPT_d02_s%02d",i+1),
				  Form("Det02 Seg%02d Phi-Theta Surface",i+1),
				  360,0,180,760,-10,370));

    segPosTP10.push_back(new GH2D(Form("posPT_d10_s%02d",i+1),
				  Form("Det10 Seg%02d Phi-Theta Surface",i+1),
				  360,0,180,760,-10,370));
  } 

  //S3 rings
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

  FillPositonVectors();
  const TVector3 incBeam = TVector3(0.0,0.0,1.0);

  TRandom* rand = new TRandom(50747227);
  const double r2d = TMath::RadToDeg();
  
  TSRIM srimB(beam_srim_name.c_str());
  TSRIM srimT(targ_srim_name.c_str());

  const double reac_en = 0.001*srimB.GetAdjustedEnergy(1000.0*beam_en,0.5*target_width,0.01);
  const double Sol2_En = KE_LAB(Theta_CM_FP(Theta_LAB_Max(reac_en),reac_en),reac_en);

  std::cout << "Correlating and histograming data..." << std::endl;
  
  Header header;
  RawData raw_data;
  while(fread(&header,header.bytes(),1,input_file)) {

    const int nE = header.evtNum;
    const int nS = header.nSdata;
    const int nT = header.nTdata;
  
    fread(&raw_data.sData,sizeof(S3Data),nS,input_file);
    fread(&raw_data.tData,sizeof(TigressData),nT,input_file);
    BuiltData data = BuildData(nE,nS,nT,raw_data);

    //S3 singles
    for(int i=0;i<data.nS3;i++) {

      int det = data.s3[i].det;
      int ring = data.s3[i].ring;
      int sec = data.s3[i].sector;
      int secr = secMap[sec-1] + 1;

      if(ring < Ri || ring > Rf)
	continue;

      double ring_en = data.s3[i].rEn;
      double sec_en = data.s3[i].sEn;

      TVector3 segPos = GetPos(det,ring,sec);
      TVector3 pos = data.s3[i].sPos;

      segPos.SetX(segPos.X() - beam_X);
      segPos.SetY(segPos.Y() - beam_Y);
      pos.SetX(pos.X() - beam_X);
      pos.SetY(pos.Y() - beam_Y);
      pos.SetZ(pos.Z() - beam_Z);
      
      if(!det) { //Upstream

	bSum->Fill(secr,sec_en);
	bSum->Fill(ring+32,ring_en);
	
	rPid0->Fill(ring,ring_en);
	sPid0->Fill(ring,sec_en);
	secD0->Fill(secr);

	if(data.s3[i].rP && data.s3[i].sP) { //Projectile
	  pSum->Fill(secr,sec_en);
	  pSum->Fill(ring+32,ring_en);

	  pThvPhUS->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	}
	if(data.s3[i].rR && data.s3[i].sR) {
	  rSum->Fill(secr,sec_en);
	  rSum->Fill(ring+32,ring_en);
	}
      }
      else { //Downstream
	
	rPid1->Fill(ring,ring_en);
	sPid1->Fill(ring,sec_en);
	secD1->Fill(secr);

	bSum->Fill(secr+64,sec_en);
	bSum->Fill(ring+96,ring_en);

	if(data.s3[i].rP && data.s3[i].sP) { //Projectile
	  sPid1_p->Fill(ring,sec_en);

	  pSum->Fill(secr+64,sec_en);
	  pSum->Fill(ring+96,ring_en);

	  //pPvPDS->Fill(segPos.Phi(),segPos.Perp());
	  
	  pThvPhDS->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	  if(ring%2) {
	    if(secr%2) {
	      pThvPhDS1->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	    else {
	      pThvPhDS2->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	  }
	  else {
	    if(secr%2) {
	      pThvPhDS2->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	    else {
	      pThvPhDS1->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	  }

	  pSvRDS->Fill(ring,secr);
	  pSecDS->Fill(secr);
	  //pPhiDS->Fill(segPos.Phi()*r2d);
	  
	  //pRingSecDS.at(ring-1)->Fill(secr);
	  pRingThvPhDS.at(ring-1)->Fill(pos.Theta()*r2d,pos.Phi()*r2d);

	  int sec_remap = secMap[sec-1];
	  
	  pSumRingDS.at(ring-1)->Fill(ring+31,ring_en);
	  pSumRingDS.at(ring-1)->Fill(sec_remap,sec_en);

	  pSumSecDS.at(sec_remap)->Fill(ring+31,ring_en);
	  pSumSecDS.at(sec_remap)->Fill(sec_remap,sec_en);
	  
	  /*
	  if(data.nS3 == 1) {
	    pSecDS_m1->Fill(sec);
	  }
	  else if(data.nS3 == 2) {
	    pSecDS_m2->Fill(sec);
	  }
	  */
	    
	}

	if(data.s3[i].rR && data.s3[i].sR) { //Recoil
	  sPid1_r->Fill(ring,sec_en);

	  rSum->Fill(secr+64,sec_en);
	  rSum->Fill(ring+96,ring_en);

	  rSec->Fill(secr);
	  //rPhi->Fill(segPos.Phi()*r2d);
	  
	  
	  //rPvP->Fill(segPos.Phi(),segPos.Perp());
	  rThvPh->Fill(pos.Theta()*r2d,pos.Phi()*r2d);

	  /*
	  rRingThvPh.at(ring-1)->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	  rRingSec.at(ring-1)->Fill(secr);

	  if(data.nS3 == 1) {
	    rSec_m1->Fill(secr);
	  }
	  else if(data.nS3 == 2) {
	    rSec_m2->Fill(secr);
	  }	  
	  */

	  int sec_remap = secMap[sec-1];
	  
	  rSumRing.at(ring-1)->Fill(ring+31,ring_en);
	  rSumRing.at(ring-1)->Fill(sec_remap,sec_en);
	  
	  rSumSec.at(sec_remap)->Fill(ring+31,ring_en);
	  rSumSec.at(sec_remap)->Fill(sec_remap,sec_en);
	  
	}
	
	if(data.nS3 == 1) {
	  rPid1m1->Fill(ring,ring_en);
	  sPid1m1->Fill(ring,sec_en);
	  
	}
	else if(data.nS3 == 2) {
	  rPid1m2->Fill(ring,ring_en);
	  sPid1m2->Fill(ring,sec_en);
	}

      }
      
    } //End S3 singles
    
    //Tigress singles
    for(int i=0;i<data.nTi;i++) {

      int det = data.tigress[i].det;
      double en = data.tigress[i].cEn;
      double core_en = rand->Gaus(en,Sigma(en));
      bool fep = data.tigress[i].fep;
      bool sup = data.tigress[i].sup;

      for(int j=0;j<data.tigress[i].nsegs;j++) {

	double exactX = data.tigress[i].x[j];
	double exactY = data.tigress[i].y[j];
	double exactZ = data.tigress[i].z[j];

	TVector3 exact_pos(exactX,exactY,exactZ);

	double phi = exact_pos.Phi();
	if(phi < 0) {
	  phi += TMath::TwoPi();
	}
	phi *= r2d;
	double theta = exact_pos.Theta()*r2d;
	
	double seg_en = data.tigress[i].sEn[j];
	int seg = data.tigress[i].segs[j];
	int num = seg + 8*(det-1);
	
	segEnergy->Fill(seg_en);
	segSum->Fill(num,seg_en);

	sPosTP->Fill(theta,phi);
	detPosTP.at(det-1)->Fill(theta,phi);

	if(det == 2)
	  segPosTP02.at(seg-1)->Fill(theta,phi);
	else if(det == 10)
	  segPosTP10.at(seg-1)->Fill(theta,phi);
	
      } //End loop over segments

      //if(BadSeg(data.tigress[i]))
      //continue;
      if(det == 6 || det == 11)
	continue;
      
      coreEnergy->Fill(core_en);
      coreSum->Fill(det,core_en);
      
      if(fep) {
	coreEn_Fep->Fill(core_en);
	FEPsum->Fill(det,core_en);
      }
      else {
	coreEn_NotFep->Fill(core_en);
      }

      if(!sup) {
	supEnergy->Fill(core_en);
	supSum->Fill(det,core_en);
      }
      
    } //End loop over tigress singles

    //Tigress Addback
    for(unsigned int i=0;i<data.tigressAB.size();i++) {

      int det = data.tigressAB[i].det;
      //int seg = data.tigressAB[i].MainSeg();
      bool sup = data.tigressAB[i].sup;
      
      double en = data.tigressAB[i].cEn;
      double coreEn = rand->Gaus(en,Sigma(en));
      
      addEnergy->Fill(coreEn);
      addSum->Fill(det,coreEn);

      if(!sup) {
	addSupEn->Fill(coreEn);
	addSupSum->Fill(det,coreEn);
      }
      
    } //end Addback Tigress loop
    
    
    //Mult2
    if(data.nTi == 2) {
      
      //if(!BadSeg(data.tigress[0]) && !BadSeg(data.tigress[1])) {

      int det1 = data.tigress[0].det;
      int det2 = data.tigress[1].det;
      double coreEn1 = data.tigress[0].cEn;
      double coreEn2 = data.tigress[1].cEn;
	
      double en1 = rand->Gaus(coreEn1,Sigma(coreEn1));
      double en2 = rand->Gaus(coreEn2,Sigma(coreEn2));

      if(det1 != 6 && det1 != 11 && det2 != 6 && det2 != 11) {
	if(en1 > 1330.0 && en1 < 1335.0)
	  gamGate->Fill(en2);
	else if(en2 > 1330.0 && en2 < 1335.0)
	  gamGate->Fill(en1);
      }
      //}

      if(data.tigress[0].fep && data.tigress[1].fep) {
		
	TVector3 posDet1 = GetPos(det1,0);
	TVector3 posDet2 = GetPos(det2,0);
	double thetaDet = posDet1.Angle(posDet2);
	
	TVector3 pos1 = data.tigress[0].MainPos();
	TVector3 pos2 = data.tigress[1].MainPos();
	double theta = pos1.Angle(pos2);

	TVector3 posEms1 = data.tigress[0].ExactMainPos();
	TVector3 posEms2 = data.tigress[1].ExactMainPos();
	double thetaEms = posEms1.Angle(posEms2);

	gammaAngDetFEP->Fill(TMath::Cos(thetaDet));
	gammaAngFEP->Fill(TMath::Cos(theta));
	gamAngExMsFEP->Fill(TMath::Cos(thetaEms));

	if(data.tigress[0].nsegs == 1 && data.tigress[1].nsegs == 1) {

	  TVector3 posE1 = data.tigress[0].ExactMainPos();
	  TVector3 posE2 = data.tigress[1].ExactMainPos();
	  double thetaE = posE1.Angle(posE2);
	  
	  gamAngExFEP->Fill(TMath::Cos(thetaE));
	}
	
      }
	
    } //End Tigress mult2
    
    //Coincidences
    if(data.nS3 > 0 && data.nTi > 0) {

      bool good =false;
      for(int i=0;i<data.nS3;i++) {

	int bDet = data.s3[i].det;
	int ring = data.s3[i].ring;
	int sector = data.s3[i].sector;
	//int sectorr = secMap[sector-1] + 1;

	if(ring < Ri || ring > Rf)
	  continue;

	double ring_en = data.s3[i].rEn;
        double sec_en = data.s3[i].sEn;

	TVector3 bPos = GetPos(bDet,ring,sector);
	bPos.SetX(bPos.X() - beam_X);
	bPos.SetY(bPos.Y() - beam_Y);
	
	if(bDet && data.s3[i].rP && data.s3[i].sP) { //Projectile DS gate

	  sPidDS->Fill(ring,sec_en);
	  rPidDS->Fill(ring,ring_en);

	  bool sol2 = false;
	  if(sec_en < Sol2_En)
	    sol2 = true;
	  
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
	  
	  for(int j=0;j<data.nTi;j++) {

	    //if(BadSeg(data.tigress[j]))
	    //continue;
	    	    
	    int det = data.tigress[j].det;
	    if(det == 6 || det == 11)
	      continue;

	    int seg = data.tigress[j].MainSeg();
	    
	    double en = data.tigress[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.tigress[j].fep;
	    bool PFEP = data.tigress[j].pfep;
	    bool SUP = data.tigress[j].sup;
	  
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);
	    
	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	    if(dopEn > 1005.0 && dopEn < 1025.0)
	      good = true;
	    
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0)
	      reac_phi += TMath::TwoPi();

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0)
	      det_phi += TMath::TwoPi();

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0)
	      planeAng += TMath::TwoPi();

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0)
	      recon_phi += TMath::TwoPi();

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0)
	      reconAng += TMath::TwoPi();
	    
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

	    if(!SUP) {
	      pSupDopDS->Fill(dopEn);
	      pSupReconEnDS->Fill(recon_en);
	    }

	    pCoreEnergyDS->Fill(coreEn);
	    pCoreSumDS->Fill(det,coreEn);
	    
	    pDopEnergyDS->Fill(dopEn);
	    pDopSumDS->Fill(det,dopEn);
	      
	    pDopvPartDS->Fill(dopEn,sec_en);

	    pCEnvRingDS->Fill(ring,coreEn);
	    pDopvRingDS->Fill(ring,dopEn);
	    pRecvRingDS->Fill(ring,recon_en);
	    
	    pThCorDS->Fill(coreEn,theta*r2d);
	    pThCrtDS->Fill(dopEn,theta*r2d);

	    pcThCrtDS->Fill(dopEn,TMath::Cos(theta));

	    pPhCorDS->Fill(coreEn,planeAng*r2d);
	    pPhCrtDS->Fill(dopEn,planeAng*r2d);

	    
	    pDetDopEnDS.at(det-1)->Fill(dopEn);
	    pDetDopSumDS.at(det-1)->Fill(seg,dopEn);

	    pDetThCorDS.at(det-1)->Fill(coreEn,theta*r2d);
	    pDetThCrtDS.at(det-1)->Fill(dopEn,theta*r2d);

	    pDetPhCorDS.at(det-1)->Fill(coreEn,planeAng*r2d);
	    pDetPhCrtDS.at(det-1)->Fill(dopEn,planeAng*r2d);
	    
	    
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
	    
	  } //End Tigress loop

	  for(unsigned int j=0;j<data.tigressAB.size();j++) {

	    int det = data.tigressAB[j].det;
	    int seg = data.tigressAB[j].MainSeg();
	    
	    double en = data.tigressAB[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool SUP = data.tigressAB[j].sup;
	 
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);
	    
	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;
	    
	    pAddDS->Fill(coreEn);
	    pAddSumDS->Fill(det,coreEn);

	    pAddDopDS->Fill(dopEn);
	    pAddDopSumDS->Fill(det,dopEn);

	    pAddReconEnDS->Fill(recon_en);

	    if(!SUP) {
	      pAddSupDopDS->Fill(dopEn);
	      pAddSupReconEnDS->Fill(recon_en);
	    }
	    
	  } //end Addback Tigress loop
	  
	} //End DS projectile gate

	else if(!bDet && data.s3[i].rP && data.s3[i].sP) { //Projectile US gate

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
	  
	  for(int j=0;j<data.nTi;j++) {

	    //if(BadSeg(data.tigress[j]))
	    //continue;
	    
	    int det = data.tigress[j].det;
	    if(det == 6 || det == 11)
	      continue;
	    
	    int seg = data.tigress[j].MainSeg();
	    
	    double en = data.tigress[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.tigress[j].fep;
	    bool PFEP = data.tigress[j].pfep;
	    bool SUP = data.tigress[j].sup;
	  
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	  
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0)
	      reac_phi += TMath::TwoPi();

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0)
	      det_phi += TMath::TwoPi();

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0)
	      planeAng += TMath::TwoPi();

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    TVector3 reconPlane = rPos.Cross(incBeam);
	    
	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0)
	      recon_phi += TMath::TwoPi();

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0)
	      reconAng += TMath::TwoPi();

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

	    if(!SUP) {
	      pSupDopUS->Fill(dopEn);
	      pSupReconUS->Fill(recon_en);
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
	  
	  } //End Tigress loop

	  for(unsigned int j=0;j<data.tigressAB.size();j++) {

	    int det = data.tigressAB[j].det;
	    int seg = data.tigressAB[j].MainSeg();
	    
	    double en = data.tigressAB[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    //bool FEP = data.tigressAB[j].fep;
	    //bool PFEP = data.tigressAB[j].pfep;
	    bool SUP = data.tigressAB[j].sup;
	 
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);
	    
	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    pAddDopUS->Fill(dopEn);
	    pAddReconUS->Fill(recon_en);
	    if(!SUP) {
	      pAddSupDopUS->Fill(dopEn);
	      pAddSupReconUS->Fill(recon_en);
	    }
	    
	  } //end Addback Tigress loop
	} ////End US projectile gate

	else if(data.s3[i].rR && data.s3[i].sR) { //Recoil gate

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

	  for(int j=0;j<data.nTi;j++) {

	    //if(BadSeg(data.tigress[j]))
	    //continue;

	    int det = data.tigress[j].det;
	    if(det == 6 || det == 11)
	      continue;
	    
	    int seg = data.tigress[j].MainSeg();
	    
	    double en = data.tigress[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.tigress[j].fep;
	    bool PFEP = data.tigress[j].pfep;
	    bool SUP = data.tigress[j].sup;
	  
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	  
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0)
	      reac_phi += TMath::TwoPi();

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0)
	      det_phi += TMath::TwoPi();

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0)
	      planeAng += TMath::TwoPi();

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;
	    if(recon_en > 1005.0 && recon_en < 1025.0)
	      good = true;

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0)
	      recon_phi += TMath::TwoPi();

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0)
	      reconAng += TMath::TwoPi();
	    
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

	    if(!SUP) {
	      rSupDopEn->Fill(dopEn);
	      rRecSupEn->Fill(recon_en);
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
	    
	  } //End Tigress loop

	  for(unsigned int j=0;j<data.tigressAB.size();j++) {

	    int det = data.tigressAB[j].det;
	    int seg = data.tigressAB[j].MainSeg();
	    
	    double en = data.tigressAB[j].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    //bool FEP = data.tigressAB[j].fep;
	    //bool PFEP = data.tigressAB[j].pfep;
	    bool SUP = data.tigressAB[j].sup;
	 
	    TVector3 sPos = GetPos(det,seg);
	    sPos.SetX(sPos.X() - beam_X);
	    sPos.SetY(sPos.Y() - beam_Y);
	    
	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    rAddDopEn->Fill(dopEn);
	    rRecAddEn->Fill(recon_en);

	    if(!SUP) {
	      rAddSupDopEn->Fill(dopEn);
	      rRecAddSupEn->Fill(recon_en);
	    }
	    
	  } //end Addback Tigress loop
	  
	} //End recoil gate
	
      } //End S3 loop

      if(good) {
	for(int i=0;i<data.nS3;i++) {

	int bDet = data.s3[i].det;
	int ring = data.s3[i].ring;

	if(ring < Ri || ring > Rf)
	  continue;

        double sec_en = data.s3[i].sEn;
	if(bDet) {

	  sPidDS_gam->Fill(ring,sec_en);
	  if(data.s3[i].rP && data.s3[i].sP)
	    psPidDS_gam->Fill(ring,sec_en);

	  if(bDet && data.s3[i].rR && data.s3[i].sR)
	    rsPidDS_gam->Fill(ring,sec_en);
	}
	
	} //End S3 loop
      } // if(good)
      
    } //End coincidences

    
  } //End while loop
  fclose(input_file);

  std::cout << "Writing histograms to file..." << std::endl;
  
  TFile* outFile = new TFile(output_filename,"RECREATE");
  
  outFile->mkdir("Tigress");
  outFile->mkdir("Tigress/Dets");
  outFile->mkdir("Tigress/Dets/Det02");
  outFile->mkdir("Tigress/Dets/Det10");
  
  outFile->mkdir("S3");
  outFile->mkdir("S3/Rings");
  outFile->mkdir("S3/Sectors");
  outFile->mkdir("S3/Rings2D");

  outFile->mkdir("Coincidence");
  outFile->mkdir("Coincidence/ProjectileDS");
  outFile->mkdir("Coincidence/ProjectileDS/Rings");
  outFile->mkdir("Coincidence/ProjectileDS/Doppler");
  outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings");
  //outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings/PhiCorr");
  //outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings/ThetaCorr");

  outFile->mkdir("Coincidence/ProjectileDS/Recon");
  outFile->mkdir("Coincidence/ProjectileDS/Recon/Rings");
  
  
  outFile->mkdir("Coincidence/ProjectileDS/TigressDets");
  outFile->mkdir("Coincidence/ProjectileDS/TigressDets/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileDS/TigressDets/ThetaCorr");
  outFile->mkdir("Coincidence/ProjectileDS/TigressDets/Summaries");
  

  outFile->mkdir("Coincidence/ProjectileUS");
  outFile->mkdir("Coincidence/ProjectileUS/Doppler");
  outFile->mkdir("Coincidence/ProjectileUS/Recon");
  //outFile->mkdir("Coincidence/ProjectileUS/Doppler/Rings");
  
  /*
  outFile->mkdir("Coincidence/ProjectileUS/S3Rings");
  outFile->mkdir("Coincidence/ProjectileUS/S3Rings/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileUS/S3Rings/ThetaCorr");
  
  outFile->mkdir("Coincidence/ProjectileUS/TigressDets");
  outFile->mkdir("Coincidence/ProjectileUS/TigressDets/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileUS/TigressDets/ThetaCorr");
  outFile->mkdir("Coincidence/ProjectileUS/TigressDets/Summaries");
  */
  
  outFile->mkdir("Coincidence/Recoil");
  outFile->mkdir("Coincidence/Recoil/Rings");
  outFile->mkdir("Coincidence/Recoil/Doppler");
  outFile->mkdir("Coincidence/Recoil/Doppler/Rings");
  outFile->mkdir("Coincidence/Recoil/Recon");
  outFile->mkdir("Coincidence/Recoil/Recon/Rings");
  
  /*
  outFile->mkdir("Coincidence/Recoil/S3Rings");
  outFile->mkdir("Coincidence/Recoil/S3Rings/PhiCorr");
  outFile->mkdir("Coincidence/Recoil/S3Rings/ThetaCorr");

  outFile->mkdir("Coincidence/Recoil/TigressDets");
  outFile->mkdir("Coincidence/Recoil/TigressDets/PhiCorr");
  outFile->mkdir("Coincidence/Recoil/TigressDets/ThetaCorr");
  outFile->mkdir("Coincidence/Recoil/TigressDets/Summaries");
  */

  outFile->cd("S3");

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

  outFile->cd("S3/Rings");
  for(auto* h : pSumRingDS)
    h->Write();

  for(auto* h : rSumRing)
    h->Write();

  outFile->cd("S3/Sectors");
  for(auto* h : pSumSecDS)
    h->Write();

  for(auto* h : rSumSec)
    h->Write();

  /*
  outFile->cd("S3/Rings");
  for(int i=0;i<24;i++) {
    pRingSecDS.at(i)->Write();
    rRingSec.at(i)->Write();
  }
  */

  outFile->cd("S3/Rings2D");
  for(int i=0;i<24;i++) {
    pRingThvPhDS.at(i)->Write();
    //rRingThvPh.at(i)->Write();
  }
  
  outFile->cd("Tigress");

  coreEnergy->Write();
  coreSum->Write();
  
  segEnergy->Write();
  segSum->Write();

  addEnergy->Write();
  addSum->Write();

  supEnergy->Write();
  supSum->Write();

  addSupEn->Write();
  addSupSum->Write();
  
  coreEn_Fep->Write();
  coreEn_NotFep->Write();
  FEPsum->Write();

  gamGate->Write();

  gammaAngDetFEP->Write();
  gammaAngFEP->Write();
  gamAngExMsFEP->Write();
  gamAngExFEP->Write();

  sPosTP->Write();

  outFile->cd("Tigress/Dets");
  for(int i=0;i<64;i++) {
    detPosTP.at(i)->Write();
  }

  outFile->cd("Tigress/Dets/Det02");
  for(int i=0;i<8;i++) {
    segPosTP02.at(i)->Write();
  }

  outFile->cd("Tigress/Dets/Det10");
  for(int i=0;i<8;i++) {
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

  pAddDS->Write();
  pAddSumDS->Write();

  pCoreEnergyDS_fep->Write();
  pCoreEnergyDS_pfep->Write();
  pCoreEnergyDS_rfep->Write();
  pCoreEnergyDS_nfep->Write();

  pCEnvRingDS->Write();
  
  outFile->cd("Coincidence/ProjectileDS/Doppler");
  
  pDopEnergyDS->Write();
  pDopSumDS->Write();

  pAddDopDS->Write();
  pAddDopSumDS->Write();

  pSupDopDS->Write();
  pAddSupDopDS->Write();

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

  pSupReconEnDS->Write();
  pAddReconEnDS->Write();
  pAddSupReconEnDS->Write();

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
  
  
  outFile->cd("Coincidence/ProjectileDS/TigressDets");
  for(int i=0;i<64;i++) {
    pDetDopEnDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/TigressDets/Summaries");
  for(int i=0;i<64;i++) {
    pDetDopSumDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/TigressDets/PhiCorr");
  for(int i=0;i<64;i++) {
    pDetPhCorDS.at(i)->Write();
    pDetPhCrtDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/TigressDets/ThetaCorr");
  for(int i=0;i<64;i++) {
    pDetThCorDS.at(i)->Write();
    pDetThCrtDS.at(i)->Write();
  }
  

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

  pSupDopUS->Write();
  pAddDopUS->Write();
  pAddSupDopUS->Write();
  
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

  pSupReconUS->Write();
  pAddReconUS->Write();
  pAddSupReconUS->Write();
  
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

  outFile->cd("Coincidence/ProjectileUS/S3Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    pRingPhCorUS.at(i)->Write();
    pRingPhCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/S3Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    pRingThCorUS.at(i)->Write();
    pRingThCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/TigressDets");
  for(int i=0;i<64;i++) {
    pDetDopEnUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/TigressDets/Summaries");
  for(int i=0;i<64;i++) {
    pDetDopSumUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/TigressDets/PhiCorr");
  for(int i=0;i<64;i++) {
    pDetPhCorUS.at(i)->Write();
    pDetPhCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/TigressDets/ThetaCorr");
  for(int i=0;i<64;i++) {
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

  rSupDopEn->Write();
  rAddDopEn->Write();
  rAddSupDopEn->Write();
  
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

  rRecSupEn->Write();
  rRecAddEn->Write();
  rRecAddSupEn->Write();
  
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
  outFile->cd("Coincidence/Recoil/S3Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    rRingPhCor.at(i)->Write();
    rRingPhCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/S3Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    rRingThCor.at(i)->Write();
    rRingThCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/TigressDets");
  for(int i=0;i<64;i++) {
    rDetDopEn.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/TigressDets/Summaries");
  for(int i=0;i<64;i++) {
    rDetDopSum.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/TigressDets/PhiCorr");
  for(int i=0;i<64;i++) {
    rDetPhCor.at(i)->Write();
    rDetPhCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/TigressDets/ThetaCorr");
  for(int i=0;i<64;i++) {
    rDetThCor.at(i)->Write();
    rDetThCrt.at(i)->Write();
  }
  */

  outFile->Close();

  std::cout << "Done!" << std::endl;

  return 0;
}
