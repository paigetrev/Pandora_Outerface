/**
 *  @file   LArRecoND/include/PandoraOuterface.h
 *
 *  @brief  Header file for PandoraOuterface.
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_OUTERFACE_H
#define PANDORA_ND_OUTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"

#ifdef USE_EDEPSIM
#include "TG4Event.h"
#endif

#include "TGeoManager.h"
#include "TGeoNode.h"

#include "LArGrid.h"
#include "LArHitInfo.h"
#include "LArRecoNDFormat.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "TFile.h"
#include "TTree.h"

namespace pandora
{
class Pandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_postreco
{

struct ParameterStruct
{
  bool runTrackFit    = false;
  bool runShowerFit   = false;
  float pixelPitch    = 0.4;
  float trackScoreCut = 0.5;

  bool applyThreshold = false;
  float thresholdVal  = 0.; // ke-

  bool voxelizeZ = false;
  float voxelZHW = 0.6; // cm
  bool useVoxelizedStartStop = false;

  int verbosity = 0;

  std::string xmlName = "";
  std::string fileName = "";
  std::string outfileName = "LArRecoND_outerface_test.root";
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Read in the file and process the post-reconstruction
 *
 *  @param  parameters the input parameters controlling aspects of post-reco
 *
 */
void ProcessPostReco(const ParameterStruct &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 *
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], ParameterStruct &parameters);

/**
 *  @brief  Parse XML file to set parameters
 *
 *  @param  parameters to take the XML file name and set the parameters
 *
 *  @return success
 */
bool ReadSettings(ParameterStruct &parameters);

/**
 *  @brief  Print the list of configurable options
 *
 *  @return false, to force abort
 */
bool PrintOptions();

/**
 * @brief Class to handle the ND Reco Output Data Model
 *
 */
class NDRecoOutputData
{
 public:
  //NDRecoOutputData(const std::string filename, const bool writeTracks); ///< default constructor
  NDRecoOutputData(const std::string filename, const bool writeTracks, const bool writeShowers);/// < added to include showers
  
  void ClearData(); ///< will reset the vectors

  void WriteToFile(); ///< write to the tree

  void CloseFile(); ///< close the TFile

  void FillBasicBranches( const std::unique_ptr<LArRecoNDFormat> &inputSpill ); ///< Fill the basic branches
  void FillTrackBranches( const std::vector<float> &startX, const std::vector<float> &startY, const std::vector<float> &startZ,
			  const std::vector<float> &dirX,  const std::vector<float> &dirY,  const std::vector<float> &dirZ,
			  const std::vector<float> &endX, const std::vector<float> &endY, const std::vector<float> &endZ,
                          const std::vector<float> &enddirX,  const std::vector<float> &enddirY,  const std::vector<float> &enddirZ,
			  const std::vector<float> &length ); ///< Fill the track fit result branches
//added line
  void FillShowerBranches( const std::vector<float> &shwrcentX, const std::vector<float> &shwrcentY,  const std::vector<float> &shwrcentZ, 
                           const std::vector<float> &shwrstartX,  const std::vector<float> &shwrstartY, const std::vector<float> &shwrstartZ,
                           const std::vector<float> &shrdirX, const std::vector<float> &shwrdirY, const std::vector<float> &shwrdirZ,
			   const std::vector<float> &shwrlength, const std::vector<int> &shwrSlice, const std::vector<int> &shwrCluster); ///< Fill the shower fit result branches
  
  void FillTrackCaloBranches( const std::vector<int> &tfSliceId, const std::vector<int> &tfPfoId,
			      const std::vector<float> &tfX, const std::vector<float> &tfY, const std::vector<float> &tfZ, const std::vector<float> &tfQ,
			      const std::vector<float> &tfRR, const std::vector<float> &tfdx, const std::vector<float> &tfdQdx );

 private:
  TFile* m_fileOut;
  TTree* m_treeOut;

  Int_t m_out_event;
  Int_t m_out_subrun;
  Int_t m_out_run;
  Int_t m_out_unixTime;
  Int_t m_out_start_t;
  Int_t m_out_end_t;
  Int_t m_out_trigger;
  std::vector<int> m_out_sliceID;
  std::vector<int> m_out_clusterID;
  std::vector<float> m_out_nuVtxX;
  std::vector<float> m_out_nuVtxY;
  std::vector<float> m_out_nuVtxZ;
  std::vector<int> m_out_n3DHits;
  std::vector<int> m_out_nUHits;
  std::vector<int> m_out_nVHits;
  std::vector<int> m_out_nWHits;
  std::vector<int> m_out_isShower;
  std::vector<float> m_out_trackScore;
  std::vector<int> m_out_recoPDG;
  std::vector<int> m_out_isRecoPrimary;
  std::vector<float> m_out_startX;
  std::vector<float> m_out_startY;
  std::vector<float> m_out_startZ;
  std::vector<float> m_out_endX;
  std::vector<float> m_out_endY;
  std::vector<float> m_out_endZ;
  std::vector<float> m_out_dirX;
  std::vector<float> m_out_dirY;
  std::vector<float> m_out_dirZ;
  std::vector<float> m_out_centroidX;
  std::vector<float> m_out_centroidY;
  std::vector<float> m_out_centroidZ;
  std::vector<float> m_out_length1;
  std::vector<float> m_out_length2;
  std::vector<float> m_out_length3;
  std::vector<float> m_out_energy;
  std::vector<int> m_out_recoHitId;
  std::vector<int> m_out_recoHitSliceId;
  std::vector<int> m_out_recoHitClusterId;
  std::vector<float> m_out_recoHitX;
  std::vector<float> m_out_recoHitY;
  std::vector<float> m_out_recoHitZ;
  std::vector<float> m_out_recoHitE;
  std::vector<int> m_out_gotMatch;
  std::vector<int> m_out_mcPDG;
  std::vector<long> m_out_mcId;
  std::vector<long> m_out_mcLocalId;
  std::vector<int> m_out_isPrimary;
  std::vector<int> m_out_nSharedHits;
  std::vector<float> m_out_completeness;
  std::vector<float> m_out_purity;
  std::vector<float> m_out_mcEnergy;
  std::vector<float> m_out_mcPx;
  std::vector<float> m_out_mcPy;
  std::vector<float> m_out_mcPz;
  std::vector<float> m_out_mcVtxX;
  std::vector<float> m_out_mcVtxY;
  std::vector<float> m_out_mcVtxZ;
  std::vector<float> m_out_mcEndX;
  std::vector<float> m_out_mcEndY;
  std::vector<float> m_out_mcEndZ;
  std::vector<int> m_out_mcNuPDG;
  std::vector<long> m_out_mcNuId;
  std::vector<int> m_out_mcNuCode;
  std::vector<float> m_out_mcNuVtxX;
  std::vector<float> m_out_mcNuVtxY;
  std::vector<float> m_out_mcNuVtxZ;
  std::vector<float> m_out_mcNuE;
  std::vector<float> m_out_mcNuPx;
  std::vector<float> m_out_mcNuPy;
  std::vector<float> m_out_mcNuPz;
  // New data products
  std::vector<float> m_out_trkfitStartX;
  std::vector<float> m_out_trkfitStartY;
  std::vector<float> m_out_trkfitStartZ;
  std::vector<float> m_out_trkfitStartDirX;
  std::vector<float> m_out_trkfitStartDirY;
  std::vector<float> m_out_trkfitStartDirZ;
  std::vector<float> m_out_trkfitEndX;
  std::vector<float> m_out_trkfitEndY;
  std::vector<float> m_out_trkfitEndZ;
  std::vector<float> m_out_trkfitEndDirX;
  std::vector<float> m_out_trkfitEndDirY;
  std::vector<float> m_out_trkfitEndDirZ;
  std::vector<float> m_out_trkfitLength;
  std::vector<int> m_out_trkfitSliceId;
  std::vector<int> m_out_trkfitPfoId;
  std::vector<float> m_out_trkfitX;
  std::vector<float> m_out_trkfitY;
  std::vector<float> m_out_trkfitZ;
  std::vector<float> m_out_trkfitQ;
  std::vector<float> m_out_trkfitRR;
  std::vector<float> m_out_trkfitdx;
  std::vector<float> m_out_trkfitdQdx;
  // added shower products
  std::vector<float> m_out_shwrfitLength;
  std::vector<float> m_out_shwrfitCentroidX;
  std::vector<float> m_out_shwrfitCentroidY;
  std::vector<float> m_out_shwrfitCentroidZ;
  std::vector<float> m_out_shwrfitStartX;
  std::vector<float> m_out_shwrfitStartY;
  std::vector<float> m_out_shwrfitStartZ;
  std::vector<float> m_out_shwrfitDirX;
  std::vector<float> m_out_shwrfitDirY;
  std::vector<float> m_out_shwrfitDirZ;
  std::vector<int> m_out_shwrSliceId;
  std::vector<int> m_out_shwrClusterId;
};

 NDRecoOutputData::NDRecoOutputData(const std::string filename, const bool writeTracks, const bool writeShowers)
 {

   m_fileOut = new TFile(filename.c_str(), "RECREATE");
   m_treeOut = new TTree("LArRecoND","LArRecoND");

   // Set the branches
   m_treeOut->Branch("event", &m_out_event);
   m_treeOut->Branch("subRun", &m_out_subrun);
   m_treeOut->Branch("run", &m_out_run);
   m_treeOut->Branch("unixTime", &m_out_unixTime);
   m_treeOut->Branch("startTime", &m_out_start_t);
   m_treeOut->Branch("endTime", &m_out_end_t);
   m_treeOut->Branch("triggers", &m_out_trigger);
   m_treeOut->Branch("sliceId", &m_out_sliceID);
   m_treeOut->Branch("clusterId", &m_out_clusterID);
   m_treeOut->Branch("nuVtxX", &m_out_nuVtxX);
   m_treeOut->Branch("nuVtxY", &m_out_nuVtxY);
   m_treeOut->Branch("nuVtxZ", &m_out_nuVtxZ);
   m_treeOut->Branch("n3DHits", &m_out_n3DHits);
   m_treeOut->Branch("nUHits", &m_out_nUHits);
   m_treeOut->Branch("nVHits", &m_out_nVHits);
   m_treeOut->Branch("nWHits", &m_out_nWHits);
   m_treeOut->Branch("isShower", &m_out_isShower);
   m_treeOut->Branch("trackScore", &m_out_trackScore);
   m_treeOut->Branch("recoPDG", &m_out_recoPDG);
   m_treeOut->Branch("isRecoPrimary", &m_out_isRecoPrimary);
   m_treeOut->Branch("startX", &m_out_startX);
   m_treeOut->Branch("startY", &m_out_startY);
   m_treeOut->Branch("startZ", &m_out_startZ);
   m_treeOut->Branch("endX", &m_out_endX);
   m_treeOut->Branch("endY", &m_out_endY);
   m_treeOut->Branch("endZ", &m_out_endZ);
   m_treeOut->Branch("dirX", &m_out_dirX);
   m_treeOut->Branch("dirY", &m_out_dirY);
   m_treeOut->Branch("dirZ", &m_out_dirZ);
   m_treeOut->Branch("centroidX", &m_out_centroidX);
   m_treeOut->Branch("centroidY", &m_out_centroidY);
   m_treeOut->Branch("centroidZ", &m_out_centroidZ);
   m_treeOut->Branch("length1", &m_out_length1);
   m_treeOut->Branch("length2", &m_out_length2);
   m_treeOut->Branch("length3", &m_out_length3);
   m_treeOut->Branch("energy", &m_out_energy);
   m_treeOut->Branch("recoHitId", &m_out_recoHitId);
   m_treeOut->Branch("recoHitSliceId", &m_out_recoHitSliceId);
   m_treeOut->Branch("recoHitClusterId", &m_out_recoHitClusterId);
   m_treeOut->Branch("recoHitX", &m_out_recoHitX);
   m_treeOut->Branch("recoHitY", &m_out_recoHitY);
   m_treeOut->Branch("recoHitZ", &m_out_recoHitZ);
   m_treeOut->Branch("recoHitE", &m_out_recoHitE);
   m_treeOut->Branch("gotMatch", &m_out_gotMatch);
   m_treeOut->Branch("mcPDG", &m_out_mcPDG);
   m_treeOut->Branch("mcId", &m_out_mcId);
   m_treeOut->Branch("mcLocalId", &m_out_mcLocalId);
   m_treeOut->Branch("isPrimary", &m_out_isPrimary);
   m_treeOut->Branch("nSharedHits", &m_out_nSharedHits);
   m_treeOut->Branch("completeness", &m_out_completeness);
   m_treeOut->Branch("purity", &m_out_purity);
   m_treeOut->Branch("mcEnergy", &m_out_mcEnergy);
   m_treeOut->Branch("mcPx", &m_out_mcPx);
   m_treeOut->Branch("mcPy", &m_out_mcPy);
   m_treeOut->Branch("mcPz", &m_out_mcPz);
   m_treeOut->Branch("mcVtxX", &m_out_mcVtxX);
   m_treeOut->Branch("mcVtxY", &m_out_mcVtxY);
   m_treeOut->Branch("mcVtxZ", &m_out_mcVtxZ);
   m_treeOut->Branch("mcEndX", &m_out_mcEndX);
   m_treeOut->Branch("mcEndY", &m_out_mcEndY);
   m_treeOut->Branch("mcEndZ", &m_out_mcEndZ);
   m_treeOut->Branch("mcNuPDG", &m_out_mcNuPDG);
   m_treeOut->Branch("mcNuId", &m_out_mcNuId);
   m_treeOut->Branch("mcNuCode", &m_out_mcNuCode);
   m_treeOut->Branch("mcNuVtxX", &m_out_mcNuVtxX);
   m_treeOut->Branch("mcNuVtxY", &m_out_mcNuVtxY);
   m_treeOut->Branch("mcNuVtxZ", &m_out_mcNuVtxZ);
   m_treeOut->Branch("mcNuE", &m_out_mcNuE);
   m_treeOut->Branch("mcNuPx", &m_out_mcNuPx);
   m_treeOut->Branch("mcNuPy", &m_out_mcNuPy);
   m_treeOut->Branch("mcNuPz", &m_out_mcNuPz);

   // And the new ones
   if ( writeTracks ) {
     m_treeOut->Branch("trkfitStartX", &m_out_trkfitStartX);
     m_treeOut->Branch("trkfitStartY", &m_out_trkfitStartY);
     m_treeOut->Branch("trkfitStartZ", &m_out_trkfitStartZ);
     m_treeOut->Branch("trkfitStartDirX", &m_out_trkfitStartDirX);
     m_treeOut->Branch("trkfitStartDirY", &m_out_trkfitStartDirY);
     m_treeOut->Branch("trkfitStartDirZ", &m_out_trkfitStartDirZ);
     m_treeOut->Branch("trkfitEndX", &m_out_trkfitEndX);
     m_treeOut->Branch("trkfitEndY", &m_out_trkfitEndY);
     m_treeOut->Branch("trkfitEndZ", &m_out_trkfitEndZ);
     m_treeOut->Branch("trkfitEndDirX", &m_out_trkfitEndDirX);
     m_treeOut->Branch("trkfitEndDirY", &m_out_trkfitEndDirY);
     m_treeOut->Branch("trkfitEndDirZ", &m_out_trkfitEndDirZ);
     m_treeOut->Branch("trkfitLength", &m_out_trkfitLength);
     m_treeOut->Branch("trkfitSliceId", &m_out_trkfitSliceId);
     m_treeOut->Branch("trkfitPfoId", &m_out_trkfitPfoId);
     m_treeOut->Branch("trkfitX", &m_out_trkfitX);
     m_treeOut->Branch("trkfitY", &m_out_trkfitY);
     m_treeOut->Branch("trkfitZ", &m_out_trkfitZ);
     m_treeOut->Branch("trkfitQ", &m_out_trkfitQ);
     m_treeOut->Branch("trkfitRR", &m_out_trkfitRR);
     m_treeOut->Branch("trkfitdx", &m_out_trkfitdx);
     m_treeOut->Branch("trkfitdQdx", &m_out_trkfitdQdx);
  }
   if( writeShowers){
     m_treeOut->Branch("shwrfitLength", &m_out_shwrfitLength);
     m_treeOut->Branch("shwrfitCentroidX", &m_out_shwrfitCentroidX);
     m_treeOut->Branch("shwrfitCentroidY", &m_out_shwrfitCentroidY);
     m_treeOut->Branch("shwrfitCentroidZ", &m_out_shwrfitCentroidZ);
     m_treeOut->Branch("shwrfitStartX", &m_out_shwrfitStartX);
     m_treeOut->Branch("shwrfitStartY", &m_out_shwrfitStartY);
     m_treeOut->Branch("shwrfitStartZ", &m_out_shwrfitStartZ);
     m_treeOut->Branch("shwrfitDirX", &m_out_shwrfitDirX);
     m_treeOut->Branch("shwrfitDirY", &m_out_shwrfitDirY);
     m_treeOut->Branch("shwrfitDirZ", &m_out_shwrfitDirZ);
     m_treeOut->Branch("shwrSliceId", &m_out_shwrSliceId);
     m_treeOut->Branch("shwrClusterId", &m_out_shwrClusterId);
   
   }
  
 }

 void NDRecoOutputData::ClearData()
 {
   m_out_event = 0;
   m_out_subrun = 0;
   m_out_run = 0;
   m_out_unixTime = 0;
   m_out_start_t = 0;
   m_out_end_t = 0;
   m_out_trigger = 0;
   m_out_sliceID.clear();
   m_out_clusterID.clear();
   m_out_nuVtxX.clear();
   m_out_nuVtxY.clear();
   m_out_nuVtxZ.clear();
   m_out_n3DHits.clear();
   m_out_nUHits.clear();
   m_out_nVHits.clear();
   m_out_nWHits.clear();
   m_out_isShower.clear();
   m_out_trackScore.clear();
   m_out_recoPDG.clear();
   m_out_isRecoPrimary.clear();
   m_out_startX.clear();
   m_out_startY.clear();
   m_out_startZ.clear();
   m_out_endX.clear();
   m_out_endY.clear();
   m_out_endZ.clear();
   m_out_dirX.clear();
   m_out_dirY.clear();
   m_out_dirZ.clear();
   m_out_centroidX.clear();
   m_out_centroidY.clear();
   m_out_centroidZ.clear();
   m_out_length1.clear();
   m_out_length2.clear();
   m_out_length3.clear();
   m_out_energy.clear();
   m_out_recoHitId.clear();
   m_out_recoHitSliceId.clear();
   m_out_recoHitClusterId.clear();
   m_out_recoHitX.clear();
   m_out_recoHitY.clear();
   m_out_recoHitZ.clear();
   m_out_recoHitE.clear();
   m_out_gotMatch.clear();
   m_out_mcPDG.clear();
   m_out_mcId.clear();
   m_out_mcLocalId.clear();
   m_out_isPrimary.clear();
   m_out_nSharedHits.clear();
   m_out_completeness.clear();
   m_out_purity.clear();
   m_out_mcEnergy.clear();
   m_out_mcPx.clear();
   m_out_mcPy.clear();
   m_out_mcPz.clear();
   m_out_mcVtxX.clear();
   m_out_mcVtxY.clear();
   m_out_mcVtxZ.clear();
   m_out_mcEndX.clear();
   m_out_mcEndY.clear();
   m_out_mcEndZ.clear();
   m_out_mcNuPDG.clear();
   m_out_mcNuId.clear();
   m_out_mcNuCode.clear();
   m_out_mcNuVtxX.clear();
   m_out_mcNuVtxY.clear();
   m_out_mcNuVtxZ.clear();
   m_out_mcNuE.clear();
   m_out_mcNuPx.clear();
   m_out_mcNuPy.clear();
   m_out_mcNuPz.clear();
   m_out_trkfitStartX.clear();
   m_out_trkfitStartY.clear();
   m_out_trkfitStartZ.clear();
   m_out_trkfitStartDirX.clear();
   m_out_trkfitStartDirY.clear();
   m_out_trkfitStartDirZ.clear();
   m_out_trkfitEndX.clear();
   m_out_trkfitEndY.clear();
   m_out_trkfitEndZ.clear();
   m_out_trkfitEndDirX.clear();
   m_out_trkfitEndDirY.clear();
   m_out_trkfitEndDirZ.clear();
   m_out_trkfitLength.clear();
   m_out_trkfitSliceId.clear();
   m_out_trkfitPfoId.clear();
   m_out_trkfitX.clear();
   m_out_trkfitY.clear();
   m_out_trkfitZ.clear();
   m_out_trkfitQ.clear();
   m_out_trkfitRR.clear();
   m_out_trkfitdx.clear();
   m_out_trkfitdQdx.clear();
//shower
   m_out_shwrfitLength.clear(); 
   m_out_shwrfitCentroidX.clear();	
   m_out_shwrfitCentroidY.clear();
   m_out_shwrfitCentroidZ.clear();
   m_out_shwrfitStartX.clear();
   m_out_shwrfitStartY.clear();
   m_out_shwrfitStartZ.clear();
   m_out_shwrfitDirX.clear();
   m_out_shwrfitDirY.clear();
   m_out_shwrfitDirZ.clear();
   m_out_shwrSliceId.clear();
   m_out_shwrClusterId.clear();

}

 void NDRecoOutputData::WriteToFile()
 {
   m_treeOut->Fill();
   ClearData();
 }

 void NDRecoOutputData::CloseFile()
 {
   m_treeOut->Write();
   m_fileOut->Close();
   std::cout << "NDRecoOutputData File has been closed." << std::endl;
 }

 void NDRecoOutputData::FillBasicBranches( const std::unique_ptr<LArRecoNDFormat> &inputSpill )
 {
   m_out_event = inputSpill->m_event;
   m_out_subrun = inputSpill->m_subrun;
   m_out_run = inputSpill->m_run;
   m_out_unixTime = inputSpill->m_unixTime;
   m_out_start_t = inputSpill->m_start_t;
   m_out_end_t = inputSpill->m_end_t;
   m_out_trigger = inputSpill->m_trigger;
   m_out_sliceID.insert( m_out_sliceID.end(), inputSpill->m_sliceID->begin(), inputSpill->m_sliceID->end());
   m_out_clusterID.insert( m_out_clusterID.end(), inputSpill->m_clusterID->begin(), inputSpill->m_clusterID->end());
   m_out_nuVtxX.insert( m_out_nuVtxX.end(), inputSpill->m_nuVtxX->begin(), inputSpill->m_nuVtxX->end());
   m_out_nuVtxY.insert( m_out_nuVtxY.end(), inputSpill->m_nuVtxY->begin(), inputSpill->m_nuVtxY->end());
   m_out_nuVtxZ.insert( m_out_nuVtxZ.end(), inputSpill->m_nuVtxZ->begin(), inputSpill->m_nuVtxZ->end());
   m_out_n3DHits.insert( m_out_n3DHits.end(), inputSpill->m_n3DHits->begin(), inputSpill->m_n3DHits->end());
   m_out_nUHits.insert( m_out_nUHits.end(), inputSpill->m_nUHits->begin(), inputSpill->m_nUHits->end());
   m_out_nVHits.insert( m_out_nVHits.end(), inputSpill->m_nVHits->begin(), inputSpill->m_nVHits->end());
   m_out_nWHits.insert( m_out_nWHits.end(), inputSpill->m_nWHits->begin(), inputSpill->m_nWHits->end());
   m_out_isShower.insert( m_out_isShower.end(), inputSpill->m_isShower->begin(), inputSpill->m_isShower->end());
   m_out_trackScore.insert( m_out_trackScore.end(), inputSpill->m_trackScore->begin(), inputSpill->m_trackScore->end());
   m_out_recoPDG.insert( m_out_recoPDG.end(), inputSpill->m_recoPDG->begin(), inputSpill->m_recoPDG->end());
   m_out_isRecoPrimary.insert( m_out_isRecoPrimary.end(), inputSpill->m_isRecoPrimary->begin(), inputSpill->m_isRecoPrimary->end());
   m_out_startX.insert( m_out_startX.end(), inputSpill->m_startX->begin(), inputSpill->m_startX->end());
   m_out_startY.insert( m_out_startY.end(), inputSpill->m_startY->begin(), inputSpill->m_startY->end());
   m_out_startZ.insert( m_out_startZ.end(), inputSpill->m_startZ->begin(), inputSpill->m_startZ->end());
   m_out_endX.insert( m_out_endX.end(), inputSpill->m_endX->begin(), inputSpill->m_endX->end());
   m_out_endY.insert( m_out_endY.end(), inputSpill->m_endY->begin(), inputSpill->m_endY->end());
   m_out_endZ.insert( m_out_endZ.end(), inputSpill->m_endZ->begin(), inputSpill->m_endZ->end());
   m_out_dirX.insert( m_out_dirX.end(), inputSpill->m_dirX->begin(), inputSpill->m_dirX->end());
   m_out_dirY.insert( m_out_dirY.end(), inputSpill->m_dirY->begin(), inputSpill->m_dirY->end());
   m_out_dirZ.insert( m_out_dirZ.end(), inputSpill->m_dirZ->begin(), inputSpill->m_dirZ->end());
   m_out_centroidX.insert( m_out_centroidX.end(), inputSpill->m_centroidX->begin(), inputSpill->m_centroidX->end());
   m_out_centroidY.insert( m_out_centroidY.end(), inputSpill->m_centroidY->begin(), inputSpill->m_centroidY->end());
   m_out_centroidZ.insert( m_out_centroidZ.end(), inputSpill->m_centroidZ->begin(), inputSpill->m_centroidZ->end());
   m_out_length1.insert( m_out_length1.end(), inputSpill->m_length1->begin(), inputSpill->m_length1->end());
   m_out_length2.insert( m_out_length2.end(), inputSpill->m_length2->begin(), inputSpill->m_length2->end());
   m_out_length3.insert( m_out_length3.end(), inputSpill->m_length3->begin(), inputSpill->m_length3->end());
   m_out_energy.insert( m_out_energy.end(), inputSpill->m_energy->begin(), inputSpill->m_energy->end());
   m_out_recoHitId.insert( m_out_recoHitId.end(), inputSpill->m_recoHitId->begin(), inputSpill->m_recoHitId->end());
   m_out_recoHitSliceId.insert( m_out_recoHitSliceId.end(), inputSpill->m_recoHitSliceId->begin(), inputSpill->m_recoHitSliceId->end());
   m_out_recoHitClusterId.insert( m_out_recoHitClusterId.end(), inputSpill->m_recoHitClusterId->begin(), inputSpill->m_recoHitClusterId->end());
   m_out_recoHitX.insert( m_out_recoHitX.end(), inputSpill->m_recoHitX->begin(), inputSpill->m_recoHitX->end());
   m_out_recoHitY.insert( m_out_recoHitY.end(), inputSpill->m_recoHitY->begin(), inputSpill->m_recoHitY->end());
   m_out_recoHitZ.insert( m_out_recoHitZ.end(), inputSpill->m_recoHitZ->begin(), inputSpill->m_recoHitZ->end());
   m_out_recoHitE.insert( m_out_recoHitE.end(), inputSpill->m_recoHitE->begin(), inputSpill->m_recoHitE->end());
   m_out_gotMatch.insert( m_out_gotMatch.end(), inputSpill->m_gotMatch->begin(), inputSpill->m_gotMatch->end());
   m_out_mcPDG.insert( m_out_mcPDG.end(), inputSpill->m_mcPDG->begin(), inputSpill->m_mcPDG->end());
   m_out_mcId.insert( m_out_mcId.end(), inputSpill->m_mcId->begin(), inputSpill->m_mcId->end());
   m_out_mcLocalId.insert( m_out_mcLocalId.end(), inputSpill->m_mcLocalId->begin(), inputSpill->m_mcLocalId->end());
   m_out_isPrimary.insert( m_out_isPrimary.end(), inputSpill->m_isPrimary->begin(), inputSpill->m_isPrimary->end());
   m_out_nSharedHits.insert( m_out_nSharedHits.end(), inputSpill->m_nSharedHits->begin(), inputSpill->m_nSharedHits->end());
   m_out_completeness.insert( m_out_completeness.end(), inputSpill->m_completeness->begin(), inputSpill->m_completeness->end());
   m_out_purity.insert( m_out_purity.end(), inputSpill->m_purity->begin(), inputSpill->m_purity->end());
   m_out_mcEnergy.insert( m_out_mcEnergy.end(), inputSpill->m_mcEnergy->begin(), inputSpill->m_mcEnergy->end());
   m_out_mcPx.insert( m_out_mcPx.end(), inputSpill->m_mcPx->begin(), inputSpill->m_mcPx->end());
   m_out_mcPy.insert( m_out_mcPy.end(), inputSpill->m_mcPy->begin(), inputSpill->m_mcPy->end());
   m_out_mcPz.insert( m_out_mcPz.end(), inputSpill->m_mcPz->begin(), inputSpill->m_mcPz->end());
   m_out_mcVtxX.insert( m_out_mcVtxX.end(), inputSpill->m_mcVtxX->begin(), inputSpill->m_mcVtxX->end());
   m_out_mcVtxY.insert( m_out_mcVtxY.end(), inputSpill->m_mcVtxY->begin(), inputSpill->m_mcVtxY->end());
   m_out_mcVtxZ.insert( m_out_mcVtxZ.end(), inputSpill->m_mcVtxZ->begin(), inputSpill->m_mcVtxZ->end());
   m_out_mcEndX.insert( m_out_mcEndX.end(), inputSpill->m_mcEndX->begin(), inputSpill->m_mcEndX->end());
   m_out_mcEndY.insert( m_out_mcEndY.end(), inputSpill->m_mcEndY->begin(), inputSpill->m_mcEndY->end());
   m_out_mcEndZ.insert( m_out_mcEndZ.end(), inputSpill->m_mcEndZ->begin(), inputSpill->m_mcEndZ->end());
   m_out_mcNuPDG.insert( m_out_mcNuPDG.end(), inputSpill->m_mcNuPDG->begin(), inputSpill->m_mcNuPDG->end());
   m_out_mcNuId.insert( m_out_mcNuId.end(), inputSpill->m_mcNuId->begin(), inputSpill->m_mcNuId->end());
   m_out_mcNuCode.insert( m_out_mcNuCode.end(), inputSpill->m_mcNuCode->begin(), inputSpill->m_mcNuCode->end());
   m_out_mcNuVtxX.insert( m_out_mcNuVtxX.end(), inputSpill->m_mcNuVtxX->begin(), inputSpill->m_mcNuVtxX->end());
   m_out_mcNuVtxY.insert( m_out_mcNuVtxY.end(), inputSpill->m_mcNuVtxY->begin(), inputSpill->m_mcNuVtxY->end());
   m_out_mcNuVtxZ.insert( m_out_mcNuVtxZ.end(), inputSpill->m_mcNuVtxZ->begin(), inputSpill->m_mcNuVtxZ->end());
   m_out_mcNuE.insert( m_out_mcNuE.end(), inputSpill->m_mcNuE->begin(), inputSpill->m_mcNuE->end());
   m_out_mcNuPx.insert( m_out_mcNuPx.end(), inputSpill->m_mcNuPx->begin(), inputSpill->m_mcNuPx->end());
   m_out_mcNuPy.insert( m_out_mcNuPy.end(), inputSpill->m_mcNuPy->begin(), inputSpill->m_mcNuPy->end());
   m_out_mcNuPz.insert( m_out_mcNuPz.end(), inputSpill->m_mcNuPz->begin(), inputSpill->m_mcNuPz->end());
 }

 void NDRecoOutputData::FillTrackBranches( const std::vector<float> &startX, const std::vector<float> &startY, const std::vector<float> &startZ,
					   const std::vector<float> &dirX,  const std::vector<float> &dirY,  const std::vector<float> &dirZ,
					   const std::vector<float> &endX, const std::vector<float> &endY, const std::vector<float> &endZ,
					   const std::vector<float> &enddirX,  const std::vector<float> &enddirY,  const std::vector<float> &enddirZ,
					   const std::vector<float> &length )
 {
   m_out_trkfitStartX.insert( m_out_trkfitStartX.end(), startX.begin(), startX.end() );
   m_out_trkfitStartY.insert( m_out_trkfitStartY.end(), startY.begin(), startY.end() );
   m_out_trkfitStartZ.insert( m_out_trkfitStartZ.end(), startZ.begin(), startZ.end() );
   m_out_trkfitStartDirX.insert( m_out_trkfitStartDirX.end(), dirX.begin(), dirX.end() );
   m_out_trkfitStartDirY.insert( m_out_trkfitStartDirY.end(), dirY.begin(), dirY.end() );
   m_out_trkfitStartDirZ.insert( m_out_trkfitStartDirZ.end(), dirZ.begin(), dirZ.end() );
   m_out_trkfitEndX.insert( m_out_trkfitEndX.end(), endX.begin(), endX.end() );
   m_out_trkfitEndY.insert( m_out_trkfitEndY.end(), endY.begin(), endY.end() );
   m_out_trkfitEndZ.insert( m_out_trkfitEndZ.end(), endZ.begin(), endZ.end() );
   m_out_trkfitEndDirX.insert( m_out_trkfitEndDirX.end(), enddirX.begin(), enddirX.end() );
   m_out_trkfitEndDirY.insert( m_out_trkfitEndDirY.end(), enddirY.begin(), enddirY.end() );
   m_out_trkfitEndDirZ.insert( m_out_trkfitEndDirZ.end(), enddirZ.begin(), enddirZ.end() );
   m_out_trkfitLength.insert( m_out_trkfitLength.end(), length.begin(), length.end() );
 }

void NDRecoOutputData::FillShowerBranches( const std::vector<float> &shwrcentX, const std::vector<float> &shwrcentY, const std::vector<float> &shwrcentZ,                                           const std::vector<float> &shwrstartX, const std::vector<float> &shwrstartY, const std::vector<float> &shwrstartZ,
                                           const std::vector<float> &shwrdirX, const std::vector<float> &shwrdirY, const std::vector<float> &shwrdirZ,                                                 const std::vector<float> &shwrlength, const std::vector<int> &shwrSlice, const std::vector<int> &shwrCluster)
{
  m_out_shwrfitCentroidX.insert(m_out_shwrfitCentroidX.end(), shwrcentX.begin(), shwrcentX.end() );
  m_out_shwrfitCentroidY.insert(m_out_shwrfitCentroidY.end(), shwrcentY.begin(), shwrcentY.end() );
  m_out_shwrfitCentroidZ.insert(m_out_shwrfitCentroidZ.end(), shwrcentZ.begin(), shwrcentZ.end() );
  m_out_shwrfitStartX.insert(m_out_shwrfitStartX.end(), shwrstartX.begin(), shwrstartX.end() );
  m_out_shwrfitStartY.insert(m_out_shwrfitStartY.end(), shwrstartY.begin(), shwrstartY.end() );
  m_out_shwrfitStartZ.insert(m_out_shwrfitStartZ.end(), shwrstartZ.begin(), shwrstartZ.end() );
  m_out_shwrfitDirX.insert(m_out_shwrfitDirX.end(), shwrdirX.begin(), shwrdirX.end() );
  m_out_shwrfitDirY.insert(m_out_shwrfitDirY.end(), shwrdirY.begin(), shwrdirY.end() );
  m_out_shwrfitDirZ.insert(m_out_shwrfitDirZ.end(), shwrdirZ.begin(), shwrdirZ.end() );
  m_out_shwrfitLength.insert(m_out_shwrfitLength.end(), shwrlength.begin(), shwrlength.end() );
  m_out_shwrSliceId.insert(m_out_shwrSliceId.end(), shwrSlice.begin(), shwrSlice.end());
  m_out_shwrClusterId.insert(m_out_shwrClusterId.end(), shwrCluster.begin(), shwrCluster.end());

}

void NDRecoOutputData::FillTrackCaloBranches( const std::vector<int> &tfSliceId, const std::vector<int> &tfPfoId,
					      const std::vector<float> &tfX, const std::vector<float> &tfY, const std::vector<float> &tfZ, const std::vector<float> &tfQ,
					      const std::vector<float> &tfRR, const std::vector<float> &tfdx, const std::vector<float> &tfdQdx )
{
  m_out_trkfitSliceId.insert( m_out_trkfitSliceId.end(), tfSliceId.begin(), tfSliceId.end() );
  m_out_trkfitPfoId.insert( m_out_trkfitPfoId.end(), tfPfoId.begin(), tfPfoId.end() );
  m_out_trkfitX.insert( m_out_trkfitX.end(), tfX.begin(), tfX.end() );
  m_out_trkfitY.insert( m_out_trkfitY.end(), tfY.begin(), tfY.end() );
  m_out_trkfitZ.insert( m_out_trkfitZ.end(), tfZ.begin(), tfZ.end() );
  m_out_trkfitQ.insert( m_out_trkfitQ.end(), tfQ.begin(), tfQ.end() );
  m_out_trkfitRR.insert( m_out_trkfitRR.end(), tfRR.begin(), tfRR.end() );
  m_out_trkfitdx.insert( m_out_trkfitdx.end(), tfdx.begin(), tfdx.end() );
  m_out_trkfitdQdx.insert( m_out_trkfitdQdx.end(), tfdQdx.begin(), tfdQdx.end() );
}

} // namespace lar_nd_postreco

#endif // #ifndef PANDORA_ND_OUTERFACE_H
