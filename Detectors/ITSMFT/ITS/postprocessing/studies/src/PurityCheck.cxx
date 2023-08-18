// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PurityCheck.cxx
/// \brief short description of the study here
/// \author Emma Yeats @cern.ch

#include "ITSStudies/PurityCheck.h"
#include "ITSStudies/ITSStudiesConfigParam.h"

#include "Framework/Task.h"
#include "ITSBase/GeometryTGeo.h"
#include "Steer/MCKinematicsReader.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "ITStracking/IOUtils.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/DCA.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "DetectorsCommonDataFormats/DetID.h"

#include <numeric>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>

namespace o2
{
namespace its
{
namespace study
{
using namespace o2::framework;
using namespace o2::globaltracking;

using namespace o2::dataformats;
using namespace o2::itsmft;
using namespace o2::its;

using GTrackID = o2::dataformats::GlobalTrackID;
using V0 = o2::dataformats::V0;
using ITSCluster = o2::BaseCluster<float>;
using mask_t = o2::dataformats::GlobalTrackID::mask_t;
using Track = o2::track::TrackParCov;
using TrackITS = o2::its::TrackITS;
using DCA = o2::dataformats::DCA;
using PID = o2::track::PID;

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

class PurityCheckStudy : public Task
{
  struct ParticleInfo {
    int event;
    int pdg;
    float pt;
    float eta;
    float phi;
    int mother;
    int first;
    unsigned short clusters = 0u;
    unsigned char isReco = 0u;
    unsigned char isFake = 0u;
    bool isPrimary = 0u;
    unsigned char storedStatus = 2; /// not stored = 2, fake = 1, good = 0
    bool canContribToVertex = false;
    std::array<std::vector<int>, 7> rofs = {
      std::vector<int>{}, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{}, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{}}; // readout frames of corresponding clusters
    o2::its::TrackITS track;
    o2::MCCompLabel lab;
  };

  struct RofInfo {
    int id = 0;
    std::vector<int> eventIds; // ID of events in rof
    std::vector<bool> usedIds; // EvtID used to calculate actual efficiency
    std::vector<bool> flaggedSimVerts;
    // Reconstructable simverts that had a
    // potential match - see access_vertex()
    std::vector<ParticleInfo> parts; // Particle usable for vertexing
    std::vector<std::vector<o2::MCCompLabel>>
      vertLabels; // Labels associated to contributors to vertex
    std::unordered_map<int, std::array<double, 3>>
      simVerts; // Map. Simulated vertices of events that can be spotted in
                // curent rof <evtId, pos[3]>
    std::vector<o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>>
      recoVerts;                     // Vertices found in current ROF
    float good_vertex_per_ROF;       // Stores the good vertices per ROF
    float fake_vertex_per_ROF;       // Stores the fake vertices per ROF
    float ambig_vertex_per_ROF;      // Stores the ambiguous vertices per ROF
    float duplicate_vertex_per_ROF;  // Duplicated good vertices per ROF
    float ratio_vertex_per_ROF;      // Ratio of good(not including duplicates) to
                                     // total reconstructed vertices, per ROF
    float ineff_vertex_per_ROF;      // Ratio of vertices not reconstructed over total
                                     // simulated, times -1, per ROF
    float ineffCount_vertex_per_ROF; // Number of simulated vertices that were not
                                     // reconstructed per ROF
    float recoeff = 0.f;             // Vertexing efficiency
    void print()
    {
      std::cout << "\n=================================== ROF " << id
                << " ============================================ \n";

      // look at all the events occuring in this ROF
      for (int i = 0; i < eventIds.size(); i++) {
        std::cout << "eventIds = " << eventIds[i] << std::endl;
      }
      std::cout << "\tsize simVerts = " << simVerts.size() << std::endl;

      // Simulated vertices
      for (auto& sV : simVerts) {
        std::cout << "\tSimulated vertex for event: " << sV.first << " vertex:"
                  << " x= " << sV.second[0] << " y= " << sV.second[1]
                  << " z= " << sV.second[2] << std::endl;
        // std::cout << "flag = " << flaggedSimVerts[sV.first] << std::endl;
        std::cout << "\t\tPotentially contributing tracks:\n";
        for (auto& part : parts) {
          if (part.lab.getEventID() == sV.first && part.canContribToVertex) {
            // organized as follows: source ID, event ID, fake(-) or correct (+) and
            // track ID. Then particle mass and PDG ID.
            std::cout << "\t\t\t" << part.lab << "\t" << part.pt << " [GeV]\t"
                      << part.pdg << std::endl;
          }
        }
        std::cout << std::endl;
      }

      // Reconstructed vertices
      for (size_t iV{0}; iV < recoVerts.size(); ++iV) {
        auto l = getMainLabel(vertLabels[iV]);
        auto eventID = l.isSet() ? l.getEventID() : -1;
        std::cout << "\tReconstructed vertex for event: " << eventID
                  << " (-1: fake):"
                  << " x= " << recoVerts[iV].getX()
                  << " y= " << recoVerts[iV].getY()
                  << " z= " << recoVerts[iV].getZ() << std::endl;
        std::cout << "\t\tContributor labels:\n";
        for (auto& contrib : vertLabels[iV]) {
          // ordered as follows: source ID, event ID, fake (-) or correct (+) and
          // track ID
          std::cout << "\t\t\t" << contrib << std::endl;
        }
      }

      // Efficiency
      if (simVerts.size() || recoVerts.size()) {
        std::cout << "\n\tEfficiency: " << recoeff * 100 << " %\n";
      }
    };
    void uniqeff()
    {
      auto c{0};
      int current{-42};
      std::sort(parts.begin(), parts.end(), [](ParticleInfo& lp, ParticleInfo& rp) {
        return lp.lab.getEventID() > rp.lab.getEventID();
      }); // sorting at this point should be harmless.
      for (auto& p : parts) {
        if (p.lab.getEventID() != current) {
          eventIds.push_back(p.lab.getEventID());
          current = p.lab.getEventID();
        }
      }

      usedIds.resize(eventIds.size(), false);
      for (size_t iV{0}; iV < vertLabels.size(); ++iV) {
        auto label = getMainLabel(vertLabels[iV]);
        for (size_t evId{0}; evId < eventIds.size(); ++evId) {
          if (eventIds[evId] == label.getEventID() && !usedIds[evId]) {
            usedIds[evId] = true;
            ++c;
          }
        }
      }
      recoeff = (float)c / (float)eventIds.size();
    }
    int rof_id()
    {
      return id;
    };
    void access_vertex(TH1F* h1_good, TH1F* h1_fake, TH1F* h1_ambig,
                       TH1F* h1_duplicate, TH1F* h1_ratio, TH1F* h1_ineff,
                       TH1F* h1_ineffCount)
    {
      // remember iD and RofInfo::id are NOT the same! iD will iterate between
      // 1,2,3,4 and RofInfo::id will be the actual rof index

      // initialize ints to keep track of the type of vertex in each ROF
      int vertex_good = 0;
      int vertex_fake = 0;
      int vertex_ambig = 0;
      int vertex_duplicate = 0;
      int vertex_total = 0;
      int sim_total = 0;
      int vertex_ineff = 0;
      // also look at the events for each vertex to help us with the duplicates
      // later
      std::vector<int> vertexEvents = {};

      flaggedSimVerts.resize(simVerts.size(), false);
      std::cout << "simVerts size = " << RofInfo::simVerts.size() << std::endl;
      // Loop over simulated vertices
      for (auto& sV : simVerts) {
        sim_total++;
        // bool flag{false};

        // Loop over reconstructed vertices
        for (size_t iV{0}; iV < recoVerts.size(); ++iV) {
          auto l = getMainLabel(vertLabels[iV]);
          auto eventID =
            l.isSet() ? l.getEventID() : -1; //-1 is just a label for an unset

          // initialize number of contributers by type
          int contrib_rec_good = 0;
          int contrib_rec_fake = 0;
          int contrib_rec_ambig = 0;

          // initialize the total number of contributors for this vertex, and the
          // amount of contributors with a DIFFERENT eventID than the vertex itself.
          int contribs_total = 0;
          int contribs_diffEventID = 0;
          int contribs_sameEventID = 0;

          // here we compare sV.second [z index] close to reco [z index] within 8mm
          int dist = abs(sV.second[2] - recoVerts[iV].getZ());

          // Only look at reco vertices within 8 mm range
          if (dist > .8)
            continue;

          // activate flag if within range (potential match)
          // flaggedSimVerts[eventID] = true;

          // Only look at reco vertices that match event number of MC vertex
          if (eventID != sV.first)
            continue;

          // initialize list (for each good reco vertex)for counting evID occurances
          // among the contributors.
          std::vector<int> evIDcontribs = {};

          // for this single reconstructed vertex, we need to iterate over the
          // contributors (tracks) and get: 	number of good contributors 	number of
          // bad contributors 	number of all contributors 	number of contributors with
          // a different eventID than the vertex eventID.
          for (auto& contrib : vertLabels[iV]) {
            contribs_total++;
            evIDcontribs.push_back(contrib.getEventID());
            if (contrib.isCorrect()) { // check if label was assigned as for
                                       // correctly identified particle
              contrib_rec_good++;
            }
            if (contrib
                  .isFake()) { // check if label was assigned as for incorrectly
                               // identified particle or not set or noise
              contrib_rec_fake++;
            }
            if (contrib.getEventID() !=
                eventID) { // check that eventID of contributor does not match the
                           // vertex event ID
              contribs_diffEventID++;
            }
            if (contrib.getEventID() ==
                eventID) { // check that eventID of contributor does match the
                           // vertex event ID
              contribs_sameEventID++;
            }
          }
          // identify fake vertices
          if (contrib_rec_good < contrib_rec_fake) {
            vertex_fake++;
            vertex_total++;
          }
          // identify ambig vertices
          if (contrib_rec_good == contrib_rec_fake ||
              contribs_diffEventID >= contribs_sameEventID) {
            vertex_ambig++;
            vertex_total++;
          } else { // must now be a good vertex. Is it first of its evID or is it a
                   // duplicate:
            // iterate through vertexEvents, checking for repeated eventIDs.
            if (std::find(vertexEvents.begin(), vertexEvents.end(), eventID) !=
                vertexEvents.end()) {
              vertex_duplicate++;
              vertex_total++;
            } else { // must be a good vertex
              vertex_good++;
              vertex_total++;

              // add this reconstructed vertex to evID
              vertexEvents.push_back(eventID);
            }
          }
        }
      }

      // start gathering the total types of vertices by ROF ()
      good_vertex_per_ROF = vertex_good + vertex_duplicate;
      fake_vertex_per_ROF = vertex_fake;
      ambig_vertex_per_ROF = vertex_ambig;
      duplicate_vertex_per_ROF = vertex_duplicate;

      // check for ROFS with no events
      if (sim_total == 0) {
        ratio_vertex_per_ROF = 0;
        ineff_vertex_per_ROF = 0;
      }
      // check for ROFs with at least one event
      else {
        // if there are good reco vertices
        // number of good/total
        ratio_vertex_per_ROF = (float)vertex_good / (float)sim_total;
        // number of remaining "not good vertices" normalized by total
        ineff_vertex_per_ROF =
          -((float)(sim_total - vertex_good) / (float)sim_total);
        ineffCount_vertex_per_ROF = sim_total - vertex_good;
      }

      // declare histograms for all the plots
      h1_good->Fill(RofInfo::id, good_vertex_per_ROF);
      h1_fake->Fill(RofInfo::id, fake_vertex_per_ROF);
      h1_ambig->Fill(RofInfo::id, ambig_vertex_per_ROF);
      h1_duplicate->Fill(RofInfo::id, duplicate_vertex_per_ROF);
      h1_ineffCount->Fill(RofInfo::id, ineffCount_vertex_per_ROF);
      h1_ratio->Fill(RofInfo::id,
                     ratio_vertex_per_ROF); // to only look at rofs with nonzero
                                            // number of simVerts.
      h1_ineff->Fill(RofInfo::id, ineff_vertex_per_ROF);
    }
  };

 public:
  PurityCheckStudy(std::shared_ptr<DataRequest> dr,
                   std::shared_ptr<o2::base::GRPGeomRequest> gr,
                   std::shared_ptr<o2::steer::MCKinematicsReader> kineReader) : mDataRequest{dr}, mGGCCDBRequest(gr), mKineReader(kineReader){};
  ~PurityCheckStudy() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;

 private:
  // Other functions
  void process(o2::globaltracking::RecoContainer&, ProcessingContext&);

  static o2::MCCompLabel getMainLabel(std::vector<o2::MCCompLabel>&); // voting algorithm, might be claener to attach to RofInfo

  // Helper functions
  void prepareOutput();
  void updateTimeDependentParams(ProcessingContext& pc);
  void plotHistograms();

  // Data
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<o2::its::GeometryTGeo> mGeomMan;
  std::shared_ptr<DataRequest> mDataRequest;
  int mDumprof;
  std::vector<int> mClusterSizes;
  gsl::span<const int> mInputITSidxs;
  std::vector<o2::MCTrack> mMCTracks;
  std::unique_ptr<o2::dataformats::MCEventHeader> mMCEventHeader;

  // Output plots
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_good{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_fake{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_ambig{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_duplicate{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_ratio{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_ineff{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_ineffCount{};
  std::unique_ptr<TH1F> mH1_ROF_Vertexing_totals{};

  std::string mOutName;
  std::shared_ptr<o2::steer::MCKinematicsReader> mKineReader;
};

void PurityCheckStudy::init(InitContext& ic)
{
  LOGP(info, "Starting purity study...");
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  o2::base::GeometryManager::loadGeometry("./");
  auto gman = o2::its::GeometryTGeo::Instance();

  prepareOutput();

  // mKineReader = std::make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");
  for (int iEvent{0}; iEvent < mKineReader->getNEvents(0); iEvent++) {
    auto mctrk = mKineReader->getTracks(0, iEvent);
    mMCTracks.insert(mMCTracks.end(), mctrk.begin(), mctrk.end());
  }

  LOGP(info, "Purity study initialized.");
}

void PurityCheckStudy::prepareOutput()
{
  auto& params = o2::its::study::ITSCheckPurityParamConfig::Instance();
  mOutName = params.outFileName;

  mH1_ROF_Vertexing_good = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_fake = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_ambig = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_duplicate = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_ratio = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_ineff = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_ineffCount = std::make_unique<TH1F>();
  mH1_ROF_Vertexing_totals = std::make_unique<TH1F>();

  mH1_ROF_Vertexing_good->SetDirectory(nullptr);
  mH1_ROF_Vertexing_fake->SetDirectory(nullptr);
  mH1_ROF_Vertexing_ambig->SetDirectory(nullptr);
  mH1_ROF_Vertexing_duplicate->SetDirectory(nullptr);
  mH1_ROF_Vertexing_ratio->SetDirectory(nullptr);
  mH1_ROF_Vertexing_ineff->SetDirectory(nullptr);
  mH1_ROF_Vertexing_ineffCount->SetDirectory(nullptr);
  mH1_ROF_Vertexing_totals->SetDirectory(nullptr);
}

void PurityCheckStudy::run(ProcessingContext& pc)
{
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  updateTimeDependentParams(pc); // Make sure this is called after recoData.collectData, which may load some conditions
  process(recoData, pc);
}

o2::MCCompLabel PurityCheckStudy::getMainLabel(std::vector<o2::MCCompLabel>& labs)
{
  o2::MCCompLabel lab;
  size_t max_count = 0;
  for (size_t i = 0; i < labs.size(); i++) {
    size_t count = 1;
    for (size_t j = i + 1; j < labs.size(); j++) {
      if (labs[i] == labs[j] && (labs[i].isSet() && labs[j].isSet()))
        count++;
    }
    if (count > max_count)
      max_count = count;
  }

  if (max_count == 1) { // pick first valid label in case of no majority
    for (size_t i = 0; i < labs.size(); i++) {
      if (labs[i].isSet())
        return labs[i];
    }
  }

  for (size_t i = 0; i < labs.size(); i++) {
    size_t count = 1;
    for (size_t j = i + 1; j < labs.size(); j++)
      if (labs[i] == labs[j])
        count++;
    if (count == max_count)
      lab = labs[i];
  }
  return lab;
}

void PurityCheckStudy::process(o2::globaltracking::RecoContainer& recoData, ProcessingContext& pc)
{
  auto& params = o2::its::study::ITSCheckPurityParamConfig::Instance();
  auto compClus = recoData.getITSClusters();
  auto compClusMC = recoData.getITSClustersMCLabels();
  LOGP(info, "Got {} clustersMC.", compClusMC->getNElements());
  auto rofrecs = recoData.getITSClustersROFRecords();
  LOGP(info, "Got {} clusters.", compClus.size());

  // auto recLabelsArr = pc.inputs().get<std::vector<o2::MCCompLabel>>("labelsVertices");
  // LOGP(info, "got {} recLabelsArr", recLabelsArr.size());
  auto recVerArr = pc.inputs().get<std::vector<Vertex>>("vertices");
  LOGP(info, "got {} recVerArr", recVerArr.size());
  auto recVerROFArr = pc.inputs().get<std::vector<ROFRecord>>("vtxROF");
  LOGP(info, "got {} recVerROFArr", recVerROFArr.size());

  auto recLabelsArr = pc.inputs().get<std::vector<o2::MCCompLabel>*>("labelsVertices");
  LOGP(info, "got {} recLabelsArr", recLabelsArr->size());

}

void PurityCheckStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this param need to be queried only once
    initOnceDone = true;
    o2::its::GeometryTGeo* geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
  }
}

void PurityCheckStudy::endOfStream(EndOfStreamContext& ec)
{
  auto& params = o2::its::study::ITSCheckPurityParamConfig::Instance();
  plotHistograms();
}

void PurityCheckStudy::plotHistograms()
{
  LOGP(info, "plotting histos");
}

void PurityCheckStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
}

DataProcessorSpec getPurityCheckStudy(std::shared_ptr<o2::steer::MCKinematicsReader> kineReader)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->requestITSClusters(true);
  std::vector<InputSpec> inputs = dataRequest->inputs;
  inputs.emplace_back("labelsVertices", "ITS", "VERTICESMCTR", 0, Lifetime::Timeframe);
  inputs.emplace_back("vertices", "ITS", "VERTICES", 0, Lifetime::Timeframe);
  inputs.emplace_back("vtxROF", "ITS", "VERTICESROF", 0, Lifetime::Timeframe);
  // branch name matching from TrackWriterSpec.cxx:L53
  // manual inputs from TrackerSpec.cxx:L375

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              false,                             // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              inputs,
                                                              true);
  return DataProcessorSpec{
    "its-study-purity-check",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<PurityCheckStudy>(dataRequest, ggRequest, kineReader)},
    Options{}};
}
} // namespace study
} // namespace its
} // namespace o2
