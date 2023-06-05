#include "MyJetAnalysis.h"

#include <trackbase_historic/SvtxTrackMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4jets/JetMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

MyJetAnalysis::MyJetAnalysis(const std::string& calojetname, const std::string& pfjetname, const std::string& truthjetname, const std::string& outputfilename)
  : SubsysReco("MyJetAnalysis_" + calojetname + " _" + pfjetname+ "_" + truthjetname)
  , m_caloJetName(calojetname)
  , m_pfJetName(pfjetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_ptThreshold(3.)
  , m_trackJetMatchingRadius(.7)
{
  m_trackdR.fill(std::numeric_limits<float>::signaling_NaN());
  m_trackpT.fill(std::numeric_limits<float>::signaling_NaN());
}

int MyJetAnalysis::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    std::cout << "MyJetAnalysis::Init - Outoput to " << m_outputFileName << std::endl;

  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  m_T = new TTree("T", "MyJetAnalysis Tree");
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("id", &m_id, "id/I");
  m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
  m_T->Branch("eta", &m_eta, "eta/F");
  m_T->Branch("phi", &m_phi, "phi/F");
  m_T->Branch("e", &m_e, "e/F");
  m_T->Branch("pt", &m_pt, "pt/F");
  m_T->Branch("truthID", &m_truthID, "truthID/I");
  m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
  m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
  m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
  m_T->Branch("truthE", &m_truthE, "truthE/F");
  m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
  m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
  m_T->Branch("trackdR", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  m_T->Branch("trackpT", m_trackpT.data(), "trackpT[nMatchedTrack]/F");

  // particle level jets
  m_T->Branch("tj_eta", &m_tj_eta, "tj_eta/F");
  m_T->Branch("tj_phi", &m_tj_phi, "tj_phi/F");
  m_T->Branch("tj_e",   &m_tj_e,   "tj_e/F");
  m_T->Branch("tj_pt",  &m_tj_pt,  "tj_pt/F");
  // detector level jets - calo
  m_T->Branch("rcj_eta",   &m_rcj_eta,  "rcj_eta/F");
  m_T->Branch("rcj_phi",   &m_rcj_phi,  "rcj_phi/F");
  m_T->Branch("rcj_e",     &m_rcj_e,    "rcj_e/F");
  m_T->Branch("rcj_pt",    &m_rcj_pt,   "rcj_pt/F");
  m_T->Branch("rcj_dist",  &m_rcj_dist, "rcj_dist/F");
  // detector level jets - particle flow
  m_T->Branch("rpj_eta",   &m_rpj_eta,  "rpj_eta/F");
  m_T->Branch("rpj_phi",   &m_rpj_phi,  "rpj_phi/F");
  m_T->Branch("rpj_e",     &m_rpj_e,    "rpj_e/F");
  m_T->Branch("rpj_pt",    &m_rpj_pt,   "rpj_pt/F");
  m_T->Branch("rpj_dist",  &m_rpj_dist, "rpj_dist/F");
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  std::cout << "MyJetAnalysis::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_T->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::InitRun(PHCompositeNode* topNode)
{

//  m_jetEvalStack = std::shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName));
//  m_jetEvalStack->get_stvx_eval_stack()->set_use_initial_vertex(initial_vertex);
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    std::cout << "MyJetAnalysis::process_event() entered" << std::endl;

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
    {
      std::cout
        << "MyJetAnalysis::process_event - Error can not find G4TruthInfo map node "
        << std::endl;
      exit(-1);
    }

  JetMap* calojets = findNode::getClass<JetMap>(topNode, m_caloJetName);
  if (!calojets)
  {
    std::cout
        << "MyJetAnalysis::process_event - Error can not find DST CaloJetMap node "
        << m_caloJetName << std::endl;
    exit(-1);
  }

  JetMap* pfjets = findNode::getClass<JetMap>(topNode, m_pfJetName);
  if (!pfjets)
  {
    std::cout
        << "MyJetAnalysis::process_event - Error can not find DST PFJetMap node "
        << m_pfJetName << std::endl;
    exit(-1);
  }


  JetMap* truthjets = findNode::getClass<JetMap>(topNode, m_truthJetName);
  if (!truthjets)
  {
    std::cout
        << "MyJetAnalysis::process_event - Error can not find DST TruthJetMap node "
        << m_truthJetName << std::endl;
    exit(-1);
  }

  ++m_event;

  // interface to tracks
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
  {
    trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
    if (!trackmap)

    {
      std::cout
          << "MyJetAnalysis::process_event - Error can not find DST trackmap node SvtxTrackMap" << std::endl;
      exit(-1);
    }
  }

  if(Verbosity() > 5)
  {
    std::cout<<"event "<<std::endl;
    int njets = truthjets->size();
    std::cout<<"njets "<<njets<<std::endl;

    for (JetMap::Iter iter = truthjets->begin(); iter != truthjets->end(); ++iter)
     {
      Jet* truthjet = iter->second;
      assert(truthjet);
      std::cout<<"jet pt "<<truthjet->get_pt()<<std::endl;   
     }
  }
      if(m_event%100==0) std::cout<<"event "<<m_event  <<std::endl;

///////////////////////
//  Truth jet
///////////////////////

    // Start with the leading jet, jets are stored in the order of increasing pt 
    JetMap::Iter iter = truthjets->end();
    --iter;  
    Jet* truthjet = iter->second;
    assert(truthjet);
    m_tj_eta= truthjet->get_eta();
    m_tj_phi= truthjet->get_phi();
    m_tj_pt = truthjet->get_pt();
    m_tj_e = truthjet->get_e();
    if (Verbosity() > 5) std::cout<<"lead jet pt "<<m_tj_pt<<std::endl;
/*
    for (Jet::ConstIter comp = truthjet->begin_comp(); comp != truthjet->end_comp(); ++comp)
      {
	PHG4Particle* truth = truthinfo->GetParticle((*comp).second); 
//	bool reject_particle = true;
//	for (int k = 0 ; k < 4;k++)
//	   {
//	     if (abs(truth->get_pid()) == chargedparticlespids[k])
//	       {
//		 reject_particle = false;
//	       }
//	   }
//	if (reject_particle) {continue;}
	float m_truthpx = truth->get_px();
	float m_truthpy = truth->get_py();
	float m_truthpz = truth->get_pz();
	float m_truthenergy = truth->get_e(); 
	float m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);
	float m_truthphi = atan2(m_truthpy , m_truthpx);
	float m_trutheta = atanh(m_truthpz / m_truthenergy); 
//	m_jc_pt.push_back(m_truthpt);
//	m_jc_phi.push_back(m_truthphi);
//	m_jc_eta.push_back(m_trutheta);
//	m_jc_index.push_back(jetnumber);
      }
*/
    if(!(passcut(m_tj_eta, m_tj_pt))) return Fun4AllReturnCodes::EVENT_OK;
   
///////////////////////
//  Calorimeter jet
///////////////////////

     JetMap::Iter iter2 = calojets->end();
     --iter2;
     Jet* calojet = iter2->second;
     assert(calojet); 
     float eta_ = calojet->get_eta(); 
     float phi_ = calojet->get_phi();
     float dphi_ = get_dphi(phi_, m_tj_phi);
     float dr_ = sqrt(pow(dphi_,2)+pow(eta_ -m_tj_eta,2));  
     m_rcj_eta = eta_;
     m_rcj_phi = phi_;
     m_rcj_dist = dr_;
     m_rcj_pt = calojet->get_pt();
     m_rcj_e = calojet->get_e();

     if(fabs(dphi_) > 2.5)
     {// if photon formed a jet
     //std::cout << "In this jet, pt: "<<m_rcj_pt << std::endl; 
     --iter2;
     calojet = iter2->second;
     assert(calojet);
     eta_ = calojet->get_eta();
     phi_ = calojet->get_phi();
     dphi_ = get_dphi(phi_, m_tj_phi);
     dr_ = sqrt(pow(dphi_,2)+pow(eta_ -m_tj_eta,2));

     if(dr_ < m_rcj_dist && calojet->get_pt() >= m_ptThreshold)
     {
     m_rcj_eta = eta_;
     m_rcj_phi = phi_;
     m_rcj_dist = dr_;
     m_rcj_pt = calojet->get_pt();
     m_rcj_e = calojet->get_e();
     }
     //std::cout << "In second jet, pt: "<< calojet->get_pt()<<" dphi: "<< dphi2_ << std::endl;

      if(Verbosity() > 5)
      {
       for (Jet::ConstIter comp = calojet->begin_comp(); comp != calojet->end_comp(); ++comp)
       {
         PHG4Particle* truth = truthinfo->GetParticle((*comp).second);
         std::cout<< "pid " <<truth->get_pid()<<std::endl;
       }
      } 

     }

///////////////////////
//  Particle Flow jet
///////////////////////

     eta_=-999; phi_=-999; dphi_=-999; dr_=-999; 
     JetMap::Iter iter1 = pfjets->end();
     --iter1;
     Jet* pfjet = iter1->second;
     assert(pfjet);
     eta_ = pfjet->get_eta();
     phi_ = pfjet->get_phi();
     dphi_ = get_dphi(phi_, m_tj_phi);
     dr_ = sqrt(pow(dphi_,2)+pow(eta_ -m_tj_eta,2));
     m_rpj_eta = eta_;
     m_rpj_phi = phi_;
     m_rpj_dist = dr_;
     m_rpj_pt = pfjet->get_pt();
     m_rpj_e = pfjet->get_e();

     if(fabs(dphi_) > 2.5)
     {// if photon formed a jet
     --iter1;
     pfjet = iter1->second;
     assert(pfjet);
     eta_ = pfjet->get_eta();
     phi_ = pfjet->get_phi();
     dphi_ = get_dphi(phi_, m_tj_phi);
     dr_ = sqrt(pow(dphi_,2)+pow(eta_ -m_tj_eta,2));

     if(dr_ < m_rpj_dist && pfjet->get_pt() >= m_ptThreshold)
     {
     m_rpj_eta = eta_;
     m_rpj_phi = phi_;
     m_rpj_dist = dr_;
     m_rpj_pt = pfjet->get_pt();
     m_rpj_e = pfjet->get_e();
     }
      
    
    }
     m_T->Fill();




/*
  for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
  {
    Jet* jet = iter->second;
    assert(jet);

    bool eta_cut = (jet->get_eta() >= m_etaRange.first) and (jet->get_eta() <= m_etaRange.second);
    bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);
    if ((not eta_cut) or (not pt_cut))
    {
      if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
      {
        std::cout << "MyJetAnalysis::process_event() - jet failed acceptance cut: ";
        std::cout << "eta cut: " << eta_cut << ", ptcut: " << pt_cut << std::endl;
        std::cout << "jet eta: " << jet->get_eta() << ", jet pt: " << jet->get_pt() << std::endl;
        jet->identify();
      }
      continue;

    }

    // fill histograms

    assert(m_hInclusiveE);
    m_hInclusiveE->Fill(jet->get_e());
    assert(m_hInclusiveEta);
    m_hInclusiveEta->Fill(jet->get_eta());
    assert(m_hInclusivePhi);
    m_hInclusivePhi->Fill(jet->get_phi());

    // fill trees - jet spectrum
    Jet* truthjet = recoeval->max_truth_jet_by_energy(jet);

    m_id = jet->get_id();
    m_nComponent = jet->size_comp();
    m_eta = jet->get_eta();
    m_phi = jet->get_phi();
    m_e = jet->get_e();
    m_pt = jet->get_pt();

    m_truthID = -1;
    m_truthNComponent = -1;
    m_truthEta = NAN;
    m_truthPhi = NAN;
    m_truthE = NAN;
    m_truthPt = NAN;

    if (truthjet)
    {
      m_truthID = truthjet->get_id();
      m_truthNComponent = truthjet->size_comp();
      m_truthEta = truthjet->get_eta();
      m_truthPhi = truthjet->get_phi();
      m_truthE = truthjet->get_e();
      m_truthPt = truthjet->get_pt();
    }

    // fill trees - jet track matching
    m_nMatchedTrack = 0;

    for (SvtxTrackMap::Iter iter = trackmap->begin();
         iter != trackmap->end();
         ++iter)
    {
      SvtxTrack* track = iter->second;

      TVector3 v(track->get_px(), track->get_py(), track->get_pz());
      const double dEta = v.Eta() - m_eta;
      const double dPhi = v.Phi() - m_phi;
      const double dR = sqrt(dEta * dEta + dPhi * dPhi);

      if (dR < m_trackJetMatchingRadius)
      {
        //matched track to jet

        assert(m_nMatchedTrack < kMaxMatchedTrack);

        m_trackdR[m_nMatchedTrack] = dR;
        m_trackpT[m_nMatchedTrack] = v.Perp();

        ++m_nMatchedTrack;
      }

      if (m_nMatchedTrack >= kMaxMatchedTrack)
      {
        std::cout << "MyJetAnalysis::process_event() - reached max track that matching a jet. Quit iterating tracks" << std::endl;
        break;
      }

    }  //    for (SvtxTrackMap::Iter iter = trackmap->begin();

    m_T->Fill();
  }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
*/

  return Fun4AllReturnCodes::EVENT_OK;
}

bool MyJetAnalysis::passcut(float eta, float pt)
{
  if (!((eta >= m_etaRange.first) and (eta <= m_etaRange.second))) return false;
  if (!((pt >= m_ptRange.first) and (pt <= m_ptRange.second))) return false;
  return true;
}

float MyJetAnalysis::get_dphi(float phi1, float phi2)
{
  float dphi = phi1 - phi2;
  if(dphi > M_PI) dphi -= 2*M_PI;
  else if(dphi < -M_PI) dphi += 2*M_PI;
  return dphi;
}

