#ifndef MACRO_FUN4ALLJETANA_C
#define MACRO_FUN4ALLJETANA_C

#include <cstdio>
#include <iostream>
using namespace std;

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <g4jets/FastJetAlgo.h>
#include <g4jets/JetReco.h>
#include <g4jets/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <caloreco/RawClusterBuilderTopo.h>
#include <particleflowreco/ParticleFlowJetInput.h>
#include <particleflowreco/ParticleFlowReco.h>

// here you need your package name (set in configure.ac)
#include <myjetanalysis/MyJetAnalysis.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libparticleflow.so)
R__LOAD_LIBRARY(libmyjetanalysis.so)


void Fun4All_JetAna(const int nevnt = 0, const int index = 0)
{

char input1[500];
sprintf(input1,"/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/tracks/nopileup/run0006/photonjet/DST_TRACKS_pythia8_PhotonJet-0000000006-%05d.root",index);
char input4[500];
sprintf(input4,"/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/trkrhit/run0006/photonjet/DST_TRUTH_pythia8_PhotonJet-0000000006-%05d.root",index);
char input3[500];
sprintf(input3,"/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/jets/nopileup/run0006/photonjet/DST_TRUTH_JET_pythia8_PhotonJet-0000000006-%05d.root",index);
char input2[500];
sprintf(input2,"/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/calocluster/run0006/photonjet/DST_CALO_CLUSTER_pythia8_PhotonJet-0000000006-%05d.root",index);

  
  gSystem->Load("libmyjetanalysis");
  gSystem->Load("libg4dst");

  int verbosity=0;
  Fun4AllServer *se = Fun4AllServer::instance();

  JetReco *towerjetreco = new JetReco();
  towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWER));
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r02");
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r03");
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r04");
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r05");
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r06");
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r07");
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r08");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);

  //Set Calo towers
  RawClusterBuilderTopo* TopoClusterEMCAL = new RawClusterBuilderTopo("EMCalRawClusterBuilderTopo");
  TopoClusterEMCAL->Verbosity(verbosity);
  TopoClusterEMCAL->set_nodename("TOPOCLUSTER_EMCAL");
  TopoClusterEMCAL->set_enable_HCal(false);
  TopoClusterEMCAL->set_enable_EMCal(true);
  TopoClusterEMCAL->set_noise(0.0025, 0.006, 0.03);
  TopoClusterEMCAL->set_significance(4.0, 2.0, 0.0);
  TopoClusterEMCAL->allow_corner_neighbor(true);
  TopoClusterEMCAL->set_do_split(true);
  TopoClusterEMCAL->set_minE_local_max(1.0, 2.0, 0.5);
  TopoClusterEMCAL->set_R_shower(0.025);
  se->registerSubsystem(TopoClusterEMCAL);

  RawClusterBuilderTopo* TopoClusterHCAL = new RawClusterBuilderTopo("HCalRawClusterBuilderTopo");
  TopoClusterHCAL->Verbosity(verbosity);
  TopoClusterHCAL->set_nodename("TOPOCLUSTER_HCAL");
  TopoClusterHCAL->set_enable_HCal(true);
  TopoClusterHCAL->set_enable_EMCal(false);
  TopoClusterHCAL->set_noise(0.0025, 0.006, 0.03);
  TopoClusterHCAL->set_significance(4.0, 2.0, 0.0);
  TopoClusterHCAL->allow_corner_neighbor(true);
  TopoClusterHCAL->set_do_split(true);
  TopoClusterHCAL->set_minE_local_max(1.0, 2.0, 0.5);
  TopoClusterHCAL->set_R_shower(0.025);
  se->registerSubsystem(TopoClusterHCAL);

  ParticleFlowReco *pf = new ParticleFlowReco();
  pf->set_energy_match_Nsigma(1.5);
  pf->Verbosity(verbosity);
  se->registerSubsystem(pf);

  JetReco *pfjetreco = new JetReco("PARTICLEFLOWJETRECO");
  pfjetreco->add_input(new ParticleFlowJetInput());
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.2), "AntiKt_ParticleFlow_r02");
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.3), "AntiKt_ParticleFlow_r03");
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_ParticleFlow_r04");
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.5), "AntiKt_ParticleFlow_r05");
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.5), "AntiKt_ParticleFlow_r06");
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.5), "AntiKt_ParticleFlow_r07");
  pfjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.5), "AntiKt_ParticleFlow_r08");
  pfjetreco->set_algo_node("ANTIKT");
  pfjetreco->set_input_node("PARTICLEFLOW");
  pfjetreco->Verbosity(verbosity);
  se->registerSubsystem(pfjetreco);


  MyJetAnalysis *myJetAnalysis = new MyJetAnalysis("AntiKt_Tower_r04","AntiKt_ParticleFlow_r04" ,"AntiKt_Truth_r04", Form("OUTPUT/photonjet_%05d.root",index));
  //  myJetAnalysis->Verbosity(0);
  // change lower pt and eta cut to make them visible using the example
  //  pythia8 file
  myJetAnalysis->setPtRange(1, 100);
  myJetAnalysis->setEtaRange(-1.1, 1.1);
  se->registerSubsystem(myJetAnalysis);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin1");
  in->fileopen(input1);
  se->registerInputManager(in);
  in = new Fun4AllDstInputManager("DSTin2");
  in->fileopen(input2);
  se->registerInputManager(in);
  in = new Fun4AllDstInputManager("DSTin3");
  in->fileopen(input3);
  se->registerInputManager(in);
  in = new Fun4AllDstInputManager("DSTin4");
  in->fileopen(input4);
  se->registerInputManager(in);

  se->run(nevnt);
  se->End();
  delete se;
  gSystem->Exit(0);
}

#endif
