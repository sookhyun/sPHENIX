void plotPerformance()
{
  TFile* fin= new TFile("OUTPUT/photonjet.root","read");
  TTree* tt = (TTree*) fin->Get("T");
  float tj_pt;
  float tj_eta;
  float tj_phi;
  float tj_e;
  float rcj_pt;
  float rcj_eta;
  float rcj_phi;
  float rcj_e;
  float rcj_dist;
  float rpj_pt;
  float rpj_eta;
  float rpj_phi;
  float rpj_e;
  float rpj_dist;


  tt->SetBranchAddress("tj_pt", &tj_pt);
  tt->SetBranchAddress("rcj_dist",&rcj_dist);
  tt->SetBranchAddress("rcj_pt", &rcj_pt);  
  tt->SetBranchAddress("rpj_dist",&rpj_dist);
  tt->SetBranchAddress("rpj_pt", &rpj_pt);



  float jetptmin=20.;
  float jetptmax=60.;
  float jetthresh=15.;

  TH1D* hangres[2]; 
  hangres[0] = new TH1D("hangres0","Angular resolution ",50,0.,0.4);
  hangres[1] = new TH1D("hangres1","Angular resolution for particle-flow-jets",50,0.,0.4);

  TH1D* hptres[2];
  hptres[0] = new TH1D("hptres0","p_{T} resolution ",100,0.5,1.5);
  hptres[1] = new TH1D("hptres1","p_{T} resolution for particle-flow-jets",100,0.5,1.5);
  const int njetptbins=5;
  float jetptbins[njetptbins+1]={jetptmin,25.,30.,35.,40.,jetptmax};
  TH1D* hjeteffnum[2];
  hjeteffnum[0] = new TH1D("hjeteffnum0","jet reconstruction efficiency",njetptbins,jetptbins); 
  hjeteffnum[1] = new TH1D("hjeteffnum1","jet reconstruction efficiency",njetptbins,jetptbins);
  TH1D* hjeteffden[2];
  hjeteffden[0] = new TH1D("hjeteffden0","jet reconstruction efficiency den",njetptbins,jetptbins);   
  
  TH1D* hpt;
  hpt=new TH1D("hpt","hpt", 20,0,100);
  TH1D* hptreco[2];
  hptreco[0] = new TH1D("hptreco0","Reconstructed p_{T} ",20,0,100);
  hptreco[1] = new TH1D("hptreco1","Reconstructed p_{T} for particle-flow-jets",20,0,100);


  unsigned int tot= tt->GetEntries();
  cout<< "tot " << tot<<endl;
  for(unsigned int ii=0; ii<tot; ++ii)
  {
    tt->GetEntry(ii); 
    hpt->Fill(tj_pt);
    if(tj_pt < jetptmin or tj_pt > jetptmax) continue;
    hjeteffden[0]->Fill(tj_pt);   
 
    if(rcj_pt > jetthresh && rcj_dist<0.4){
      hjeteffnum[0]->Fill(tj_pt);
      float res_calo= (rcj_pt)/tj_pt;
      hptres[0]->Fill(res_calo);
    }
    if(rcj_dist<0.4)  hptreco[0]->Fill(rpj_pt);

    
    if(rpj_pt > jetthresh && rpj_dist<0.4){
      hjeteffnum[1]->Fill(tj_pt);
      float res_pf = (rpj_pt)/tj_pt;
      hptres[1]->Fill(res_pf);
    }
    if(rpj_dist<0.4)  hptreco[1]->Fill(rpj_pt);

    
    hangres[0]->Fill(rcj_dist);
    hangres[1]->Fill(rpj_dist);
    
  }

  gStyle->SetOptStat(0);
  TLatex tl= TLatex();
  TCanvas* c1=new TCanvas("c1","c1",1200,400);
  c1->Divide(3);
  c1->cd(1);
  hptres[0]->SetXTitle("p_{T}(reco)/p_{T}(truth)");
  hptres[0]->Draw();
  hptres[0]->GetXaxis()->SetTitleOffset(1.4);
  hptres[0]->SetMaximum(hptres[0]->GetMaximum()*1.4);
  tl.SetTextColor(kBlue);
  tl.DrawText(.6,hptres[0]->GetMaximum()*0.9,Form("Calo-jet Mean: %2.2f, RMS: %2.2f",hptres[0]->GetMean(), hptres[0]->GetRMS()));
  tl.SetTextColor(kRed);
  tl.DrawText(.6,hptres[0]->GetMaximum()*0.8,Form("PF-jet Mean: %2.2f, RMS: %2.2f", hptres[1]->GetMean(), hptres[1]->GetRMS()));
  hptres[1]->SetLineColor(kRed);
  hptres[1]->Draw("same");
  c1->cd(2);
  hangres[0]->SetXTitle("dR_{reco-truth} = #sqrt{d#phi^{2} + d#eta^{2}}");
  hangres[0]->Draw();
  hangres[1]->SetLineColor(kRed);
  hangres[1]->Draw("same");
  c1->cd(3);
  TEfficiency* tjeteff0 =  new TEfficiency(*hjeteffnum[0], *hjeteffden[0]);
  tjeteff0->SetTitle("Reconstruction efficiency ; p_{T} [GeV/c];");
  //tjeteff0->SetMinimum(0.5);
  tjeteff0->SetMarkerStyle(20);
  tjeteff0->SetMarkerSize(2);
  tjeteff0->SetLineColor(kBlue);
  tjeteff0->SetMarkerColor(kBlue);
  tjeteff0->Draw();
  gPad->Update(); 
  auto graph = tjeteff0->GetPaintedGraph(); 
  graph->SetMinimum(0.9);
//graph->SetMaximum(1); 
  gPad->Update(); 
  TEfficiency* tjeteff1 =  new TEfficiency(*hjeteffnum[1], *hjeteffden[0]);
  tjeteff1->SetTitle("Reconstruction efficiency for particle-flow-jets; p_{T} [GeV/c];");
  tjeteff1->SetLineColor(kRed);
  tjeteff1->SetMarkerStyle(20);
  tjeteff1->SetMarkerSize(2);
  tjeteff1->SetMarkerColor(kRed);
  tjeteff1->Draw("same");

  tl.SetTextColor(kBlack);
  tl.DrawLatex(25, 0.95, Form("jet p_{T} threshold: %2.1f GeV/c", jetthresh));

  TCanvas* c2=new TCanvas("c2","c2",800,400);
  c2->Divide(2);
  c2->cd(1);
  hpt->Draw();  
  hptreco[0]->Draw("same");
  hptreco[1]->SetLineColor(kRed);
  hptreco[1]->Draw("same");

}
