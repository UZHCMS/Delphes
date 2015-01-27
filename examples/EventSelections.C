//          ./EventSelections delphes_output.root outfile.root



//remember to add your plots here
struct MyPlots
{
  TH1 *hJetBDeltaR;
  TH1 *hJetTrackDeltaR;
  TH1 *hJetEta;
  
  TH1 *hBJetEta;
  TH1 *hBJetPt;

  TH1 *hMuonEta;
  TH1 *hMuonPt;
 
  TH1 *hEleEta;
  TH1 *hElePt;
  
  TH1 *hPhoEta;
  TH1 *hPhoPt;

  TH1 *hMET;
  
  TH1 *hParticleEta;
  TH1 *hBQuarkEta;
  
  TH1 *hBJetIP;
  TH1 *hLightJetIP; 
  TH1 *hCJetIP;  
  TH1 *h_bJetTrackPt;
  TH1 *h_bJetLeadingTrackPt;
  TH1 *h_bJetTrackMultiplicity;
  TH1 *hBJetIP_first;
  TH1 *hBJetIP_second;
  TH1 *hBJetIP_third;
  
  TH1 *hCJetIP_first;
  TH1 *hCJetIP_second;
  TH1 *hCJetIP_third;
  
  TH1 *hLightJetIP_first;
  TH1 *hLightJetIP_second;
  TH1 *hLightJetIP_third;
  
  TH2 *h_IPversusTrackPt;
  
  
};

//------------------------------------------------------------------------------
class ExRootResult;
class ExRootTreeReader;
//------------------------------------------------------------------------------

//Function to compute dR
double DeltaR(TVector3 candidate1, TVector3 candidate2) 
{
  Double_t deta = candidate1.Eta()-candidate2.Eta();
  Double_t dphi = candidate1.Phi()-candidate2.Phi();
  return sqrt( deta*deta+dphi*dphi );
}
//------------------------------------------------------------------------------

//Function to compute the impact parameter
double ComputeIP(Track *LeadingTrack) {
  
  Double_t dxy;
  TVector3 vLeadingTrack;
  vLeadingTrack.SetPtEtaPhi(LeadingTrack->PT, LeadingTrack->Eta, LeadingTrack->Phi);

   Int_t Bz = 3.8; //magnetic field
   const Double_t c_light = 2.99792458E8; //c
   
   
   Double_t r;
   Double_t x_c, y_c, r_c, phi_c, phi_0, phi;
   Double_t rcu, rc2, xd, yd;
   Double_t pT = (LeadingTrack->PT);
   Double_t X = LeadingTrack->X * 1.0E-3; //position in meters
   Double_t Y = LeadingTrack->Y * 1.0E-3; //position in meters
   Double_t Z = LeadingTrack->Z * 1.0E-3; //position in meters
   
   r = pT / (LeadingTrack->Charge * Bz) * 1.0E9/c_light; // helix radius in [m]
   phi_0 = TMath::ATan2(vLeadingTrack.Py(), vLeadingTrack.Px()); // [rad] in [-pi, pi]
  
   // 2. helix axis coordinates
   x_c = X + r*TMath::Sin(phi_0);
   y_c = Y - r*TMath::Cos(phi_0);
   r_c = TMath::Hypot(x_c, y_c);
   phi_c = TMath::ATan2(y_c, x_c);
   phi = phi_c;
   if(x_c < 0.0) phi += TMath::Pi();
   
   rcu = TMath::Abs(r);
   rc2 = r_c*r_c;
   
   // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
   xd = x_c*x_c*x_c - x_c*rcu*r_c + x_c*y_c*y_c;
   xd = (rc2 > 0.0) ? xd / rc2 : -999;
   yd = y_c*(-rcu*r_c + rc2);
   yd = (rc2 > 0.0) ? yd / rc2 : -999;
//   //     zd = z + (TMath::Sqrt(xd*xd + yd*yd) - TMath::Sqrt(x*x + y*y))*pz/pt;
  
  // calculate impact paramater
  dxy = (xd*vLeadingTrack.Py() - yd*vLeadingTrack.Px())/LeadingTrack->PT;
  dxy = dxy * 1.0E2; //return value in cm
  
//  dxy = (- LeadingTrack->X *vLeadingTrack.Py() + LeadingTrack->Y  *vLeadingTrack.Px())/LeadingTrack->PT;
  
  return dxy;
}

//------------------------------------------------------------------------------

//Function to compute the median

Double_t GetMedian(TH1 *h1) { 
   //compute the median for 1-d histogram h1 
   Int_t nbins = h1->GetXaxis()->GetNbins(); 
   Double_t *x = new Double_t[nbins]; 
   Double_t *y = new Double_t[nbins]; 
   for (Int_t i=0;i<nbins;i++) {
      x[i] = h1->GetXaxis()->GetBinCenter(i+1); 
      y[i] = h1->GetBinContent(i+1); 
   } 
   Double_t median = TMath::Median(nbins,x,y); 
   delete [] x; 
   delete [] y; 
   return median; 
}
//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, MyPlots *plots)
{
  
  THStack *stack;
  THStack *stack2;
  THStack *stack3;
  TLegend *legend;
  TLegend *legend2;
// TPaveText *comment;
  
  // AddHist1D(name, title, xlabel, ylabel, nxbins, xmin, xmax
  plots->hJetBDeltaR = result->AddHist1D(
    "hJetBDeltaR", "#DeltaR_{jet, b quark}",
    "#DeltaR_{jet, quark}", "number of jets",
      100, 0., 0.6);
  
  plots->hJetTrackDeltaR = result->AddHist1D(
    "hJetTrackDeltaR", "#DeltaR_{jet, track}",
    "R^{jet} - R^{track})", "number of jets",
     100, 0., 0.6);
  
  plots->hJetEta = result->AddHist1D(
    "jetEta", "Jet Eta",
    "Jet eta", "number of jets",
    100, -2.5, 2.5);
  
  plots->hBJetPt = result->AddHist1D(
    "hBJetPt", "b-jet p_{T}",
    "p_{T}", "number of jets",
    100, 20., 1000.);
    
      plots->hBJetEta = result->AddHist1D(
    "hBJetEta", "b-jet #eta",
    "#eta", "number of jets",
    100, -2.5, 2.5);
    
    plots->hBQuarkEta = result->AddHist1D(
    "hBQuarkEta", "b-quark #eta",
    "#eta", "number of b quarks",
    100, -2.5, 2.5);
  
  plots->hParticleEta = result->AddHist1D(
    "hParticleEta", "Gen Particle Eta",
    "Particle eta","number of jets",
    100, -2.5, 2.5);
  
  plots->hBJetIP = result->AddHist1D(
    "hBJetIP", "b-jet track IP",
    "b-jet track IP","number of tracks",
    400, -0.2, 0.2);
  
  plots->hCJetIP = result->AddHist1D(
    "hCJetIP", "c-jet track IP",
    "c-jet track IP","number of tracks",
  400, -0.2, 0.2);  
  plots->hLightJetIP = result->AddHist1D(
    "hLightJetIP", "Light jet track IP",
    "Light jet track IP","number of tracks",
  400, -0.2, 0.2);    

plots->h_bJetTrackPt = result->AddHist1D(
    "h_bJetTrackPt", "b-jet track p_{T}",
    "Track p_{T}","number of tracks",
    60, 0., 60.);
    
    plots->h_bJetLeadingTrackPt = result->AddHist1D(
    "h_bJetLeadingTrackPt", "b-jet leading track p_{T}",
    "Leading track p_{T}","number of tracks",
    300, 0., 60.);   
    
    plots->h_bJetTrackMultiplicity = result->AddHist1D(
    "h_bJetTrackMultiplicity", "b-jet track multiplicity",
    "Track multiplicity","number of jets",
    300, 0., 25.); 
    
        plots->hLightJetIP_first = result->AddHist1D(
    "hLightJetIP_first", "Light-jet first track IP",
    "track IP","number of tracks",
      400, -0.2, 0.2);
      
          plots->hLightJetIP_second = result->AddHist1D(
    "hLightJetIP_second", "Light-jet second track IP",
    "track IP","number of tracks",
      400, -0.2, 0.2);
  
    plots->hLightJetIP_third = result->AddHist1D(
    "hLightJetIP_third", "Light-jet third track IP",
    "track IP","number of tracks",
     400, -0.2, 0.2);
      
          plots->hCJetIP_first = result->AddHist1D(
    "hCJetIP_first", "c-jet first track IP",
    "track IP","number of tracks",
      400, -0.2, 0.2);
      
          plots->hCJetIP_second = result->AddHist1D(
    "hCJetIP_second", "c-jet second track IP",
    "track IP","number of tracks",
      400, -0.2, 0.2);
  
    plots->hCJetIP_third = result->AddHist1D(
    "hCJetIP_third", "c-jet third track IP",
    "track IP","number of tracks",
     400, -0.2, 0.2);
      
    plots->hBJetIP_first = result->AddHist1D(
    "hBJetIP_first", "b-jet first track IP",
    "track IP","number of tracks",
      400, -0.2, 0.2);
    
    plots->hBJetIP_second = result->AddHist1D(
    "hBJetIP_second", "b-jet second track IP",
    "track IP","number of tracks",
      400, -0.2, 0.2);
  
    plots->hBJetIP_third = result->AddHist1D(
    "hBJetIP_third", "b-jet third track IP",
    "track IP","number of tracks",
     400, -0.2, 0.2);
  
  plots->hBJetIP->SetLineColor(kRed);
  plots->hLightJetIP->SetLineColor(kBlue);
  plots->hLightJetIP->SetLineColor(kGreen);
  
  // book 1 stack of 2 histograms
  
  stack = result->AddHistStack("IP all tracks b, c, light", "IP for b, c and light");
  stack->Add(plots->hBJetIP);
  stack->Add(plots->hLightJetIP);
  stack->Add(plots->hCJetIP);
  
  plots->hBJetIP_first->SetLineColor(kRed);
  plots->hBJetIP_second->SetLineColor(kBlue);
  plots->hBJetIP_third->SetLineColor(kGreen);
  
  plots->hCJetIP_first->SetLineColor(kBlue);
  plots->hLightJetIP_first->SetLineColor(kGreen);
  
  stack2 = result->AddHistStack("b-jet IP", "b-jet first, second and third highest IP track");
  stack2->Add(plots->hBJetIP_first);
  stack2->Add(plots->hBJetIP_second);
  stack2->Add(plots->hBJetIP_third);
  
  stack3 = result->AddHistStack("Leading track IP", "Leading track IP for b, c and light ");
  stack3->Add(plots->hBJetIP_first);
  stack3->Add(plots->hCJetIP_first);
  stack3->Add(plots->hLightJetIP_first);
  
  legend2 = result->AddLegend(0.72, 0.86, 0.98, 0.98);
  legend2->AddEntry(plots->hBJetIP_first, "First", "f");
  legend2->AddEntry(plots->hBJetIP_second, "Second", "f");
  legend2->AddEntry(plots->hBJetIP_third, "Third", "f");
  legend2->Draw();
  result->Attach(stack2, legend2);
  
  // book legend for stack of 2 histograms
  
  legend = result->AddLegend(0.72, 0.86, 0.98, 0.98);
  legend->AddEntry(plots->hBJetIP, "b", "l");
  legend->AddEntry(plots->hCJetIP, "c ", "l");
  legend->AddEntry(plots->hLightJetIP, "light ", "l");
  
  
  // attach legend to stack (legend will be printed over stack in .eps file)
  
  result->Attach(stack, legend);
  
}

//------------------------------------------------------------------------------
void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
{
  
  
  
  //Branches
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  //   TClonesArray *branchEFlowTrack = treeReader->UseBranch("Track");
  
  //Objects 
  
  GenParticle *particle;
  
  Jet *jet;
  Jet *bjet1;
  Jet *bjet2;
  Jet *ljet1;
  Jet *ljet2;

  Electron *ele;
  Electron *ele1;
  Electron *ele2;

  Photon *pho;
  Muon *mu;
  Muon *mu1;
  Muon *mu2;

  MissingET *MET;
  Track *track;
  Track *LeadingTrack;
  
  TVector3 vParticle; 
  TVector3 vJet;
  TLorentzVector vBJet1;
  TLorentzVector vBJet2;
  TLorentzVector vLJet1;
  TLorentzVector vLJet2;
  TVector3 vTrack;

  TVector3 vElectron;
  TLorentzVector vElectron1;
  TLorentzVector vElectron2;

  TVector3 vPhoton;
  TVector3 vMuon;
  TLorentzVector vMuon1;
  TLorentzVector vMuon2;


  //Cuts
  Double_t pTMin  = 20;
  Double_t pTMin_ele  = 20;
  Double_t etaMin = 0.0001;
  Double_t etaMax = 2.5;
  Double_t dRmin  = 0.5;
  Double_t PartondRmin  = 0.3;
  
  Double_t PartonPTMin = 1.0;
  Double_t PartonEtaMax  = 4.0;   
  
 
  Int_t i, j, k, je, jm, jg, jme, pdgCode, pdgCodeMax, trackpT, trackpTMax; 
  Bool_t Epass, Gpass;
  
  Double_t dxy; //Impact Parameter
  
  Long64_t allEntries = treeReader->GetEntries();
  Long64_t entry;  
  
  cout << "** Reading " << allEntries << " events" << endl;
  
  
  // Loop over events
  //     for(entry = 0; entry < 1; ++entry)
  for(entry = 0; entry < allEntries; ++entry)
  {
    jet = 0;
    bjet1 = 0;
    bjet2 = 0;
    ljet1 = 0;
    ljet2 = 0;
    

    ele = 0;
    ele1 = 0;
    ele2 = 0;

    pho = 0;
    mu = 0;
    mu1 = 0;
    mu2 = 0;

    MET = 0;
    track = 0;
    LeadingTrack = 0;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::cout<< branchElectron->GetEntriesFast()<<std::endl;
    std::cout<< branchMuon->GetEntriesFast()<<std::endl;
    std::cout<< branchPhoton->GetEntriesFast()<<std::endl;
    // Loop over Electrons
    std::cout<<"loop electrons"<<std::endl;
    for(je = 0; je < branchElectron->GetEntriesFast(); ++je) {
      ele = (Electron*)branchElectron->At(je);

      if(!(ele->PT>10)) continue;
      if(!(TMath::Abs(ele->Eta) < 2.5)) continue;
      if(TMath::Abs(ele->Eta) > 1.44 && TMath::Abs(ele->Eta)< 1.56) continue;      
      Epass = kTRUE;
      if(TMath::Abs(ele->Eta)<1.44){
      	if(!(
  	     ele->EhadOverEem < 0.01
  	     && ele->IsolationVar < 0.15
  	     //&& TMath::Abs(ele->d0) < 0.04
  	     //&& TMath::Abs(ele->dz) <0.2
  	     )) Epass = kFALSE;
      }
      else{
      	if(!(
  	     ele->IsolationVar < 0.15
  	     //&& TMath::Abs(ele->d0) < 0.04
  	     //&& TMath::Abs(ele->dz) <0.2
  	     )) Epass = kFALSE;
      }
      if(!(Epass)) continue;
      // plots->hEleEta->Fill(ele->Eta);
      // plots->hElePt->Fill(ele->PT);

      // Store leading pT electron as ele1 and next leading as ele2
    	if (!ele1 && ele->PT > pTMin_ele){
    	  ele1=ele;
    	} else if (ele->PT > ele1->PT){
    	  ele2 = ele1;
    	  ele1 = ele;
    	} else if (!ele2){
    	  ele2 = ele;
    	} else if (ele->PT > ele2->PT){
    	  ele2 = ele;
    	}

    }
    
    // Loop over Photons
    for(jg = 0; jg < branchPhoton->GetEntriesFast(); ++jg){
      pho = (Photon*)branchPhoton->At(jg);
      std::cout<<"in photon loop"<<std::endl;
      if(!(pho->PT>20)) continue;
      
      // Constrain photons to barrel and end caps
      if(!(TMath::Abs(pho->Eta)) < 2.5) continue;
      if(TMath::Abs(pho->Eta) > 1.44 && TMath::Abs(pho->Eta)< 1.56) continue;
      
      Gpass = kTRUE;
      if(TMath::Abs(pho->Eta)<1.44){
      	if(!(
  	     pho->EhadOverEem < 0.05
  	     && (pho->IsolationVar*pho->PT) < 7.4 + 0.045 * pho->PT
  	     )) Gpass = kFALSE;
      }
      else{
      	if(!(
  	     pho->EhadOverEem < 0.05
  	     && (pho->IsolationVar * pho->PT) < 5.2 + 0.04*pho->PT
  	     )) Gpass = kFALSE;
      }
      if(!(Gpass)) continue;

      // plots->hPhoEta->Fill(pho->Eta);
      // plots->hPhoPt->Fill(pho->PT);
    }
    std::cout<<"finished photons"<<std::endl;
    
    std::cout<<"looping muons"<<std::endl;
    // Loop over Muons
    for(jm = 0; jm < branchMuon->GetEntriesFast(); ++jm){
      
      mu = (Muon*)branchMuon->At(jm);
	 
  	 if(!(mu->PT > 10)) continue;
  	 if(!(TMath::Abs(mu->Eta) <2.4)) continue;
	 
  	 if(mu->IsolationVar > 0.15) continue;
  	 //if(TMath::Abs(mu->dz) >0.5) continue;
  	 //if(mu->nPixHits < 1) continue;
  	 //if(mu->nTkHis < 11) continue;
  	 //if(TMath::Abs(dxy) < 2 mm) continue;
	 
	 std::cout<<"1 muons"<<std::endl;
  	 // Store leading pT muon as mu1 and next leading as mu2
    // for HH->bbZZ
      
    if((!mu1) /*&& (mu->PT > 20)*/){ //FIXME
      std::cout<<"1.5 muons"<<std::endl;
      mu1 = mu;
    } else if (mu->PT > mu1->PT){
      std::cout<<"1.6 muons"<<std::endl;
      mu2=mu1;
      mu1=mu;
    } else if (!mu2){
      std::cout<<"1.7 muons"<<std::endl;
      mu2=mu;
    } else if (mu->PT > mu2->PT){
      std::cout<<"1.8 muons"<<std::endl;
      mu2=mu;
    }
     std::cout<<"2 muons"<<std::endl;
	 
     // plots->hMuonEta->Fill(mu->Eta);
     // plots->hMuonPt->Fill(mu->PT);
     std::cout<<"3 muons"<<std::endl;
	 
   }
   std::cout<<"finished muons"<<std::endl;
   
   std::cout<<"looping MET"<<std::endl;
   
   // Missing Energy
   for(jme = 0; jme < branchMissingET->GetEntriesFast(); ++jme){
      MET = (MissingET*) branchMissingET->At(jme);
      if(MET->MET > 25) continue;
   }
   std::cout<<"finished MET"<<std::endl;
   
   std::cout<<"looping jets"<<std::endl;
    
    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      Double_t PartonJetdR = -1;
      
      pdgCodeMax = -1;
      
      jet = (Jet*) branchJet->At(i);
      
      //Apply cuts
      if(jet->PT < pTMin || TMath::Abs(jet->Eta) > etaMax || TMath::Abs(jet->Eta) < etaMin ) continue;
      
      //Define a jet 3-vector
      vJet.SetPtEtaPhi(jet->PT, jet->Eta, jet->Phi);
      
      plots->hJetEta->Fill(jet->Eta);
      
      
      // Loop over GenParticles
      for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
      {
      	particle = (GenParticle*)branchParticle->At(j);
	
      	if(particle->PT < PartonPTMin || TMath::Abs(particle->Eta) > PartonEtaMax) continue;
	
      	vParticle.SetPtEtaPhi(particle->PT, particle->Eta, particle->Phi);
	
      	plots->hParticleEta->Fill(particle->Eta);
	
      	// Get jet flavor
      	pdgCode = TMath::Abs(particle->PID);
      	if(pdgCode == 5) plots->hBQuarkEta->Fill(particle->Eta);
	
      	if(pdgCode != 21 && pdgCode > 5) continue;
      	if(pdgCode == 21) pdgCode = 0;	  
	
      	Double_t dR = DeltaR(vParticle, vJet);
      	if(dR < PartondRmin)
      	{
      	  if(pdgCodeMax < pdgCode){
      	    PartonJetdR = dR;
      	    pdgCodeMax = pdgCode;
      	  }
      	}
      } // end gen particle loop 
      
      if(pdgCodeMax == 0) pdgCodeMax = 21;
      if(pdgCodeMax == -1) pdgCodeMax = 0;
      
      if( (pdgCodeMax != 21 && pdgCodeMax > 5) || pdgCodeMax == 0 ) continue;
      
      // If b jet
      std::cout<<"if b"<<std::endl;
      if(pdgCodeMax == 5)
      {	  
      	plots->hJetBDeltaR->Fill(PartonJetdR);
      	plots->hBJetEta->Fill(jet->Eta);
      	plots->hBJetPt->Fill(jet->PT);
      }
      
      trackpTMax = -1;
      Int_t NrTracks = 0;
      
      std::vector<double> IP;
      
      // Loop over tracks
      std::cout<<"loop tracks"<<std::endl;
      for(k = 0; k < branchEFlowTrack->GetEntriesFast(); ++k)
      {
      	track = (Track*)branchEFlowTrack->At(k);
      	if(track->PT < 1.0) continue;
	
      	vTrack.SetPtEtaPhi(track->PT, track->Eta, track->Phi);
      	Double_t dR = DeltaR(vTrack, vJet);
	
      	// Get highest pT track in jet cone
      	if(dR > dRmin)continue;
	
      	Double_t dxyAllTracks = ComputeIP(track);
      	IP.push_back(dxyAllTracks);
	
      	NrTracks +=1;
      	trackpT = track->PT;
      	plots->hJetTrackDeltaR->Fill(dR);
	
      	if(pdgCodeMax == 5){
      		plots->h_bJetTrackPt->Fill(trackpT);
    		}
		
    	  if(trackpTMax < trackpT){
    	    trackpTMax = trackpT;
    	    LeadingTrack = track;
    	  }
      } // end track loop
      std::cout<<"sort IP"<<std::endl;
      std::sort(IP.begin(),IP.end());
      std::reverse(IP.begin(),IP.end());
      
      std::cout<<"dxy"<<std::endl;
      // calculate Impact parameter
      dxy = -10;
      if (LeadingTrack)
        dxy = ComputeIP(LeadingTrack);
      
      // If b jet
      std::cout<<"if b2"<<std::endl;
      if(pdgCodeMax == 5)
      {	  
        std::cout<<"if b2.1"<<std::endl;
      	// Store leading pT bjet as bjet1 and next leading as bjet2
	      if (!bjet1 /*&& jet->PT > 30*/){
          std::cout<<"if b2.2"<<std::endl;
      		bjet1=jet;
	      } else if (jet->PT > bjet1->PT){
          std::cout<<"if b2.3"<<std::endl;
      		bjet2=bjet1;
      		bjet1=jet;
	      } else if (!bjet2 && jet->PT > 30){
          std::cout<<"if b2.4"<<std::endl;
          bjet2=jet;
	      } else if (jet->PT > bjet2->PT){
          std::cout<<"if b2.5"<<std::endl;
          bjet2=jet;
	      }

       //         plots->hBJetIP->Fill(dxy);
       //         plots->h_bJetLeadingTrackPt->Fill(LeadingTrack->PT);
       //         plots->h_bJetTrackMultiplicity->Fill(NrTracks);
       //
       // if(IP.size()>0)plots->hBJetIP_first->Fill(IP[0]);
       //        if(IP.size()>1)plots->hBJetIP_second->Fill(IP[1]);
       //        if(IP.size()>2)plots->hBJetIP_third->Fill(IP[2]);
      
       IP.clear();
      
      } // end if b jet
      // If c jet
      else if(pdgCodeMax == 4)
      {	  
        std::cout<<"if c2"<<std::endl;
        plots->hCJetIP->Fill(dxy);
	
        if(IP.size()>0)plots->hCJetIP_first->Fill(IP[0]);  
        if(IP.size()>1)plots->hCJetIP_second->Fill(IP[1]);  
        if(IP.size()>2)plots->hCJetIP_third->Fill(IP[2]);  
      
        IP.clear();
      
      } // end if c jet
      // If light jet
      else if( (pdgCodeMax > 0) && (pdgCodeMax <=3 || pdgCodeMax ==21) )
      {	  
        std::cout<<"if light2"<<std::endl;
        plots->hLightJetIP->Fill(dxy);
	
        if(IP.size()>0)plots->hLightJetIP_first->Fill(IP[0]);  
        if(IP.size()>1)plots->hLightJetIP_second->Fill(IP[1]);  
        if(IP.size()>2)plots->hLightJetIP_third->Fill(IP[2]);  
      
        IP.clear();
	      if(jet->PT >20){
          // Store leading pt light jet to ljet1 and subleading to ljet2
          if (!ljet1){
            ljet1 = jet;
          } else if (jet->PT > ljet1->PT){
            ljet2=ljet1;
            ljet1=jet;
          } else if (!ljet2){
            ljet2=jet;
          } else if (jet->PT > ljet2->PT){
            ljet2=jet;
          }
	      }
      
      } // end if light jet
      
    }	// end jet loop
    
    std::cout<<"finished jets"<<std::endl;


    //  HHH->bbZZ ->bblljj
    if(bjet1 && bjet2){  
      vBJet1.SetPtEtaPhiM(bjet1->PT, bjet1->Eta, bjet1->Phi, bjet1->Mass);
      vBJet2.SetPtEtaPhiM(bjet2->PT, bjet2->Eta, bjet2->Phi, bjet2->Mass);
      // bbWW restricts dR_{bb} < 3.1  60 < M_bb < 160  
      if( ((vBJet1+vBJet2).M() <90) || ((vBJet1+vBJet2).M() >150) ) continue;

      if(ljet1 &&ljet2){
        vLJet1.SetPtEtaPhiM(ljet1->PT, ljet1->Eta, ljet1->Phi, ljet1->Mass);
        vLJet2.SetPtEtaPhiM(ljet2->PT, ljet2->Eta, ljet2->Phi, ljet2->Mass);
  
        if( ((vLJet1+vLJet2).M() <60) || ((vLJet1+vLJet2).M() >120) ) continue;
  
        if(ele1 && ele2){
          vElectron1.SetPtEtaPhiM(ele1->PT, ele1->Eta, ele1->Phi, 0.000511);
          vElectron2.SetPtEtaPhiM(ele2->PT, ele2->Eta, ele2->Phi, 0.000511);
  
          if( ((vElectron1+vElectron2).M() <80) || ((vElectron1+vElectron2).M() > 100) ) continue;
        }
        else if(mu1 && mu2){
          vMuon1.SetPtEtaPhiM(mu1->PT, mu1->Eta, mu1->Phi, 0.105);
          vMuon2.SetPtEtaPhiM(mu2->PT, mu2->Eta, mu2->Phi, 0.105);
    
          if( ((vMuon1+vMuon2).M() <80) || ((vMuon1+vMuon2).M() > 100) ) continue;
        }
      }
    } 
    std::cout<<"finished ele and mu"<<std::endl;

  } // end event loop
  
  std::cout<<"finished event loop"<<std::endl;
  
  plots->hBJetIP->Scale(1./plots->hBJetIP->Integral());
  plots->hCJetIP->Scale(1./plots->hCJetIP->Integral());
  plots->hLightJetIP->Scale(1./plots->hLightJetIP->Integral());
  
  plots->hLightJetIP_first->Scale(1./plots->hLightJetIP_first->Integral());
  plots->hCJetIP_first->Scale(1./plots->hCJetIP_first->Integral());
  plots->hBJetIP_first->Scale(1./plots->hBJetIP_first->Integral());
  
	std::cout <<plots->hLightJetIP->ClassName()<<std::endl;
  Double_t median = GetMedian(plots->h_bJetTrackPt);
  
  std::cout<<"Median b-jetTrackPt"<<median<<std::endl;
  
  
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  
  
  result->Print("png");
}

//------------------------------------------------------------------------------


void EventSelections (const char *inputFile, const char *outputFile)
{
  gSystem->Load("libDelphes");
  
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();
  
  MyPlots *plots = new MyPlots;
  
  BookHistograms(result, plots);
  
  AnalyseEvents(treeReader, plots);
  
  //  PrintHistograms(result, plots);
  
  result->Write(outputFile);
  
  cout << "** Exiting..." << endl;
  
  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}
