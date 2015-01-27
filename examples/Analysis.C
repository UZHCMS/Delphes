//          ./Analysis delphes_output.root outfile.root

extern double PTmedian_first;
double PTmedian_second;
double PTmedian_third;
double PTmedian_fourth;

//remember to add your plots here
struct MyPlots
{
  TH1 *hJetBDeltaR;
  TH1 *hJetTrackDeltaR;
  TH1 *hJetEta;
  
  TH1 *hBJetEta;
  TH1 *hBJetPt;
  
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
  
  TH1 *hBJetEta_first;
  TH1 *hBJetEta_second;
  TH1 *hBJetEta_third;
  TH1 *hBJetEta_fourth;

  TH1 *hBJetPT_Etafirst;
  TH1 *hBJetPT_Etasecond;
  TH1 *hBJetPT_Etathird;
  TH1 *hBJetPT_Etafourth;

  
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
  THStack *stack_Bjets_eta;
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


  plots->hBJetEta_first = result->AddHist1D(
					       "hBJetEta_first", "B-jet first eta",
					       "eta","number of tracks",
					       100, -2.5, 2.5);
      
  plots->hBJetEta_second = result->AddHist1D(
						"hBJetEta_second", "B-jet second eta",
						"eta","number of tracks",
						100, -2.5, 2.5);
  
  plots->hBJetEta_third = result->AddHist1D(
					       "hBJetEta_third", "B-jet third eta",
					       "eta","number of tracks",
					       100, -2.5, 2.5);

  plots->hBJetEta_fourth = result->AddHist1D(
					       "hBJetEta_fourth", "B-jet third eta",
					       "eta","number of tracks",
					       100, -2.5, 2.5);
 
  


  plots->hBJetPT_Etafirst = result->AddHist1D(
					       "hBJetPT_Etafirst", "B-jet Pt first eta",
					       "eta","number of tracks",
					       80, 20., 1000.);
      
  plots->hBJetPT_Etasecond = result->AddHist1D(
						"hBJetPT_Etasecond", "B-jet Pt second eta",
						"eta","number of tracks",
						80, 20., 1000.);
  
  plots->hBJetPT_Etathird = result->AddHist1D(
					       "hBJetPT_Etathird", "B-jet Pt third eta",
					       "eta","number of tracks",
					       80, 20., 1000.);

  plots->hBJetPT_Etafourth = result->AddHist1D(
					       "hBJetPT_Etafourth", "B-jet Pt third eta",
					       "eta","number of tracks",
					       80, 20., 1000.);
 
 
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


  plots->hBJetEta_first->SetLineColor(kRed);
  plots->hBJetEta_second->SetLineColor(kBlue);
  plots->hBJetEta_third->SetLineColor(kGreen);
  plots->hBJetEta_fourth->SetLineColor(kMagenta);

  plots->hBJetPT_Etafirst->SetLineColor(kRed);
  plots->hBJetPT_Etasecond->SetLineColor(kBlue);
  plots->hBJetPT_Etathird->SetLineColor(kGreen);
  plots->hBJetPT_Etafourth->SetLineColor(kMagenta);
  
  stack_Bjets_eta = result->AddHistStack("Leading B-jet in Eta", "Leading b-jet eta ");
  stack_Bjets_eta->Add(plots->hBJetEta_first);
  stack_Bjets_eta->Add(plots->hBJetEta_second);
  stack_Bjets_eta->Add(plots->hBJetEta_third);
  stack_Bjets_eta->Add(plots->hBJetEta_fourth);

  
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
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  //   TClonesArray *branchEFlowTrack = treeReader->UseBranch("Track");
  
  //Objects 
  
  GenParticle *particle;
  
  Jet *jet;
  
  Track *track;
  Track *LeadingTrack;
  
  TVector3 vParticle; 
  TVector3 vJet;
  TVector3 vTrack;
  
  //Cuts
  Double_t pTMin  = 20;
  Double_t etaMin = 0.01;
  Double_t etaMax = 1.6;
  Double_t dRmin  = 0.5;
  Double_t PartondRmin  = 0.3;
  
  Double_t PartonPTMin = 1.0;
  Double_t PartonEtaMax  = 4.0;   
  
  
  Int_t i, j, k, pdgCode, pdgCodeMax, trackpT, trackpTMax; 
  
  Double_t dxy; //Impact Parameter
  
  Long64_t allEntries = treeReader->GetEntries();
  Long64_t entry;  
  
  
  cout << "** Reading " << allEntries << " events" << endl;
  
  
  
  // Loop over events
  //     for(entry = 0; entry < 1; ++entry)
  for(entry = 0; entry < allEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      
      std::vector<std::pair<double,double>> eta_bjets;	       
      
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
	    } 
	  
	  if(pdgCodeMax == 0) pdgCodeMax = 21;
	  if(pdgCodeMax == -1) pdgCodeMax = 0;
	  
	  if( (pdgCodeMax != 21 && pdgCodeMax > 5) || pdgCodeMax == 0 ) continue;
	  
	  // If b jet
	  if(pdgCodeMax == 5)
	    {	  
	      plots->hJetBDeltaR->Fill(PartonJetdR);
	      plots->hBJetEta->Fill(jet->Eta);
	      plots->hBJetPt->Fill(jet->PT);
	      
	      std::pair <double,double> eta_pt = std::make_pair(jet->Eta, jet->PT);
	      //std::cout << "foo: " << eta_pt.first << ", " << eta_pt.second << '\n';
	      eta_bjets.push_back(eta_pt);
	    }
	  
	  trackpTMax = -1;
	  Int_t NrTracks = 0;
	  
	  std::vector<double> IP;
	  
	  // Loop over tracks
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
	    }
      
	  std::sort(IP.begin(),IP.end());
	  std::reverse(IP.begin(),IP.end());
	  
	  
	  // calculate Impact parameter
	  dxy = ComputeIP(LeadingTrack);
	  
	  // If b jet
	  if(pdgCodeMax == 5)
	    {	  
	      plots->hBJetIP->Fill(dxy);  
	      plots->h_bJetLeadingTrackPt->Fill(LeadingTrack->PT);
	      plots->h_bJetTrackMultiplicity->Fill(NrTracks);
	      
	      if(IP.size()>0)plots->hBJetIP_first->Fill(IP[0]);  
	      if(IP.size()>1)plots->hBJetIP_second->Fill(IP[1]);  
	      if(IP.size()>2)plots->hBJetIP_third->Fill(IP[2]);  
	      
	      IP.clear();
	      
	    }
	  
	  // If c jet
	  else if(pdgCodeMax == 4)
	    {	  
	      plots->hCJetIP->Fill(dxy);
	      
	      if(IP.size()>0)plots->hCJetIP_first->Fill(IP[0]);  
	      if(IP.size()>1)plots->hCJetIP_second->Fill(IP[1]);  
	      if(IP.size()>2)plots->hCJetIP_third->Fill(IP[2]);  
	      
	      IP.clear();
	      
	    }
      
	  // If light jet
	  else if( (pdgCodeMax > 0) && pdgCodeMax <=3 || pdgCodeMax ==21 )
	    {	  
	      plots->hLightJetIP->Fill(dxy);
	      
	      if(IP.size()>0)plots->hLightJetIP_first->Fill(IP[0]);  
	      if(IP.size()>1)plots->hLightJetIP_second->Fill(IP[1]);  
	      if(IP.size()>2)plots->hLightJetIP_third->Fill(IP[2]);  
	      
	      IP.clear();
	      
	    }
	  
	}

      std::sort(eta_bjets.begin(),eta_bjets.end());
      std::reverse(eta_bjets.begin(),eta_bjets.end());

      //std::cout << "number of bjets " << eta_bjets.size() << std:: endl;

      //if (eta_bjets.size() > 0) std::cout << "eta bjet 1 " << eta_bjets[0].first << std::endl;
      //if (eta_bjets.size() > 0) std::cout << "eta bjet 1 " << eta_bjets[0].second << std::endl;
    
      //if (eta_bjets.size() > 1) std::cout << "eta bjet 2 " << eta_bjets[1].first << std::endl;
      //if (eta_bjets.size() > 1) std::cout << "eta bjet 2 " << eta_bjets[1].second << std::endl;
      
      //if (eta_bjets.size() > 2) std::cout << "eta bjet 3 " << eta_bjets[2].first << std::endl;
      //if (eta_bjets.size() > 2) std::cout << "eta bjet 3 " << eta_bjets[2].second << std::endl;
    
      //if (eta_bjets.size() > 3) std::cout << "eta bjet 4 " << eta_bjets[3].first << std::endl;
      //if (eta_bjets.size() > 3) std::cout << "eta bjet 4 " << eta_bjets[3].second << std::endl;
      
     	  
      if(eta_bjets.size()>0)plots->hBJetEta_first->Fill(eta_bjets[0].first);  
      if(eta_bjets.size()>1)plots->hBJetEta_second->Fill(eta_bjets[1].first);  
      if(eta_bjets.size()>2)plots->hBJetEta_third->Fill(eta_bjets[2].first);
      if(eta_bjets.size()>3)plots->hBJetEta_fourth->Fill(eta_bjets[3].first);

      if(eta_bjets.size()>0)plots->hBJetPT_Etafirst->Fill(eta_bjets[0].second);  
      if(eta_bjets.size()>1)plots->hBJetPT_Etasecond->Fill(eta_bjets[1].second);  
      if(eta_bjets.size()>2)plots->hBJetPT_Etathird->Fill(eta_bjets[2].second);
      if(eta_bjets.size()>3)plots->hBJetPT_Etafourth->Fill(eta_bjets[3].second);

    

    }

  plots->hBJetIP->Scale(1./plots->hBJetIP->Integral());
  plots->hCJetIP->Scale(1./plots->hCJetIP->Integral());
  plots->hLightJetIP->Scale(1./plots->hLightJetIP->Integral());
  
  plots->hLightJetIP_first->Scale(1./plots->hLightJetIP_first->Integral());
  plots->hCJetIP_first->Scale(1./plots->hCJetIP_first->Integral());
  plots->hBJetIP_first->Scale(1./plots->hBJetIP_first->Integral());
  
  std::cout <<plots->hLightJetIP->ClassName()<<std::endl;
  Double_t median = GetMedian(plots->h_bJetTrackPt);
  
  std::cout<<"Median b-jetTrackPt"<<median<<std::endl;
  
  PTmedian_first  =GetMedian(plots->hBJetPT_Etafirst);
  PTmedian_second =GetMedian(plots->hBJetPT_Etasecond);
  PTmedian_third  =GetMedian(plots->hBJetPT_Etathird);
  PTmedian_fourth =GetMedian(plots->hBJetPT_Etafourth);


  std::cout << "TEST2 " <<  GetMedian(plots->hBJetPT_Etafirst) << std::endl; 
  std::cout << "TEST2 " <<  GetMedian(plots->hBJetPT_Etasecond) << std::endl; 
  std::cout << "TEST2 " <<  GetMedian(plots->hBJetPT_Etathird) << std::endl; 
  std::cout << "TEST2 " <<  GetMedian(plots->hBJetPT_Etafourth) << std::endl; 


  
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  
  
  result->Print("png");
}

//------------------------------------------------------------------------------


void Analysis (const char *inputFile, const char *outputFile)
{
  gSystem->Load("libDelphes");
  
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();
  
  MyPlots *plots = new MyPlots;
  
  BookHistograms(result, plots);
  
  AnalyseEvents(treeReader, plots);
  
  //PrintHistograms(result, plots);

  TPaveText *pt;
  pt = result->AddComment(0.14,0.70,0.33,0.86);
  pt->AddText(Form("Median p_{T}=%.5f", PTmedian_first));
  result->Attach(plots->hBJetEta_first, pt);
  
  result->Write(outputFile);

  cout << "** Exiting..." << endl;
  
  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

