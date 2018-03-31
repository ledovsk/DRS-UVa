#define pulse_cxx
#include "pulse.hh"
#include "RecoWaveForm.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraph.h>


using namespace std;

void pulse::Loop()
{
//   In a ROOT session, you can do:
//      root> .L pulse.C
//      root> pulse t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}




void pulse::SimpleCheck(int channelID = 10)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  TH1D *hist001 = new TH1D("hist001", " reco amplitude " , 1000, 0., 10000.);
  TH1D *hist002 = new TH1D("hist002", " reco time " , 1000, 0., 200.);
  TH1D *hist003 = new TH1D("hist003", " maximum sample " , 1000, 0., 3000.);
  TH1D *hist004 = new TH1D("hist004", " pedestal mean " , 1000, -100., 100.);
  TH1D *hist005 = new TH1D("hist005", " pedestal rms " , 1000, 0., 100.);
  TH1D *hist006 = new TH1D("hist006", " n Ratios " , 30, -0.5, 29.5);

  TH1D *hist010 = new TH1D("hist010", " ntracks " , 10, -0.5, 9.5);
  TH1D *hist011 = new TH1D("hist011", " log10(chi2) " , 1000, -7, 4);
  TH2D *hist012 = new TH2D("hist012", " track Y vs X " , 1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist013 = new TH2D("hist013", " track Y vs X with LYSO hit" , 1000, -10000, 40000, 1000, -10000, 40000);
  
  
  // Setup DRS channel for reconstruction with Ratios. It is assumed
  // to be a SiPM attached to LYSO (not MCP)

  int    dN = 3;
  double dT = 0.2;
  double parT[5] = { -1.20339, 3.94932, -4.77487, 5.69944, -0.630818 };
  double parA[5] = { 0, 0.0330576, 0.292408, 0.651366, 0.000751415 };
  ChannelDRS *ch = new ChannelDRS(5, parT, 5, parA, 0.1, 0.4, 0.1, dT, dN );

  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry % 100 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double waveForm[1024];
    double timeAxis[1024];
    int group = channelID / 9;
    for(int is=0; is<1024; is++){
      waveForm[is] = channel[channelID][is];
      timeAxis[is] = time[group][is];
    }
    
    RecoWaveForm *reco = new RecoWaveForm(waveForm, timeAxis, ch);

    hist001->Fill( fabs(reco->aReco()) );
    hist002->Fill( reco->tReco() );
    hist003->Fill( fabs(reco->amp(reco->imin())) );
    hist004->Fill( reco->pedMean() );
    hist005->Fill( reco->pedRMS() );
    hist006->Fill( reco->nRatios() );

    hist010->Fill( ntracks );
    if(ntracks>0 && chi2>1e-7){
      hist011->Fill( log10(chi2) );
      hist012->Fill( xIntercept, yIntercept );
      //      if( fabs(reco->aReco())>10.0*reco->pedRMS() ){
      if( fabs(reco->aReco())> 300.0  && reco->pedRMS()<10.0 ){
	hist013->Fill( xIntercept, yIntercept );
      }
    }
    
    
  }

  TFile *fout = new TFile("simpleCheck.root","recreate");
  hist001->Write();
  hist002->Write();
  hist003->Write();
  hist004->Write();
  hist005->Write();
  hist006->Write();
  hist010->Write();
  hist011->Write();
  hist012->Write();
  hist013->Write();
  fout->Close();
}


void pulse::CalibrateOneChannel(int channelID=10)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  int    dN = 3;
  double dT = 0.2;
  double parT[5] = { -1.20339, 3.94932, -4.77487, 5.69944, -0.630818 };
  double parA[5] = { 0, 0.0330576, 0.292408, 0.651366, 0.000751415 };
  ChannelDRS *ch = new ChannelDRS(5, parT, 5, parA, 0.1, 0.4, 0.1, dT, dN );
  

  // First iteration: find better parameters for time
  // reconstruction. This can be done on clipped waveforms as
  // well. Limit this step to 1000 good waveforms

  TGraphErrors *grRatios = new TGraphErrors();
  grRatios->SetName("grRatios");
  grRatios->SetTitle("T(R) before 1st iteration");
  grRatios->GetYaxis()->SetName("Ratio_{i}");
  grRatios->GetXaxis()->SetName("T_{i} - T_{0} (ns)");
  int ngr = 0;
  int nWFs = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry % 100 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double waveForm[1024];
    double timeAxis[1024];
    int group = channelID / 9;
    for(int is=0; is<1024; is++){
      waveForm[is] = channel[channelID][is];
      timeAxis[is] = time[group][is];
    }
    
    RecoWaveForm *reco = new RecoWaveForm(waveForm, timeAxis, ch);

    if(reco->amp(reco->imin())<-200.
       && reco->pedRMS()<5.0 && reco->tReco() < 200.0){
      if(nWFs>2000) break;
      nWFs++;
      for(int ip=0; ip<reco->nRatios(); ip++){
	grRatios->SetPoint(ngr, reco->ratioValue(ip), reco->timeValue(ip) - reco->tReco());
	grRatios->SetPointError(ngr, reco->ratioError(ip), 0.025);
	ngr++;
      }
    }
      
  }

  TF1 *ffit = new TF1("ffit","pol4",0.0, 0.8);
  for(int i=0; i<5; i++){
    ffit->SetParameter(i,parT[i]);
  }
  grRatios->Fit("ffit","W","",0.0,0.8);
  for(int i=0; i<5; i++){
    parT[i] = ffit->GetParameter(i);
  }
  cout << " new calibration for T(R): " << endl;
  cout << " parT[5] = { ";
  for(int i=0; i<4; i++)
    cout << parT[i] << ", ";
  cout << parT[4] << " };" << endl;

  // Apply new calibration

  for(int i=0; i<5; i++){
    ch->SetParTvsR(i, parT[i]);
  }

  // Second iteration: average pulse shape A(T). We need unclipped
  // waveforms for this

  TProfile *hAT = new TProfile("hAT", "average A(T)", 500, -50, 50);
  TH1D *hpedRMS = new TH1D("hpedRMS", "noise RMS; rms (adc)", 100, 0., 10.);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry % 100 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double waveForm[1024];
    double timeAxis[1024];
    int group = channelID / 9;
    for(int is=0; is<1024; is++){
      waveForm[is] = channel[channelID][is];
      timeAxis[is] = time[group][is];
    }
    
    RecoWaveForm *reco = new RecoWaveForm(waveForm, timeAxis, ch);

    if(reco->amp(reco->imin())>-1800. && reco->amp(reco->imin())<-200.
       && reco->pedRMS()<5.0 && reco->tReco() < 200.0){
      hpedRMS->Fill( reco->pedRMS() );
      for(int ip=0; ip<1024; ip++){
	hAT->Fill( reco->tStart() + ip * ch->timeStep() - reco->tReco(), reco->amp(ip) / reco->amp(reco->imin()) ); 
      }
    }
      
  }

  // Analyze average pulse shape A(T)

  double ampMax = -1e+9;
  double timeMax = 0;
  for(int ib=1; ib<=hAT->GetNbinsX(); ib++){
    if(hAT->GetBinContent(ib)>ampMax){
      ampMax = hAT->GetBinContent(ib);
      timeMax = hAT->GetBinCenter(ib);
    }
  }


  TGraph *grAvsR = new TGraph();
  TGraph *grTvsR = new TGraph();
  ngr = 0;
  double tCurrent = -10.0;
  while(tCurrent<20.0){
    double a1 = hAT->Interpolate(tCurrent);
    double a2 = hAT->Interpolate(tCurrent + ch->ratioStep() * ch->timeStep());
    if( a1>1e-3 && a2>1e-3 && tCurrent<timeMax && a2<ampMax*0.8){
      double r = a1 / a2;
      grTvsR->SetPoint(ngr, r, tCurrent);
      grAvsR->SetPoint(ngr, r, a1 / ampMax);
      ngr++;
    }
    tCurrent += 0.1;
  }

  TF1 *ffitT = new TF1("ffitT", "pol4", 0., 1.);
  grTvsR->Fit("ffitT","W","",0.0, 0.8);
  grTvsR->SetName("grTvsR");
  grTvsR->SetTitle("T(R) after 1st iteration");
  grTvsR->SetMarkerStyle(20);
  
  TF1 *ffitA = new TF1("ffitA", "pol4", 0., 1.);
  ffitA->FixParameter(0,0);
  grAvsR->Fit("ffitA","W","",0.0, 0.8);
  grAvsR->SetName("grAvsR");
  grAvsR->SetTitle("A(R) after 1st iteration");
  grAvsR->SetMarkerStyle(20);

  for(int i=0; i<5; i++){
    parT[i] = ffitT->GetParameter(i);
    parA[i] = ffitA->GetParameter(i);
  }

  cout << " new calibration for T(R): " << endl;
  cout << " parT[5] = { ";
  for(int i=0; i<4; i++)
    cout << parT[i] << ", ";
  cout << parT[4] << " };" << endl;

  cout << " new calibration for A(R): " << endl;
  cout << " parA[5] = { ";
  for(int i=0; i<4; i++)
    cout << parA[i] << ", ";
  cout << parA[4] << " };" << endl;


  TFile *fout = new TFile("calibrate.root","recreate");
  grRatios->Write();
  hpedRMS->Write();
  hAT->Write();
  grTvsR->Write();
  grAvsR->Write();
  fout->Close();
}



void pulse::CheckAllDRSChannels()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  TH1D *hist097 = new TH1D("hist097", " ped mean ch7" , 100, -100, 100);
  TH1D *hist098 = new TH1D("hist098", " ped RMS ch7" , 100, 0, 100);
  TH1D *hist099 = new TH1D("hist099", " ped RMS ch7" , 100, 0, 100);

  TH2D *hist101 = new TH2D("hist101", " ped mean vs channel" ,   36, -0.5, 35.5, 100, -2500, 2500);
  TH2D *hist102 = new TH2D("hist102", " ped RMS vs channel" ,    36, -0.5, 35.5, 100, 0., 100);
  TH2D *hist103 = new TH2D("hist103", " min sample vs channel" , 36, -0.5, 35.5, 1024, -0.5, 1023.5);
  TH2D *hist104 = new TH2D("hist104", " max amp vs channel" ,    36, -0.5, 35.5, 100, -2500, 2500);
  TH2D *hist105 = new TH2D("hist105", " min amp vs channel" ,    36, -0.5, 35.5, 100, -2500, 2500);
  TH2D *hist106 = new TH2D("hist106", " good amp vs channel" ,   36, -0.5, 35.5, 2510 ,-2500, 10);
  TH2D *hist107 = new TH2D("hist107", " channel pairing" ,       36, -0.5, 35.5, 36 , -0.5, 35.5);
  TH1D *hist108 = new TH1D("hist108", " MCP min amplitude" ,    1000, 0.0, 1000.);

  TH2D *hist111 = new TH2D("hist111", " nbars vs ntracks" ,    6, -0.5, 5.5, 6, -0.5, 5.5);
  TH1D *hist112 = new TH1D("hist112", " log10(chi2)" ,    100, -2, 6);
  TH2D *hist120 = new TH2D("hist120", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist121 = new TH2D("hist121", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist122 = new TH2D("hist122", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist123 = new TH2D("hist123", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist124 = new TH2D("hist124", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist125 = new TH2D("hist125", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);

  TH2D *hist126 = new TH2D("hist126", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist127 = new TH2D("hist127", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist128 = new TH2D("hist128", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);
  TH2D *hist129 = new TH2D("hist129", " Y vs X" ,    1000, -10000, 40000, 1000, -10000, 40000);

  TH1D *hist130 = new TH1D("hist130", " ch 0" , 1000, 100, -900);
  TH2D *hist131 = new TH2D("hist131", " ch 2 vs ch 1" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist132 = new TH2D("hist132", " ch 4 vs ch 3" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist133 = new TH2D("hist133", " ch 6 vs ch 5" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist134 = new TH2D("hist134", " ch 11 vs ch 10" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist135 = new TH2D("hist135", " ch 13 vs ch 12" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist136 = new TH2D("hist136", " ch 29 vs ch 28" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist137 = new TH2D("hist137", " ch 31 vs ch 30" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist138 = new TH2D("hist138", " ch 33 vs ch 32" , 200, 100,-3000, 200, 100,-3000);
  TH2D *hist139 = new TH2D("hist139", " ch 15 vs ch 14" , 200, 100,-3000, 200, 100,-3000);

  // First Ten WaveForms
  TGraph *grWF[36][10];
  for(int ich=0; ich<36; ich++){
    for(int igr=0; igr<10; igr++){
      grWF[ich][igr] = new TGraph();
      char grName[120];
      sprintf(grName,"grWF_%d_%d", ich, igr);
      grWF[ich][igr]->SetName(grName);
    }
  }

  int ngr = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry % 100 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double minAmplitude[36];
    double maxAmplitude[36];
    int    minSample[36];
    double pedMean[36];
    double pedRMS[36];

    bool isInteresting = false;

    // ch7 is empty. It should have nice pedestal and RMS
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for(int is=10; is<1000; is++){
      if(fabs(channel[7][is])>1e-6){
	sum0 += 1.;
	sum1 += channel[7][is];
	sum2 += channel[7][is] * channel[7][is];
      }
    }
    if(sum0<0.5) continue;
    pedMean[7] = sum1 / sum0;
    pedRMS[7]  = sqrt( sum2 / sum0 - pedMean[7] * pedMean[7] );
    hist097->Fill( pedMean[7] );
    hist098->Fill( pedRMS[7] );
      
    if(fabs(pedMean[7])>10) continue;


    
    hist099->Fill( pedRMS[7] );
    if(pedRMS[7]>10) continue;
    
    for(int ich=0; ich<36; ich++){

      minAmplitude[ich] =  1e+9;
      maxAmplitude[ich] = -1e+9;
      minSample[ich]    = -1;
      pedMean[ich]      =  1e+9;
      pedRMS[ich]       = -1e+9;

      sum0 = 0;
      sum1 = 0;
      sum2 = 0;
      for(int is=10; is<1024; is++){
	if(is>10 && is<110 && fabs(channel[ich][is])>1e-9){
	  sum0 += 1.;
	  sum1 += channel[ich][is];
	  sum2 += channel[ich][is] * channel[ich][is];
	}
	if(channel[ich][is] <= minAmplitude[ich]){
	  minAmplitude[ich] = channel[ich][is];
	  minSample[ich]    = is;
	}
	if(channel[ich][is] >= maxAmplitude[ich]){
	  maxAmplitude[ich] = channel[ich][is];
	}
      }
      if(sum0>0){
	pedMean[ich] = sum1 / sum0;
	pedRMS[ich]  = sqrt( sum2 / sum0 - pedMean[ich] * pedMean[ich] );
      }
      
      if( ich==1 || ich==2 ||
	  ich==3 || ich==4 ||
	  ich==5 || ich==6 ||
	  ich==10 || ich==11 ||
	  ich==12 || ich==13 ){
	if( minAmplitude[ich]<-500) isInteresting = true;
      }else if( ich==28 || ich==29 ||
		ich==30 || ich==31 ){
	//	if( minAmplitude[ich]<-200) isInteresting = true;
      }

      hist101->Fill( ich, pedMean[ich] );
      hist102->Fill( ich, pedRMS[ich] );
      hist103->Fill( ich, minSample[ich] );
      hist104->Fill( ich, maxAmplitude[ich] );
      hist105->Fill( ich, minAmplitude[ich] );
      hist106->Fill( ich, minAmplitude[ich] );
    }

    hist130->Fill(minAmplitude[0]);
    hist131->Fill(minAmplitude[1], minAmplitude[2]);
    hist132->Fill(minAmplitude[3], minAmplitude[4]);
    hist133->Fill(minAmplitude[5], minAmplitude[6]);
    hist134->Fill(minAmplitude[10], minAmplitude[11]);
    hist135->Fill(minAmplitude[12], minAmplitude[13]);
    hist136->Fill(minAmplitude[28], minAmplitude[29]);
    hist137->Fill(minAmplitude[30], minAmplitude[31]);
    hist138->Fill(minAmplitude[32], minAmplitude[33]);
    hist139->Fill(minAmplitude[14], minAmplitude[15]);

    // How many bars fired?
    int nbars = 0;
    int idL[5] = { 1, 3, 5, 10, 12 };
    int idR[5] = { 2, 4, 6, 11, 13 };
    for(int ib=0; ib<5; ib++){
      if(minAmplitude[idL[ib]] < -500 && minAmplitude[idR[ib]] < -500) nbars++;
    }
    
    hist111->Fill(ntracks, nbars);
    if(ntracks==1 && chi2>0){
      hist112->Fill( log10(chi2) );
      if(chi2<1e+5){
	hist120->Fill( xIntercept, yIntercept );
	if(nbars==1){
	  if( minAmplitude[1]<-500 && minAmplitude[2]<-500  )
	    hist121->Fill( xIntercept, yIntercept );
	  if( minAmplitude[3]<-500 && minAmplitude[4]<-500 )
	    hist122->Fill( xIntercept, yIntercept );
	  if( minAmplitude[5]<-500 && minAmplitude[6]<-500 )
	    hist123->Fill( xIntercept, yIntercept );
	  if( minAmplitude[10]<-500 && minAmplitude[11]<-500 )
	    hist124->Fill( xIntercept, yIntercept );
	  if( minAmplitude[12]<-500 && minAmplitude[13]<-500 )
	    hist125->Fill( xIntercept, yIntercept );
	}
	if(nbars<=1 &&
	   ((minAmplitude[28]<-200 && minAmplitude[29]<-200) || (minAmplitude[30]<-200 && minAmplitude[31]<-200))){
	  hist126->Fill( xIntercept, yIntercept );
	}
	if(minAmplitude[0]<-20){
	  hist127->Fill( xIntercept, yIntercept );
	}
	if(minAmplitude[14]<-1500){
	  hist128->Fill( xIntercept, yIntercept );
	}
	if(minAmplitude[15]<-1500){
	  hist129->Fill( xIntercept, yIntercept );
	}
	
	if(isInteresting){
	  hist108->Fill(fabs(minAmplitude[0]));
    }

      }
    }

    //    isInteresting = false;
    //    if(minAmplitude[1]<-1500 && minAmplitude[2]<-1500)
    //      isInteresting = true;
    

    
    // find pairing
    for(int ich=0; ich<36; ich++){
      if(ich==0 || ich==8 || ich==9 || ich==17 || ich==18 || ich==26 || ich==27 || ich==35) continue;
      if( minAmplitude[ich]<-1000 ){
	int npairs = 0;
	for(int jch=0; jch<36; jch++){
	  if(jch==0 || jch==8 || jch==9 || jch==17 || jch==18 || jch==26 || jch==27 || jch==35) continue;
	  if(jch==ich) continue;
	  if( minAmplitude[jch]<-500 ){
	     npairs++;
	  }
	}
	if(npairs<5){
	  for(int jch=0; jch<36; jch++){
	    if(jch==0 || jch==8 || jch==9 || jch==17 || jch==18 || jch==26 || jch==27 || jch==35) continue;
	    if(jch==ich) continue;
	    if( minAmplitude[jch]<-500 ){
	      hist107->Fill( ich, jch );
	    }
	  }
	}
      }
    }
    
    // save first 10 interesting waveforms
    if(isInteresting && ngr<10){
      for(int ich=0; ich<36; ich++){
	std::cout << ngr << "  " << ich << "  " << minAmplitude[ich] << std::endl;
	int group = ich / 9;
	for(int is=0; is<1024; is++){
	  grWF[ich][ngr]->SetPoint(is, time[group][is], channel[ich][is]);
	}
      }
      ngr++;
    }
  }

  TFile *fout = new TFile("checkAllDRS.root","recreate");
  hist097->Write();
  hist098->Write();
  hist099->Write();
  hist101->Write();
  hist102->Write();
  hist103->Write();
  hist104->Write();
  hist105->Write();
  hist106->Write();
  hist107->Write();
  hist108->Write();
  hist111->Write();
  hist112->Write();
  hist120->Write();
  hist121->Write();
  hist122->Write();
  hist123->Write();
  hist124->Write();
  hist125->Write();
  hist126->Write();
  hist127->Write();
  hist128->Write();
  hist129->Write();

  hist130->Write();
  hist131->Write();
  hist132->Write();
  hist133->Write();
  hist134->Write();
  hist135->Write();
  hist136->Write();
  hist137->Write();
  hist138->Write();
  hist139->Write();

  for(int ich=0; ich<36; ich++){
    for(int igr=0; igr<std::min(10,ngr); igr++){
      grWF[ich][igr]->Write();
    }
  }
  fout->Close();
}





