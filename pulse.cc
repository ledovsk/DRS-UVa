#define pulse_cxx
#include "pulse.hh"
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



void pulse::CheckAllDRSChannels()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  TH2D *h001 = new TH2D("h001", " ped mean vs channel" ,   36, -0.5, 35.5, 100, -4000, 4000);
  TH2D *h002 = new TH2D("h002", " ped RMS vs channel" ,    36, -0.5, 35.5, 100, -100,  100);
  TH2D *h003 = new TH2D("h003", " amax vs channel" ,       36, -0.5, 35.5, 100, -100,  4000);
  TH2D *h004 = new TH2D("h004", " 4max vs channel" ,       36, -0.5, 35.5, 200,  0,    200);
  TH1D *h005 = new TH1D("h005", " ntracks ",               10, -0.5,  9.5);
  TH1D *h006 = new TH1D("h006", " log10(chi2) ",          200, -6.0, 4.0);
  TH2D *h007 = new TH2D("h007", " track Y vs X ",         100, -10,  50,   100, -10,   50);

  // Few events with WaveForms
  const int nWFmax = 10;
  TGraph *grWF[36][nWFmax];
  for(int ich=0; ich<36; ich++){
    for(int igr=0; igr<nWFmax; igr++){
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

    bool isInteresting = false;

    double pedMean[36];
    double pedRMS[36];

    for(int ic=0; ic<36; ic++){
      pedMean[ic] = 0;
      pedRMS[ic]  = 0;
      
      double sum0 = 0;
      double sum1 = 0;
      double sum2 = 0;
      
      for(int it=10; it<60; it++){
	if(fabs(channel[ic][it])>1e-6){  // skip zeros, they are not real
	  sum0 += 1.;
	  sum1 += channel[ic][it];
	  sum2 += channel[ic][it] * channel[ic][it];
	}
      }
      if(sum0>2){
	pedMean[ic] = sum1 / sum0;
	pedRMS[ic]  = sqrt( sum2 / sum0 - pedMean[ic] * pedMean[ic] );
      }
    }

    double ampMax[36];
    double timeMax[36];
    for(int ic=0; ic<36; ic++){
      ampMax[ic] = 1e+9;
      timeMax[ic] = 0;
      for(int it=10; it<1010; it++){
	if(channel[ic][it]<=ampMax[ic]){
	  ampMax[ic]  = channel[ic][it];
	  timeMax[ic] = time[ic][it];
	}
      }
      ampMax[ic] *= -1.;
    }

    for(int ic=0; ic<36; ic++){
      h001->Fill( ic, pedMean[ic] );
      h002->Fill( ic, pedRMS[ic] );
      h003->Fill( ic, ampMax[ic] );
      h004->Fill( ic, timeMax[ic] );
    }

    h005->Fill( ntracks );
    if(ntracks==1){
      if(chi2>1e-6){
	h006->Fill( log10(chi2) );
      }
      h007->Fill( xIntercept, yIntercept );
    }
    
    isInteresting = true;

    // save interesting waveforms
    if(isInteresting && ngr<nWFmax){
      for(int ic=0; ic<36; ic++){
	int group = ic / 9;
	for(int it=0; it<1024; it++){
	  grWF[ic][ngr]->SetPoint(it, time[group][it], channel[ic][it]);
	}
      }
      ngr++;
    }
  }

  TFile *fout = new TFile("output.root","recreate");
  for(int ich=0; ich<36; ich++){
    for(int igr=0; igr<std::min(nWFmax,ngr); igr++){
      grWF[ich][igr]->Write();
    }
  }
  h001->Write();
  h002->Write();
  h003->Write();
  h004->Write();
  h005->Write();
  h006->Write();
  h007->Write();
  fout->Close();
}





