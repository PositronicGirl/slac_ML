#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <TROOT.h>
#include <TFile.h>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "vector"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void AnalyseEvents(ExRootTreeReader *treeReader)
{
  TClonesArray *branchCaloTower = treeReader->UseBranch("CaloTower");
  //TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
  TClonesArray *branchJetAK8 = treeReader->UseBranch("JetAK8");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  //Electron *electron;
  //Photon *pho
  Tower *tower;

  Jet *jet;
  TObject *object;

  //TLorentzVector momentum;

  //  Float_t eem, ehad;
  //Bool_t skip;

  Long64_t entry;

  Int_t i, j, k, pdgcode;
  Int_t etaBinEdge, etaBins = 31; //half of total number of bins
  Float_t etaBinStep = .0714;
  Float_t etaLowEdge = -2.2134;
  
  const Float_t pi = 3.14159;
  Int_t phiBinEdge, phiBins=44; //half of total number of  bins
  Float_t phiBinStep = pi/phiBins;
  Int_t count_jet=0;
  Int_t count_tower=0;
  
  TTree *outputTree = new TTree("jetgrids","jetgrids");
  outputfile= new TFile("delphes-tree-grids.root","recreate");
  vector<float> CaloTower_ET;
  vector<float> JetAK8_ETA;
  vector<float> JetAK8_PHI;
  vector<float> CaloEdgeEtaMin;
  vector<float> CaloEdgeEtaMax;
  vector<float> CaloEdgePhiMax;
  vector<float> CaloEdgePhiMin;
  vector<float> CaloTower_ETA;
  vector<float> CaloTower_PHI;
  

  outputTree->Branch("CaloTower_ET",&CaloTower_ET);
  outputTree->Branch("JetAK8_ETA",&JetAK8_ETA);
  outputTree->Branch("JetAK8_PHI",&JetAK8_PHI);
  outputTree->Branch("JetAK8_ET",&JetAK8_ETA);
  outputTree->Branch("CaloEdgeEtaMin",&CaloEdgeEtaMin);
  outputTree->Branch("CaloEdgeEtaMax",&CaloEdgeEtaMax);
  outputTree->Branch("CaloEdgePhiMin",&CaloEdgePhiMin);
  outputTree->Branch("CaloEdgePhiMax",&CaloEdgePhiMax);
  outputTree->Branch("CaloTower_ETA",&CaloTower_ETA);
  outputTree->Branch("CaloTower_PHI",&CaloTower_PHI);
  
  for(entry = 0; entry < allEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    JetAK8_ETA.clear();
    JetAK8_PHI.clear();
    CaloTower_ET.clear();
    CaloEdgeEtaMin.clear();
    CaloEdgeEtaMax.clear();
    CaloEdgePhiMin.clear();
    CaloEdgePhiMax.clear();
    CaloTower_ETA.clear();
    CaloTower_PHI.clear();

  
    for(i = 0; i < branchJetAK8->GetEntriesFast(); ++i){
      jet = (Jet*) branchJetAK8->At(i);
     
      //	  momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
      
      //Selecting Jet PT
      if (jet->PT >= 200 and jet->PT <= 300){
	
	
	// Find jet center
	double etacenter = jet->Eta;
	double phicenter = jet->Phi;
	

	// toss jets too close to the eta edges
	if(abs(etacenter) > (abs(etaLowEdge)-12*etaBinStep)) continue; 
	
	cout<< "Jet#:"<<i<<" etacenter:"<< etacenter<<" phicenter:"<< phicenter<<endl;
	count_jet+=1;

	JetAK8_ETA.push_back(etacenter);
	JetAK8_PHI.push_back(phicenter);


	// Find distance from lowest eta.
	double etaToLowEdge = etacenter - etaLowEdge;
	double phiToLowEdge = phicenter;  
	
	// Find the bin indices -- putting this as "int" should truncate/round to the closest
	int etaIndex = etaToLowEdge/etaBinStep;
	int phiIndex = phiToLowEdge/phiBinStep;

	// For phi we can wrap, so find the right high and low phi values
	// For simplicity, let's always use the same edge of the bin (the left edge)
	int phiIndexDn = phiIndex-12;
	int phiIndexUp = phiIndex+13;
	if(phiIndexDn < -32){
	  // phiIndexDn = 44 - (12-phiIndex);
	  // phiIndexUp = 12+phiIndex;
	  phiIndexDn = (phiIndex+32)+44; 
	  phiIndexUp = (13+phiIndex);
	}else if(phiIndexUp > 32){
	  //phiIndexDn = phiIndex-12;
          //phiIndexUp = 44-phiIndex;

	  phiIndexDn = (phiIndex-12);
	  phiIndexUp = (phiIndex-44)-32;
	}
	double phiLowEdge = phiIndexDn*phiBinStep;
	double phiHighEdge = phiIndexUp*phiBinStep;
	//
	//count_tower=0;
	// Loop over towers to find the towers in this region
	for(k = 0; k < branchCaloTower->GetEntriesFast(); ++k){
	  tower = (Tower*) branchCaloTower->At(k);
	  
	  // skip towers that are lower than low eta or higher than high eta
	  if(tower->Edges[0] < (etaLowEdge + (etaIndex-12)*etaBinStep) || tower->Edges[0] > (etaLowEdge + (etaIndex+13)*etaBinStep)) continue;
	  
	  
	  // skip towers that are in between low phi and high phi (since it wraps)
	  if(tower->Edges[2] < phiLowEdge || tower->Edges[2] > phiHighEdge) continue;
	  count_tower+=1;
	  // cout<<"PhiLow:"<<phiLowEdge <<" Phihigh:"<<phiHighEdge<<" toweredge:"<<tower->Edges[2]<<endl;
	  CaloTower_ET.push_back(tower->ET);
	  CaloEdgeEtaMin.push_back(tower->Edges[0]);
	  CaloEdgeEtaMax.push_back(tower->Edges[1]);
	  CaloEdgePhiMin.push_back(tower->Edges[2]);
	  CaloEdgePhiMax.push_back(tower->Edges[3]);
	  CaloTower_PHI.push_back(tower->Phi);
	  CaloTower_ETA.push_back(tower->Eta);
	  
	  cout << "LowEdgeEta:"<<tower->Edges[0] << " HighEdgeEta:" << tower->Edges[1]<< " Lowedge Phi:" << tower->Edges[2] << " ET:" << tower->ET<<endl;
	  cout<< "calo Phi: "<<tower->Phi<<" calo Eta:"<<tower->Eta<<endl; //" ET Phi:"<<tower->Phi<<" ET Eta:"<<tower->Eta <<endl;//<< " Tower#:" <<count_tower<< endl;
	}
      }
    }
    outputTree->Fill(); 
  }
  cout<<"Jet numbers:  "<<count_jet<<" Tower numbers: "<<count_tower<<endl;
  outputTree->Write();
}

	      
// 	      // Constraints on Eta
// 	      for(etaBinEdge = -etaBins+12 ; etaBinEdge <etaBins-12; ++etaBinEdge)
// 		{
// 		  if ( (jet->Eta <= etaBinStep*(etaBinEdge+1)) and jet->Eta >= etaBinStep*etaBinEdge)
// 		    {
// 		      //Constraints on Phi
// 		      for (phiBinEdge = -phiBins; phiBinEdge < pi; ++phiBinEdge)
// 			{ 
// 			  if ((jet->Phi <= phiBinStep*(phiBinEdge+1)) and jet->Phi >= phiBinStep*phiBinEdge)
// 			    {
// 			      count_jet+=1;
			      
// 			      //Looping through CaloTowers
// 			      for (k = 0 ; k <branchCaloTower->GetEntriesFast(); ++k)
// 				{
// 				  tower= (Tower*) branchCaloTower->At(k);
				  
// 				  //Checking if Jet is within edges of tower
// 				  if( (jet->Eta >= tower->Edges[0] and jet->Eta<= tower->Edges[1]) and (jet->Phi>=tower->Edges[2] and jet->Phi<=tower->Edges[3]))
// 				    {
// 				      count_tower+=1;
// 				      cout <<"CaloEdges: "<<tower->Edges[0]<<", "<<tower->Edges[2] <<endl;
// 				      cout <<"JetEta: "<<jet->Eta <<" JetPhi: "<<jet->Phi<<endl;
// 				      cout <<" CaloE: "<<tower->ET<<endl;
// 				    }
				  
// 				  // cout <<"tower: "<<branchCaloTower->GetEntriesFast()<<endl;
// 				}//cout<<"JEt:"<<branchJetAK8->GetEntriesFast()<<" ,Cal: "<<branchCaloTower->GetEntriesFast()<<endl;
			      
			      
			      
// 			      //			   cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;
// 			      //  cout<< " etabins: "<<etaBinStep*etaBinEdge<<", "<< etaBinStep*(etaBinEdge+1) <<" phibins: "<< phiBinStep*phiBinEdge <<", "<<phiBinStep*(phiBinEdge+1) <<endl;
// 			    }
// 			}
// 		    }
// 		}
// 	    }
	  
	  
// 	  /*
	    
// 	  // Loop over all jet's constituents
// 	  for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
// 	  {
// 	  object = jet->Constituents.At(j);
	  
// 	  // Check if the constituent is accessible
// 	  if(object == 0) continue;
	  
// 	  if(object->IsA() == GenParticle::Class())
// 	  {
// 	  particle = (GenParticle*) object;
// 	  // cout << "    GenPart pt: " << particle->PT << ", eta: " << particl	\
// 	  particle->Eta << ", phi: " << particle->Phi << endl;
// 	  momentum += particle->P4();
// 	  }
// 	  else if(object->IsA() == Tower::Class())
// 	  {
// 	  tower = (Tower*) object;
// 	  //		   cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta <\< ", phi: " << tower->Phi << endl;
// 	  momentum += tower->P4();
// 	  }
// 	  }*/
// 	  //plots->fJetDeltaPT->Fill((jet->PT - momentum.Pt())/jet->PT);
// 	}
//     }
//   cout<<"Jet numbers:  "<<count_jet<<" Tower numbers: "<<count_tower<<endl;
// }






void MakeJetGrids_jh( const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  AnalyseEvents(treeReader);
}