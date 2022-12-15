#define VLLAna_Trigger_cxx
// The class definition in VLLAna_Trigger.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("VLLAna_Trigger.C")
// root> T->Process("VLLAna_Trigger.C","some options")
// root> T->Process("VLLAna_Trigger.C+")
//


#include "VLLAna_Trigger.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <iostream>
#include <iomanip>

void VLLAna_Trigger::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   nEvtTotal          =0;
   nEvtTrigger        =0;
   nEvtRan            =0;
   //Create the histogram file
   _HstFile = new TFile(_HstFileName,"recreate");
   //Call the function to book the histograms we declared in Hists.
   BookHistograms();
}

void VLLAna_Trigger::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void VLLAna_Trigger::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  _HstFile->Write();
  _HstFile->Close();
  
  //Output to screen
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total triggered events ran = "<<nEvtTrigger<<endl;
  cout<<"Total events  = "<<nEvtTotal<<endl;

}

void VLLAna_Trigger::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

Bool_t VLLAna_Trigger::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);
  
  if(_verbosity==0 && nEvtTotal%1000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  
  nEvtTotal++;
  h.nevt[0]->Fill(0);

  h.metfilter[0]->Fill(*Flag_goodVertices);
  h.metfilter[1]->Fill(*Flag_globalSuperTightHalo2016Filter);
  h.metfilter[2]->Fill(*Flag_HBHENoiseFilter);
  h.metfilter[3]->Fill(*Flag_HBHENoiseIsoFilter);
  h.metfilter[4]->Fill(*Flag_EcalDeadCellTriggerPrimitiveFilter);
  h.metfilter[5]->Fill(*Flag_BadPFMuonFilter);
  h.metfilter[6]->Fill(*Flag_eeBadScFilter);
  
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  h.metfilter[7]->Fill(GoodEvt2016);
  h.metfilter[8]->Fill(GoodEvt2017);
  h.metfilter[9]->Fill(GoodEvt2018);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;

  trigger2018 = (_year==2018 ? (_lep==1 ? *HLT_IsoMu24==1 : _lep==0 && *HLT_Ele32_WPTight_Gsf) : 1);
  trigger2017 = (_year==2017 ? (_lep==1 ? *HLT_IsoMu27==1 : _lep==0 && (*HLT_Ele32_WPTight_Gsf)) : 1);
  trigger2016 = (_year==2016 ? (_lep==1 ? (*HLT_IsoMu24==1) : _lep==0 && *HLT_Ele27_WPTight_Gsf) : 1);
  triggerRes = trigger2018 && trigger2017 && trigger2016;
  h.metfilter[10]->Fill(GoodEvt);
  
  if(GoodEvt){
    nEvtRan++;
    if(triggerRes){
      nEvtTrigger++;
      
      /****************************************************************
       *                     Trigger Object                           *
       ****************************************************************/  
      trigObj.clear();
      h.nevt[1]->Fill(*nTrigObj);
      for(unsigned int i =0; i<(*nTrigObj); i++){
	h.trigObject[0]->Fill(TrigObj_pt[i]);
	h.trigObject[1]->Fill(TrigObj_eta[i]);
	h.trigObject[2]->Fill(TrigObj_phi[i]);
	h.trigObject[3]->Fill(TrigObj_id[i]);
	if(fabs(TrigObj_id[i])==13){
	  h.nevt[2]->Fill(*nTrigObj);
	  h.trigObject[4]->Fill(TrigObj_pt[i]);
	  h.trigObject[5]->Fill(TrigObj_eta[i]);
	  h.trigObject[6]->Fill(TrigObj_phi[i]);
	}
	if(fabs(TrigObj_id[i])==11){
	  h.nevt[3]->Fill(*nTrigObj);
	  h.trigObject[7]->Fill(TrigObj_pt[i]);
	  h.trigObject[8]->Fill(TrigObj_eta[i]);
	  h.trigObject[9]->Fill(TrigObj_phi[i]);
	}
      }
      
      
      /****************************************************************
       *                          Muons                               *
       ****************************************************************/
  
      tagMu.clear();
      probeMu.clear();
      goodMu.clear();
      
      for(unsigned int i=0; i<(*nMuon); i++){
	Lepton temp;
	temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105);
	temp.id = -13*Muon_charge[i];
	temp.ind = i;
	temp.charge = Muon_charge[i];
	temp.dz=Muon_dz[i];
	temp.dxy=Muon_dxy[i];
	
	//Creating array for the good muon
	bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_pfRelIso04_all[i]<0.15 && Muon_mediumId[i] && fabs(temp.dxy)<0.05 && fabs(temp.dz)<0.1;
	if(passCuts){
	  goodMu.push_back(temp);
	}
	SortPt(0);
      }
      
      //Muon Array
      h.goodmuons[0]->Fill((int)goodMu.size());
      for(int i=0; i<(int)goodMu.size(); i++){
	h.goodmuons[1]->Fill(goodMu.at(i).v.Pt());
	h.goodmuons[2]->Fill(goodMu.at(i).v.Eta());
	h.goodmuons[3]->Fill(goodMu.at(i).v.Phi());
      }
      
      float dRmin_lead = 999.0;
      float dRmin_sublead = 999.0;
      float dR_leadMuon_subleadMuon = 999.0;
      unsigned int index_1 = 1;
      unsigned int index_2 = 1;
      
      //dRMin between Leading muon and trigger object
      if((int)goodMu.size()==2 && goodMu.at(0).v.Pt()>24){
	for(unsigned int j =0; j<(*nTrigObj); j++){
	  if(TrigObj_id[j]==13){
	    float dR_leadMuon_trigObject = sqrt(pow(delta_phi(TrigObj_phi[j],goodMu.at(0).v.Phi()),2)+pow(delta_Eta(TrigObj_eta[j],goodMu.at(0).v.Eta()),2));
	    h.dR[2]->Fill(dR_leadMuon_trigObject);
	    if(dR_leadMuon_trigObject<dRmin_lead){
	      dRmin_lead = dR_leadMuon_trigObject;
	      index_1=j;
	    }
	  }
	}
	h.dR[0]->Fill(dRmin_lead);
      }
      
      
      //dRMin between sub leading muon and trigger object
      if((int)goodMu.size()==2){
	for(unsigned int j =0; j<(*nTrigObj); j++){
	  if(TrigObj_id[j]==13 && j!= index_1){
	    float dR_subleadMuon_trigObject = sqrt(pow(delta_phi(TrigObj_phi[j],goodMu.at(1).v.Phi()),2)+pow(delta_Eta(TrigObj_eta[j],goodMu.at(1).v.Eta()),2));
	    h.dR[3]->Fill(dR_subleadMuon_trigObject);
	    if(dR_subleadMuon_trigObject<dRmin_sublead){
	      dRmin_sublead = dR_subleadMuon_trigObject;
	      index_2=j;
	    }
	  }
	}
	h.dR[1]->Fill(dRmin_sublead);
      }
      
      //dR matching between leading and subleading muon
      bool isanyclose = false;
      if((int)goodMu.size()==2){
	if((goodMu.at(0).v.DeltaR(goodMu.at(1).v))<0.4)
	  isanyclose = true;
      }

      //Defining two regions barrel and endcap
      bool isbarrel = false;
      bool isendcap = false;
      if((int)goodMu.size()==2){
	for(int i = 0; i<(int)goodMu.size(); i++){
	  if(fabs(goodMu.at(i).v.Eta()) < 1.479)
	    isbarrel = true;
	  if(fabs(goodMu.at(i).v.Eta()) > 1.479 && fabs(goodMu.at(i).v.Eta()) < 2.4)
	    isendcap = true;
	}
      }

      //Defining Numerator and denominator histograms for MC sample
      if(_data == 0 && (int)goodMu.size() == 2){
	float invMass_leadMu_subleadMu = (goodMu.at(0).v+goodMu.at(1).v).M();
	if(invMass_leadMu_subleadMu > 76 && invMass_leadMu_subleadMu < 105){
	  h.test[0]->Fill(invMass_leadMu_subleadMu);
	  if(dRmin_lead<0.2 && goodMu.at(0).v.Pt()>24 && !isanyclose && goodMu.at(0).id*goodMu.at(1).id == -169 && isbarrel){
	    h.goodmuons[4]->Fill(goodMu.at(1).v.Pt());
	    h.goodmuons[12]->Fill(goodMu.at(1).v.Pt());
	    if(dRmin_sublead<0.2 && TrigObj_pt[index_2]>24){
	      h.goodmuons[5]->Fill(goodMu.at(1).v.Pt());
	      h.goodmuons[13]->Fill(goodMu.at(1).v.Pt());
	    }
	  }
	  if(dRmin_lead<0.2 && goodMu.at(0).v.Pt()>24 && !isanyclose && goodMu.at(0).id*goodMu.at(1).id == -169 && isendcap){
	    h.goodmuons[6]->Fill(goodMu.at(1).v.Pt());
	    h.goodmuons[14]->Fill(goodMu.at(1).v.Pt());
	    if(dRmin_sublead<0.2 && TrigObj_pt[index_2]>24){
	      h.goodmuons[7]->Fill(goodMu.at(1).v.Pt());
	      h.goodmuons[15]->Fill(goodMu.at(1).v.Pt());
	    }
	  } 
	}
      }
      
      if(_data == 1 && (int)goodMu.size() == 2){
	float invMass_leadMu_subleadMu_data = (goodMu.at(0).v+goodMu.at(1).v).M();
	float dphi_leadMu_subleadMu = delta_phi(goodMu.at(0).v.Phi(),goodMu.at(1).v.Phi());
	h.test[2]->Fill(dphi_leadMu_subleadMu);
	if(invMass_leadMu_subleadMu_data>76 && invMass_leadMu_subleadMu_data<105){
	  h.test[1]->Fill(invMass_leadMu_subleadMu_data);
	  if(dRmin_lead<0.2 && goodMu.at(0).v.Pt()>24 && !isanyclose && goodMu.at(0).id*goodMu.at(1).id == -169 && isbarrel){
	    h.goodmuons[8]->Fill(goodMu.at(1).v.Pt());
	    h.goodmuons[16]->Fill(goodMu.at(1).v.Pt());
	    if(dRmin_sublead<0.2 && TrigObj_pt[index_2]>24){
	      h.goodmuons[9]->Fill(goodMu.at(1).v.Pt());
	      h.goodmuons[17]->Fill(goodMu.at(1).v.Pt());
	    }
	  }
	  if(dRmin_lead<0.2 && goodMu.at(0).v.Pt()>24 && !isanyclose && goodMu.at(0).id*goodMu.at(1).id == -169 && isendcap){
	    h.goodmuons[10]->Fill(goodMu.at(1).v.Pt());
	    h.goodmuons[18]->Fill(goodMu.at(1).v.Pt());
	    if(dRmin_sublead<0.2 && TrigObj_pt[index_2]>24){
	      h.goodmuons[11]->Fill(goodMu.at(1).v.Pt());
	      h.goodmuons[19]->Fill(goodMu.at(1).v.Pt());
	    }
	  }
	}
      }

      /****************************************************************
       *                          Electrons                           *
       ****************************************************************/
  
      tagEle.clear();
      probeEle.clear();
      goodEle.clear();
      
      for(unsigned int i=0; i<(*nElectron); i++){
	Lepton temp;
	temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511);
	temp.id = -11*Electron_charge[i];
	temp.ind = i;
	temp.charge = Electron_charge[i];
	temp.dz=Electron_dz[i];
	temp.dxy=Electron_dxy[i];

	bool isprompt = false;
	if(fabs(temp.v.Eta())<=1.479)
	  if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	    isprompt = true;      
	if(fabs(temp.v.Eta())>1.479)
	  if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	    isprompt = true;
      	
	//Creating array for the good electron
	bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Electron_cutBased[i]>2;
	passCuts = passCuts && isprompt;
	if(passCuts){
	  goodEle.push_back(temp);
	}
	SortPt(1);
      }
      
      //Electron Array
      h.goodelectrons[0]->Fill((int)goodEle.size());
      for(int i=0; i<(int)goodEle.size(); i++){
	h.goodelectrons[1]->Fill(goodEle.at(i).v.Pt());
	h.goodelectrons[2]->Fill(goodEle.at(i).v.Eta());
	h.goodelectrons[3]->Fill(goodEle.at(i).v.Phi());
      }
      
      float dRmin_leadele = 999.0;
      float dRmin_subleadele = 999.0;
      float dR_leadele_subleadele = 999.0;
      unsigned int index_1e = 1;
      unsigned int index_2e = 1;
      
      //dRMin between Leading electron and trigger object
      if((int)goodEle.size()==2 && goodEle.at(0).v.Pt()>32){
	for(unsigned int j =0; j<(*nTrigObj); j++){
	  if(TrigObj_id[j]==11){
	    float dR_leadele_trigObject = sqrt(pow(delta_phi(TrigObj_phi[j],goodEle.at(0).v.Phi()),2)+pow(delta_Eta(TrigObj_eta[j],goodEle.at(0).v.Eta()),2));
	    h.dR[6]->Fill(dR_leadele_trigObject);
	    if(dR_leadele_trigObject<dRmin_leadele){
	      dRmin_leadele = dR_leadele_trigObject;
	      index_1e=j;
	    }
	  }
	}
	h.dR[4]->Fill(dRmin_leadele);
      }
      
      
      //dRMin between sub leading electron and trigger object
      if((int)goodEle.size()==2){
	for(unsigned int j =0; j<(*nTrigObj); j++){
	  if(TrigObj_id[j]==11 && j!= index_1e){
	    float dR_subleadele_trigObject = sqrt(pow(delta_phi(TrigObj_phi[j],goodEle.at(1).v.Phi()),2)+pow(delta_Eta(TrigObj_eta[j],goodEle.at(1).v.Eta()),2));
	    h.dR[7]->Fill(dR_subleadele_trigObject);
	    if(dR_subleadele_trigObject<dRmin_subleadele){
	      dRmin_subleadele = dR_subleadele_trigObject;
	      index_2e=j;
	    }
	  }
	}
	h.dR[5]->Fill(dRmin_subleadele);
      }

      //dR matching between leading and subleading electron
      bool isanyeleclose = false;
      if((int)goodEle.size()==2){
	if((goodEle.at(0).v.DeltaR(goodEle.at(1).v))<0.4)
	  isanyeleclose = true;
      }

      //Defining two regions barrel and endcap
      bool isbarrelele = false;
      bool isendcapele = false;
      if((int)goodEle.size()==2){
	for(int i = 0; i<(int)goodEle.size(); i++){
	  if(fabs(goodEle.at(i).v.Eta()) < 1.479)
	    isbarrelele = true;
	  if(fabs(goodEle.at(i).v.Eta()) > 1.479 && fabs(goodEle.at(i).v.Eta()) < 2.4)
	    isendcapele = true;
	}
      }

      //Defining Numerator and denominator histograms for MC sample
      if(_data == 0 && (int)goodEle.size() == 2){
	float invMass_leadele_subleadele = (goodEle.at(0).v+goodEle.at(1).v).M();
	if(invMass_leadele_subleadele > 76 && invMass_leadele_subleadele < 105){
	  h.test[3]->Fill(invMass_leadele_subleadele);
	  if(dRmin_leadele<0.2 && goodEle.at(0).v.Pt()>32 && !isanyeleclose && goodEle.at(0).id*goodEle.at(1).id == -121 && isbarrelele){
	    h.goodelectrons[4]->Fill(goodEle.at(1).v.Pt());
	    h.goodelectrons[12]->Fill(goodEle.at(1).v.Pt());
	    if(dRmin_subleadele<0.2 && TrigObj_pt[index_2e]>32){
	      h.goodelectrons[5]->Fill(goodEle.at(1).v.Pt());
	      h.goodelectrons[13]->Fill(goodEle.at(1).v.Pt());
	    }
	  }
	  if(dRmin_leadele<0.2 && goodEle.at(0).v.Pt()>32 && !isanyeleclose && goodEle.at(0).id*goodEle.at(1).id == -121 && isendcapele){
	    h.goodelectrons[6]->Fill(goodEle.at(1).v.Pt());
	    h.goodelectrons[14]->Fill(goodEle.at(1).v.Pt());
	    if(dRmin_subleadele<0.2 && TrigObj_pt[index_2e]>32){
	      h.goodelectrons[7]->Fill(goodEle.at(1).v.Pt());
	      h.goodelectrons[15]->Fill(goodEle.at(1).v.Pt());
	    }
	  } 
	}
      }
      
      if(_data == 1 && (int)goodEle.size() == 2){
	float invMass_leadele_subleadele_data = (goodEle.at(0).v+goodEle.at(1).v).M();
	float dphi_leadele_subleadele = delta_phi(goodEle.at(0).v.Phi(),goodEle.at(1).v.Phi());
	h.test[5]->Fill(dphi_leadele_subleadele);
	if(invMass_leadele_subleadele_data>76 && invMass_leadele_subleadele_data<105){
	  h.test[4]->Fill(invMass_leadele_subleadele_data);
	  if(dRmin_leadele<0.2 && goodEle.at(0).v.Pt()>32 && !isanyeleclose && goodEle.at(0).id*goodEle.at(1).id == -121 && isbarrelele){
	    h.goodelectrons[8]->Fill(goodEle.at(1).v.Pt());
	    h.goodelectrons[16]->Fill(goodEle.at(1).v.Pt());
	    if(dRmin_subleadele<0.2 && TrigObj_pt[index_2e]>32){
	      h.goodelectrons[9]->Fill(goodEle.at(1).v.Pt());
	      h.goodelectrons[17]->Fill(goodEle.at(1).v.Pt());
	    }
	  }
	  if(dRmin_leadele<0.2 && goodEle.at(0).v.Pt()>32 && !isanyeleclose && goodEle.at(0).id*goodEle.at(1).id == -121 && isendcapele){
	    h.goodelectrons[10]->Fill(goodEle.at(1).v.Pt());
	    h.goodelectrons[18]->Fill(goodEle.at(1).v.Pt());
	    if(dRmin_subleadele<0.2 && TrigObj_pt[index_2e]>32){
	      h.goodelectrons[11]->Fill(goodEle.at(1).v.Pt());
	      h.goodelectrons[19]->Fill(goodEle.at(1).v.Pt());
	    }
	  }
	}
      }



      
      
      
    }
  }
  
  return kTRUE;
  
}
void VLLAna_Trigger::BookHistograms()
{
  //Event Hists
  h.nevt[0]       = new TH1F("nEvents","Total events",4,0,4);
  h.nevt[1]       = new TH1F("nTrigObj","nTrigObj",10,0,10);
  h.nevt[2]       = new TH1F("nTrigObj_Id13","nTrigObj_Id13",10,0,10);
  h.nevt[3]       = new TH1F("nTrigObj_Id11","nTrigObj_Id11",10,0,10);
  for(int i =0; i<4; i++){
    h.nevt[i]->Sumw2();
  }
  
  //MET Filter
  h.metfilter[0]  = new TH1F("METfilter_goodVertices","METfilter_goodVertices",5,-1,4);
  h.metfilter[1]  = new TH1F("METfilter_globalSuperTightHalo2016Filter","METfilter_globalSuperTightHalo2016Filter",5,-1,4);
  h.metfilter[2]  = new TH1F("METfilter_HBHENoiseFilter","METfilter_HBHENoiseFilter",5,-1,4);
  h.metfilter[3]  = new TH1F("METfilter_HBHENoiseIsoFilter","METfilter_HBHENoiseIsoFilter",5,-1,4);
  h.metfilter[4]  = new TH1F("METfilter_EcalDeadCellTriggerPrimitiveFilter","METfilter_EcalDeadCellTriggerPrimitiveFilter",5,-1,4);
  h.metfilter[5]  = new TH1F("METfilter_BadPFMuonFilter","METfilter_BadPFMuonFilter",5,-1,4);
  h.metfilter[6]  = new TH1F("METfilter_eeBadScFilter","METfilter_eeBadScFilter",5,-1,4);
  h.metfilter[7]  = new TH1F("METfilter_GoodEvt2016","METfilter_GoodEvt2016",5,-1,4);
  h.metfilter[8]  = new TH1F("METfilter_GoodEvt2017","METfilter_GoodEvt2017",5,-1,4);
  h.metfilter[9]  = new TH1F("METfilter_GoodEvt2018","METfilter_GoodEvt2018",5,-1,4);
  h.metfilter[10] = new TH1F("METfilter_GoodEvt","METfilter_GoodEvt",5,-1,4);
  for(int i =0; i<11; i++){
    h.metfilter[i]->Sumw2();
  }
  //Trigger Object
  h.trigObject[0]  = new TH1F("TrigObject_pT","TrigObject_pT",1000,0,1000);
  h.trigObject[1]  = new TH1F("TrigObject_eta","TrigObject_eta",60,-3,3);
  h.trigObject[2]  = new TH1F("TrigObject_phi","TrigObject_phi",60,-3,3);
  h.trigObject[3]  = new TH1F("TrigObject_id","TrigObject_id",20,0,20);
  h.trigObject[4]  = new TH1F("TrigObject_pT_id13","TrigObject_pT_id13",1000,0,1000);
  h.trigObject[5]  = new TH1F("TrigObject_eta_id13","TrigObject_eta_id13",60,-3,3);
  h.trigObject[6]  = new TH1F("TrigObject_phi_id13","TrigObject_phi_id13",60,-3,3);
  h.trigObject[7]  = new TH1F("TrigObject_pT_id11","TrigObject_pT_id11",1000,0,1000);
  h.trigObject[8]  = new TH1F("TrigObject_eta_id11","TrigObject_eta_id11",60,-3,3);
  h.trigObject[9]  = new TH1F("TrigObject_phi_id11","TrigObject_phi_id11",60,-3,3);
  for(int i =0; i<10; i++){
    h.trigObject[i]->Sumw2();
  }
  //Muons
  h.goodmuons[0]       = new TH1F("ngoodMuons","ngoodMuons", 10, 0, 10);
  h.goodmuons[1]       = new TH1F("goodMuon_pT","goodMuon_pT",1000,0,1000);
  h.goodmuons[2]       = new TH1F("goodMuon_Eta","goodMuon_Eta",60,-3,3);
  h.goodmuons[3]       = new TH1F("goodMuon_Phi","goodMuon_Phi",60,-3,3);

  //Custom Binning
  float x_bin_pT[11] = {10,20,22,24,26,28,30,40,60,100,250};
  h.goodmuons[4]       = new TH1F("AllprobeMuon_barrelpT","AllprobeMuon_barrelpT",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);
  h.goodmuons[5]       = new TH1F("AllprobeMuon_barrelpT_trigObjectMatched","AllprobeMuon_barrelpT_trigObjectMatched",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);
  h.goodmuons[6]       = new TH1F("AllprobeMuon_endcappT","AllprobeMuon_endcappT",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);
  h.goodmuons[7]       = new TH1F("AllprobeMuon_endcappT_trigObjectMatched","AllprobeMuon_endcappT_trigObjectMatched",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);
  h.goodmuons[8]       = new TH1F("AllprobeMuon_barrelpT_data","AllprobeMuon_barrelpT_data",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);
  h.goodmuons[9]       = new TH1F("AllprobeMuon_barrelpT_trigObjectMatched_data","AllprobeMuon_barrelpT_trigObjectMatched_data",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);

  h.goodmuons[10]       = new TH1F("AllprobeMuon_endcappT_data","AllprobeMuon_endcappT_data",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);
  h.goodmuons[11]       = new TH1F("AllprobeMuon_endcappT_trigObjectMatched_data","AllprobeMuon_endcappT_trigObjectMatched_data",(sizeof(x_bin_pT)/ sizeof(x_bin_pT[0])-1),x_bin_pT);

  //Normal Binning
  h.goodmuons[12]       = new TH1F("AllprobeMuon_barrelpT_N","AllprobeMuon_barrelpT_N",500,0,500);
  h.goodmuons[13]       = new TH1F("AllprobeMuon_barrelpT_trigObjectMatched_N","AllprobeMuon_barrelpT_trigObjectMatched_N",500,0,500);
  h.goodmuons[14]       = new TH1F("AllprobeMuon_endcappT_N","AllprobeMuon_endcappT_N",500,0,500);
  h.goodmuons[15]       = new TH1F("AllprobeMuon_endcappT_trigObjectMatched_N","AllprobeMuon_endcappT_trigObjectMatched_N",500,0,500);
  h.goodmuons[16]       = new TH1F("AllprobeMuon_barrelpT_data_N","AllprobeMuon_barrelpT_data_N",500,0,500);
  h.goodmuons[17]       = new TH1F("AllprobeMuon_barrelpT_trigObjectMatched_data_N","AllprobeMuon_barrelpT_trigObjectMatched_data_N",500,0,500);
  h.goodmuons[18]       = new TH1F("AllprobeMuon_endcappT_data_N","AllprobeMuon_endcappT_data_N",500,0,500);
  h.goodmuons[19]       = new TH1F("AllprobeMuon_endcappT_trigObjectMatched_data_N","AllprobeMuon_endcappT_trigObjectMatched_data",500,0,500);

  for(int i =0; i<20; i++){
    h.goodmuons[i]->Sumw2();
  }

  //Electrons
  h.goodelectrons[0]       = new TH1F("ngoodElectrons","ngoodElectrons", 10, 0, 10);
  h.goodelectrons[1]       = new TH1F("goodElectron_pT","goodElectron_pT",1000,0,1000);
  h.goodelectrons[2]       = new TH1F("goodElectron_Eta","goodElectron_Eta",60,-3,3);
  h.goodelectrons[3]       = new TH1F("goodElectron_Phi","goodElectron_Phi",60,-3,3);

  //Custom Binning
  float x_bin_elepT[14] = {10,20,24,28,30,32,34,36,38,40,60,100,150,250};
  h.goodelectrons[4]       = new TH1F("AllprobeElectron_barrelpT","AllprobeElectron_barrelpT",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[5]       = new TH1F("AllprobeElectron_barrelpT_trigObjectMatched","AllprobeElectron_barrelpT_trigObjectMatched",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[6]       = new TH1F("AllprobeElectron_endcappT","AllprobeElectron_endcappT",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[7]       = new TH1F("AllprobeElectron_endcappT_trigObjectMatched","AllprobeElectron_endcappT_trigObjectMatched",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[8]       = new TH1F("AllprobeElectron_barrelpT_data","AllprobeElectron_barrelpT_data",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[9]       = new TH1F("AllprobeElectron_barrelpT_trigObjectMatched_data","AllprobeElectron_barrelpT_trigObjectMatched_data",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[10]       = new TH1F("AllprobeElectron_endcappT_data","AllprobeElectron_endcappT_data",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);
  h.goodelectrons[11]       = new TH1F("AllprobeElectron_endcappT_trigObjectMatched_data","AllprobeElectron_endcappT_trigObjectMatched_data",(sizeof(x_bin_elepT)/ sizeof(x_bin_elepT[0])-1),x_bin_elepT);

  //Normal Binning
  h.goodelectrons[12]       = new TH1F("AllprobeElectron_barrelpT_N","AllprobeElectron_barrelpT_N",500,0,500);
  h.goodelectrons[13]       = new TH1F("AllprobeElectron_barrelpT_trigObjectMatched_N","AllprobeElectron_barrelpT_trigObjectMatched_N",500,0,500);
  h.goodelectrons[14]       = new TH1F("AllprobeElectron_endcappT_N","AllprobeElectron_endcappT_N",500,0,500);
  h.goodelectrons[15]       = new TH1F("AllprobeElectron_endcappT_trigObjectMatched_N","AllprobeElectron_endcappT_trigObjectMatched_N",500,0,500);
  h.goodelectrons[16]       = new TH1F("AllprobeElectron_barrelpT_data_N","AllprobeElectron_barrelpT_data_N",500,0,500);
  h.goodelectrons[17]       = new TH1F("AllprobeElectron_barrelpT_trigObjectMatched_data_N","AllprobeElectron_barrelpT_trigObjectMatched_data_N",500,0,500);
  h.goodelectrons[18]       = new TH1F("AllprobeElectron_endcappT_data_N","AllprobeElectron_endcappT_data_N",500,0,500);
  h.goodelectrons[19]       = new TH1F("AllprobeElectron_endcappT_trigObjectMatched_data_N","AllprobeElectron_endcappT_trigObjectMatched_data",500,0,500);

  for(int i =0; i<20; i++){
    h.goodelectrons[i]->Sumw2();
  }
  
  //dR
  h.dR[0]        = new TH1F("dRmin_LeadMuon_trigObject","dRmin_LeadMuon_trigObject",120,0,0.3);
  h.dR[1]        = new TH1F("dRmin_SubLeadMuon_trigObject","dRmin_SubLeadMuon_trigObject",120,0,0.3);
  h.dR[2]        = new TH1F("dR_LeadMuon_trigObject","dR_LeadMuon_trigObject",120,0,0.3);
  h.dR[3]        = new TH1F("dR_SubLeadMuon_trigObject","dR_SubLeadMuon_trigObject",120,0,0.3);
  h.dR[4]        = new TH1F("dRmin_Leadele_trigObject","dRmin_Leadele_trigObject",120,0,0.3);
  h.dR[5]        = new TH1F("dRmin_SubLeadele_trigObject","dRmin_SubLeadele_trigObject",120,0,0.3);
  h.dR[6]        = new TH1F("dR_Leadele_trigObject","dR_Leadele_trigObject",120,0,0.3);
  h.dR[7]        = new TH1F("dR_SubLeadele_trigObject","dR_SubLeadele_trigObject",120,0,0.3);
  
  for(int i =0; i<8; i++){
    h.dR[i]->Sumw2();
  }
  //Test Plots
  h.test[0]      = new TH1F("InvMass_Muons_MC","InvMass_Muons_MC",40,70,110);
  h.test[1]      = new TH1F("InvMass_Muons_data","InvMass_Muons_data",200,0,200);
  h.test[2]      = new TH1F("dPhi_Muons_data","dPhi_Muons_data",80,-4,4);
  h.test[3]      = new TH1F("InvMass_Electrons_MC","InvMass_Electrons_MC",40,70,110);
  h.test[4]      = new TH1F("InvMass_Electrons_data","InvMass_Electrons_data",200,0,200);
  h.test[5]      = new TH1F("dPhi_Electrons_data","dPhi_Electrons_data",80,-4,4);
  for(int i =0; i<6; i++){
    h.test[i]->Sumw2();
  }
}


float VLLAna_Trigger::delta_phi(float phi1, float phi2)
{
  //Calculate the correct deltaPhi=phi1-phi2
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float VLLAna_Trigger::delta_Eta(float eta1, float eta2)
{
  //Calculate the correct deltaPhi=phi1-phi2
  float deta = fabs(fabs(eta1) - fabs(eta2));
  return deta;
}

void VLLAna_Trigger::SortPt(int opt){

  if(opt==0){
    for(int i=0; i<(int)goodMu.size()-1; i++){
      for(int j=i+1; j<(int)goodMu.size(); j++){
	if( goodMu[i].v.Pt() < goodMu[j].v.Pt() ) swap(goodMu.at(i),goodMu.at(j));
      }
    }
  }
  if(opt==1){
    for(int i=0; i<(int)goodEle.size()-1; i++){
      for(int j=i+1; j<(int)goodEle.size(); j++){
	if( goodEle[i].v.Pt() < goodEle[j].v.Pt() ) swap(goodEle.at(i),goodEle.at(j));
      }
    }
  }
  
}

