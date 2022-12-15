#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
using std::cout;
using std::endl;

void decorate(TGraph *h, TString gtitle, float gmax, float gmin, 
	      int markercolor, int markerstyle, int linecolor, int linewidth);
void decorateHist(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill);
void decorate(TLegend *g, float textSize, TString legendheader);
float get_nevents(TH1F *hst, float bin_lo, float bin_hi);
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi);
void SetOverflowBin(TH1F *histo);

void overlayMuon()
{
  //TFile *fdy   = new TFile("hst_DY.root");
  //TFile *fdata = new TFile("hst_SingleMuon.root");
  TFile *fdy   = new TFile("/home/ykumar/Yash/Work/VLLAnalysis_2022/Trigger_Efficiency/hstfiles/Muons/VLL2018_1L2JAnalysis_TriggerEff_Dec12_v1/VLL2018_1L2JAnalysis_TriggerEff_Dec12_v1_DYJetsToLL_M50_sample.root");
  TFile *fdata = new TFile("/home/ykumar/Yash/Work/VLLAnalysis_2022/Trigger_Efficiency/hstfiles/Muons/VLL2018_1L2JAnalysis_TriggerEff_Dec12_v1/hst_SingleMuon.root");
  
  // Determine the bin-edges; efficiency will be calculated in these bins
  float lowed[] =  {10,20,22,24,26,28,30,40,60,100};
  float hied[]  =  {20,22,24,26,28,30,40,60,100,250};

  int nbin = 10; 
  TString plots_mc[4] = {"AllprobeMuon_barrelpT",
			 "AllprobeMuon_barrelpT_trigObjectMatched",
			 "AllprobeMuon_endcappT",
			 "AllprobeMuon_endcappT_trigObjectMatched",
  };
  TString plots_data[4] = {"AllprobeMuon_barrelpT_data",
			   "AllprobeMuon_barrelpT_trigObjectMatched_data",
			   "AllprobeMuon_endcappT_data",
			   "AllprobeMuon_endcappT_trigObjectMatched_data",
  };

  TString plotname1, plotname2;
  TH1F *hmc[50],*hdata[50];
  vector<float> binwiseError_hmc[5], binwiseIntegral_hmc[5], binwiseError_hdata[5], binwiseIntegral_hdata[5], multiplicativeErrorpropagator[5];

  for(int i=0; i<4; i++){
    //Cleaning arrays
    binwiseError_hmc[i].clear();
    binwiseIntegral_hmc[i].clear();
    binwiseError_hdata[i].clear();
    binwiseIntegral_hdata[i].clear();
    //Defining plotnames
    plotname1 = plots_mc[i];
    plotname2 = plots_data[i];
    //Now open the respective histograms from the file
    hmc[i] = (TH1F*)fdy ->Get(plotname1);
    hdata[i] = (TH1F*)fdata ->Get(plotname2);
    //Setting overflow and ubnderflow bin
    SetOverflowBin(hmc[i]);
    SetOverflowBin(hdata[i]);
    //Calculating bin content and bin error
    for(int j =0; j<(hmc[i]->GetNbinsX()); j++){
      binwiseIntegral_hmc[i].push_back(hmc[i]->GetBinContent(j));
      binwiseError_hmc[i].push_back(hmc[i]->GetBinError(j));
    }
    for(int k =0; k<(hdata[i]->GetNbinsX()); k++){
      binwiseIntegral_hdata[i].push_back(hdata[i]->GetBinContent(k));
      binwiseError_hdata[i].push_back(hdata[i]->GetBinError(k));
    }
  }

  for(int k = 0; k<(float)binwiseError_hmc[0].size(); k++){
    multiplicativeErrorpropagator[0].push_back(sqrt((pow(binwiseError_hmc[0][k],2)/pow(binwiseIntegral_hmc[0][k],2))+(pow(binwiseError_hmc[1][k],2)/pow(binwiseIntegral_hmc[1][k],2))));
    multiplicativeErrorpropagator[1].push_back(sqrt((pow(binwiseError_hmc[2][k],2)/pow(binwiseIntegral_hmc[2][k],2))+(pow(binwiseError_hmc[3][k],2)/pow(binwiseIntegral_hmc[3][k],2))));
  }
    
  for(int k = 0; k<(float)binwiseError_hdata[0].size(); k++){
    multiplicativeErrorpropagator[2].push_back(sqrt((pow(binwiseError_hdata[0][k],2)/pow(binwiseIntegral_hdata[0][k],2))+(pow(binwiseError_hdata[1][k],2)/pow(binwiseIntegral_hdata[1][k],2))));
    multiplicativeErrorpropagator[3].push_back(sqrt((pow(binwiseError_hdata[2][k],2)/pow(binwiseIntegral_hdata[2][k],2))+(pow(binwiseError_hdata[3][k],2)/pow(binwiseIntegral_hdata[3][k],2))));
  }
  
  //Declare some arrays to store the efficiency
  int n = 0;
  float x_b[150],y_b_endcap[150],y_b_barrel[150],y_b_endcap_data[150],y_b_barrel_data[150],ex[150],ey_barrel[150],ey_endcap[150],ey_barrel_data[150],ey_endcap_data[150];
  double et_lo,et_hi;
  float et_step=0;
  for(int i=0; i<150; i++){
    x_b[i]=y_b_barrel[i]=y_b_endcap[i]=y_b_barrel_data[i]=y_b_endcap_data[i]=ex[i]=ey_barrel[i]=ey_endcap[i]=ey_barrel_data[i]=ey_endcap_data[i]=0; 
  }

  for(int i=0; i<nbin; i++){
    et_lo = lowed[i];
    et_hi = hied[i];
    et_step = et_hi - et_lo;
    
    float denValue_barrel = get_nevents(hmc[0],et_lo,et_hi);
    float numValue_barrel = get_nevents(hmc[1],et_lo,et_hi);
    float denValue_endcap = get_nevents(hmc[2],et_lo,et_hi);
    float numValue_endcap = get_nevents(hmc[3],et_lo,et_hi);
    float denValue_barrel_data = get_nevents(hdata[0],et_lo,et_hi);
    float numValue_barrel_data = get_nevents(hdata[1],et_lo,et_hi);
    float denValue_endcap_data = get_nevents(hdata[2],et_lo,et_hi);
    float numValue_endcap_data = get_nevents(hdata[3],et_lo,et_hi);
    
    x_b[n] = (et_lo+et_hi)/2;
    ex[n] = (et_hi - et_lo)/2.;
        
    if(numValue_barrel>0 && denValue_barrel>0){
      y_b_barrel[n]=numValue_barrel/denValue_barrel;
      ey_barrel[n] = (numValue_barrel/denValue_barrel)*multiplicativeErrorpropagator[0][i];
    }
    if(numValue_endcap>0 && denValue_endcap>0){
      y_b_endcap[n]=numValue_endcap/denValue_endcap;
      ey_endcap[n] = (numValue_endcap/denValue_endcap)*multiplicativeErrorpropagator[1][i];
    }
    if(numValue_barrel_data>0 && denValue_barrel_data>0){
      y_b_barrel_data[n]=numValue_barrel_data/denValue_barrel_data;
      ey_barrel_data[n] = (numValue_barrel_data/denValue_barrel_data)*multiplicativeErrorpropagator[2][i];
    }
    if(numValue_endcap_data>0 && denValue_endcap_data>0){
      y_b_endcap_data[n]=numValue_endcap_data/denValue_endcap_data;
      ey_endcap_data[n] = (numValue_endcap_data/denValue_endcap_data)*multiplicativeErrorpropagator[3][i];
    }
    else
      y_b_barrel_data[n] = y_b_endcap_data[n] = 0.0;
    
    n++;
  }

  // Open a canvas
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  //Declare the graph
  TGraphErrors *g1 = new TGraphErrors(nbin,x_b,y_b_barrel,ex,ey_barrel);
  TGraphErrors *g2 = new TGraphErrors(nbin,x_b,y_b_barrel_data,ex,ey_barrel_data);
  //decorate the graph
  decorate(g1,";Probe p_{T} [GeV]; Efficiency",0.5,1.1,kRed,21,kRed,2);
  decorate(g2,";Probe p_{T} [GeV]; Efficiency",0.5,1.1,kBlue-7,21,kBlue-7,2);
  g1->SetTitle("TriggerEfficiency_Barrel");
  auto legend = new TLegend(0.85,0.60,0.55,0.70);
  legend->AddEntry(g1,"MC_DY","lep");
  legend->AddEntry(g2,"data","lep");
  legend->SetLineColor(17);
  legend->SetFillStyle(2);
  g1->Draw("ap");
  g2->Draw("p same");
  legend->Draw();
  //c1->SetLogx(1);
  g1->GetYaxis()->SetRangeUser(0.0,1.1);
  c1->SaveAs("TriggerEfficiency_Barrel.png");

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  //Declare the graph
  TGraphErrors *g3 = new TGraphErrors(nbin,x_b,y_b_endcap,ex,ey_endcap);
  TGraphErrors *g4 = new TGraphErrors(nbin,x_b,y_b_endcap_data,ex,ey_endcap_data);
  //decorate the graph
  decorate(g3,";Probe p_{T} [GeV]; Efficiency",0.5,1.1,kRed,21,kRed,2);
  decorate(g4,";Probe p_{T} [GeV]; Efficiency",0.5,1.1,kBlue-7,21,kBlue-7,2);
  g3->SetTitle("TriggerEfficiency_Endcap");
  legend->SetFillStyle(2);
  g3->Draw("ap");
  g4->Draw("p same");
  legend->Draw();
  //c2->SetLogx(1);
  g3->GetYaxis()->SetRangeUser(0.0,1.1);
  c2->SaveAs("TriggerEfficiency_Endcap.png");
  
  
}


float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
    int bin_width = hst->GetBinWidth(1);
    int ibin_begin = 1;
    float nevents = 0.;
    while ( hst->GetBinCenter(ibin_begin) < bin_lo )
        ibin_begin++;
    int ibin_end = ibin_begin;
    while ( hst->GetBinCenter(ibin_end) < bin_hi )
        ibin_end++;
    for ( int i=ibin_begin; i<ibin_end; i++ )
        nevents += hst->GetBinContent(i);

    return nevents;
}
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
    int bin_width = hst->GetBinWidth(1);
    int ibin_begin = 1;
    Double_t nevents = 0.;
    while ( hst->GetBinCenter(ibin_begin) < bin_lo )
      ibin_begin++;
    int ibin_end = ibin_begin;
    while ( hst->GetBinCenter(ibin_end) < bin_hi )
      ibin_end++;
    for ( int i=ibin_begin; i<ibin_end; i++ )
      nevents += pow(hst->GetBinError(i),2);
    nevents = sqrt(nevents);
    return nevents;
}
void decorate(TGraph *h, TString gtitle, float gmax, float gmin, 
	      int markercolor, int markerstyle, int linecolor, int linewidth)
{
  h->SetTitle(gtitle);
  h->SetMarkerColor(markercolor); h->SetMarkerStyle(markerstyle);
  h->SetLineColor(linecolor); h->SetLineWidth(linewidth);
  h->SetMaximum(gmax); h->SetMinimum(gmin);
}
void decorateHist(TH1*h,const char* xtitle, const char* ytitle, const char* title,
		  int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {
  
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  
  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);

  h->SetMarkerSize(1.1);
  h->SetTitle(title);
}

void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(2);
  g->SetFillColor(4);
  g->SetBorderSize(1);
  g->SetLineColor(0);
  //Usually legends should not have headers
  //but if you want one, uncomment the next line.
  //g->SetHeader(legendheader);
}

void SetOverflowBin(TH1F *histo)
{
  int nbins;
  nbins = histo->GetNbinsX();
  histo->SetBinContent(nbins, histo->GetBinContent(nbins) + histo->GetBinContent(nbins+100));//Overflow
  histo->SetBinContent(1, histo->GetBinContent(1)+ histo->GetBinContent(-100));//Underflow
}
