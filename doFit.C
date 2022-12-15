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

double f_line(double *x, double *par)
{
  double fitval = par[0]*x[0] + par[1];
  return fitval;
}
double f_expo(double *x, double *par)
{
  double fitval = par[1]*exp(par[0]*x[0]);
  return fitval;
}
double f_con(double *x, double *pad)
{
  double fitval = pad[0];
  return fitval;
}
double f_erf(double *x, double *par)
{
  double arg = ( x[0]-par[0] )/( 2.*par[1] ) ;
  double fitval = 0.5*par[2]*(1.+TMath::Erf(arg));
  return fitval;
}
double f_turn(double *x, double *par)
{
  double fitval = par[1] - par[0]/x[0];
  return fitval;
}

void doFit()
{

  //First we input the data
  int npoints = 14;  //Number of data points
  //Input the x and y values as arrays in the next two lines.
  double x[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5};
  double y[] = {83,82.5,80,79,78,77,76,75,74,73,72,71,70,69.5};

  //Input the uncertainties on the x and y values either as arrays
  // double ex[] = {val,val,val,....};

  // or in the following, the x-error (ex) is set to 0 for each point
  double ex[14] = {0.};
  // and the y-error (ey) is set to 0.5% of the y value for that point
  double ey[14] = {0.};
  for(int i=0; i<14; i++) ey[i] = 0.005*y[i];


  // Declare the fit relevant conditions
  float chi2,ndof;
  // The syntax is TF1("some_name",NameOfFunction,Range_low,Range_high,NumberOfParameters)
  // The NameOfFunction should be from the functions listed below, which are
  // f_line, f_expo, f_con, f_erf, f_turn
  // Here we fit the exponential fit, from 0.1 to 7.0, and it takes 2 parameters
  TF1 *f1 = new TF1("fit_exp",f_expo,0.1,7,2);
  // We provide initial values for the fit
  // Declare an array with two values, and provide them
  double par[2]; par[0]=21.; par[1]=-15.; 
  // and set those parameters for the function f1
  f1->SetParameters(par);
  // Set some visual properties of the fit.
  f1->SetLineColor(kRed); f1->SetLineWidth(2);

  // Open a canvas, with the title c1, with width=800, height=600 pixels
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  // Declare a graph, with x along x-axis, y along y-axis with corresponding errors
  TGraphErrors *g1 = new TGraphErrors(npoints,x,y,ex,ey);
  // Set some visual properties of the graph
  g1->SetMarkerColor(kBlack); g1->SetMarkerStyle(20); g1->SetMarkerSize(0.9);
  g1->SetLineColor(kBlack); g1->SetLineWidth(2);
  // Set the title in this format: "Title;X-axis Title;Y-axis Title"
  g1->SetTitle(";X;Y");
  //gr->GetXaxis()->SetTitle("#tau_{candidate} p_{T} [GeV]");
  //gr->GetYaxis()->SetTitle("Fake Tau Rate");
  g1->Draw("ap");

  // Now fit the above declared function to graph g1
  g1->Fit("fit_exp","WEMR");
  f1->Draw("same");
  // Get and display the parameters of the fit
  // and their uncertainties
  chi2 = f1->GetChisquare(); ndof = f1->GetNDF();
  f1->GetParameters(par);
  double *epar = (double*)f1->GetParErrors();
  double e_par[4];
  for ( int i=0; i<3; i++ )e_par[i] = *(epar+i);
  for(int j=0; j<3; j++)
    cout<<"par["<<j<<"] = "<<par[j]<<" +/- "<<e_par[j]<<endl;     
  cout<<"chi2/ndof = "<<chi2<<"/"<<ndof<<endl;

  

}
