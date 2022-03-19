#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace std;
using namespace RooFit ;

// functions
double eff_calc(int);

int main(){
  double masses[4] = {110,125,130,140};
  double effs[4] = {eff_calc(110),eff_calc(125),eff_calc(130),eff_calc(140)};

  TCanvas* c3 = new TCanvas("c3","c3",600, 400);
  c3->cd();

  TGraph* gr1 = new TGraph(4,masses,effs);
  gr1->SetName("Single_fit");
  gr1->SetTitle("Single_fit");
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerColor(2);
  gr1->SetDrawOption("AP");
  gr1->SetLineColor(2);
  gr1->SetLineWidth(4);
  //gr1->SetFillStyle(0);

  gr1->SetTitle("Efficiencies(m_{h})");
  gr1->GetXaxis()->SetTitle("m_{yy} [GeV]");
  gr1->GetYaxis()->SetTitle("Eff");
  gr1->Draw("ALP");

  c3->SaveAs("Plots/efficiencies.pdf");
  return 0;
}


double eff_calc(int MASS){
  // usefull variables
  string name = to_string(MASS);
  int intervall = 30;
  // read variables
  float weight;
  int cutFlow;
  // defining objetcs + datas
  const char *repo = "../Data/PowhegPy8_NNLOPS_ggH";
  const char *end = ".root";
  string path = repo + name + end;
  TFile * input = new TFile(path.c_str(),"read"); // reading the 140 GeV file;
  // datas
  TTree *datatree = (TTree*)input->Get("myCatNtuple;1");
  datatree->SetBranchAddress("cutFlow", &cutFlow);
  datatree->SetBranchAddress("EventWeight", &weight);
  // counters
  float k_all = 0;
  float k_cut = 0;

  for(int i=0;i<datatree->GetEntries();i++){
    datatree->GetEvent(i);
    k_all+=weight;
    // photons selection with cutFlow > 13
    if(cutFlow > 13){
      k_cut+=weight;
    }
  }
  double eff = k_cut/k_all;
  return eff;
}
