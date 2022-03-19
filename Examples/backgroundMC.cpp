#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;

void backgroundMC(){
  // reading file
  float m_yy, weight;
  int cutFlow;
  TFile *input = new TFile("../Data/Sherpa2_myy_90_175.root","read"); // reading the 140 GeV file;
  // datas
  TTree *datatree = (TTree*)input->Get("myCatNtuple;1");
  datatree->SetBranchAddress("m_yy", &m_yy);
  datatree->SetBranchAddress("EventWeight", &weight);
  datatree->SetBranchAddress("cutFlow", &cutFlow);

  RooRealVar m("m_yy","mass",110,170);
  RooRealVar w("weight","weight",-1000,1000);

  RooDataSet * dataset = new RooDataSet("dataset","dataset",RooArgSet(m,w),WeightVar(w));
  for(int i=0;i<datatree->GetEntries();i++){
    datatree->GetEvent(i);
    m.setVal(m_yy*pow(10,-3));
    // photons selection with cutFlow > 13
    if(cutFlow > 13){
      if(m_yy*pow(10,-3)>=110){
        dataset->add(RooArgSet(m),weight*138965.16);
      }
    }
  }

  // preparing the exponential fit
  RooRealVar c("c","c",-0.001,-0.025,-0.0001);
  RooExponential * Efit = new RooExponential("Efit","Efit",m,c);

  /*// preparing the exponential fit 2° poly
  RooRealVar a("a","a",1,0.0001,5);
  RooRealVar b("b","b",1,0.0001,5);
  RooPolynomial * Pfit = new RooPolynomial("Pfit","Pfit",m,RooArgList(a,b),1);*/

  // fitting
  Efit->fitTo(*dataset, AsymptoticError(true));
  //Pfit->fitTo(*dataset, AsymptoticError(true));
  // plotting
  TCanvas* c1 = new TCanvas("c1","c1",600,400) ;

  TPaveText *text = new TPaveText(0.9, 0.75, 0.66, 0.9,"brNDC");
  text->AddText("Exp(x*c):");
  string line = "c = "+to_string(c.getVal())+"+-"+to_string(c.getAsymErrorHi());
  text->AddText(line.c_str());
  //text->AddText(line1.c_str());
  //text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);

  RooPlot *frame = m.frame();
  frame->SetTitle("backgroundMC m_yy distribution + exp(c*x) fit"); // 140 GeV file
  // settings the axis names
  frame->SetXTitle("m_{h} [GeV]");
  frame->SetYTitle("Events");

  // ploting the datas and fits
  dataset->plotOn(frame,Name("Data"));
  Efit->plotOn(frame,Name("Fit"),LineColor(2));
  //Pfit->plotOn(frame,Name("PFit"),LineColor(2));

  // adding the legend
  TLegend leg(0.15, 0.8, 0.48, 0.9);
  leg.AddEntry(frame->findObject("Data"), "Data", "lep");
  leg.AddEntry(frame->findObject("Fit"), "Fit", "L");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.025);

  //drawing
  frame->Draw();
  leg.DrawClone();
  text->Draw();

  // saving plot
  c1->SaveAs("Plots/backgroundMC_dist.pdf");





  /*// preparing the fit
  TF1 * Efit = new TF1("Exp fit","exp([0]+x*[1])",110,160);
  Efit->SetLineColor(kRed);
  Efit->SetLineStyle(2);
  Efit->SetParNames("A","B");

  // createHistogram
  TH1* h_fromdata = dataset->createHistogram("m_yy",m);

  TCanvas* c1 = new TCanvas("c1","c1",600,400) ;
  c1->cd();
  h_fromdata->Draw();
  h_fromdata->Fit(Efit);*/
}
