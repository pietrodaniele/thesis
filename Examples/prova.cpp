#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TLatex.h"

using namespace RooFit;

void prova(){
   TCanvas *c = new TCanvas("c", "c", 800,800);
   c->Draw();
   TPad *p1 = new TPad("p1","p1",0.,0.2,1.,1.);
   p1->SetFillColor(kRed);
   p1->Draw();
   c->cd(0);
   TPad *p2 = new TPad("p2","p2",0.,0.,1.,0.2);
   p2->SetFillColor(kBlue);
   p2->Draw();

}
