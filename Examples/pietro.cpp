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
void pietro(){
  int x[2] = {1,2};
  int y[2] = {1,2};
  int z[2] = {1,1};

  TGraph* gr1 = new TGraph(2,x,y);
  TGraph* gr2 = new TGraph(2,x,z);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1,"lp");
  mg->Add(gr2,"cp");
  mg->Draw("a");
}
