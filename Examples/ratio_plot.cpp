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

void ratio_plot(int argc, char **argv){
  // errors check
  if(argc!=1){
    cout << "Must insert the energy value" << endl;
    return 1;
  }
  // reading the file .root
  TString path("../Data/PowhegPy8_NNLOPS_ggH110.root");
  TFile *input = new TFile(path,"read"); // reading the 140 GeV file;
}
