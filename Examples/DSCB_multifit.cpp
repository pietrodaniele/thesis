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
#include "function.h"

using namespace RooFit;

// #############################################################################
// ||                                                                         ||
// || * fit a DSCB with 3 consecutive fits                                    ||
// ||   - tails fixed, fit on mu and sigma                                    ||
// ||   - fit on tails, mu and sigma fixed                                    ||
// ||   - fit on all params                                                   ||
// #############################################################################

// functions
RooDataSet *read_data(int);
void multifit_and_plot(int);

// main function
void DSCB_multifit(){

  multifit_and_plot(110);
  multifit_and_plot(125);
  multifit_and_plot(130);
  multifit_and_plot(140);

  return 0;
}

RooDataSet *read_data(int MASS){
  // usefull variables
  string name = to_string(MASS);
  int intervall = 30;
  // read variables
  float m_yy, weight;
  int cutFlow;
  // defining objetcs + datas
  const char *repo = "../Data/PowhegPy8_NNLOPS_ggH";
  const char *end = ".root";
  string path = repo + name + end;
  TFile *input = new TFile(path.c_str(),"read"); // reading the 140 GeV file;
  // datas
  TTree *datatree = (TTree*)input->Get("myCatNtuple;1");
  datatree->SetBranchAddress("m_yy", &m_yy);
  datatree->SetBranchAddress("EventWeight", &weight);
  datatree->SetBranchAddress("cutFlow", &cutFlow);

  RooRealVar m("m_yy","mass",MASS-intervall,MASS+intervall);
  RooRealVar w("weight","weight",-1000,1000);

  RooDataSet * dataset = new RooDataSet("dataset","dataset",RooArgSet(m,w),WeightVar(w));
  for(int i=0;i<datatree->GetEntries();i++){
    datatree->GetEvent(i);
    m.setVal(m_yy*pow(10,-3));
    // photons selection with cutFlow > 13
    if(cutFlow > 13){
      dataset->add(RooArgSet(m),weight*128965.16);
    }
  }
  return dataset;
}

void multifit_and_plot(int MASS){
  // dataset
  RooDataSet * dataset = read_data(MASS);
  // usefull variables
  int intervall = 30;
  RooRealVar m("m_yy","mass",MASS-intervall,MASS+intervall);

  //############################################################################
  // fit on mu and sigma
  RooRealVar mean("mu","mu",MASS,MASS-intervall/2.,MASS+intervall/2.);
  RooRealVar width("width","width",2.5,0.,30.);
  RooRealVar a1("a1","a1",1);
  RooRealVar p1("p1","p1",1);
  RooRealVar a2("a2","a2",1);
  RooRealVar p2("p2","p2",1);

  RooCrystalBall * dcbPdf_mu_sigma = new RooCrystalBall("dcbPdf","DoubleSidedCB",m,mean,width,a1,p1,a2,p2);
  // fit to datas
  dcbPdf_mu_sigma->fitTo(*dataset, Range("signal"), AsymptoticError(true));

  cout << "--------------------------------------------" << endl;
  cout << "DCB fit:" << endl;
  cout << "mu = " << mean.getVal() << " +- " << mean.getAsymErrorHi() << endl;
  cout << "width = " << width.getVal() << " +- " << width.getAsymErrorHi() << endl;
  cout << "p1 = " << p1.getVal() << " +- " << p1.getAsymErrorHi() << endl;
  cout << "a1 = " << a1.getVal() << " +- " << a1.getAsymErrorHi() << endl;
  cout << "a2 = " << a2.getVal() << " +- " << a2.getAsymErrorHi() << endl;
  cout << "p2 = " << p2.getVal() << " +- " << p2.getAsymErrorHi() << endl;
  cout << "############################################" << endl;
  cout << endl;

  //############################################################################
  // fit on tails
  RooRealVar mean_tails("mu","mu",mean.getVal());
  RooRealVar width_tails("width","width",width.getVal());
  RooRealVar a1_tails("a1","a1",1,0.001,100.);
  RooRealVar p1_tails("p1","p1",1,0.001,100.);
  RooRealVar a2_tails("a2","a2",1,0.001,100.);
  RooRealVar p2_tails("p2","p2",1,0.001,100.);

  RooCrystalBall * dcbPdf_tails = new RooCrystalBall("dcbPdf","DoubleSidedCB",m,mean_tails,width_tails,a1_tails,p1_tails,a2_tails,p2_tails);
  // fit to datas
  dcbPdf_tails->fitTo(*dataset, Range("signal"), AsymptoticError(true));

  cout << "--------------------------------------------" << endl;
  cout << "DCB fit:" << endl;
  cout << "mu = " << mean_tails.getVal() << " +- " << mean_tails.getAsymErrorHi() << endl;
  cout << "width = " << width_tails.getVal() << " +- " << width_tails.getAsymErrorHi() << endl;
  cout << "p1 = " << p1_tails.getVal() << " +- " << p1_tails.getAsymErrorHi() << endl;
  cout << "a1 = " << a1_tails.getVal() << " +- " << a1_tails.getAsymErrorHi() << endl;
  cout << "a2 = " << a2_tails.getVal() << " +- " << a2_tails.getAsymErrorHi() << endl;
  cout << "p2 = " << p2_tails.getVal() << " +- " << p2_tails.getAsymErrorHi() << endl;
  cout << "############################################" << endl;
  cout << endl;

  //############################################################################
  // fit on all
  RooRealVar mean_all("mu","mu",mean_tails.getVal(),MASS-intervall/2.,MASS+intervall/2.);
  RooRealVar width_all("width","width",width_tails.getVal(),0.,30.);
  RooRealVar a1_all("a1","a1",a1_tails.getVal(),0.001,100.);
  RooRealVar p1_all("p1","p1",p1_tails.getVal(),0.001,100.);
  RooRealVar a2_all("a2","a2",a2_tails.getVal(),0.001,100.);
  RooRealVar p2_all("p2","p2",p2_tails.getVal(),0.001,100.);

  RooCrystalBall * dcbPdf_all = new RooCrystalBall("dcbPdf","DoubleSidedCB",m,mean_all,width_all,a1_all,p1_all,a2_all,p2_all);
  // fit to datas
  dcbPdf_all->fitTo(*dataset, Range("signal"), AsymptoticError(true));

  cout << "--------------------------------------------" << endl;
  cout << "DCB fit:" << endl;
  cout << "mu = " << mean_all.getVal() << " +- " << mean_all.getAsymErrorHi() << endl;
  cout << "width = " << width_all.getVal() << " +- " << width_all.getAsymErrorHi() << endl;
  cout << "p1 = " << p1_all.getVal() << " +- " << p1_all.getAsymErrorHi() << endl;
  cout << "a1 = " << a1_all.getVal() << " +- " << a1_all.getAsymErrorHi() << endl;
  cout << "a2 = " << a2_all.getVal() << " +- " << a2_all.getAsymErrorHi() << endl;
  cout << "p2 = " << p2_all.getVal() << " +- " << p2_all.getAsymErrorHi() << endl;
  cout << "############################################" << endl;
  cout << endl;

  // adding text
  TPaveText *text_cut  = new TPaveText(0.05,0.88,0.77,0.98,"NDC");
  text_cut->AddText("Ph selection cutFlow > 13");
  text_cut->SetBorderSize(0);
  text_cut->SetFillStyle(0);
  text_cut->SetTextAlign(12);
  text_cut->SetTextFont(42);
  text_cut->SetTextSize(0.02);

  string intro[6] = {"#mu = ","#sigma = ","a_{1} = ","p_{1} = ","a_{2} = ","p_{2} = "};
  string values[6] = {to_string(mean_all.getVal()),to_string(width_all.getVal()),to_string(a1_all.getVal()),to_string(p1_all.getVal()),to_string(a2_all.getVal()),to_string(p2_all.getVal())};
  string errors[6] = {to_string(mean_all.getAsymErrorHi()),to_string(width_all.getAsymErrorHi()),to_string(a1_all.getAsymErrorHi()),to_string(p1_all.getAsymErrorHi()),to_string(a2_all.getAsymErrorHi()),to_string(p2_all.getAsymErrorHi())};
  string pm = " +- ";


  TPaveText *text_fit  = new TPaveText(0.9, 0.65, 0.7, 0.9,"NDC");
  text_fit->SetFillStyle(0);
  text_fit->SetTextAlign(12);
  text_fit->SetTextFont(42);
  text_fit->SetTextSize(0.025);
  text_fit->AddText("DCB fit:");
  string line;
  for(int j=0; j<6; j++){
    line = intro[j]+values[j]+pm+errors[j];
    text_fit->AddText(line.c_str());
  }

  // plotting
  TCanvas *c2 = new TCanvas("c2","m_yy",400,100,1200,700);
  c2->cd();
  RooPlot *frame = m.frame();
  // SetTitle
  const char *title_init = "m_{yy} distribution + DCB fit [";
  const char *title_end = " GeV] multifit";
  string title_plot = title_init + to_string(MASS) + title_end;
  frame->SetTitle(title_plot.c_str()); // 140 GeV file
  // settings the axis names
  frame->SetXTitle("m_yy [GeV]");
  frame->SetYTitle("Events");

  // ploting the datas and fits
  dataset->plotOn(frame,Name("Data Hist"));
  dcbPdf_all->plotOn(frame,Name("DCB Fit"),LineColor(2));

  // legend
  TLegend leg(0.15, 0.8, 0.48, 0.9);
  leg.AddEntry(frame->findObject("Data Hist"), "Data Hist", "lep");
  leg.AddEntry(frame->findObject("DCB Fit"), "DCB Fit", "L");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.025);

  // drawing
  frame->Draw();
  text_fit->Draw();
  text_cut->Draw();
  leg.DrawClone();

  // saving plot
  const char *save_init = "Plots/myy_";
  const char *save_end = "GeV_DCBfit_cutFlow_multifit.pdf";
  string path_plot = save_init + to_string(MASS) + save_end;
  c2->SaveAs(path_plot.c_str());

  return;
}
