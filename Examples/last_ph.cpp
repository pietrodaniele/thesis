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

void last_ph(){
  // write the results file
  ofstream write; //ofstream is the class for fstream package
  write.open("Fit_Results/fit_res_last_ph.txt"); //open is the method of ofstream
  // TCanvas
  TCanvas *c1 = new TCanvas("C1","m_yy",400,100,1200,700);
  // mass
  int mass = 140;
  int intervall = 30;
  float m_yy, weight;

  // write on fit results the mass
  write << "############################################" << endl;
  write << "Fit results of "+to_string(mass)+" GeV distribution:" << endl;
  // defining objetcs + datas
  TFile *input = new TFile("../Data/PowhegPy8_NNLOPS_ggH140.root","read"); // reading the 140 GeV file;
    // datas
  TTree *datatree = (TTree*)input->Get("myCatNtuple;1");
  datatree->SetBranchAddress("m_yy", &m_yy);
  datatree->SetBranchAddress("EventWeight", &weight);


  RooRealVar m("m_yy","mass",mass-intervall,mass+intervall);
  RooRealVar w("weight","weight",-1000,1000);

  RooDataSet * dataset = new RooDataSet("dataset","dataset",RooArgSet(m,w),WeightVar(w));
  for(int i=0;i<datatree->GetEntries();i++){
    datatree->GetEvent(i);
    m.setVal(m_yy*pow(10,-3));
    dataset->add(RooArgSet(m),weight*128965.16);
  }
  dataset->Print();


  // RooFit

  // RooFit RooCrystalBall_fit
  //DCB parameters
  RooRealVar mu("mu","mu",mass,mass-15,mass+15);
  RooRealVar width("width","width",2.5,0.,30.);
  RooRealVar a1("a1","a1",1,0.,100.);
  RooRealVar p1("p1","p1",1,0.,100.);
  RooRealVar a2("a2","a2",1,0.,100.);
  RooRealVar p2("p2","p2",1,0.,100.);
  RooCrystalBall * dcbPdf = new RooCrystalBall("dcbPdf","DoubleSidedCB",m,mu,width,a1,p1,a2,p2);
  // fit to datas
  dcbPdf->fitTo(*dataset, Range("signal"), SumW2Error(true));

  write << "--------------------------------------------" << endl;
  write << "DCB fit:" << endl;
  write << "mu = " << mu.getVal() << " +- " << mu.getAsymErrorHi()<< endl;
  write << "width = " << width.getVal() << " +- " << width.getAsymErrorHi()<< endl;
  write << "a1 = " << a1.getVal() << " +- " << a1.getAsymErrorHi()<< endl;
  write << "p1 = " << a1.getVal() << " +- " << p1.getAsymErrorHi()<< endl;
  write << "a2 = " << a2.getVal() << " +- " << a2.getAsymErrorHi()<< endl;
  write << "p2 = " << a2.getVal() << " +- " << p2.getAsymErrorHi()<< endl;
  write << "############################################" << endl;
  write << endl;
  // / drawing datas
  TPad *pad1 = new TPad("data+fit","data+fit",0.,0.25,1.,1.);
  pad1->Draw();
  pad1->cd();
  // printing results on plot
  string intro[6] = {"mu = ","width = ","a1 = ","p1 = ","a2 = ","p2 = "};
  string values[6] = {to_string(mu.getVal()),to_string(width.getVal()),to_string(a1.getVal()),to_string(p1.getVal()),to_string(a2.getVal()),to_string(p2.getVal())};
  string errors[6] = {to_string(mu.getAsymErrorHi()),to_string(width.getAsymErrorHi()),to_string(a1.getAsymErrorHi()),to_string(p1.getAsymErrorHi()),to_string(a2.getAsymErrorHi()),to_string(p2.getAsymErrorHi())};
  string pm = " +- ";


  TPaveText *text = new TPaveText(0.9, 0.65, 0.7, 0.9,"brNDC");
  text->AddText("DCB fit:");
  for(int j=0; j<6; j++){
    string line = intro[j]+values[j]+pm+errors[j];
    text->AddText(line.c_str());
  }
  //text->AddText(line1.c_str());
  //text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);
  //
  RooPlot *frame = m.frame();
  frame->SetTitle("m_yy distribution + DCB fit [140 GeV]"); // 140 GeV file
  // settings the axis names
  frame->SetXTitle("m_yy [GeV]");
  frame->SetYTitle("Events");

  // ploting the datas and fits
  dataset->plotOn(frame,Name("Data Hist"));
  dcbPdf->plotOn(frame,Name("DCB Fit"),LineColor(2));

  // adding the legend
  TLegend leg(0.15, 0.8, 0.48, 0.9);
  leg.AddEntry(frame->findObject("Data Hist"), "Data Hist", "lep");
  leg.AddEntry(frame->findObject("DCB Fit"), "DCB Fit", "L");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.025);

  //drawing
  frame->Draw();
  leg.DrawClone();
  text->Draw();
  // Save the canvas
  //c1->SaveAs("Plots/myy_140GeV_DCBfit.pdf");
  write.close();

  // ratio plots
  c1->cd(0);
  TPad *pad2 = new TPad("data/fit","data/fit",0.,0.,1.,0.25);
  pad2->Draw();
  pad2->cd();
  RooDataSet * datafit = dcbPdf->generate(RooArgSet(m),datatree->GetEntries());
  TH1* h_fromdata = dataset->createHistogram("m_yy_data",m);
  h_fromdata->Scale(1./h_fromdata->Integral(), "width");
  TH1* h_fromfit = datafit->createHistogram("m_yy_data",m);
  h_fromfit->Scale(1./h_fromfit->Integral(), "width");

  h_fromdata->Divide(h_fromdata,h_fromfit);
  h_fromdata->SetTitle("");
  h_fromdata->GetXaxis()->SetTitle("");
  //h_fromdata->SetTitleFont(0.3,"X");
  h_fromdata->SetLabelSize(0.08,"X");
  h_fromdata->GetYaxis()->SetTitle("data/fit");
  h_fromdata->GetYaxis()->SetRangeUser(0., 2.);
  h_fromdata->SetLabelSize(0.08,"Y");
  h_fromdata->SetTitleSize(0.08,"Y");
  // drawing
  h_fromdata->Draw();
  // plot a y=1 line
  TLine * line = new TLine(mass-intervall,1,mass+intervall,1);
  line->SetLineColor(kGreen);
  line->Draw();
  // saving plot
  c1->SaveAs("Plots/myy_140GeV_DCBfit.pdf");
  return;
}
