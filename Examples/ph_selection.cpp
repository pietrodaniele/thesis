// libraries
#include <iostream>
#include <fstream>

// root libraries
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

using namespace RooFit;

void ph_selection(){
    // write the results file
    ofstream write; //ofstream is the class for fstream package
    write.open("Fit_Results/fit_res_ph_sel.txt"); //open is the method of ofstream
    write << "ph selection: cutFlow > 13" << endl;
    // TCanvas
    TCanvas *c1 = new TCanvas("C1","m_yy",200,10,1200,700);
    c1->Divide(2,2);
    TCanvas *c2 = new TCanvas("C2","m_yy+fit",200,10,1200,700);
    c2->Divide(2,2);
    // text selection
    // adding text
    TPaveText *text = new TPaveText(0.05,0.88,0.77,0.98,"brNDC");
    text->AddText("Ph selection cutFlow > 13");
    text->SetBorderSize(0);
    text->SetFillStyle(0);
    text->SetTextAlign(12);
    text->SetTextFont(42);
    text->SetTextSize(0.04);
    // masses
    int masses[4] = {110,125,130,140}; // 125 GeV non riesco a scaricarlo
    for(int m=0; m<4; m++){
      // write on fit results the mass
      write << "############################################" << endl;
      write << "Fit results of "+to_string(masses[m])+" GeV distribution:" << endl;
      // defining objetcs
      TH1F *hist;
      // due to TFile gives me some problems, I must run this awfull code
      TFile *input;
      if(m==0){
        input = new TFile("../Data/PowhegPy8_NNLOPS_ggH110.root","read"); // reading the 110 GeV file
        hist = new TH1F("m_yy","m_yy distribution [110 GeV]",100,masses[m]-30,masses[m]+30);
      }
      if(m==1){
        input = new TFile("../Data/PowhegPy8_NNLOPS_ggH125.root","read"); // reading the 125 GeV file
        hist = new TH1F("m_yy","m_yy distribution [125 GeV]",100,masses[m]-30,masses[m]+30);
      }
      if(m==2){
        input = new TFile("../Data/PowhegPy8_NNLOPS_ggH130.root","read"); // reading the 130 GeV file
        hist = new TH1F("m_yy","m_yy distribution [130 GeV]",100,masses[m]-30,masses[m]+30);

      }
      if(m==3){
        input = new TFile("../Data/PowhegPy8_NNLOPS_ggH140.root","read"); // reading the 140 GeV file
        hist = new TH1F("m_yy","m_yy distribution [140 GeV]",100,masses[m]-30,masses[m]+30);
      }
      TTree *tree = (TTree*)input->Get("myCatNtuple;1"); // loading the TTree file
      int entries = tree->GetEntries(); // number of rows

      // loading the columns
      float m_yy, weight;
      int cutFlow;

      tree->SetBranchAddress("m_yy", &m_yy);
      tree->SetBranchAddress("EventWeight", &weight);
      tree->SetBranchAddress("cutFlow", &cutFlow);

      for(int i=0; i<entries; i++){
        // getting the colums values
        tree->GetEntry(i);
        // photons selection with cutFlow > 13
        if(cutFlow > 13){
          hist->Fill(m_yy*pow(10,-3),weight);
          // cout << cutFlow << endl;
        }
      }

      // RooFit
      RooRealVar x("x","m_yy",masses[m]-30,masses[m]+30);
      RooDataHist dh("dh", "dh", x, Import(*hist));

      // RooFit gaussian_fit
      // settings
      RooRealVar mean("mean","mean of gaussian",masses[m],masses[m]-30,masses[m]+30);
      RooRealVar sigma("sigma","width of gaussian",2.5,0.,10.);
      RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma); // buil guassian pdf in terms of x, mean and sigma
      // fit pdf to datas
      gauss.fitTo(dh);
      // printing results
      write << "--------------------------------------------" << endl;
      write << "Guassian fit:" << endl;
      write << "mean = " << mean.getVal() << " +- " << mean.getAsymErrorHi()<< endl;
      write << "sigma = " << sigma.getVal() << " +- " << sigma.getAsymErrorHi()<< endl;


      // RooFit RooCrystalBall_fit
      //DCB parameters
      RooRealVar mu("mu","mu",masses[m],masses[m]-30,masses[m]+30);
      RooRealVar width("width","width",2.5,0.,30.);
      RooRealVar a1("a1","a1",1,0.,100.);
      RooRealVar p1("p1","p1",1,0.,100.);
      RooRealVar a2("a2","a2",1,0.,100.);
      RooRealVar p2("p2","p2",1,0.,100.);
      RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB",x,mu,width,a1,p1,a2,p2);
      // fit to datas
      dcbPdf.fitTo(dh);
      // printing results
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

      // drawing datas
      c1->cd(m+1);
      hist->SetXTitle("m_yy [GeV]");
      hist->SetYTitle("Events");
      hist->Draw();
      text->Draw();

      // RooFit plots
      c2->cd(m+1);
      // creating the 110, 130, 140 GeV plots
      RooPlot *frame = x.frame();
      if(m==0){
        frame->SetTitle("m_yy distribution + fits [110 GeV]"); // 110 GeV file
      }
      if(m==1){
        frame->SetTitle("m_yy distribution + fits [125 GeV]"); // 125 GeV file
      }
      if(m==2){
        frame->SetTitle("m_yy distribution + fits [130 GeV]"); // 130 GeV file
      }
      if(m==3){
      frame->SetTitle("m_yy distribution + fits [140 GeV]"); // 140 GeV file
      }
      // settings the axis names
      frame->SetXTitle("m_yy [GeV]");
      frame->SetYTitle("Events");
      // ploting the datas and fits
      dh.plotOn(frame,Name("Data Hist"));
      gauss.plotOn(frame,Name("Guassian Fit"),LineColor(3));
      dcbPdf.plotOn(frame,Name("DCB Fit"),LineColor(2));
      // adding the legend
      TLegend leg(0.11, 0.5, 0.5, 0.8);
      leg.AddEntry(frame->findObject("Data Hist"), "Data Hist", "L");
      leg.AddEntry(frame->findObject("Guassian Fit"), "Guassian Fit", "L");
      leg.AddEntry(frame->findObject("DCB Fit"), "DCB Fit", "L");
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      //drawing
      frame->Draw();
      leg.DrawClone();
      text->Draw();
     }
     // Save the canvas
      c1->SaveAs("Plots/myy_dist_ph_sel.pdf");
      c2->SaveAs("Plots/myy+fit_ph_sel.pdf");
      // closing results file
      write.close();
     return;
}
