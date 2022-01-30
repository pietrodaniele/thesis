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

void ph_all(){
  // write the results file
  ofstream write; //ofstream is the class for fstream package
  write.open("Fit_Results/fit_res.txt"); //open is the method of ofstream
  // TCanvas
  TCanvas *c1 = new TCanvas("C1","m_yy",200,10,1200,700);
  c1->Divide(2,2);
  TCanvas *c2 = new TCanvas("C2","m_yy+fit",200,10,1200,700);
  c2->Divide(2,2);

  // masses
  int masses[4] = {110,125,130,140};
  int intervall = 30;
  for(int m=0; m<4; m++){
    // write on fit results the mass
    write << "############################################" << endl;
    write << "Fit results of "+to_string(masses[m])+" GeV distribution:" << endl;
    // defining objetcs
    // histogram
    string hist_title = "m_yy distribution ["+to_string(masses[m])+" GeV]";
    TH1F *hist = new TH1F("m_yy",hist_title.c_str(),100,masses[m]-intervall,masses[m]+intervall);
    // input file
    string path = "../Data/PowhegPy8_NNLOPS_ggH"+to_string(masses[m])+".root";
    TFile *input = new TFile(path.c_str(),"read"); // reading the ... GeV file
    // TTree
    TTree *tree = (TTree*)input->Get("myCatNtuple;1"); // loading the TTree file

    // loading the columns
    float m_yy, weight;

    tree->SetBranchAddress("m_yy", &m_yy);
    tree->SetBranchAddress("EventWeight", &weight);

    RooRealVar Mass("m_yy","mass",masses[m]-intervall,masses[m]+intervall);
    RooRealVar w("weight","weight",-1000,1000);

    RooDataSet * dataset = new RooDataSet("dataset","dataset",RooArgSet(Mass,w),WeightVar(w));
    for(int i=0;i<tree->GetEntries();i++){
      tree->GetEvent(i);
      Mass.setVal(m_yy*pow(10,-3));
      hist->Fill(m_yy*pow(10,-3),weight*128965.16);
      dataset->add(RooArgSet(Mass),weight*128965.16);
    }
    dataset->Print();

    // RooFit
    RooDataHist dh("dh", "dh", Mass, Import(*hist));

    // RooFit gaussian_fit
    // settings
    RooRealVar mean("mean","mean of gaussian",masses[m],masses[m]-intervall/2,masses[m]+intervall/2);
    RooRealVar sigma("sigma","width of gaussian",2.5,0.,10.);
    RooGaussian gauss("gauss","gaussian PDF",Mass,mean,sigma); // buil guassian pdf in terms of x, mean and sigma
    // fit pdf to datas
    gauss.fitTo(*dataset, Range("signal"), SumW2Error(true));
    // printing results
    write << "--------------------------------------------" << endl;
    write << "Guassian fit:" << endl;
    write << "mean = " << mean.getVal() << " +- " << mean.getAsymErrorHi()<< endl;
    write << "sigma = " << sigma.getVal() << " +- " << sigma.getAsymErrorHi()<< endl;


    // RooFit RooCrystalBall_fit
    //DCB parameters
    RooRealVar mu("mu","mu",masses[m],masses[m]-15,masses[m]+15);
    RooRealVar width("width","width",2.5,0.,30.);
    RooRealVar a1("a1","a1",1,0.,100.);
    RooRealVar p1("p1","p1",1,0.,100.);
    RooRealVar a2("a2","a2",1,0.,100.);
    RooRealVar p2("p2","p2",1,0.,100.);
    RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB",Mass,mu,width,a1,p1,a2,p2);
    // fit to datas
    dcbPdf.fitTo(*dataset, Range("signal"), SumW2Error(true));
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

    // RooFit plots
    c2->cd(m+1);
    // printing results on plot
    TPaveText *text = new TPaveText(0.9, 0.45, 0.6, 0.9,"brNDC");
    string pm = " +- ";
    // gaussian fit
    string intro_gauss[2] = {"mean = ","sigma = "};
    string values_gauss[2] = {to_string(mean.getVal()),to_string(sigma.getVal())};
    string errors_gauss[2] = {to_string(mu.getAsymErrorHi()),to_string(width.getAsymErrorHi())};

    text->AddText("Gauss fit:");
    for(int j=0; j<2; j++){
      string line = intro_gauss[j]+values_gauss[j]+pm+errors_gauss[j];
      text->AddText(line.c_str());
    }
    // separator
    text->AddText("###################");

    // DCB fit
    string intro_dcb[6] = {"mu = ","width = ","a1 = ","p1 = ","a2 = ","p2 = "};
    string values_dcb[6] = {to_string(mu.getVal()),to_string(width.getVal()),to_string(a1.getVal()),to_string(p1.getVal()),to_string(a2.getVal()),to_string(p2.getVal())};
    string errors_dcb[6] = {to_string(mu.getAsymErrorHi()),to_string(width.getAsymErrorHi()),to_string(a1.getAsymErrorHi()),to_string(p1.getAsymErrorHi()),to_string(a2.getAsymErrorHi()),to_string(p2.getAsymErrorHi())};

    text->AddText("DCB fit:");
    for(int j=0; j<6; j++){
      string line = intro_dcb[j]+values_dcb[j]+pm+errors_dcb[j];
      text->AddText(line.c_str());
    }
    // text settings
    text->SetFillStyle(0);
    text->SetTextAlign(12);
    text->SetTextFont(42);
    text->SetTextSize(0.04);

    // creating the 110, 130, 140 GeV plots
    RooPlot *frame = Mass.frame();
    string title = "m_yy distribution + fits ["+to_string(masses[m])+" GeV]";
    frame->SetTitle(title.c_str()); // 110 GeV file
    // settings the axis names
    frame->SetXTitle("m_yy [GeV]");
    frame->SetYTitle("Events");
    // ploting the datas and fits
    dataset->plotOn(frame,Name("Data Hist"));
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
  c1->SaveAs("Plots/myy_dist.pdf");
  c2->SaveAs("Plots/myy+fit.pdf");
  // closing results file
  write.close();

  return;
}
