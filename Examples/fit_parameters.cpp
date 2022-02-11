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

// fitting all DCB parameters with a linear fit varing the mass
// cutFlow selection

void fit_parameters(){
  // all parameters obtained from previous fit
  double masses[4] = {110,125,130,140};
  double masses_error[4] = {0,0,0,0};
  double mu[4] = {110.122045, 125.120605, 130.119203, 140.128386};
  double mu_error[4] = {0.031419, 0.032085, 0.034132, 0.041798};
  double sigma[4] = {1.683887, 1.788779, 1.823120, 1.906895};
  double sigma_error[4] = {0.037484, 0.037591, 0.040010, 0.048537};
  double a1[4] = {1.737295, 1.723697, 1.739125, 1.748865};
  double a1_error[4] = {0.123130, 0.116492, 0.122213, 0.144540};
  double p1[4] = {4.825106, 4.912847, 4.751721, 4.606459};
  double p1_error[4] = {1.078196, 1.066697, 1.061026, 1.199767};
  double a2[4] = {1.602378, 1.563756, 1.511050, 1.628135};
  double a2_error[4] = {0.145021, 0.130597, 0.126748, 0.178149};
  double p2[4] = {15.929337, 21.855431, 30.468083, 18.320550};
  double p2_error[4] = {9.688985, 16.271775, 28.915310, 15.547188};

  // preparing the fit
  TF1 * Lfit = new TF1("Linear fit","[0]+x*[1]",100,150);
  Lfit->SetLineColor(kRed);
  Lfit->SetLineStyle(2);
  Lfit->SetParNames("A","B");

  //############################################################################
  // mu fit
  TCanvas* c1 = new TCanvas("c1","mu fit",200,10,700,500);
  TGraphErrors* gr = new TGraphErrors(4,masses,mu,masses_error,mu_error);
  gr->SetTitle("Mu distribution with a linear fit");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetTitle("m_yy [GeV]");
  gr->GetYaxis()->SetTitle("m_yy fit");
  gr->Draw("ALP");
  // fitting
  gr->Fit(Lfit);
  cout << Lfit->GetChisquare() << endl;
  // TLegend
  TLegend* leg = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
  leg->AddEntry(gr,"Exp data");
  leg->AddEntry(Lfit,"Fit");
  leg->Draw("SAME");
  // fit results
  TPaveText *text = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
  string line = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
  text->AddText(line.c_str());
  string line_chi = "ChiQuadro = "+to_string(Lfit->GetChisquare());
  text->AddText(line_chi.c_str());
  text->Draw();
  // saving the plot
  c1->SaveAs("Plots/mu_fit.pdf");
  //############################################################################
  // sigma fit
  TCanvas* c2 = new TCanvas("c2","sigma fit",200,10,700,500);
  TGraphErrors* gr2 = new TGraphErrors(4,masses,sigma,masses_error,sigma_error);
  gr2->SetTitle("Sigma distribution with a linear fit");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);
  gr2->GetXaxis()->SetTitle("m_yy [GeV]");
  gr2->GetYaxis()->SetTitle("sigma fit");
  gr2->Draw("ALP");
  // fitting
  gr2->Fit(Lfit);
  cout << Lfit->GetChisquare() << endl;
  // TLegend
  TLegend* leg2 = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
  leg2->AddEntry(gr2,"Exp data");
  leg2->AddEntry(Lfit,"Fit");
  leg2->Draw("SAME");
  // fit results
  TPaveText *text2 = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
  string line2 = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
  text2->AddText(line2.c_str());
  string line_chi2 = "ChiQuadro = "+to_string(Lfit->GetChisquare());
  text2->AddText(line_chi2.c_str());
  text2->Draw();
  // saving the plot
  c2->SaveAs("Plots/sigma_fit.pdf");
  //############################################################################
  // a1 fit
  TCanvas* c3 = new TCanvas("c3","a1 fit",200,10,700,500);
  TGraphErrors* gr3 = new TGraphErrors(4,masses,a1,masses_error,a1_error);
  gr3->SetTitle("a1 distribution with a linear fit");
  gr3->SetMarkerColor(4);
  gr3->SetMarkerStyle(21);
  gr3->GetXaxis()->SetTitle("m_yy [GeV]");
  gr3->GetYaxis()->SetTitle("a1 fit");
  gr3->Draw("ALP");
  // fitting
  gr3->Fit(Lfit);
  cout << Lfit->GetChisquare() << endl;
  // TLegend
  TLegend* leg3 = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
  leg3->AddEntry(gr3,"Exp data");
  leg3->AddEntry(Lfit,"Fit");
  leg3->Draw("SAME");
  // fit results
  TPaveText *text3 = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
  string line3 = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
  text3->AddText(line3.c_str());
  string line_chi3 = "ChiQuadro = "+to_string(Lfit->GetChisquare());
  text3->AddText(line_chi3.c_str());
  text3->Draw();
  // saving the plot
  c3->SaveAs("Plots/a1_fit.pdf");
  //############################################################################
  // p1 fit
  TCanvas* c4 = new TCanvas("c4","p1 fit",200,10,700,500);
  TGraphErrors* gr4 = new TGraphErrors(4,masses,p1,masses_error,p1_error);
  gr4->SetTitle("p1 distribution with a linear fit");
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(21);
  gr4->GetXaxis()->SetTitle("m_yy [GeV]");
  gr4->GetYaxis()->SetTitle("p1 fit");
  gr4->Draw("ALP");
  // fitting
  gr4->Fit(Lfit);
  cout << Lfit->GetChisquare() << endl;
  // TLegend
  TLegend* leg4 = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
  leg4->AddEntry(gr4,"Exp data");
  leg4->AddEntry(Lfit,"Fit");
  leg4->Draw("SAME");
  // fit results
  TPaveText *text4 = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
  string line4 = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
  text4->AddText(line4.c_str());
  string line_chi4 = "ChiQuadro = "+to_string(Lfit->GetChisquare());
  text4->AddText(line_chi4.c_str());
  text4->Draw();
  // saving the plot
  c4->SaveAs("Plots/p1_fit.pdf");
  //############################################################################
  // a2 fit
  TCanvas* c5 = new TCanvas("c5","a2 fit",200,10,700,500);
  TGraphErrors* gr5 = new TGraphErrors(4,masses,a2,masses_error,a2_error);
  gr5->SetTitle("a2 distribution with a linear fit");
  gr5->SetMarkerColor(4);
  gr5->SetMarkerStyle(21);
  gr5->GetXaxis()->SetTitle("m_yy [GeV]");
  gr5->GetYaxis()->SetTitle("a2 fit");
  gr5->Draw("ALP");
  // fitting
  gr5->Fit(Lfit);
  cout << Lfit->GetChisquare() << endl;
  // TLegend
  TLegend* leg5 = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
  leg5->AddEntry(gr5,"Exp data");
  leg5->AddEntry(Lfit,"Fit");
  leg5->Draw("SAME");
  // fit results
  TPaveText *text5 = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
  string line5 = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
  text5->AddText(line5.c_str());
  string line_chi5 = "ChiQuadro = "+to_string(Lfit->GetChisquare());
  text5->AddText(line_chi5.c_str());
  text5->Draw();
  // saving the plot
  c5->SaveAs("Plots/a2_fit.pdf");
  //############################################################################
  // a2 fit
  TCanvas* c6 = new TCanvas("c6","p2 fit",200,10,700,500);
  TGraphErrors* gr6 = new TGraphErrors(4,masses,p2,masses_error,p2_error);
  gr6->SetTitle("p2 distribution with a linear fit");
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(21);
  gr6->GetXaxis()->SetTitle("m_yy [GeV]");
  gr6->GetYaxis()->SetTitle("p2 fit");
  gr6->Draw("ALP");
  // fitting
  gr6->Fit(Lfit);
  cout << Lfit->GetChisquare() << endl;
  // TLegend
  TLegend* leg6 = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
  leg6->AddEntry(gr6,"Exp data");
  leg6->AddEntry(Lfit,"Fit");
  leg6->Draw("SAME");
  // fit results
  TPaveText *text6 = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
  string line6 = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
  text6->AddText(line6.c_str());
  string line_chi6 = "ChiQuadro = "+to_string(Lfit->GetChisquare());
  text6->AddText(line_chi6.c_str());
  text6->Draw();
  // saving the plot
  c6->SaveAs("Plots/p2_fit.pdf");
  return 0;
}
