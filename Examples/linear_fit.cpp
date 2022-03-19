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

// #############################################################################
// ||                                                                         ||
// || * Uploading the single fit results                                      ||
// || * Applying to these fit results a linear fit                            ||
// || * Saving the lin fit results (A,B) in the "linear_fit.txt" file in      ||
// ||   Fit_Results directory:                                                ||
// ||   <fit_parameter> A B                                                   ||
// ||   ...                                                                   ||
// ||                                                                         ||
// #############################################################################

// functions
vector<double> param_lin(double *, double *, double *, double *, string);

void linear_fit(){
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

  // linear fit results ==> (a,b) from y = a + b*x
  vector<double> mu_par = param_lin(masses,masses_error,mu,mu_error,"mu");
  vector<double> sigma_par = param_lin(masses,masses_error,sigma,sigma_error,"sigma");
  vector<double> a1_par = param_lin(masses,masses_error,a1,a1_error,"a1");
  vector<double> p1_par = param_lin(masses,masses_error,p1,p1_error,"p1");
  vector<double> a2_par = param_lin(masses,masses_error,a2,a2_error,"a2");
  vector<double> p2_par = param_lin(masses,masses_error,p2,p2_error,"p2");

  // writing fit results
  ofstream write;
  write.open("Fit_Results/linear_fit.txt"); //open is the method of ofstream
  write << "mu " << mu_par[0] << " " << mu_par[1] << endl;
  write << "sigma " << sigma_par[0] << " " << sigma_par[1] << endl;
  write << "a1 " << a1_par[0] << " " << a1_par[1] << endl;
  write << "p1 " << p1_par[0] << " " << p1_par[1] << endl;
  write << "a2 " << a2_par[0] << " " << a2_par[1] << endl;
  write << "p2 " << p2_par[0] << " " << p2_par[1] << endl;
  write.close();

}


vector<double> param_lin(double *M, double *M_err, double *par, double *par_err, string name){
    vector<double> a_b;

    // preparing the fit
    TF1 * Lfit = new TF1("Linear fit","[0]+x*[1]",100,150);
    Lfit->SetLineColor(kRed);
    Lfit->SetLineStyle(2);
    Lfit->SetParNames("A","B");

    // fit
    TCanvas* c1 = new TCanvas("c1","mu fit",200,10,700,500);
    TGraphErrors* gr = new TGraphErrors(4,M,par,M_err,par_err);
    const char *dist = " distribution with a linear fit";
    string title = name + dist;
    gr->SetTitle(title.c_str());
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetTitle("m_{h} [GeV]");
    const char *ytit = " fit";
    string ytitle = name + ytit;
    gr->GetYaxis()->SetTitle(ytitle.c_str());
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
    const char *repo = "Plots/";
    const char *end = "_fit.pdf";
    string path = repo + name + end;

    c1->SaveAs(path.c_str());

    // saving linear fit parameters
    a_b.push_back(Lfit->GetParameter(0));
    a_b.push_back(Lfit->GetParameter(1));

    return a_b;
}
