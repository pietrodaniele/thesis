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
// || * fitting all DCB parameters with a linear fit varing the mass;         ||
// || * cutFlow selection;                                                    ||
// || * after the linear fit, the results are used for a global fit:          ||
// ||   - the DCB params are:                                                 ||
// ||     + mu(m_h)=a + b*m_h; and sigma(m_h)=a + b*m_h;                      ||
// ||     + all param = a + b*m_h;                                            ||
// ||     + then the two methods are compared;                                ||
// ||                                                                         ||
// #############################################################################

// functions
vector<double> param_lin(double *, double *, double *, double *, string);
RooDataSet *read_data(int);
void fit_and_plot(vector<double>, vector<double>, int);
void fit_and_plot_all(vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, int);


// main function
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

  // linear fit results ==> (a,b) from y = a + b*x
  vector<double> mu_par = param_lin(masses,masses_error,mu,mu_error,"mu");
  vector<double> sigma_par = param_lin(masses,masses_error,sigma,sigma_error,"sigma");
  vector<double> a1_par = param_lin(masses,masses_error,a1,a1_error,"a1");
  vector<double> p1_par = param_lin(masses,masses_error,p1,p1_error,"p1");
  vector<double> a2_par = param_lin(masses,masses_error,a2,a2_error,"a2");
  vector<double> p2_par = param_lin(masses,masses_error,p2,p2_error,"p2");

  // fit and plot the distribution with mu(m_h) and sigma(m_h)
  fit_and_plot(mu_par,sigma_par,110);
  fit_and_plot(mu_par,sigma_par,125);
  fit_and_plot(mu_par,sigma_par,130);
  fit_and_plot(mu_par,sigma_par,140);

  // fit and plot the distribution with all params(m_h)
  fit_and_plot_all(mu_par,sigma_par,a1_par,p1_par,a2_par,p2_par,110);
  fit_and_plot_all(mu_par,sigma_par,a1_par,p1_par,a2_par,p2_par,125);
  fit_and_plot_all(mu_par,sigma_par,a1_par,p1_par,a2_par,p2_par,130);
  fit_and_plot_all(mu_par,sigma_par,a1_par,p1_par,a2_par,p2_par,140);

  return 0;
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
    gr->GetXaxis()->SetTitle("m_yy [GeV]");
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

void fit_and_plot(vector<double> mu_par, vector<double> sigma_par, int MASS){
  // dataset
  RooDataSet * dataset_110 = read_data(MASS);
  // global fit
  // 110 GeV
  // int MASS = 110;
  int intervall = 30;
  RooRealVar m("m_yy","mass",MASS-intervall,MASS+intervall);
  // fit_parameters
  // ---------------------------------------------------------------------------
  RooRealVar m_h("m_h", "Higgs mass", MASS);
  // ---------------------------------------------------------------------------
  RooRealVar A_mu("A_mu","Bias mu",mu_par[0]);
  RooRealVar B_mu("B_mu","Ang coeff mu",mu_par[1]);
  RooFormulaVar mean("mu","@0 + @1*m_h",RooArgList(A_mu,B_mu,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar A_width("A_width","Bias width",sigma_par[0]);
  RooRealVar B_width("B_width","Ang coeff width",sigma_par[1]);
  RooFormulaVar width_110("width","@0 + @1*m_h",RooArgList(A_width,B_width,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar a1_110("a1","a1",1,0.0001,100.);
  RooRealVar p1_110("p1","p1",1,0.0001,100.);
  RooRealVar a2_110("a2","a2",1,0.0001,100.);
  RooRealVar p2_110("p2","p2",1,0.0001,100.);
  RooCrystalBall * dcbPdf = new RooCrystalBall("dcbPdf","DoubleSidedCB",m,mean,width_110,a1_110,p1_110,a2_110,p2_110);
  // fit to datas
  dcbPdf->fitTo(*dataset_110, Range("signal"), AsymptoticError(true));

  cout << "--------------------------------------------" << endl;
  cout << "DCB fit:" << endl;
  cout << "mu = " << mean.getVal() << endl;
  cout << "width = " << width_110.getVal() << endl;
  cout << "p1 = " << p1_110.getVal() << " +- " << p1_110.getAsymErrorHi() << endl;
  cout << "a1 = " << a1_110.getVal() << " +- " << a1_110.getAsymErrorHi() << endl;
  cout << "a2 = " << a2_110.getVal() << " +- " << a2_110.getAsymErrorHi() << endl;
  cout << "p2 = " << p2_110.getVal() << " +- " << p2_110.getAsymErrorHi() << endl;
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
  string values[6] = {to_string(mean.getVal()),to_string(width_110.getVal()),to_string(a1_110.getVal()),to_string(p1_110.getVal()),to_string(a2_110.getVal()),to_string(p2_110.getVal())};
  string errors[6] = {" "," ",to_string(a1_110.getAsymErrorHi()),to_string(p1_110.getAsymErrorHi()),to_string(a2_110.getAsymErrorHi()),to_string(p2_110.getAsymErrorHi())};
  string pm = " +- ";


  TPaveText *text_fit  = new TPaveText(0.9, 0.65, 0.7, 0.9,"NDC");
  text_fit->SetFillStyle(0);
  text_fit->SetTextAlign(12);
  text_fit->SetTextFont(42);
  text_fit->SetTextSize(0.025);
  text_fit->AddText("DCB fit:");
  string line;
  for(int j=0; j<6; j++){
    if (j<2){
      line = intro[j]+values[j];
    }else{
      line = intro[j]+values[j]+pm+errors[j];
    }
    text_fit->AddText(line.c_str());
  }

  // plotting
  TCanvas *c2 = new TCanvas("c2","m_yy",400,100,1200,700);
  c2->cd();
  RooPlot *frame = m.frame();
  // SetTitle
  const char *title_init = "m_{yy} distribution + DCB fit [";
  const char *title_end = " GeV] with #mu(m_{h}) and #sigma(m_{h})";
  string title_plot = title_init + to_string(MASS) + title_end;
  frame->SetTitle(title_plot.c_str()); // 140 GeV file
  // settings the axis names
  frame->SetXTitle("m_yy [GeV]");
  frame->SetYTitle("Events");

  // ploting the datas and fits
  dataset_110->plotOn(frame,Name("Data Hist"));
  dcbPdf->plotOn(frame,Name("DCB Fit"),LineColor(2));

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
  const char *save_end = "GeV_DCBfit_cutFlow_mu_sigma_mh.pdf";
  string path_plot = save_init + to_string(MASS) + save_end;
  c2->SaveAs(path_plot.c_str());

  return;
}

void fit_and_plot_all(vector<double> mu_par, vector<double> sigma_par, vector<double> a1_par, vector<double> p1_par, vector<double> a2_par, vector<double> p2_par, int MASS){
  // dataset
  RooDataSet * dataset = read_data(MASS);
  // global fit
  // 110 GeV
  // int MASS = 110;
  int intervall = 30;
  RooRealVar m("m_yy","mass",MASS-intervall,MASS+intervall);
  // fit_parameters
  // ---------------------------------------------------------------------------
  RooRealVar m_h("m_h", "Higgs mass", MASS);
  // ---------------------------------------------------------------------------
  RooRealVar A_mu("A_mu","Bias mu",mu_par[0]);
  RooRealVar B_mu("B_mu","Ang coeff mu",mu_par[1]);
  RooFormulaVar mean("mu","@0 + @1*m_h",RooArgList(A_mu,B_mu,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar A_width("A_width","Bias width",sigma_par[0]);
  RooRealVar B_width("B_width","Ang coeff width",sigma_par[1]);
  RooFormulaVar width("width","@0 + @1*m_h",RooArgList(A_width,B_width,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar A_a1("A_a1","Bias a1",a1_par[0]);
  RooRealVar B_a1("B_a1","Ang coeff a1",a1_par[1]);
  RooFormulaVar a1("a1","@0 + @1*m_h",RooArgList(A_a1,B_a1,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar A_p1("A_p1","Bias p1",p1_par[0]);
  RooRealVar B_p1("B_p1","Ang coeff p1",p1_par[1]);
  RooFormulaVar p1("p1","@0 + @1*m_h",RooArgList(A_p1,B_p1,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar A_a2("A_a2","Bias a2",a2_par[0]);
  RooRealVar B_a2("B_a2","Ang coeff a2",a2_par[1]);
  RooFormulaVar a2("a2","@0 + @1*m_h",RooArgList(A_a2,B_a2,m_h));
  // ---------------------------------------------------------------------------
  RooRealVar A_p2("A_p2","Bias p2",p2_par[0]);
  RooRealVar B_p2("B_p2","Ang coeff p2",p2_par[1]);
  RooFormulaVar p2("p2","@0 + @1*m_h",RooArgList(A_p2,B_p2,m_h));
  // ---------------------------------------------------------------------------
  RooCrystalBall * dcbPdf = new RooCrystalBall("dcbPdf","DoubleSidedCB",m,mean,width,a1,p1,a2,p2);
  // fit to datas
  dcbPdf->fitTo(*dataset, Range("signal"), AsymptoticError(true));

  cout << "--------------------------------------------" << endl;
  cout << "DCB fit:" << endl;
  cout << "mu = " << mean.getVal() << endl;
  cout << "width = " << width.getVal() << endl;
  cout << "p1 = " << p1.getVal() << endl;
  cout << "a1 = " << a1.getVal() << endl;
  cout << "a2 = " << a2.getVal() << endl;
  cout << "p2 = " << p2.getVal() << endl;
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
  string values[6] = {to_string(mean.getVal()),to_string(width.getVal()),to_string(a1.getVal()),to_string(p1.getVal()),to_string(a2.getVal()),to_string(p2.getVal())};
  // string errors[6] = {" "," ",to_string(a1_110.getAsymErrorHi()),to_string(p1_110.getAsymErrorHi()),to_string(a2_110.getAsymErrorHi()),to_string(p2_110.getAsymErrorHi())};
  // string pm = " +- ";


  TPaveText *text_fit  = new TPaveText(0.9, 0.65, 0.7, 0.9,"NDC");
  text_fit->SetFillStyle(0);
  text_fit->SetTextAlign(12);
  text_fit->SetTextFont(42);
  text_fit->SetTextSize(0.025);
  text_fit->AddText("DCB fit:");
  string line;
  for(int j=0; j<6; j++){
    string line = intro[j]+values[j];
    text_fit->AddText(line.c_str());
  }

  // plotting
  TCanvas *c2 = new TCanvas("c2","m_yy",400,100,1200,700);
  c2->cd();
  RooPlot *frame = m.frame();
  // SetTitle
  const char *title_init = "m_{yy} distribution + DCB fit [";
  const char *title_end = " GeV] with all params(m_{h})";
  string title_plot = title_init + to_string(MASS) + title_end;
  frame->SetTitle(title_plot.c_str()); // 140 GeV file
  // settings the axis names
  frame->SetXTitle("m_yy [GeV]");
  frame->SetYTitle("Events");

  // ploting the datas and fits
  dataset->plotOn(frame,Name("Data Hist"));
  dcbPdf->plotOn(frame,Name("DCB Fit"),LineColor(2));

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
  const char *save_end = "GeV_DCBfit_cutFlow_all_mh.pdf";
  string path_plot = save_init + to_string(MASS) + save_end;
  c2->SaveAs(path_plot.c_str());

  return;
}
