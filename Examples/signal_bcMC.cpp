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

// functions
RooDataSet *read_data(int, RooRealVar&);
vector<double> param_lin(double *, double *, double *, double *, string);
void plot_param(double [4], double [4], double [4], string);
void plot_results(string, RooFormulaVar, RooFormulaVar, RooFormulaVar, RooFormulaVar, RooFormulaVar, RooFormulaVar,
                  RooRealVar, RooDataSet, RooSimultaneous, RooCategory);


// all params are function of m_h
// cutFlow > 13 applied


void signal_bcMC(){
  RooWorkspace w("w");
  // -------------------------------------------------------------------------------------------------------------------------
  // creating the signal model and the background model
  // *************************************************************************************************************************
  // background model
  w.factory("Exponential:fit_bkg(x[105,160], c[-0.001,-0.025,-0.0001])");
  // *************************************************************************************************************************
  // signal model
  w.factory("Gaussian:fit_110(x,mu_110[110.081076],sigma_110[2.099786])");
  /*w.factory("Gaussian:fit_125(x,mu_125[125.076653],sigma_125[2.222825])");
  w.factory("Gaussian:fit_130(x,mu_130[130.084487],sigma_130[2.270348])");
  w.factory("Gaussian:fit_140(x,mu_140[140.068206],sigma_140[2.363402])");
  // Create discrete observable to label channels
  w.factory("index[110,125,130,140]");

  // Create joint pdf (RooSimultaneous)
  w.factory("SIMUL:jointModel(index,110=fit_110,125=fit_125,130=fit_130,140=fit_140)");*/
  // *************************************************************************************************************************
  // joining two models
  w.factory("SUM:model(fit_110,fit_bkg)");  // for extended model




  // signal model
  /*w.factory("CrystalBall:fit_110(fit_110,DoubleSidedCB,x,MEAN_110[mean_110.getVal()],\
            WIDTH_110[width_110.getVal()],A1_110[a1_110.getVal()],P1_110[p1_110.getVal()],\
            A2_110[a2_110.getVal()],P2_110[P2_110.getVal()])");

  w.factory("CrystalBall:fit_125(x,MEAN_125[mean_125.getVal()],\
            WIDTH_125[width_125.getVal()],A1_125[a1_125.getVal()],P1_125[p1_125.getVal()],\
            A2_125[a2_125.getVal()],P2_125[P2_125.getVal()])");

  w.factory("CrystalBall:fit_130(x,MEAN_130[mean_130.getVal()],\
            WIDTH_130[width_130.getVal()],A1_130[a1_130.getVal()],P1_130[p1_130.getVal()],\
            A2_130[a2_130.getVal()],P2_130[P2_130.getVal()])");

  w.factory("CrystalBall:fit_140(x,MEAN_140[mean_140.getVal()],\
            WIDTH_140[width_140.getVal()],A1_140[a1_140.getVal()],P1_140[p1_140.getVal()],\
            A2_140[a2_140.getVal()],P2_140[P2_140.getVal()])");*/
















  return;
}

/*
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

// creating a single mass range
RooRealVar m("m","mass",80,170);
// -------------------------------------------------------------

// Create observables
int intervall = 30;
int mass110 = 110;
int mass125 = 125;
int mass130 = 130;
int mass140 = 140;
RooRealVar A_mu("A_mu","Bias mu",mu_par[0]);
RooRealVar B_mu("B_mu","Ang coeff mu",mu_par[1]);
RooRealVar A_width("A_width","Bias width",sigma_par[0]);
RooRealVar B_width("B_width","Ang coeff width",sigma_par[1]);
RooRealVar A_a1("A_a1","Bias a1",a1_par[0]);
RooRealVar B_a1("B_a1","Ang coeff a1",a1_par[1]);
RooRealVar A_p1("A_p1","Bias p1",p1_par[0]);
RooRealVar B_p1("B_p1","Ang coeff p1",p1_par[1]);
RooRealVar A_a2("A_a2","Bias a2",a2_par[0]);
RooRealVar B_a2("B_a2","Ang coeff a2",a2_par[1]);
RooRealVar A_p2("A_p2","Bias p2",p2_par[0]);
RooRealVar B_p2("B_p2","Ang coeff p2",p2_par[1]);

// dataset
RooDataSet * dataset110 = read_data(110,m);
RooDataSet * dataset125 = read_data(125,m);
RooDataSet * dataset130 = read_data(130,m);
RooDataSet * dataset140 = read_data(140,m);
//  110 model
// variables 110 GeV
// fit_parameters
// ---------------------------------------------------------------------------
RooRealVar m_h_110("m_h_110", "Higgs mass", mass110);
// ---------------------------------------------------------------------------
RooFormulaVar mean_110("mu_110","@0 + @1*m_h_110",RooArgList(A_mu,B_mu,m_h_110));
// ---------------------------------------------------------------------------
RooFormulaVar width_110("width_110","@0 + @1*m_h_110",RooArgList(A_width,B_width,m_h_110));
// ---------------------------------------------------------------------------
RooFormulaVar a1_110("a1_110","@0 + @1*m_h_110",RooArgList(A_a1,B_a1,m_h_110));
// ---------------------------------------------------------------------------
RooFormulaVar p1_110("p1_110","@0 + @1*m_h_110",RooArgList(A_p1,B_p1,m_h_110));
// ---------------------------------------------------------------------------
RooFormulaVar a2_110("a2_110","@0 + @1*m_h_110",RooArgList(A_a2,B_a2,m_h_110));
// ---------------------------------------------------------------------------
RooFormulaVar p2_110("p2_110","@0 + @1*m_h_110",RooArgList(A_p2,B_p2,m_h_110));
// ---------------------------------------------------------------------------
RooCrystalBall dcbPdf_110("dcbPdf_110","DoubleSidedCB",m,mean_110,width_110,a1_110,p1_110,a2_110,p2_110);

//  125 model
// variables 125 GeV
// fit_parameters
// ---------------------------------------------------------------------------
RooRealVar m_h_125("m_h_125", "Higgs mass", mass125);
// ---------------------------------------------------------------------------
RooFormulaVar mean_125("mu_125","@0 + @1*m_h_125",RooArgList(A_mu,B_mu,m_h_125));
// ---------------------------------------------------------------------------
RooFormulaVar width_125("width_125","@0 + @1*m_h_125",RooArgList(A_width,B_width,m_h_125));
// ---------------------------------------------------------------------------
RooFormulaVar a1_125("a1_125","@0 + @1*m_h_125",RooArgList(A_a1,B_a1,m_h_125));
// ---------------------------------------------------------------------------
RooFormulaVar p1_125("p1_125","@0 + @1*m_h_125",RooArgList(A_p1,B_p1,m_h_125));
// ---------------------------------------------------------------------------
RooFormulaVar a2_125("a2_125","@0 + @1*m_h_125",RooArgList(A_a2,B_a2,m_h_125));
// ---------------------------------------------------------------------------
RooFormulaVar p2_125("p2_125","@0 + @1*m_h_125",RooArgList(A_p2,B_p2,m_h_125));
// ---------------------------------------------------------------------------
RooCrystalBall dcbPdf_125("dcbPdf_125","DoubleSidedCB",m,mean_125,width_125,a1_125,p1_125,a2_125,p2_125);

//  130 model
// variables 130 GeV
// fit_parameters
// ---------------------------------------------------------------------------
RooRealVar m_h_130("m_h_130", "Higgs mass", mass130);
// ---------------------------------------------------------------------------
RooFormulaVar mean_130("mu_130","@0 + @1*m_h_130",RooArgList(A_mu,B_mu,m_h_130));
// ---------------------------------------------------------------------------
RooFormulaVar width_130("width","@0 + @1*m_h_130",RooArgList(A_width,B_width,m_h_130));
// ---------------------------------------------------------------------------
RooFormulaVar a1_130("a1_130","@0 + @1*m_h_130",RooArgList(A_a1,B_a1,m_h_130));
// ---------------------------------------------------------------------------
RooFormulaVar p1_130("p1_130","@0 + @1*m_h_130",RooArgList(A_p1,B_p1,m_h_130));
// ---------------------------------------------------------------------------
RooFormulaVar a2_130("a2_130","@0 + @1*m_h_130",RooArgList(A_a2,B_a2,m_h_130));
// ---------------------------------------------------------------------------
RooFormulaVar p2_130("p2_130","@0 + @1*m_h_130",RooArgList(A_p2,B_p2,m_h_130));
// ---------------------------------------------------------------------------
RooCrystalBall dcbPdf_130("dcbPdf_130","DoubleSidedCB",m,mean_130,width_130,a1_130,p1_130,a2_130,p2_130);

//  140 model
// variables 140 GeV
// fit_parameters
// ---------------------------------------------------------------------------
RooRealVar m_h_140("m_h_140", "Higgs mass", mass140);
// ---------------------------------------------------------------------------
RooFormulaVar mean_140("mu_140","@0 + @1*m_h_140",RooArgList(A_mu,B_mu,m_h_140));
// ---------------------------------------------------------------------------
RooFormulaVar width_140("width","@0 + @1*m_h_140",RooArgList(A_width,B_width,m_h_140));
// ---------------------------------------------------------------------------
RooFormulaVar a1_140("a1_140","@0 + @1*m_h_140",RooArgList(A_a1,B_a1,m_h_140));
// ---------------------------------------------------------------------------
RooFormulaVar p1_140("p1_140","@0 + @1*m_h_140",RooArgList(A_p1,B_p1,m_h_140));
// ---------------------------------------------------------------------------
RooFormulaVar a2_140("a2_140","@0 + @1*m_h_140",RooArgList(A_a2,B_a2,m_h_140));
// ---------------------------------------------------------------------------
RooFormulaVar p2_140("p2_140","@0 + @1*m_h_140",RooArgList(A_p2,B_p2,m_h_140));
// ---------------------------------------------------------------------------
RooCrystalBall dcbPdf_140("dcbPdf_140","DoubleSidedCB",m,mean_140,width_140,a1_140,p1_140,a2_140,p2_140);

// C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
// ---------------------------------------------------------------------------

// Define category to distinguish 110,125,130 and 140 GeV samples events
RooCategory sample("sample","sample") ;
sample.defineType("myy_110") ;
sample.defineType("myy_125") ;
sample.defineType("myy_130") ;
sample.defineType("myy_140") ;

// Construct combined dataset in (x,sample)
RooDataSet combData("combData","combined data",m,Index(sample),Import("myy_110",*dataset110),
                    Import("myy_125",*dataset125),Import("myy_130",*dataset130),Import("myy_140",*dataset140)) ;



// C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
// -----------------------------------------------------------------------------------

// Construct a simultaneous pdf using category sample as index
RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

// Associate model with the physics state and model_ctl with the control state
simPdf.addPdf(dcbPdf_110,"myy_110") ;
simPdf.addPdf(dcbPdf_125,"myy_125") ;
simPdf.addPdf(dcbPdf_130,"myy_130") ;
simPdf.addPdf(dcbPdf_140,"myy_140") ;



// P e r f o r m   a   s i m u l t a n e o u s   f i t
// ---------------------------------------------------

// Perform simultaneous fit of model to data and model_ctl to data_ctl
simPdf.fitTo(combData) ;
*/

void plot_results(string name, RooFormulaVar mean, RooFormulaVar width, RooFormulaVar a1, RooFormulaVar p1, RooFormulaVar a2,RooFormulaVar p2,
                  RooRealVar m, RooDataSet combData, RooSimultaneous simPdf, RooCategory sample){
  // plotting
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
  string errors[6] = {" "," "," "," "," "," "};
  string pm = " +- ";

  TPaveText *text_fit  = new TPaveText(0.9, 0.65, 0.7, 0.9,"NDC");
  text_fit->SetFillStyle(0);
  text_fit->SetTextAlign(12);
  text_fit->SetTextFont(42);
  text_fit->SetTextSize(0.025);
  text_fit->AddText("DCB fit:");
  string line;
  for(int j=0; j<6; j++){
    line = intro[j]+values[j];
    text_fit->AddText(line.c_str());
  }

  TCanvas* c1 = new TCanvas("c1","c1",600,400) ;
  c1->cd();
  const char *title_init = "m_{yy} distribution + DCB global fit [";
  const char *title_end = " GeV] with all params(m_{h})";
  string title_plot = title_init + name + title_end;
  RooPlot* frame = m.frame(Title(title_plot.c_str()));
  const char * name_comb = "sample==sample::myy_";
  string path_comb = name_comb + name;
  combData.plotOn(frame,Name("Data"),Cut(path_comb.c_str())) ;
  const char * name_sim = "myy_";
  string path_sim = name_sim + name;
  simPdf.plotOn(frame,Name("Fit"),Slice(sample,path_sim.c_str()),ProjWData(sample,combData)) ;
  // settings the axis names
  frame->SetXTitle("m_{h} [GeV]");
  frame->SetYTitle("Events");

  // legend
  TLegend leg(0.15, 0.8, 0.48, 0.9);
  leg.AddEntry(frame->findObject("Data"), "Data Hist", "lep");
  leg.AddEntry(frame->findObject("Fit"), "DCB Fit", "L");
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
  const char *save_end = "GeV_DCBfit_cutFlow_global_all.pdf";
  string path_plot = save_init + name + save_end;
  c1->SaveAs(path_plot.c_str());

  return;
}

void plot_param(double P_sing[4], double P_all[4], double Err_all[4], string name){

  TCanvas* c3 = new TCanvas("c3","c3",600, 400);
  double masses[4] = {110,125,130,140};
  TGraph* gr1 = new TGraph(4,masses,P_sing);
  gr1->SetName("Single_fit");
  gr1->SetTitle("Single_fit");
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerColor(2);
  gr1->SetDrawOption("APL");
  gr1->SetLineColor(2);
  gr1->SetLineWidth(4);
  gr1->SetFillStyle(0);

  TGraph* gr2 = new TGraphErrors(4,masses,P_all,Err_all);
  gr2->SetName("Sim_fit");
  gr2->SetTitle("Sim_fit");
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerColor(3);
  gr2->SetDrawOption("APL");
  gr2->SetLineColor(3);
  gr2->SetLineWidth(4);
  gr2->SetFillStyle(0);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr2);
  const char *title = "Single and global fit comp: ";
  string plot_title = title + name;
  mg->SetTitle(plot_title.c_str());
  mg->GetXaxis()->SetTitle("m_{h} [GeV]");
  mg->GetYaxis()->SetTitle("Params");
  mg->Draw("ALP");

  c3->BuildLegend();

  const char *repo = "Plots/";
  const char *end = "_param_all.pdf";
  string path = repo + name + end;

  c3->SaveAs(path.c_str());
}

RooDataSet *read_data(int MASS, RooRealVar& m){
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

  // RooRealVar m("m_yy","mass",MASS-intervall,MASS+intervall);
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
